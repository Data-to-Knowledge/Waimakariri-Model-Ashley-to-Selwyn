"""
Author: matth
Date Created: 22/06/2017 3:53 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import os
from warnings import warn
import pickle
import flopy
from matplotlib.colors import from_levels_and_colors
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# todo set up exceptions is required varibles are needed
# todo look though new modeling and add any useful tools to this class e.g convergence check or something around bc overlap
# todo document the shit out of this
class ModelTools(object):
    ulx = None
    uly = None
    grid_space = None
    layers = None
    rows = None
    cols = None
    sdp = None
    _no_flow_calc = None
    temp_file_dir = None
    model_version_name = None
    base_mod_path = None
    _elv_calculator = None

    def __init__(self, model_version_name, sdp=None, ulx=None, uly=None, rotation=0,layers=None, rows=None, cols=None,
                 grid_space=None, no_flow_calc=None, temp_file_dir=None, elv_calculator=None, base_mod_path=None,base_map_path=None):
        """
        rotation and differential grid spacing not supported
        :param model_version_name: required the model version name
        :param sdp: the path to the folder containing supporting data (for ease of reference)
        :param ulx: upper left x coordinate for the model grid
        :param uly:  the upper left y coordinate for the model grid
        :param rotation: the rotation of the grid
        :param layers: the number of layers in the model
        :param rows: the number of rows in the model
        :param cols: the number of columns in the model
        :param grid_space: the grid spacing (must be equally spaced at present)
        :param no_flow_calc: the function which is used to return the ibound
        :param temp_file_dir: a path to a directory to dump temporary files (mostly for shapefile conversion)
        :param elv_calculator: the function whcih is used to return the elevation data base (top,bottoms) shape (k,i,j)
        :param base_mod_path: the path to the base model (depreciated
        :param base_map_path: the path to a geotiff to allow underlaying it below plots
        """
        #todo add a flopy.utils.SpatialReference object as an attriubte to this class with try/except
        self.rotation = rotation
        self.ulx = ulx
        self.uly = uly
        self.grid_space = grid_space
        self.layers = layers
        self.rows = rows
        self.cols = cols
        self.sdp = sdp
        if not os.path.exists(sdp):
            os.makedirs(sdp)
        self._no_flow_calc = no_flow_calc
        self.temp_file_dir = temp_file_dir
        self.model_version_name = model_version_name
        self.pickle_dir = '{}/pickled_files'.format(self.sdp)
        if not os.path.exists(self.pickle_dir):
            os.makedirs(self.pickle_dir)
        self.base_mod_path = base_mod_path
        self._elv_calculator = elv_calculator
        self.model_array_shape = (self.layers, self.rows, self.cols)
        self.base_map_path = base_map_path

    def _get_xlim_ylim(self):
        """
        returns the x and y lim of the grid
        :return:
        """
        x_min, y_min = self.ulx, self.uly - self.grid_space * self.rows
        x_max = x_min + self.grid_space*self.cols
        y_max = self.uly
        return (x_min,x_max), (y_min,y_max)

    def get_model_x_y(self,grid=True):
        """
        returns the x and y coordinates
        :param grid: if True return the grid, else return the 1d arrays
        :return:
        """
        pixelWidth = pixelHeight = self.grid_space

        x_min, y_min = self.ulx, self.uly - self.grid_space * self.rows

        # mid points of the array
        x = np.linspace(x_min + 100, (x_min + 100 + pixelWidth * (self.cols - 1)), self.cols)
        y = np.linspace((y_min + 100 + pixelHeight * (self.rows - 1)), y_min + 100, self.rows)
        if not grid:
            return x,y
        model_xs, model_ys = np.meshgrid(x, y)
        return model_xs, model_ys

    def get_empty_model_grid(self, _3d=False):
        """
        returns a numpy zeros array of the model grid
        :param _3d: Boolean if True return a grid of (k,i,j) else (i,j)
        :return:
        """
        if _3d:
            return np.zeros((self.layers, self.rows, self.cols))
        else:
            return np.zeros((self.rows, self.cols))

    def add_mxmy_to_df(self, df):  # todo make faster
        """
        adds the x and y centroid of the cells to a dataframe with either row,col or i,j in the keys
        :param df:
        :return:
        """
        df = deepcopy(df)
        if 'mx' in df.keys() or 'my' in df.keys():
            warn('mx or my present in dataframe, this will be overwritten')

        try:
            for i in df.index:
                row, col = df.loc[i, ['row', 'col']]
                x, y = self.convert_matrix_to_coords(row, col)
                df.loc[i, 'mx'] = x
                df.loc[i, 'my'] = y
        except KeyError:
            for i in df.index:
                row, col = df.loc[i, ['i', 'j']]
                x, y = self.convert_matrix_to_coords(row, col)
                df.loc[i, 'mx'] = x
                df.loc[i, 'my'] = y
        return df

    def model_where(self, condition):
        """
        quick function to return list of model indexes from a np.where style condition assumes that the array in the
        condition has the shape (k,i,j,) or (i,j)
        :param condition:
        :return: return tuples of i,j or kij
        """
        idx = np.where(condition)
        idx2 = [i for i in zip(*idx)]
        return idx2

    def df_to_array(self, dataframe, value_to_map, _3d=False):
        """
        maps a dataframe to a model array
        :param dataframe: dataframe with at least keys (k),i,j and value
        :param value_to_map: key to the value to map onto the array
        :param _3d: boolean if true map to a 3d array with kIJ
        :return: array
        """
        if _3d:
            layval, rowval, colval = 'k', 'i', 'j'
            if value_to_map not in dataframe.keys():
                raise ValueError('value to map: {}, not found'.format(value_to_map))
            if not np.in1d(['k', 'i', 'j'], dataframe.keys()).all():
                if np.in1d(['layer', 'row', 'col'], dataframe.keys()).all():
                    layval, rowval, colval = 'layer', 'row', 'col'
                else:
                    raise ValueError('dataframe missing keys')
            outdata = np.zeros((self.layers, self.rows, self.cols)) * np.nan
            for i in dataframe.index:
                layer, row, col = dataframe.loc[i, [layval, rowval, colval]].astype(int)
                outdata[int(layer), int(row), int(col)] = dataframe.loc[i, value_to_map]
        else:
            rowval, colval = 'i', 'j'
            if not np.in1d(['i', 'j', value_to_map], dataframe.keys()).all():
                if np.in1d(['row', 'col', value_to_map], dataframe.keys()).all():
                    rowval, colval = 'row', 'col'
                else:
                    raise ValueError('dataframe missing keys')
            outdata = np.zeros((self.rows, self.cols)) * np.nan
            for i in dataframe.index:
                row, col = dataframe.loc[i, [rowval, colval]]
                outdata[int(row), int(col)] = dataframe.loc[i, value_to_map]

        return outdata

    def convert_matrix_to_coords(self, row, col, layer=None, elv_db=None):
        """
        convert from matix indexing to real world coordinates
        :param row:
        :param col:
        :param layer:
        :param elv_db: 3d numpy array of bottom elevation the 0th index is the top layer of the surface of the model
        :return: lon, lat if layer is not specified or lon, lat, elv
        """
        if elv_db is None and layer is not None:
            # below is a pickled nparray of bottoms (with layer 1 top on the top of the stack
            elv_db = self.calc_elv_db()

        col = int(col)
        row = int(row)
        model_xs, model_ys = self.get_model_x_y()
        lon = model_xs[0, col]
        lat = model_ys[row, 0]

        if layer is None:
            return lon, lat
        else:
            elv = elv_db[layer:layer + 2, row, col].mean()
            return lon, lat, elv

    def convert_coords_to_matix(self, lon, lat, elv=None, elv_db=None,
                                return_AE=False):  # todo can this be made quicker vectorize this shit...
        """
        convert from real world coordinates to model indexes by comparing value to center point of array.  Where the value
        is on an edge defaults to top left corner (e.g. row 1, col 1, layer 1)
        :param lon: longitude NZTM
        :param lat: latitude NZTM
        :param elv: elevation Lyttleton 43 MSL #not sure
        :param elv_db: 3d numpy array of bottom elevation the 0th index is the top layer of the surface of the model
        :param return_AE: if true return (layer,z0), (row, y0), (col, x0) where x0 ect are the preportion along the cell
                          for use with modpath particles
        :return: (row, col) if elv is not presented or (layer, row, col)
        """
        model_xs, model_ys = self.get_model_x_y()
        if lon > model_xs[0, :].max() or lon < model_xs[0, :].min():
            raise AssertionError('coords not in model domain')

        if lat > model_ys[:, 0].max() or lat < model_ys[:, :].min():
            raise AssertionError('coords not in model domain')

        layer, row, col = None, None, None
        if elv_db is None and elv is not None:
            # below is a pickled np.array of bottoms (with layer 1 top on the top of the stack)
            elv_db = self.calc_elv_db()
        # find the index of the closest middle point
        # convert lon to col
        col = np.abs(model_xs[0, :] - lon).argmin()  # todo see if this is identical with the gdal process

        # convert lat to row
        row = np.abs(model_ys[:, 0] - lat).argmin()

        if return_AE:
            x0 = (self.grid_space/2 - (model_xs[0, :] - lon)[col]) / self.grid_space
            y0 = (self.grid_space/2 - (model_ys[:, 0] - lat)[row]) / self.grid_space

        if elv is None:
            if return_AE:
                return (row, y0), (col, x0)
            else:
                return row, col
        else:
            # convert elv to layer
            top = elv_db[:-1, row, col]
            bottom = elv_db[1:, row, col]
            if elv > top.max():
                warn('above top of database setting layer to zero')
                layer = 0
                z0 = 0
            else:
                layer = np.where((top > elv) & (bottom <= elv))
                if len(layer[0]) > 1:
                    raise ValueError(
                        'returns multiple indexes for layer')
                elif len(layer[0]) == 0:
                    raise AssertionError('elevation: {} outside of bounds'.format(elv))
                layer = layer[0][0]
                if return_AE:
                    z0 = (top[layer] - elv) / (top[layer] - bottom[layer])
            if return_AE:
                return (layer, z0), (row, y0), (col, x0)
            else:
                return layer, row, col

    def convert_elv_to_k(self, row, col, elv, elv_db=None, return_AE=False):
        """
        as per convert coords to matrix, but row, col already known
        :param row:
        :param col:
        :param elv:
        :param elv_db:
        :param return_AE:
        :return:
        """
        if elv_db is None:
            # below is a pickled np.array of bottoms (with layer 1 top on the top of the stack)
            elv_db = self.calc_elv_db()
        top = elv_db[:-1, row, col]
        bottom = elv_db[1:, row, col]
        if elv > top.max():
            warn('above top of database setting layer to zero')
            layer = 0
            z0 = 0
        else:
            layer = np.where((top > elv) & (bottom <= elv))
            if len(layer[0]) > 1:
                raise ValueError(
                    'returns multiple indexes for layer')
            elif len(layer[0]) == 0:
                raise AssertionError('elevation: {} outside of bounds'.format(elv))
            layer = layer[0][0]
            if return_AE:
                z0 = (top[layer] - elv) / (top[layer] - bottom[layer])
        if return_AE:
            return layer, z0
        else:
            return layer

    def get_all_adjacent_cells(self, cell_locs, _3d=False, return_loc=True):
        """
        returns indexs for all adjacent cells in list of cell location
        :param cell_locs: list of cell locations(k,i,j) or (i,j) see below
        :param _3d: if true requires (k,i,j) and returns all in 3d else returns all in 2d requires input of (i,j)
        :param return_loc: bool if True include the indicies for the cell
        :return:
        """

        cell_locs = np.atleast_2d(cell_locs)

        out_locs = []
        for loc in cell_locs:
            if return_loc:
                out_locs.append(loc)

            if _3d:
                k = loc[0]
                i = loc[1]
                j = loc[2]
                if k != self.layers - 1:
                    out_locs.append((k + 1, i, j))
                if k != 0:
                    out_locs.append((k - 1, i, j))
                if i != self.rows - 1:
                    out_locs.append((k, i + 1, j))
                if i != 0:
                    out_locs.append((k, i - 1, j))
                if j != self.cols - 1:
                    out_locs.append((k, i, j + 1))
                if j != 0:
                    out_locs.append((k, i, j - 1))
            else:
                i = loc[0]
                j = loc[1]
                if i != self.rows - 1:
                    out_locs.append((i + 1, j))
                if i != 0:
                    out_locs.append((i - 1, j))
                if j != self.cols - 1:
                    out_locs.append((i, j + 1))
                if j != 0:
                    out_locs.append((i, j - 1))

        return out_locs

    def array_to_raster(self, path, array, no_flow_layer=None):
        """
        saves a 2d model array as a raster
        :param path: path to save the raster
        :param array: array to save (must be 2d model array)
        :param no_flow_layer: None or int the no-flow layer to mask the array on
        :return:
        """
        from osgeo import osr, gdal
        if array.shape != (self.rows,self.cols):
            raise ValueError('array must match model grid shape')
        if no_flow_layer is not None:
            no_flow = self.get_no_flow(no_flow_layer).astype(bool)
            array[~no_flow] = -99
        output_raster = gdal.GetDriverByName('GTiff').Create(path, array.shape[1], array.shape[0], 1,
                                                             gdal.GDT_Float32)  # Open the file
        geotransform = (self.ulx, self.grid_space, 0, self.uly, 0, -self.grid_space)
        output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
        srs = osr.SpatialReference()  # Establish its coordinate encoding
        srs.ImportFromEPSG(2193)  # This one specifies NZTM.
        output_raster.SetProjection(srs.ExportToWkt())  # Exports the coordinate system
        # to the file
        band = output_raster.GetRasterBand(1)
        band.WriteArray(array)  # Writes my array to the raster
        band.FlushCache()
        band.SetNoDataValue(-99)

    def plt_matrix(self, array, vmin=None, vmax=None, title=None, no_flow_layer=0, ax=None, color_bar=True,
                   base_map=False, plt_background=True, cbar_lab=None, **kwargs):
        """

        :param array: they array to plot (of i,j)
        :param vmin: vmin to pass to the plot
        :param vmax: vmax to pass to the plot
        :param title: a title to include on the plot
        :param no_flow_layer: the layer to plot the no flow boundry on
        :param ax: None or a matplotlib ax to plot on top of, not frequently used (execpt in subplots)
        :param color_bar: Boolean if true plot a color bar scale
        :param base_map: Boolean if True underlie the plot with greyscale basemap of the region the map is defined by
                         self.base_map_pathsets default alpha to 0.5
        :param plt_background: Boolean if True plot the noflow boundary as black
        :param cbar_lab: string to lable the cbar
        :param kwargs: other kwargs passed to matplotlib.pcolor alpha is the only kwarg that is also passed to the background
        :return: fig, ax
        """
        alpha = 1
        array = deepcopy(array.astype(float))
        if vmax is None:
            vmax = np.nanmax(array)
        if vmin is None:
            vmin = np.nanmin(array)

        if ax is None:
            fig, (ax) = plt.subplots(figsize=(18.5, 9.5))
            if title is not None:
                ax.set_title(title)
            ax.set_aspect('equal')
        else:
            fig = ax.figure
            ax.set_aspect('equal')
        model_xs, model_ys = self.get_model_x_y()
        if base_map:
            if 'alpha' not in kwargs:
                alpha = 0.5
            if self.base_map_path is None:
                raise ValueError('in order to use base_map self.base_map_path must be defined')
            from osgeo.gdal import Open

            ds = Open(self.base_map_path)
            width = ds.RasterXSize
            height = ds.RasterYSize
            gt = ds.GetGeoTransform()
            minx = gt[0]
            miny = gt[3] + width * gt[4] + height * gt[5]
            maxx = gt[0] + width * gt[1] + height * gt[2]
            maxy = gt[3]

            image = ds.ReadAsArray()
            if image.ndim == 3: # if a rgb image then plot as greyscale
                image = image.mean(axis=0)
            ll = (minx, miny)
            ur = (maxx, maxy)

            ax.imshow(image, extent=[ll[0], ur[0], ll[1], ur[1]], cmap='gray') #
            xlim, ylim = self._get_xlim_ylim()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
            kwargs.pop('alpha')
        if alpha == 1:
            linewidth = None
            edgecolors = 'face'
        else:
            edgecolors = 'face'
            linewidth = 0
        if self._no_flow_calc is not None and no_flow_layer is not None:
            no_flow = self.get_no_flow(no_flow_layer).astype(bool)
            ax.contour(model_xs, model_ys, no_flow)
            temp = np.zeros((self.rows, self.cols))
            temp[no_flow] = np.nan
            if plt_background:
                cmap, norm = from_levels_and_colors([-1, -0.5, 1, 2], ['black', 'black', 'black']) #todo could pcolor be changed to imshow
                ax.pcolor(model_xs, model_ys, np.ma.masked_invalid(temp), cmap=cmap, alpha=alpha, edgecolors=edgecolors)
            array[~no_flow] = np.nan
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
            kwargs.pop('cmap')
        else:
            cmap = 'plasma'
        pcm = ax.pcolor(model_xs, model_ys, np.ma.masked_invalid(array),
                            cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha, edgecolors=edgecolors, linewidth=linewidth, **kwargs)
        if color_bar:
            cbar = fig.colorbar(pcm, ax=ax, extend='max')
            if cbar_lab is not None:
                cbar.ax.set_ylabel(cbar_lab, rotation=270)
            cbar.set_alpha(1)
            cbar.draw_all()


        return fig, ax

    def geodb_to_model_array(self, path, shape_name, attribute, alltouched=False,
                             area_statistics=False, fine_spacing=10, resample_method='average'
                             ):
        """
        create model array from a shapefile in a geodatabase
        :param path: path to geodatabase
        :param shape_name: name of the shape file
        :param attribute: name of the attribute to convert to matrix
        :param area_statistics: boolean if true first burn the raster at a finer detail and do some statistics to
                                group up
        :param fine_spacing: int m of spacing to burn the raster to (limits the polygon overlay),
                             only used if area_statistics is True
        :param resample_method: key to define resampling options
                                near:
                                    nearest neighbour resampling (default, fastest algorithm, worst interpolation
                                    quality).
                                bilinear:
                                    bilinear resampling.
                                cubic:
                                    cubic resampling.
                                cubicspline:
                                    cubic spline resampling.
                                lanczos:
                                    Lanczos windowed sinc resampling.
                                average:
                                    average resampling, computes the average of all non-NODATA contributing pixels.
                                     (GDAL >= 1.10.0)
                                mode:
                                    mode resampling, selects the value which appears most often of all the
                                    sampled points. (GDAL >= 1.10.0)
                                max:
                                    maximum resampling, selects the maximum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                min:
                                    minimum resampling, selects the minimum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                med:
                                    median resampling, selects the median value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q1:
                                    first quartile resampling, selects the first quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q3:
                                    third quartile resampling, selects the third quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)

        :return: np array of the matrix
        """
        from osgeo import ogr
        driver = ogr.GetDriverByName("OpenFileGDB")
        ds = driver.Open(path, 0)
        source_layer = ds.GetLayer(shape_name)

        outdata = self._layer_to_model_array(source_layer, attribute, alltouched=alltouched,
                                             area_statistics=area_statistics, fine_spacing=fine_spacing,
                                             resample_method=resample_method)
        return outdata

    def shape_file_to_model_array(self, path, attribute, alltouched=False,
                                  area_statistics=False, fine_spacing=10, resample_method='average'):
        """
        shape file to vistas array
        :param path: path to shapefile
        :param attribute: attribute name to convert
        :param area_statistics: boolean if true first burn the raster at a finer detail and do some statistics to
                                group up
        :param fine_spacing: int m of spacing to burn the raster to (limits the polygon overlay),
                             only used if area_statistics is True
        :param resample_method: key to define resampling options
                                near:
                                    nearest neighbour resampling (default, fastest algorithm, worst interpolation
                                    quality).
                                bilinear:
                                    bilinear resampling.
                                cubic:
                                    cubic resampling.
                                cubicspline:
                                    cubic spline resampling.
                                lanczos:
                                    Lanczos windowed sinc resampling.
                                average:
                                    average resampling, computes the average of all non-NODATA contributing pixels.
                                     (GDAL >= 1.10.0)
                                mode:
                                    mode resampling, selects the value which appears most often of all the
                                    sampled points. (GDAL >= 1.10.0)
                                max:
                                    maximum resampling, selects the maximum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                min:
                                    minimum resampling, selects the minimum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                med:
                                    median resampling, selects the median value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q1:
                                    first quartile resampling, selects the first quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q3:
                                    third quartile resampling, selects the third quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)

        :return:  np.array of the data in model format
        """
        # open shape file
        from osgeo import ogr
        source_ds = ogr.Open(path)
        source_layer = source_ds.GetLayer()

        outdata = self._layer_to_model_array(source_layer, attribute, alltouched=alltouched,
                                             area_statistics=area_statistics, fine_spacing=fine_spacing,
                                             resample_method=resample_method)
        return outdata

    def _layer_to_model_array(self, source_layer, attribute, alltouched=False,
                              area_statistics=False, fine_spacing=10, resample_method='average'): #todo propogate these to upper level
        """
        hidden function to convert a source layer to a rasterized np array
        :param source_layer: from either function above
        :param attribute: attribute to convert.
        :param area_statistics: boolean if true first burn the raster at a finer detail and do some statistics to
                                group up
        :param fine_spacing: int m of spacing to burn the raster to (limits the polygon overlay),
                             only used if area_statistics is True
        :param resample_method: key to define resampling options
                                near:
                                    nearest neighbour resampling (default, fastest algorithm, worst interpolation
                                    quality).
                                bilinear:
                                    bilinear resampling.
                                cubic:
                                    cubic resampling.
                                cubicspline:
                                    cubic spline resampling.
                                lanczos:
                                    Lanczos windowed sinc resampling.
                                average:
                                    average resampling, computes the average of all non-NODATA contributing pixels.
                                     (GDAL >= 1.10.0)
                                mode:
                                    mode resampling, selects the value which appears most often of all the
                                    sampled points. (GDAL >= 1.10.0)
                                max:
                                    maximum resampling, selects the maximum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                min:
                                    minimum resampling, selects the minimum value from all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                med:
                                    median resampling, selects the median value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q1:
                                    first quartile resampling, selects the first quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
                                q3:
                                    third quartile resampling, selects the third quartile value of all non-NODATA
                                    contributing pixels. (GDAL >= 2.0.0)
        :return:  np array of rasterized data
        """
        from osgeo import gdal, osr #todo see if I can use tempfile for this
        if not os.path.exists(self.temp_file_dir):
            os.makedirs(self.temp_file_dir)

        if area_statistics:
            assert self.grid_space % fine_spacing == 0, 'fine spacing {} must be an ' \
                                                        'even factor of gridspacing: {}'.format(fine_spacing,
                                                                                                self.grid_space)
            trans = self.grid_space/fine_spacing
            cols = int(self.cols * trans)
            rows = int(self.rows * trans)
            pixelWidth = pixelHeight = fine_spacing  # depending how fine you want your raster
        else:
            cols = self.cols
            rows = self.rows
            pixelWidth = pixelHeight = self.grid_space  # depending how fine you want your raster

        x_min, y_min = self.ulx, self.uly - self.grid_space * self.rows

        target_ds = gdal.GetDriverByName('GTiff').Create('{}/temp.tif'.format(self.temp_file_dir), cols, rows,
                                                         1,
                                                         gdal.GDT_Float64)
        target_ds.SetGeoTransform((x_min, pixelWidth, 0, y_min, 0, pixelHeight))
        band = target_ds.GetRasterBand(1)
        NoData_value = -999999
        band.SetNoDataValue(NoData_value)
        band.FlushCache()
        if alltouched:
            opt = ["ALL_TOUCHED=TRUE", "ATTRIBUTE={}".format(attribute)]
        else:
            opt = ["ATTRIBUTE={}".format(attribute)]
        gdal.RasterizeLayer(target_ds, [1], source_layer, options=opt)

        target_dsSRS = osr.SpatialReference()
        target_dsSRS.ImportFromEPSG(2193)
        target_ds.SetProjection(target_dsSRS.ExportToWkt())
        target_ds = None
        temp_file = '{}/temp.tif'.format(self.temp_file_dir)
        if area_statistics:
            gdal.Warp('{}/temp2.tif'.format(self.temp_file_dir), temp_file,
                      xRes=self.grid_space, yRes=self.grid_space, resampleAlg=resample_method)
            os.remove(temp_file)
            temp_file = '{}/temp2.tif'.format(self.temp_file_dir)

        outdata = gdal.Open(temp_file).ReadAsArray()
        if not area_statistics:
            outdata = np.flipud(outdata)
        outdata[np.isclose(outdata, -999999)] = np.nan
        os.remove(temp_file)

        return outdata

    def get_well_postions(self, well_nums, screen_handling='middle', raise_exct=True, error_log_path=None,
                          one_val_per_well=False):
        """
        gets model indexs from a well position
        :param well_nums: value or iterable of well numbers
        :param screen_handling: how to calculate elevation from wells with screen details
                                middle: returns all unique layers for the set
                                all: returns all unique layers intersected by screens(e.g. from top and bottom elevation)
        :param raise_exct: raise exceptions or continue, see error_log_path for handeling if False
        :param error_log_path: path to error log where values are written if raise_exct is True.  if None, exceptions passed
                                as warnings
        :param one_val_per_well: if True, lists of layer values will be condensed and the most frequent layer value will be
                                 returned ensuring there is only one elevation per well. Where there is a tie the lowest
                                 index layer (e.g. the smallest number) of the tied values will be returned
        :return:
        """
        from core.ecan_io import rd_sql, sql_db
        import pandas as pd

        if error_log_path is not None:
            f = open(error_log_path, 'w')

        # check inputs
        if screen_handling not in ['middle', 'all']:
            raise ValueError('screen_handling not appropriately specified')

        if isinstance(well_nums, str):
            return_value = True
        else:
            return_value = False

        # load supporting information
        # below is a pickled nparray of bottoms (with layer 1 top on the top of the stack
        elv_db = self.calc_elv_db()
        well_details = rd_sql(**sql_db.wells_db.well_details)
        well_details = well_details.set_index('WELL_NO')
        screen_details = rd_sql(**sql_db.wells_db.screen_details)

        well_nums = np.atleast_1d(well_nums)

        if well_nums.ndim != 1:
            raise ValueError('unexpected input dimensions')

        # initalize output arrays
        layer = np.zeros(well_nums.shape, dtype=object) * np.nan
        row = np.zeros(well_nums.shape, dtype=object) * np.nan
        col = np.zeros(well_nums.shape, dtype=object) * np.nan
        for i, well in enumerate(well_nums):
            try:
                # first get row and column
                ref_level = well_details.loc[well, 'REFERENCE_RL']
                lat = well_details.loc[well, 'NZTMY']
                lon = well_details.loc[well, 'NZTMX']
                ground_level = well_details.loc[well, 'GROUND_RL'] + ref_level
                if pd.isnull(ground_level):
                    rt, ct = self.convert_coords_to_matix(lon, lat)
                    ground_level = elv_db[0, rt, ct]  # take as top of model

                screen_num = well_details.loc[well, 'Screens']
                depth = ground_level - well_details.loc[well, 'DEPTH']  # this is from ground level

                # set elevation
                elv = None
                if screen_num == 0:
                    elv = depth  # also from ground level
                    layer_temp, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, elv, elv_db)
                elif screen_num == 1:
                    top = ground_level - screen_details['TOP_SCREEN'][screen_details['WELL_NO'] == well]
                    bottom = ground_level - screen_details['BOTTOM_SCREEN'][screen_details['WELL_NO'] == well]
                    if len(top) != 1 or len(bottom) != 1:
                        raise ValueError('well: {} incorrect number of screen values'.format(well))
                    if screen_handling == 'middle':
                        elv = (top.iloc[0] + bottom.iloc[0]) / 2  # also from ground level
                        layer_temp, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, elv, elv_db)
                    elif screen_handling == 'all':
                        layer_temp = []
                        lt1, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, top, elv_db)
                        lt2, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, bottom, elv_db)
                        layer_temp.extend(range(lt1, lt2 + 1))  # assume the screen is continous
                else:
                    layer_temp = []
                    for j in range(1, screen_num + 1):
                        top = ground_level - screen_details['TOP_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                                          (screen_details['SCREEN_NO'] == j)]
                        bottom = ground_level - screen_details['BOTTOM_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                                                (screen_details['SCREEN_NO'] == j)]
                        if len(top) != 1 or len(bottom) != 1:
                            raise ValueError('well: {} incorrect number of screen values for screen {}'.format(well, j))
                        if screen_handling == 'middle':
                            elv = (top.iloc[0] + bottom.iloc[0]) / 2  # also from ground level
                            lt, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, elv, elv_db)
                            layer_temp.append(lt)
                        elif screen_handling == 'all':
                            lt1, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, top, elv_db)
                            lt2, row_temp, col_temp = self.convert_coords_to_matix(lon, lat, bottom, elv_db)
                            layer_temp.extend(range(lt1, lt2 + 1))  # assume the screen is continous

                if pd.isnull(np.atleast_1d(layer_temp)).sum() != 0:
                    raise ValueError('well: {} does not have screen, reference level, '
                                     'ground level, or depth data'.format(well))

                if isinstance(layer_temp, list):
                    if one_val_per_well:
                        layer_temp = max(set(layer_temp),
                                         key=layer_temp.count)  # todo it may be preferable to calculate a mid screen
                    else:
                        layer_temp = list(set(layer_temp))
                layer[i] = layer_temp
                row[i] = row_temp
                col[i] = col_temp
            except Exception as val:
                if raise_exct:
                    raise Exception('{} {}'.format(well, val))
                elif error_log_path is None:
                    warn('{} {}'.format(well, val))
                else:
                    f.write('{} {}\n'.format(well, val))

        if raise_exct:
            # check for more missing data
            if pd.isnull(layer).sum() != 0:
                raise ValueError('wells {} do not have data in layer'.format(well_nums[pd.isnull(layer)]))
            elif pd.isnull(row).sum() != 0:
                raise ValueError('wells {} do not have data in row'.format(well_nums[pd.isnull(row)]))
            elif pd.isnull(col).sum() != 0:
                raise ValueError('wells {} do not have data in col'.format(well_nums[pd.isnull(col)]))

        if return_value:
            return layer[0], row[0], col[0]
        else:
            return layer, row, col

    def calc_elv_db(self, recalc=False):  # todo make this stored in the class so that it isn't loaded multple times
        """
        calculates the elevation database (bottoms with the top of layer one on top or loads pickel
        :param recalc: force recalculates the elv_db from the gns data even if it is not present
        :return: elv_db
        """
        pickle_path = '{}/elv_dp.p'.format(self.pickle_dir)
        if os.path.exists(pickle_path) and not recalc:
            elv_db = pickle.load(open(pickle_path))
        else:
            elv_db = self._elv_calculator()
            pickle.dump(elv_db, open(pickle_path, 'w'))
        return elv_db

    def get_no_flow(self, layer=None, recalc=False):
        pickle_path = '{}/no_flow.p'.format(self.pickle_dir)
        if os.path.exists(pickle_path) and not recalc:
            no_flow = pickle.load(open(pickle_path))
        else:
            no_flow = self._no_flow_calc()
            pickle.dump(no_flow, open(pickle_path, 'w'))
        if layer is None:
            return no_flow
        else:
            return no_flow[layer]

    def get_base_model(self, recalc=False):  # todo raise eception if model path not passed perhaps delete
        """
        This has not been used much perhaps depreciate
        :param recalc:
        :return:
        """
        import flopy
        m = flopy.modflow.Modflow.load(self.base_mod_path, model_ws=os.path.dirname(self.base_mod_path),
                                       forgive=False)
        return m

    def get_package_spd(self, package, recalc=False):  # todo delete
        """
        this had not been used much perhaps depreciate
        :param package:
        :param recalc:
        :return:
        """
        import pickle

        picklepath = '{}/base_{}.p'.format(self.pickle_dir, package)

        if os.path.exists(picklepath) and not recalc:
            base_data = pickle.load(open(picklepath))
            return base_data

        org_m = self.get_base_model()
        if package.lower in ['str', 'sfr', 'drn', 'wel']:
            base_data = org_m.get_package().stress_period_data.data[0]
        elif package.lower == 'rch':
            raise NotImplementedError()  # todo
        else:
            raise ValueError('unexpected package {}, may not be implemented'.format(package))
        pickle.dump(base_data, open(picklepath, 'w'))

        return base_data

        # #todo think about adding get stress period data... input the package names and pickle the file if needed so you don't always have to load the thing

    # todo eventually handle the layer,col,row vs kij keys
    def convert_well_data_to_stresspd(self, well_data_in):
        """
        take a dataframe aggregates the wells and returns the well stress period data for the well package
        :param well_data_in: pd.Dataframe needs columns layer, row, col, flux
        :return:
        """
        # convert a dataframe of well features (x,y,z,flux,type, etc.) to well standard stress period data
        # do something similar for concentration data?

        # do groupby statistics to sum everything that is in the same layer, col, row
        g = well_data_in.groupby(['layer', 'row', 'col'])
        well_ag = g.aggregate({'flux': np.sum}).reset_index()

        outdata = well_ag.loc[:, ['layer', 'row', 'col', 'flux']]
        outdata = outdata.rename(columns={'layer': 'k', 'row': 'i', 'col': 'j'}).to_records(False)
        outdata = outdata.astype(flopy.modflow.ModflowWel.get_default_dtype())
        return outdata

    # todo think about adding default model_maps


    # todo add base map and cross section
    # todo add elvation check
    # todo add recalc_pickles

    def check_layer_overlap(self, required_overlap=0.50, top=None, bot=None, use_elv_db=False, layer=None,
                            return_min=False):
        """
        check that there is at least a certain ammount of overlap between cells in a layer \
        (must match model discritation) users must specifiy either (top and bot) or (use_elv_db and layer)
        :param required_overlap: fraction of over lap to return if retuning boolean array
        :param top: top array
        :param bot: bottom array
        :param use_elv_db: if true use the elevation database
        :param layer: the layer to assess if using the elvivation database
        :param return_min: if true retun the minimum overlap
        :return: either a boolean array or teh array of the minimum fraction overlap
        """
        if (top is None or bot is None) and not (use_elv_db and layer is not None):
            raise ValueError('more arguments need to be specified: requiers either use_elv_db and layer or top and bot')

        if use_elv_db:
            elv = self.calc_elv_db()
            top = elv[layer]
            bot = elv[layer + 1]

        if top.shape != (self.rows, self.cols) or bot.shape != (self.rows, self.cols):
            raise ValueError('top or bottom shape does not match model discritisation')
        outdata = np.zeros((self.rows, self.cols))
        for i in range(self.rows):
            for j in range(self.cols):
                cell_top = top[i, j]
                cell_bot = bot[i, j]
                thickness = cell_top - cell_bot
                other_cells = self.get_all_adjacent_cells((i, j), return_loc=False)
                other_tops = np.array([min(top[tuple(e)], cell_top) for e in other_cells])
                other_bots = np.array([max(bot[tuple(e)], cell_bot) for e in other_cells])
                other_thick_frac = (other_tops - other_bots) / thickness
                outdata[i, j] = other_thick_frac.min()

        if not return_min:
            outdata = outdata < return_min

        return outdata
