"""
 Author: Matt Hanson
 Created: 13/04/2020 10:40 AM

 this script is used to create a GIS directory so that many of the inputs can be easily viewed
 """

import numpy as np
import pandas as pd
import netCDF4 as nc
import geopandas as gpd
from waimak_extended_boundry.extended_boundry_model_tools import smt
import os
from waimak_extended_boundry.model_run_tools.model_setup.realisation_id import get_stocastic_set
import matplotlib.pyplot as plt


def _save_to_array(temp_data, outdir, var, crop_by_noflow, inpath):
    no_flow_layer = None
    if temp_data.shape == (smt.rows, smt.cols):
        if crop_by_noflow:
            no_flow_layer = 0

        smt.array_to_raster(os.path.join(outdir, '{}.tif'.format(var)), temp_data, no_flow_layer=no_flow_layer)

    elif temp_data.shape == (smt.layers, smt.rows, smt.cols):
        for l in range(smt.layers):
            if crop_by_noflow:
                no_flow_layer = l
            smt.array_to_raster(os.path.join(outdir, '{}_layer_{}.tif'.format(var, l)),
                                temp_data[l],
                                no_flow_layer=no_flow_layer)

    else:
        print('skipping unexpected shape for {}: {} in {}'.format(var, temp_data.shape, inpath))


def _save_plot(temp_data, outdir, var, crop_by_noflow, inpath):
    temp_data = temp_data.astype(float)
    no_flow_layer = None
    if temp_data.shape == (smt.rows, smt.cols):
        if crop_by_noflow:
            no_flow_layer = 0
            ibnd = smt.get_no_flow(0)
            temp_data[ibnd == 0] = np.nan

        fig, ax = smt.plt_matrix(temp_data, no_flow_layer=no_flow_layer, base_map=True, title=var)
        fig.savefig(os.path.join(outdir, '{}.png'.format(var)))
        plt.close()

    elif temp_data.shape == (smt.layers, smt.rows, smt.cols):
        for l in range(smt.layers):
            if crop_by_noflow:
                no_flow_layer = l
                ibnd = smt.get_no_flow(l)
                temp_data[l][ibnd == 0] = np.nan
            fig, ax = smt.plt_matrix(temp_data[l],
                                     no_flow_layer=no_flow_layer, base_map=True, title='{}_layer_{}'.format(var, l))
            fig.savefig(os.path.join(outdir, '{}_layer_{}.png'.format(var, l)), )
            plt.close(fig)

    else:
        print('skipping unexpected shape for {}: {} in {}'.format(var, temp_data.shape, inpath))


def small_nc_to_array(inpath, base_outdir, variables=None, crop_by_noflow=False):
    """
    for arrays that can be held in memory
    :param inpath:
    :param base_outdir:
    :param variables:
    :return:
    """
    dataset = nc.Dataset(inpath)

    outdir = os.path.join(base_outdir, os.path.basename(inpath).replace('.nc', ''))
    outdir_plot = os.path.join(base_outdir, 'plots_{}'.format(os.path.basename(inpath).replace('.nc', '')))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir_plot):
        os.makedirs(outdir_plot)

    if variables is None:
        variables = set(dataset.variables.keys())

    with open(os.path.join(outdir, 'README.txt'), 'w') as f:
        f.write(repr(dataset))
        for v in variables:
            f.write('\n\n')
            f.write(repr(dataset.variables[v]))

    with open(os.path.join(outdir_plot, 'README.txt'), 'w') as f:
        f.write(repr(dataset))
        for v in variables:
            f.write('\n\n')
            f.write(repr(dataset.variables[v]))

    for var in variables:
        temp_data = np.array(dataset.variables[var])
        _save_to_array(temp_data, outdir, var, crop_by_noflow, inpath)
        _save_plot(temp_data, outdir_plot, var, crop_by_noflow, inpath)


def big_nc_to_array(inpath, base_outdir, variables=None, crop_by_noflow=False, nsmc_nums_nm='stochastic',
                    functions=(np.nanmedian, np.nanstd), functions_nm=('nanmedian', 'nanstd')):
    """
    for arrays that cannot be held in memory (iterate over layers)
    :param inpath:
    :param base_outdir:
    :param variables:
    :param crop_by_noflow:
    :param nsmc_nums_nm:
    :param functions:
    :param functions_nm:
    :return:
    """
    dataset = nc.Dataset(inpath)

    if nsmc_nums_nm == 'all':
        nsmc_nums = np.array(dataset.variables['nsmc_num'][:])
    elif nsmc_nums_nm == 'stochastic':
        nsmc_nums = get_stocastic_set(False)
    elif nsmc_nums_nm == 'optimised':
        nsmc_nums = [-1]
    else:
        raise ValueError('unexpected value for nsmc_nums_nm: {}'.format(nsmc_nums_nm))

    nidx = np.in1d(np.array(dataset.variables['nsmc_num'][:]), nsmc_nums)
    if not nidx.any():
        print('skipping {} {}'.format(nsmc_nums_nm, inpath))
        return None

    nsmc_nums = np.atleast_1d(nsmc_nums)
    functions = np.atleast_1d(functions)

    if len(nsmc_nums) == 1:
        functions = [None]
        functions_nm = ['none']

    for f, fn in zip(functions, functions_nm):
        outdir = os.path.join(base_outdir, '{}_{}_amalg_{}'.format(os.path.basename(inpath).replace('.nc', ''),
                                                                   nsmc_nums_nm,
                                                                   fn
                                                                   ))
        outdir_plots = os.path.join(base_outdir,
                                    'plots_{}_{}_amalg_{}'.format(os.path.basename(inpath).replace('.nc', ''),
                                                                  nsmc_nums_nm,
                                                                  fn
                                                                  ))

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(outdir_plots):
            os.makedirs(outdir_plots)

        if variables is None:
            variables = set(dataset.variables.keys())

        with open(os.path.join(outdir, 'README.txt'), 'w') as rdme:
            rdme.write(repr(dataset))
            for v in variables:
                rdme.write('\n\n')
                rdme.write(repr(dataset.variables[v]))
        with open(os.path.join(outdir_plots, 'README.txt'), 'w') as rdme:
            rdme.write(repr(dataset))
            for v in variables:
                rdme.write('\n\n')
                rdme.write(repr(dataset.variables[v]))

        for var in variables:
            print(var, fn, os.path.basename(inpath))
            if not dataset.variables[var].ndim == 3 and not dataset.variables[var].ndim == 4:
                print('skipping {} as weird dimentions {} in {}'.format(var, dataset.variables[var].dimensions, inpath))
                continue

            _3d = True
            if dataset.variables[var].ndim == 3:
                _3d = False
            temp_data = smt.get_empty_model_grid(_3d) * np.nan

            if not _3d:
                temp = np.array(dataset.variables[var][nidx])
                if f is None:
                    temp_data = temp
                else:
                    temp_data = f(temp, axis=0)

            for l in range(smt.layers):
                temp = np.array(dataset.variables[var][nidx, l])
                if f is None:
                    temp_data[l] = f
                else:
                    temp_data[l] = f(temp, axis=0)

            _save_to_array(temp_data, outdir, '{}'.format(var), crop_by_noflow, inpath)
            _save_plot(temp_data, outdir_plots, '{}'.format(var), crop_by_noflow, inpath)


def create_nc_datasets(outdir, indir=r"D:\Waimakariri_model_input_data"):  # todo test!!!!
    """
    handle the netcdfs big and small
    :param outdir:
    :param indir:
    :return:
    """

    small_ncs = [
        'required/chbs.nc',
        'required/elv_db.nc',
        'required/no_flow.nc',
        'required/non_head_target_locations.nc',
        'required/recharge_arrays.nc',

    ]
    small_nc_vars = [
        None,
        None,
        None,
        None,
        ['opt_rch']]

    for nc, var in zip(small_ncs, small_nc_vars):

        if 'no_flow' in nc:
            crop = False
        else:
            crop = True
        small_nc_to_array(os.path.join(indir, nc), outdir, var, crop)

    stocastic_ncs = [
        'recommended/lifestyle_fraction.nc',
        'recommended/sbda_lowdrain_fraction.nc',
        'recommended/sbda_highdrain_fraction.nc',
        'recommended/dairy_highdrain_fraction.nc',
        'recommended/dairy_lowdrain_fraction.nc',
        'recommended/GMP_mednload_ucn.nc',

    ]

    for f in stocastic_ncs:
        big_nc_to_array(os.path.join(indir, f), outdir, variables=None, crop_by_noflow=True, nsmc_nums_nm='stochastic',
                        functions=(np.nanmedian, np.nanstd), functions_nm=('nanmedian', 'nanstd'))

        big_nc_to_array(os.path.join(indir, f), outdir, variables=None, crop_by_noflow=True, nsmc_nums_nm='optimised')

    post_filter1_ncs = [
        'recommended/post_filter1_mednload_unc.nc',
        'recommended/post_filter1_hydraulic_properties.nc',
        'recommended/post_filter1_emma_unc_riv.nc',
        'recommended/post_filter1_hds.nc',
        'recommended/post_filter1_cell_budgets.nc',

    ]

    for f in post_filter1_ncs:
        big_nc_to_array(os.path.join(indir, f), outdir, variables=None, crop_by_noflow=True, nsmc_nums_nm='stochastic',
                        functions=(np.nanmedian, np.nanstd), functions_nm=('nanmedian', 'nanstd'))

        big_nc_to_array(os.path.join(indir, f), outdir, variables=None, crop_by_noflow=True, nsmc_nums_nm='all',
                        functions=(np.nanmedian, np.nanstd), functions_nm=('nanmedian', 'nanstd'))

        big_nc_to_array(os.path.join(indir, f), outdir, variables=None, crop_by_noflow=True, nsmc_nums_nm='optimised')


def create_from_hdfs(outdir, indir=r"D:\Waimakariri_model_input_data"):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outdir_plots = os.path.join(os.path.dirname(outdir),
                                'plots_{}'.format(os.path.basename(outdir)))
    if not os.path.exists(outdir_plots):
        os.makedirs(outdir_plots)

    rdme = open(os.path.join(outdir,'README.txt'), 'w')
    rdme.write('the following are the metadata for the shapefiles:\n')
    hdf_paths = [
        'required/base_well_data.hdf',
        'required/sfr_data.hdf',
        'required/all_wells_row_col_layer.hdf',
        'required/base_drn_data.hdf',

    ]

    for hdf in hdf_paths:
        hdf = os.path.join(indir,hdf)
        rdme.write('\n\n{}\n'.format(os.path.basename(hdf)))

        store = pd.HDFStore(hdf)

        keys = [e.strip('/') for e in store.keys()]
        metadata = {e: store.get_storer(e).attrs.metadata for e in keys}
        for k, v in metadata.items():
            rdme.write('   {}: {}\n'.format(k,v))

        store.close()

        for key in keys:
            print(os.path.basename(hdf), key)
            temp_data = pd.read_hdf(hdf, key)
            temp_data = temp_data.apply(pd.to_numeric, errors='ignore')

            try:
                temp_data.drop('consent',1, inplace=True)
            except:
                pass

            try:
                temp_data = gpd.GeoDataFrame(temp_data,
                                             geometry=gpd.points_from_xy(temp_data['mx'], temp_data['my']),
                                             crs='EPSG:2193')
                temp_data.to_file(os.path.join(outdir, '{}_{}.shp'.format(os.path.splitext(os.path.basename(hdf))[0],
                                                                         key)))
                fig, (ax_map, ax_text) = plt.subplots(ncols=2, figsize=(15, 10), constrained_layout=True)
                # make shapefile plot
                smt.plt_matrix(smt.get_empty_model_grid()*np.nan,base_map=True, ax=ax_map, color_bar=False)
                temp_data.plot()
                temp_data.plot(ax=ax_map, color='r')

                # plot the attribute data
                attributes = np.array(list(temp_data.keys()))
                num = attributes.shape[0] // 2 * 2
                attributes = attributes[0:num].reshape((num // 2, 2))

                ax_text.table(cellText=attributes, loc='upper center')
                ax_text.text(0.02,0.02,'metadata: {}'.format(metadata[key]['description']), wrap=True)

                ax_map.set_title('{}_{}'.format(os.path.splitext(os.path.basename(hdf))[0],
                                                                         key))
                ax_text.set_title('attributes')
                fig.set_constrained_layout(True)
                fig.savefig(os.path.join(outdir_plots, '{}_{}.png'.format(os.path.splitext(os.path.basename(hdf))[0],
                                                                         key)))
                plt.close(fig)

            except KeyError:
                continue




if __name__ == '__main__':
    import subprocess

    subprocess.call('powercfg -change standby-timeout-ac 0')  # stop it from sleeping...

    create_from_hdfs(r"D:\Waimakariri_model_input_data\gis_database\from_hdfs")
    print('creating from ncs')
    create_nc_datasets(r"D:\Waimakariri_model_input_data\gis_database\from_ncs")

    subprocess.call('powercfg -change standby-timeout-ac 30')  # and start it sleeping again
