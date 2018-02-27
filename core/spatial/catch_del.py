# -*- coding: utf-8 -*-
"""
Functions to delineate catchments.
"""
from pandas import DataFrame, read_csv, concat
from geopandas import read_file, GeoDataFrame, GeoSeries
from numpy import in1d, append, isnan
from core.ecan_io import rd_sql
from core.spatial.network import find_upstream_rec, extract_rec_catch, agg_rec_catch
from core.spatial.vector import closest_line_to_pts
from core.misc import select_sites


def catch_net(catch_sites_csv, catch_sites_col=['GRIDCODE', 'SITE']):
    """
    Function to create a dataframe of all the upstream catchments from specific locations. Input should be the output from ArcGIS script (2 columns of catchments and sites).
    """

    ## Read in data
    catch_sites_names=['catch', 'site']
    catch_sites = read_csv(catch_sites_csv)[catch_sites_col]
    catch_sites.columns = catch_sites_names

    ## Reorganize and select intial catchments
    catch_sites1 = catch_sites[catch_sites.site != 0]
    catch_sites2 = catch_sites1[catch_sites1.catch != catch_sites1.site]

    catch_unique = catch_sites2.catch.unique()

    singles = catch_sites1.catch[~catch_sites1.catch.duplicated(keep=False)]

    ## Network processing
    df = catch_sites2
    index1 = catch_unique
    catch_set2 = []
    for i in index1:
        catch1 = df.loc[df.catch == i, 'site'].values
        catch_set1 = catch1
        check1 = in1d(df.catch, catch1)
        while sum(check1) >= 1:
#            if sum(check1) > len(catch1):
#                print('Index numbering is wrong!')
            catch2 = df[check1].site.values.flatten()
            catch3 = catch2[~isnan(catch2)]
            catch_set1 = append(catch_set1, catch3)
            check1 = in1d(df.catch, catch3)
            catch1 = catch3
        catch_set2.append(catch_set1.tolist())

    df2 = DataFrame(catch_set2, index=index1)
    return([df2, singles.values])


def agg_catch(catch_del_shp, catch_sites_csv, catch_sites_col=['GRIDCODE', 'SITE'], catch_col='GRIDCODE'):
    """
    Function to take the output of the ArcGIS catchment delineation polygon shapefile and cathcment sites csv and return a shapefile with appropriately delineated polygons.
    """

    ## Catchment areas shp
    catch = read_file(catch_del_shp)[[catch_col, 'geometry']]

    ## dissolve the polygon
    catch3 = catch.dissolve(catch_col)

    ## Determine upstream catchments
    catch_df, singles_df = catch_net(catch_sites_csv, catch_sites_col)

    base1 = catch3[in1d(catch3.index, singles_df)].geometry
    for i in catch_df.index:
            t1 = append(catch_df.loc[i, :].dropna().values, i)
            t2 = GeoSeries(catch3[in1d(catch3.index, t1)].unary_union, index=[i])
            base1 = GeoSeries(concat([base1, t2]))

    ## Convert to GeoDataFrame (so that all functions can be applied to it)
    base2 = GeoDataFrame(base1.index, geometry=base1.geometry.values, crs=catch.crs)
    base2.columns = ['site', 'geometry']
    return(base2)


def rec_catch_del(sites_shp, sites_col='site', catch_output=None):
    """
    Catchment delineation using the REC streams and catchments.

    sites_shp -- Points shapfile of the sites along the streams.\n
    sites_col -- The column name of the site numbers in the sites_shp.\n
    catch_output -- The output polygon shapefile path of the catchment delineation.
    """

    ### Parameters
    server = 'SQL2012PROD05'
    db = 'GIS'
    streams_table = 'MFE_NZTM_REC'
    streams_cols = ['NZREACH', 'NZFNODE', 'NZTNODE']
    catch_table = 'MFE_NZTM_RECWATERSHEDCANTERBURY'
    catch_cols = ['NZREACH']

    ### Modifications {NZREACH: {NZTNODE/NZFNODE: node # to change}}
    mods = {13053151: {'NZTNODE': 13055874}, 13048353: {'NZTNODE': 13048851}, 13048498: {'NZTNODE': 13048851}}

    ### Load data
    rec_streams = rd_sql(server, db, streams_table, streams_cols, geo_col=True)
    rec_catch = rd_sql(server, db, catch_table, catch_cols, geo_col=True)
    pts = select_sites(sites_shp)

    ### make mods
    for i in mods:
        rec_streams.loc[rec_streams['NZREACH'] == i, mods[i].keys()] = mods[i].values()

    ### Find closest REC segment to points
    pts_seg = closest_line_to_pts(pts, rec_streams, line_site_col='NZREACH', dis=400)
    nzreach = pts_seg.copy().NZREACH.unique()

    ### Find all upstream reaches
    reaches = find_upstream_rec(nzreach, rec_shp=rec_streams)

    ### Extract associated catchments
    rec_catch = extract_rec_catch(reaches, rec_catch_shp=rec_catch)

    ### Aggregate individual catchments
    rec_shed = agg_rec_catch(rec_catch)
    rec_shed.columns = ['NZREACH', 'geometry', 'area']
    rec_shed1 = rec_shed.merge(pts_seg.drop('geometry', axis=1), on='NZREACH')

    ### Export and return
    rec_shed1.to_file(catch_output)
    return(rec_shed1)










