# -*- coding: utf-8 -*-
"""
Network processing and analysis functions.
"""

############################################
#### Network processing and analysis


def nx_shp(shp_pts, shp_lines, site_col='site'):
    """
    Function to import shapefiles into nx networks. The lines become the edges and the points are used as the node (names). The points shapefile should have a site name column that is precisely at the intersection of the shp_lines.
    """
    from networkx import read_shp, relabel_nodes
    from numpy import where, in1d

    ## Read in shapefiles
    t2 = read_shp(shp_pts)
    t3 = read_shp(shp_lines)

    ## extract the site names
    sites = [i[1][site_col] for i in t2.nodes(data=True)]

    ## Set up rename dictionaries
    rename1 = {(i[0], i[1]): (round(i[0]), round(i[1])) for i in t3.nodes()}
    rename2 = {(round(i[0][0]), round(i[0][1])): i[1]['site'] for i in t2.nodes(data=True)}

    ## Rename nodes
    g2 = relabel_nodes(t3, rename1)
    g3 = relabel_nodes(g2, rename2)

    ## Remove unnecessary nodes
    remove1 = [g3.nodes()[i] for i in where(~in1d(g3.nodes(), sites))[0]]
    g3.remove_nodes_from(remove1)

    return(g3)


def str_paths(nx1):
    """
    Function to process the stream network to determine the nodes and edges to downstream gauging sites. The input is from the output of nx_shp. The output is a two dictorionary object list of site nodes and site edges.
    """
    from networkx import all_pairs_shortest_path, all_pairs_dijkstra_path_length

    def iter1(g1, d2, site):
        from numpy import argmin

        keys1 = g1.keys()
        sites2 = [i for i in keys1 if ((i != site) & (i < 10000000))]
        if not sites2:
            output = [site]
        else:
            len1 = [d2[site][i] for i in sites2]
            down_site = sites2[argmin(len1)]
            output = g1[down_site]
        return(output)

    ## Determine all paths
    p1 = all_pairs_shortest_path(nx1)
    d1 = all_pairs_dijkstra_path_length(nx1, None, 'len')

    ## Make list of all sites
    sites = [i for i in nx1.nodes() if (i < 10000000)]

    ## Extract the paths for all sites (into a dict)
    p2 = {i: p1[i] for i in sites}
    d2 = {i: d1[i] for i in sites}

    site_nodes = {i: iter1(p2[i], d2, i) for i in p2}
    site_paths = {i: [j[2] for j in nx1.out_edges(site_nodes[i], data='num')][0:-1] for i in site_nodes}
    return([site_nodes, site_paths])


#####################################################
#### MFE REC streams network and catchments


def find_upstream_rec(nzreach, rec_shp=r'\\fs02\ManagedShares2\Data\Surface Water\shared\GIS_base\vector\streams\rec-canterbury-2010.shp'):
    """
    Function to estimate all of the reaches (and nodes) upstream of specific reaches. Input is a list/array/Series of NZREACH IDs.
    """
    from core.ecan_io import rd_sql
    from pandas import concat
    from geopandas import read_file, GeoDataFrame

    ### Parameters
#    server = 'SQL2012PROD05'
#    db = 'GIS'
#    table = 'MFE_NZTM_REC'
#    cols = ['NZREACH', 'NZFNODE', 'NZTNODE']
#
#    ### Load data
#    rec = rd_sql(server, db, table, cols)
    if isinstance(rec_shp, GeoDataFrame):
        rec = rec_shp.copy().drop('geometry', axis=1)
    elif isinstance(rec_shp, str):
        if rec_shp.endswith('shp'):
            rec = read_file(rec_shp).drop('geometry', axis=1)
    else:
        raise TypeError('rec_shp must be either a GeoDataFrame or a shapefile')

    ### Run through all nzreaches
    reaches_lst = []
    for i in nzreach:
        reach1 = rec[rec.NZREACH == i]
        up1 = rec[rec.NZTNODE.isin(reach1.NZFNODE)]
        while not up1.empty:
            reach1 = concat([reach1, up1])
            up1 = rec[rec.NZTNODE.isin(up1.NZFNODE)]
        reach1.loc[:, 'start'] = i
        reaches_lst.append(reach1)

    reaches = concat(reaches_lst)
    reaches.set_index('start', inplace=True)
    return(reaches)


def extract_rec_catch(reaches, rec_catch_shp=r'\\fs02\ManagedShares2\Data\Surface Water\shared\GIS_base\vector\catchments\river-environment-classification-watershed-canterbury-2010.shp'):
    """
    Function to extract the catchment polygons from the rec catchments layer. Appends to reaches layer.
    """
    from core.ecan_io import rd_sql
    from pandas import concat
    from geopandas import read_file, GeoDataFrame

    ### Parameters
#    server = 'SQL2012PROD05'
#    db = 'GIS'
#    table = 'MFE_NZTM_RECWATERSHEDCANTERBURY'
#    cols = ['NZREACH']
#
    sites = reaches.NZREACH.unique().astype('int32').tolist()
#
#    ### Extract reaches from SQL
#    catch1 = rd_sql(server, db, table, cols, where_col='NZREACH', where_val=sites, geo_col=True)
#    catch2 = catch1.dissolve('NZREACH')
    if isinstance(rec_catch_shp, GeoDataFrame):
        catch0 = rec_catch_shp.copy()
    elif isinstance(rec_catch_shp, str):
        if rec_catch_shp.endswith('shp'):
            catch0 = read_file(rec_catch_shp)
    else:
        raise TypeError('rec_shp must be either a GeoDataFrame or a shapefile')

    catch1 = catch0[catch0.NZREACH.isin(sites)]
    catch2 = catch1.dissolve('NZREACH').reset_index()[['NZREACH', 'geometry']]

    ### Combine with original sites
    catch3 = catch2.merge(reaches.reset_index(), on='NZREACH')
    catch3.crs = catch0.crs

    return(catch3)


def agg_rec_catch(rec_catch):
    """
    Simple function to aggregate rec catchments.
    """
    rec_shed = rec_catch[['start', 'geometry']].dissolve('start')
    rec_shed.index = rec_shed.index.astype('int32')
    rec_shed['area'] = rec_shed.area
    rec_shed.crs = rec_catch.crs
    return(rec_shed.reset_index())




















