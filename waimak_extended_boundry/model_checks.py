# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/07/2017 1:28 PM
"""

from __future__ import division
from waimak_extended_boundry.m_packages.sfr2_packages import _get_reach_data, _get_segment_data
from waimak_extended_boundry.m_packages.wel_packages import get_wel_spd
from waimak_extended_boundry.m_packages.drn_packages import _get_drn_spd
from waimak_extended_boundry.m_packages.rch_packages import _get_rch
from waimak_extended_boundry import smt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
import os



#todo a lot of this could go into model tools (maybe)

def check_no_overlapping_features():
    no_flow = smt.get_no_flow()
    no_flow[no_flow<0]=0
    no_flow = pd.DataFrame(smt.model_where(~no_flow.astype(bool)), columns=['k','i','j'])
    no_flow['bc_type'] = 'no_flow'

    sfr_data = pd.DataFrame(_get_reach_data(smt.reach_v))
    sfr_data['bc_type'] = 'sfr'

    well_data = get_wel_spd(3,True).rename(columns={'row':'i','col':'j', 'layer':'k'})
    well_data['bc_type'] = 'well'

    drn_data = _get_drn_spd(smt.reach_v,smt.wel_version)
    drn_data['bc_type'] = 'drn'
    drn_data = drn_data.loc[~drn_data.group.str.contains('carpet')] # don't worry about carpet drains

    all_data = pd.concat((no_flow,sfr_data,well_data,drn_data))

    temp = all_data.loc[all_data.duplicated(['i','j','k'],keep=False)]
    types = ['well','sfr','drn','no_flow']
    colors = ['red','green','yellow','pink']
    sizes = [80,60,40,20]
    for layer in range(smt.layers):
        fig,ax = plt.subplots(1,1)
        ax.set_aspect('equal')
        for t,c,s in zip(types,colors,sizes):
            print(t)
            temp2 = temp.loc[(temp.bc_type==t) & (temp.k==layer)]
            ax.scatter(temp2.j,temp2.i*-1,c=c,s=s,label=t)
            ax.set_title('layer {}'.format(layer))
        plt.legend()
        plt.show(fig)

    if any(all_data.duplicated(['i','j','k'],keep=False)): #I did not fix the well drain overlap s of waimak because we don't care
        #there are some boundry flux drain overlaps, but is shoudn't affect stream depletion assessments (everything is relavitve)
        raise ValueError ('There are duplicate boundry conditions')

def check_layer_overlap():
    for i in range(smt.layers):
        fig, ax = smt.plt_matrix(smt.check_layer_overlap(use_elv_db=True,layer=i,required_overlap=0.50),title='layer {}'.format(i),vmax=2,vmin=-2,no_flow_layer=i)
        plt.show(fig)

def check_elv_db():
    no_flow = smt.get_no_flow(0)
    elv = smt.calc_elv_db()
    elv[~(np.repeat(no_flow[np.newaxis,:,:].astype(bool),smt.layers,axis=0))] = np.nan
    tops = elv[0:-1]
    bots = elv[1:]
    if any((bots>tops).flatten()):
        raise ValueError('some bottoms are higher than tops')

def check_noflow_overlap():
    no_flow = np.abs(smt.get_no_flow())
    nabove = no_flow[0:-1]
    nbelow = no_flow[1:]
    if any((nbelow>nabove).flatten()):
        raise ValueError('there exist pockets of flow under no_flow')

def check_null_spd():
    sfr_data = pd.DataFrame(_get_reach_data(smt.reach_v))
    if any(np.array(sfr_data.isnull()).flatten()):
        raise ValueError('null data in SFR spd')

    well_data = get_wel_spd(smt.wel_version)
    if any(np.array(well_data.loc[:,[u'col', u'flux', u'layer', u'row']].isnull()).flatten()):
        raise ValueError('null data in well spd')

    drn_data = _get_drn_spd(smt.reach_v, smt.wel_version)
    if any(np.array(drn_data.loc[:,['k','i','j','elev','cond']].isnull()).flatten()):
        raise ValueError('null data in drain spd')

    rch = _get_rch()
    if any(np.isnan(rch).flatten()):
        raise ValueError('null data in recharge spd')


def create_digital_appendix(root_dir,dpi):
    base_dir = '{}/individual_pdfs'.format(root_dir)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    paths = []
    paths.extend(check_well_loc_discharge(base_dir,dpi))
    paths.extend(check_drain_locations(base_dir,dpi))
    paths.extend(check_sfr_locations(base_dir,dpi))
    paths.extend(check_rch_spatially(base_dir,dpi))
    paths.extend(check_elevations_spatially(base_dir,dpi))



# #locations of wells
def check_well_loc_discharge(base_dir,dpi):
    paths = []
    org_well = get_wel_spd(smt.wel_version)
    wells = org_well.loc[org_well.type=='well']
    races = org_well.loc[org_well.type=='race']
    rivers = org_well.loc[org_well.type=='river']
    boundry_flux = org_well.loc[org_well.type.str.contains('boundry_flux')]

    # waimak = 4, chch_wm = 7, selwyn=8
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    # 1 = ashley, 2 = eyre, 3 = cust
    sub_zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/Sub_cwms_zones.shp".format(smt.sdp),'id_num')
    flux_stats = pd.DataFrame(index=['Ashley sub-zone', 'Cust sub-zone', 'Eyre sub-zone','Waimakariri CWMS',
                                     'CHCH-WM CWMS','Selwyn CWMS'],
                              columns=['pumping wells', 'WILS race losses', 'Losing Reaches', 'Boundary Flux'])
    for flux, name in zip([wells, races, rivers, boundry_flux], ['pumping wells', 'WILS race losses', 'Losing Reaches', 'Boundary Flux']):
        flux_array = smt.df_to_array(flux,'flux')
        for zone_id, zone_name in zip([1, 3, 2,
                                       4, 7, 8],
                                      ['Ashley sub-zone', 'Cust sub-zone', 'Eyre sub-zone',
                                       'Waimakariri CWMS', 'CHCH-WM CWMS','Selwyn CWMS']):
            if zone_id > 3:
                id_matrix = zones
            else:
                id_matrix = sub_zones

            flux_stats.loc[zone_name,name] = np.nansum(flux_array[np.isclose(id_matrix,zone_id)])
    flux_stats = flux_stats/86400
    fig,ax = plt.subplots(1,1,figsize=(14, 4.5))
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=flux_stats.values,
              rowLabels=flux_stats.index,
              colLabels=flux_stats.columns,
              cellLoc='center', rowLoc='center',
              loc='center')
    ax.set_title('fluxes by zone negative numbers equate to flux out of the model')
    path = '{}/well_table.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)


    # plots
    cmap, norm = from_levels_and_colors([-1, 0, 1,2], ['blue', 'black', 'white'])
    x,y =np.meshgrid(np.arange(0,smt.cols),np.arange(0,-smt.rows,-1))
    for layer in range(smt.layers):
        temp = wells.loc[wells.layer==layer]
        no_flow = smt.get_no_flow(layer)
        fig,ax = plt.subplots(1,1,figsize=(18.5, 9.5))
        ax.set_aspect('equal')
        ax.pcolormesh(x,y,no_flow,cmap=cmap,norm=norm)
        ax.scatter(temp.col+0.5,-1*temp.row-0.5,c='r')
        ax.set_title('pumping wells in layer {}'.format(layer+1))
        path = '{}/pumping_wells{}.png'.format(base_dir,layer)
        fig.savefig(path,dpi=dpi)
        paths.append(path)
        plt.close(fig)

    fig,ax = plt.subplots(1,1,figsize=(18.5, 9.5))
    no_flow = smt.get_no_flow(0)
    ax.set_aspect('equal')
    ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm)
    ax.scatter(races.col+0.5,-1*races.row-0.5, c='y',label='WILLS Races')
    ax.scatter(boundry_flux.col+0.5,-1*boundry_flux.row-0.5, c='g', label='N boundary flux')
    ax.scatter(rivers.col+0.5,-1*rivers.row-0.5, c='purple',label='losing rivers')
    ax.set_title('other well boundary conditions')
    ax.legend()
    path = '{}/well_types.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)


    return paths

def check_drain_locations(base_dir,dpi):
    paths = []
    drn_data = _get_drn_spd(smt.reach_v,smt.wel_version)
    cmap, norm = from_levels_and_colors([-1, 0, 1,2], ['blue', 'black', 'white'])
    cmap2, norm2 = from_levels_and_colors([1, 2], ['gold'])
    x, y = smt.get_model_x_y()
    for group in set(drn_data.group):
        temp = smt.df_to_array(drn_data.loc[drn_data.group == group],'k')
        temp[np.isfinite(temp)] = 1
        fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
        no_flow = smt.get_no_flow(0)
        ax.set_aspect('equal')
        ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm)
        ax.pcolormesh(x, y, np.ma.masked_invalid(temp), cmap=cmap2,norm=norm2)
        ax.set_title('{} drains'.format(group))
        path = '{}/{}drains.png'.format(base_dir, group)
        fig.savefig(path,dpi=dpi)
        paths.append(path)
        plt.close(fig)

    return paths

def check_sfr_locations(base_dir,dpi):
    paths = []
    str_data = pd.DataFrame(_get_reach_data(smt.reach_v))
    seg_data = pd.DataFrame(_get_segment_data(smt.seg_v)).set_index('nseg')
    str_data.loc[:, 'flow'] = 0
    for i in str_data.index:
        if str_data.loc[i,'ireach'] == 1:
            str_data.loc[i,'flow'] = seg_data.loc[str_data.loc[i, 'iseg'],'flow']

    str_temp=smt.df_to_array(str_data,'flow')


    cmap, norm = from_levels_and_colors([-1, 0, 1,2], ['blue', 'black', 'white'])
    cmap2, norm2 = from_levels_and_colors([0, 1,99999999], ['green', 'fuchsia'])
    x,y =np.meshgrid(np.arange(0,smt.cols),np.arange(0,-smt.rows,-1))
    fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
    no_flow = smt.get_no_flow(0)
    ax.set_aspect('equal')
    ax.pcolormesh(x[:190,:], y[:190,:], no_flow[:190,:], cmap=cmap, norm=norm)
    ax.pcolormesh(x[:190,:],y[:190,:], np.ma.masked_invalid(str_temp)[:190,:],cmap=cmap2, norm=norm2)
    ax.set_title('stream reaches, pink represents surface water influx')
    path = '{}/streams.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)


    # belop has not been saved
    # plot up/downstream segments
    str_segs = smt.df_to_array(str_data,'iseg')
    outseg_data = seg_data.reset_index().set_index('outseg')


    cmap3, norm3 = from_levels_and_colors([-1, 2], ['green'])
    cmap4, norm4 = from_levels_and_colors([-1, 2], ['fuchsia'])

    for seg in set(str_data.iseg):
        try:
            outsegs = np.atleast_1d(outseg_data.loc[seg,'nseg'])
        except KeyError:
            outsegs = None

        fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
        ax.set_aspect('equal')
        ax.pcolormesh(x[:190, :], y[:190, :], no_flow[:190, :], cmap=cmap, norm=norm)
        current_stream_idx = np.where(str_segs==seg)
        current_stream = np.zeros((smt.rows,smt.cols))*np.nan
        current_stream[current_stream_idx] = 1

        ax.pcolormesh(x[:190, :], y[:190, :], np.ma.masked_invalid(current_stream)[:190, :], cmap=cmap3, norm=norm3)
        if outsegs is not None:
            tribs = np.zeros((smt.rows,smt.cols))*np.nan
            for outseg in outsegs:
                tribs_idx = np.where(np.isclose(str_segs,outseg))
                tribs[tribs_idx] = 1
            ax.pcolormesh(x[:190, :], y[:190, :], np.ma.masked_invalid(tribs)[:190, :], cmap=cmap4, norm=norm4)

        ax.set_title('stream reaches, outflow stream segement is green tributaries are pink')

    return paths

#rch spatially
def check_rch_spatially(base_dir,dpi):
    paths = []
    # plot rch spatially
    rch = _get_rch()
    cmap, norm = from_levels_and_colors([-1, 0, 1,2], ['blue', 'black', 'white'])
    x,y =np.meshgrid(np.arange(0,smt.cols),np.arange(0,-smt.rows,-1))
    fig,ax = plt.subplots(1,1,figsize=(18.5, 9.5))
    no_flow = smt.get_no_flow(0)
    ax.set_aspect('equal')
    ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm)
    rch[no_flow<=0]=np.nan
    pcm = ax.pcolormesh(x,y,np.ma.masked_invalid(rch)*1000,cmap='plasma')
    ax.set_title('recharge, mm/day')
    fig.colorbar(pcm, ax=ax, extend='max')
    path = '{}/rch.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)

    return paths

#check elevations
def check_elevations_spatially(base_dir,dpi):
    paths = []
    # plot 5m contours of base of each layer in the database
    elv = smt.calc_elv_db()
    cmap, norm = from_levels_and_colors([-1, 0, 1,2], ['blue', 'black', 'white'])
    x,y =np.meshgrid(np.arange(0,smt.cols),np.arange(0,-smt.rows,-1))
    for i in range(smt.layers+1):
        if i ==0:
            layer = 0
            title = 'top'
            step = 20
        else:
            layer = i-1
            title = 'bottom of layer {}'.format(i)
            step = 10
        if i==11:
            step=20
        fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
        no_flow = smt.get_no_flow(layer)
        ax.set_aspect('equal')
        levels = np.arange(np.nanmin(elv[i]), np.nanmax(elv[i]), step)
        levels = [int(np.round(e)) for e in levels]
        ct = ax.contour(x,y,elv[i], colors='k', levels=levels)
        ax.set_title(title)
        ax.clabel(ct, fontsize=9, inline=1, fmt='%3.0f')
        ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm, alpha=0.7)
        path = '{}/elv{}.png'.format(base_dir,i)
        fig.savefig(path,dpi=dpi)
        paths.append(path)
        plt.close(fig)


    # plot 5m contours of thickness for layer 0 and 10
    fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
    no_flow = smt.get_no_flow(0)
    ax.set_aspect('equal')
    thick = elv[0]-elv[1]
    levels = np.arange(np.nanmin(thick), np.nanmax(thick), 10)
    levels = [int(np.round(e)) for e in levels]
    ct = ax.contour(x, y, thick, colors='k', levels=levels)
    ax.set_title('layer 1 thickness')
    ax.clabel(ct, fontsize=9, inline=1, fmt='%3.0f')
    ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm, alpha=0.7)
    path = '{}/layer1 thickness.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)


    fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
    no_flow = smt.get_no_flow(10)
    ax.set_aspect('equal')
    thick = elv[10]-elv[11]
    levels = np.arange(np.nanmin(thick), np.nanmax(thick), 20)
    levels = [int(np.round(e)) for e in levels]
    ct = ax.contour(x, y, thick, colors='k', levels=levels)
    ax.set_title('layer 11 thickness')
    ax.clabel(ct, fontsize=9, inline=1, fmt='%3.0f')
    ax.pcolormesh(x, y, no_flow, cmap=cmap, norm=norm, alpha=0.7)
    path = '{}/layer11thick.png'.format(base_dir)
    fig.savefig(path,dpi=dpi)
    paths.append(path)
    plt.close(fig)

    return paths

if __name__ == '__main__':
    check_no_overlapping_features() #passed # this has failed for now we are just ignoring it due to cost benifit
    check_layer_overlap() #passed
    check_elv_db() #passed
    check_noflow_overlap() #passed
    check_null_spd() #passed
    create_digital_appendix(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\digital_appendix_06-09-2017",None)
    print('all passed')
