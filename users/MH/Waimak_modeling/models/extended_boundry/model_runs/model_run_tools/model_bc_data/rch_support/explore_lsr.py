# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/09/2017 2:17 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.rch_packages import _get_rch_comparison
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.LSR_arrays import \
    get_lsrm_base_array, get_lsr_base_period_inputs, get_ird_base_array
import itertools
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
import matplotlib.patches as mpatches

def plt_cc_rch(naturalised, pc5, title):
    """
    plot up the recharge as lines for the chch data
    :param naturalised: boolean
    :param pc5: Boolean
    :param title: title to pass to the plot
    :return: fig, ax
    """
    periods = range(2010, 2100, 20)
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    rcps_lines = {'RCP4.5': '-', 'RCP8.5': '--'}
    rcm_colors = {'BCC-CSM1.1': 'k', 'CESM1-CAM5': 'r', 'GFDL-CM3': 'g', 'GISS-EL-R': 'b', 'HadGEM2-ES': 'orange',
                  'NorESM1-M': 'purple'}
    sen = 'current'
    if naturalised:
        sen = 'nat'
    if pc5:
        sen = 'pc5'

    fig, ax = plt.subplots(figsize=(18.5, 9.5))
    ax.set_title(title)
    ax.set_xlabel('year')
    ax.set_ylabel('mean LSR / mean opt LSR')
    at = 'period_mean'
    # get data
    for rcp in rcps:
        linestl = rcps_lines[rcp]
        for rcm in rcms:
            color = rcm_colors[rcm]
            data = []
            for period in periods:
                temp = get_lsrm_base_array(sen, rcp, rcm, period, at)
                base_rch = get_lsrm_base_array(*get_lsr_base_period_inputs(sen, rcp, rcm, period, at))
                data.append(np.nanmean(temp / base_rch))
            data = np.array(data)
            ax.plot(periods, data, linestyle=linestl, color=color, label=rcp)
    ax.legend()
    return fig, ax


def plt_vcsn():
    raise NotImplementedError
    # check the vcsn data but have already done a bit of this


def comp_amalg_types(rcm, rcp, naturalised=False, pc5=False):
    """
    plot up and check the differnt amalgimation types
    :param rcm: the rcm
    :param rcp: the rcp
    :param naturalised: boolean
    :param pc5: boolean
    :return: fig,ax
    """
    periods = range(2010, 2100, 20)
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    des = 'current_practice'
    sen = 'current'
    if naturalised:
        des = 'naturalised'
        sen = 'nat'
    if pc5:
        des = 'pc5'
        sen = 'pc5'

    fig, ax = plt.subplots(figsize=(18.5, 9.5))
    ax.set_title('{}_{}_{}'.format(rcm, rcp, des))
    ax.set_xlabel('year')
    ax.set_ylabel('mean LSR / mean opt LSR')

    for at, color in zip(amalg_types, ['k', 'r', 'b']):
        data = []
        for period in periods:
            temp = get_lsrm_base_array(sen, rcp, rcm, period, at)
            base_rch = get_lsrm_base_array(*get_lsr_base_period_inputs(sen, rcp, rcm, period, at))
            data.append(np.nanmean(temp / base_rch))
        data = np.array(data)
        ax.plot(periods, data, color=color, label=at)
    ax.legend()
    return fig, ax


def plt_change_maps(base_dir):
    """
    plot maps of the changes in LSR to a base period
    :param base_dir: the directory to save the plots in
    :return:
    """
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    sens = ['current', 'pc5', 'nat']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    per = 2090
    no_flow = smt.get_no_flow(1)
    no_flow[no_flow < 0] = 0
    no_flow = no_flow.astype(bool)
    for rcp, rcm, sen, at in itertools.product(rcps, rcms, sens, amalg_types):
        mod_rch = get_lsrm_base_array(sen, rcp, rcm, per, at)
        base_rch = get_lsrm_base_array(*get_lsr_base_period_inputs(sen, rcp, rcm, per, at))
        plt_data = mod_rch / base_rch
        plt_data[~no_flow] = np.nan
        ttl = 'cc_vs_rcppast_{}-{}-{}-{}-{}'.format(sen, rcp, rcm, per, at)
        mapfig, mapax = smt.plt_matrix(plt_data, title=ttl, cmap='RdBu',vmin=0.2,vmax=1.8)
        outdir = os.path.join(base_dir, at)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        mapfig.savefig(os.path.join(outdir, '{}.png'.format(ttl)))
        plt.close(mapfig)


def plt_ts_at_comp(base_dir):
    """
    plot all three of the time series
    :param base_dir: the dir to save the plots in
    :return:
    """
    # 100%
    fig, ax = plt_cc_rch(False, True, 'pc5')
    fig.savefig(os.path.join(base_dir, ax.title._text + '.png'))
    # 80%
    fig2, ax2 = plt_cc_rch(False, False, 'no_change')
    fig2.savefig(os.path.join(base_dir, ax2.title._text + '.png'))

    # Nat
    fig3, ax3 = plt_cc_rch(True, False, 'naturalised')
    fig3.savefig(os.path.join(base_dir, ax3.title._text + '.png'))

    # amalg types
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    for rcm, rcp in itertools.product(rcms, rcps):
        fig, ax = comp_amalg_types(rcm, rcp)
        fig.savefig(os.path.join(base_dir, ax.title._text + '.png'))
    for rcm, rcp in itertools.product(rcms, rcps):
        fig, ax = comp_amalg_types(rcm, rcp, True)
        fig.savefig(os.path.join(base_dir, ax.title._text + '.png'))


def plt_rch_budget(base_dir):
    """
    plot budgets of recharge for different sensarios
    :param base_dir: dir to save the plots in
    :return:
    """
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    sens = ['current', 'pc5', 'nat']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    rcps_lines = {'RCP4.5': '-', 'RCP8.5': '--'}
    rcm_colors = {'BCC-CSM1.1': 'k', 'CESM1-CAM5': 'r', 'GFDL-CM3': 'g', 'GISS-EL-R': 'b', 'HadGEM2-ES': 'orange',
                  'NorESM1-M': 'purple'}
    periods = range(2010, 2100, 20)

    for at, (zone, name) in itertools.product(amalg_types,
                                              zip(['chch', 'waimak', 'selwyn', ['chch', 'waimak', 'selwyn']],
                                                  ['chch', 'waimak', 'selwyn', 'total'])):
        zidx = get_zone_array_index(zone)
        fig, axes = plt.subplots(nrows=3, figsize=(18.5, 9.5))
        ttl = 'forward_{}-{}'.format(at, name)
        fig.suptitle(ttl)
        mins = []
        maxs = []
        for ax, sen in zip(axes.flatten(), sens):
            for rcp in rcps:
                for rcm in rcms:
                    temp_data = []
                    for per in periods:
                        temp_rch = get_lsrm_base_array(sen, rcp, rcm, per, at) * 200 * 200
                        temp_data.append(temp_rch[zidx].sum() / 86400)
                    ax.plot(periods, temp_data, linestyle=rcps_lines[rcp], color=rcm_colors[rcm], label=rcp)
            ax.set_ylabel('zonal rch m3/s')
            ax.set_title(sen)
            mins.append(ax.get_ylim()[0])
            maxs.append(ax.get_ylim()[1])
        for ax in axes.flatten():
            ax.set_ylim([np.min(mins), np.max(maxs)])
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles,labels)
        fig.tight_layout()
        fig.savefig(os.path.join(base_dir, '{}.png'.format(ttl)))
        plt.close(fig)

    # rcppast
    fig,axs = plt.subplots(nrows=3,ncols=4,figsize=(18.5, 9.5))
    fig.suptitle('rcp_past')
    mins = []
    maxs = []
    for i,sen in enumerate(sens):
        for j, (zone,name) in enumerate(zip(['chch', 'waimak', 'selwyn', ['chch', 'waimak', 'selwyn']],
                                              ['chch', 'waimak', 'selwyn', 'total'])):
            ax = axs[i,j]
            ax.set_title('{}|{}'.format(sen,name))
            names = []
            data = []
            zidx = get_zone_array_index(zone)
            for k, at in enumerate(amalg_types):
                names.append(name)
                temp = []
                for rcm in rcms:
                    rch = get_lsrm_base_array(sen,'RCPpast',rcm,1980,at)*200*200
                    temp.append(rch[zidx].sum()/86400)
                data.append(temp)
            ax.boxplot(data)
            ax.set_xticklabels(amalg_types)
            mins.append(ax.get_ylim()[0])
            maxs.append(ax.get_ylim()[1])
    for ax in axs.flatten():
        ax.set_ylim([np.min(mins), np.max(maxs)])
    fig.tight_layout()
    fig.savefig(os.path.join(base_dir,'rcppast.png'))
    plt.close(fig)


    # model_per and vcsn
    model_per = _get_rch_comparison()/86400
    model_per.index = ['mod_per_{}'.format(e) for e in model_per.index]

    # vcsn
    temp = ['vcsn_{}'.format(e) for e in sens]
    vcsn = pd.DataFrame(index=temp, columns=['chch', 'waimak', 'selwyn', 'total'])
    for sen, (zone, name) in itertools.product(sens,zip(['chch', 'waimak', 'selwyn', ['chch', 'waimak', 'selwyn']],
                                              ['chch', 'waimak', 'selwyn', 'total'])):
        zidx = get_zone_array_index(zone)
        rch = get_lsrm_base_array(sen,None,None,None,'mean')*200*200
        vcsn.loc['vcsn_{}'.format(sen),name] = rch[zidx].sum()/86400
    outdata = pd.concat((model_per,vcsn))
    outdata.to_csv(os.path.join(base_dir,'mod_per_vcsn_budgets.csv'))

def boxplot_cc_projections(base_dir):
    """
    plot up boxplots of the senarios of recharge budgets
    :param base_dir: dir to save the plots
    :return:
    """
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    sens = ['current', 'pc5', 'nat']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    rcp_colors = {'RCP4.5': 'b', 'RCP8.5': 'r'}
    periods = range(2010, 2100, 20)


    for at, (zone, name) in itertools.product(amalg_types,
                                              zip(['chch', 'waimak', 'selwyn',
                                                   ['chch', 'waimak', 'selwyn'], 'inland_waimak', 'coastal_waimak'],
                                                  ['chch', 'waimak', 'selwyn',
                                                   'total', 'inland_waimak', 'coastal_waimak'])):
        zidx = get_zone_array_index(zone)
        fig, axes = plt.subplots(nrows=3, figsize=(18.5, 9.5))
        ttl = 'boxplot_forward_{}-{}'.format(at, name)
        fig.suptitle(ttl)
        mins = []
        maxs = []
        for ax, sen in zip(axes.flatten(), sens):
            for rcp in rcps:
                plt_data = []
                temp_data = []
                for rcm in rcms:
                    temp_rch = get_lsrm_base_array(sen,'RCPpast',rcm,1980,at) *200 *200
                    temp_data.append(temp_rch[zidx].sum() / 86400)
                plt_data.append(temp_data)
                for per in periods:
                    temp_data = []
                    for rcm in rcms:
                        temp_rch = get_lsrm_base_array(sen, rcp, rcm, per, at) * 200 * 200
                        temp_data.append(temp_rch[zidx].sum() / 86400)

                    plt_data.append(temp_data)
                if rcp == 'RCP8.5':
                    add_val = 0.1
                else:
                    add_val = 0
                bp = ax.boxplot(plt_data,positions=np.array(range(len(periods)+1)) + add_val)
                for box in bp['boxes']:
                    # change outline color
                    box.set(color=rcp_colors[rcp], linewidth=2, alpha=0.5)

            ax.set_ylabel('zonal rch m3/s')
            ax.set_title(sen)
            ax.set_xticklabels([1980] + periods)
            mins.append(ax.get_ylim()[0])
            maxs.append(ax.get_ylim()[1])
        for ax in axes.flatten():
            ax.set_ylim([np.min(mins), np.max(maxs)])
        handles = [mpatches.Patch(color='r', label='RCP8.5'),mpatches.Patch(color='b', label='RCP4.5')]
        labels = ['RCP8.5', 'RCP4.5']
        fig.legend(handles,labels)
        fig.tight_layout()
        fig.savefig(os.path.join(base_dir, '{}.png'.format(ttl)))
        plt.close(fig)

def boxplot_cc_projections_ird(base_dir):
    """
    plot boxplots of the irrigation demand budgets
    :param base_dir: dir to save the plots in
    :return:
    """
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    sens = ['current', 'pc5', 'nat']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    rcp_colors = {'RCP4.5': 'b', 'RCP8.5': 'r'}
    periods = range(2010, 2100, 20)


    for at, (zone, name) in itertools.product(amalg_types,
                                              zip(['chch', 'waimak', 'selwyn',
                                                   ['chch', 'waimak', 'selwyn'], 'inland_waimak', 'coastal_waimak'],
                                                  ['chch', 'waimak', 'selwyn',
                                                   'total', 'inland_waimak', 'coastal_waimak'])):
        zidx = get_zone_array_index(zone)
        fig, axes = plt.subplots(nrows=3, figsize=(18.5, 9.5))
        ttl = 'ird_boxplot_forward_{}-{}'.format(at, name)
        fig.suptitle(ttl)
        mins = []
        maxs = []
        for ax, sen in zip(axes.flatten(), sens):
            for rcp in rcps:
                plt_data = []
                temp_data = []
                for rcm in rcms:
                    temp_rch = get_ird_base_array(sen,'RCPpast',rcm,1980,at) *200 *200
                    temp_data.append(np.nansum(temp_rch[zidx]) / 86400)
                plt_data.append(temp_data)
                for per in periods:
                    temp_data = []
                    for rcm in rcms:
                        temp_rch = get_ird_base_array(sen, rcp, rcm, per, at) * 200 * 200
                        temp_data.append(np.nansum(temp_rch[zidx]) / 86400)

                    plt_data.append(temp_data)
                if rcp == 'RCP8.5':
                    add_val = 0.1
                else:
                    add_val = 0
                bp = ax.boxplot(plt_data,positions=np.array(range(len(periods)+1)) + add_val)
                for box in bp['boxes']:
                    # change outline color
                    box.set(color=rcp_colors[rcp], linewidth=2, alpha=0.5)

            ax.set_ylabel('zonal rch m3/s')
            ax.set_title(sen)
            ax.set_xticklabels([1980] + periods)
            mins.append(ax.get_ylim()[0])
            maxs.append(ax.get_ylim()[1])
        for ax in axes.flatten():
            ax.set_ylim([np.min(mins), np.max(maxs)])
        handles = [mpatches.Patch(color='r', label='RCP8.5'),mpatches.Patch(color='b', label='RCP4.5')]
        labels = ['RCP8.5', 'RCP4.5']
        fig.legend(handles,labels)
        fig.tight_layout()
        fig.savefig(os.path.join(base_dir, '{}.png'.format(ttl)))
        plt.close(fig)



if __name__ == '__main__':
    # runs
    base_dir_all = env.sci(
        "Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\lsr_checks")
    # note this means RCP past is not tested, and vcsn only minimally tested
    base_dir_comp = os.path.join(base_dir_all, 'comparisons')
    if not os.path.exists(base_dir_comp):
        os.makedirs(base_dir_comp)
    base_dir_maps = os.path.join(base_dir_all, 'maps')
    if not os.path.exists(base_dir_maps):
        os.makedirs(base_dir_maps)
    base_dir_budget = os.path.join(base_dir_all, 'budgets')
    if not os.path.exists(base_dir_budget):
        os.makedirs(base_dir_budget)
    test_type = [4]

    if 0 in test_type:
        plt_ts_at_comp(base_dir_comp)
    if 1 in test_type:
        plt_change_maps(base_dir_maps)
    if 2 in test_type:
        plt_rch_budget(base_dir_budget)
    if 3 in test_type:
        boxplot_cc_projections(base_dir_budget)
    if 4 in test_type:
        boxplot_cc_projections_ird(base_dir_budget)
