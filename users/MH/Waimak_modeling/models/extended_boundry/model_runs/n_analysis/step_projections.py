# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/07/2018 11:53 AM
"""

from __future__ import division
from core import env
import numpy as np
from copy import deepcopy
import pandas as pd
import itertools
from percentage_reduction_maps import gen_stream_targets, gen_well_targets, get_current_pathway_n, gen_interzone_targets
from number_of_steps import lclasses, get_average_n, get_fractions


def con_at_steps(step_reductions, pa_00, num_steps=10):
    """
    calculate the number of steps to reach a selected target scheme
    :param step_reductions: a dictionary of the percentage reduction to apply at each step e.g.
                            {'dairy_lowdrain' :20,
                            'sbda_lowdrain':10,
                            'dairy_highdrain':20,
                            'sbda_highdrain':5,
                            'lifestyle':0}
    :param pa_00: boolean if True, set teh pa rules to 0,0
    :param num_steps: number of steps to project
    :return:
    """
    assert isinstance(step_reductions, dict)
    assert set(step_reductions.keys()) == lclasses
    use_step_reductions = {}

    for k, v in step_reductions.items():
        use_step_reductions[k] = 1 - v / 100

    average_n = get_average_n()
    mode = {k: '50%' for k in average_n.keys()}
    current_paths = get_current_pathway_n(mode=mode, conservative_zones='use_mix', from_mt3d_runs=True,
                                          mt3d_add_pa=not pa_00, inc_interzone=True)

    fractions = get_fractions()
    outdata = pd.DataFrame(index=pd.MultiIndex.from_product((current_paths.keys(), range(1, num_steps + 1)),
                                                            names=['site', 'step']),
                           columns=[e + '_reduction' for e in lclasses] + ['n_con']
                           )
    for site in current_paths.keys():
        # calculate the average n concentration for the missing components (check this)
        n_proj = deepcopy(current_paths[site])
        n_proj_temp = deepcopy(current_paths[site])
        for lclass in lclasses:
            n_proj_temp += - average_n[site][lclass] * fractions[site][lclass]

        other_n_con = n_proj_temp / fractions[site]['other']

        # apply reduction, calculate new concentration
        for step in range(1, num_steps + 1):
            temp_use_reductions = {}
            for lclass in lclasses:  # assign the reduction to use and save it
                temp = step_reductions[lclass] * step/100
                if temp >= 1:
                    temp = 1
                temp_use_reductions[lclass] = temp
                outdata.loc[(site, step), '{}_reduction'.format(lclass)] = temp
            new_n = other_n_con * fractions[site]['other']
            for lclass in lclasses:
                new_n += average_n[site][lclass] * (1-temp_use_reductions[lclass]) * fractions[site][lclass]
            outdata.loc[(site, step), 'n_con'] = new_n

    return outdata

def run_scenarios(num_steps=10):
# todo make a wrapper that adds the outdata to a bigger dataframe by scen name

    pas = [True, False]

    # make reduction options
    dairy_options = [10, 25]
    sheep_options = [0, 5]
    lifestyle_options = [0]
    include_ld = [True, False]

    reduction_options = {}
    for d, s, ly, ld in itertools.product(dairy_options, sheep_options, lifestyle_options, include_ld):
        ld_name = 'exc'
        if ld:
            ld_name = 'inc'

        name = 'd{}_s{}_L{}_{}ld'.format(d, s, ly, ld_name)
        temp = {}
        temp['dairy_highdrain'] = d
        temp['sbda_highdrain'] = s
        temp['lifestyle'] = ly

        if ld:
            temp['dairy_lowdrain'] = d
            temp['sbda_lowdrain'] = s
        else:
            temp['dairy_lowdrain'] = 0
            temp['sbda_lowdrain'] = 0

        reduction_options[name] = temp
    scens = {}
    for red, pc5pa00 in itertools.product(reduction_options.keys(), pas):
        pa_name = 'without'
        if pc5pa00:
            pa_name = 'with'
        outname = '{}_{}_pc5pa00'.format(red, pa_name)
        scens[outname] = {'step_reductions': red,
                          'pa_00': pc5pa00}


    average_n = get_average_n()
    outdata = pd.DataFrame(index=pd.MultiIndex.from_product((average_n.keys(), range(1, num_steps + 1)),
                                                            names=['site', 'step']),
                           columns=pd.MultiIndex.from_product((scens.keys(),
                                                               [e + '_reduction' for e in lclasses] + ['n_con']),
                                                              names=['scenario', 'data'])
                           )
    for name, scen in scens.items():
        temp = con_at_steps(num_steps=num_steps, **scen)
        outdata[name] = temp #todo check

def plot_scenarios(outdir, data):
    # todo make a plotting function with two sub plots 1 of reduction and one of concentration by scenario
    raise NotImplementedError



if __name__ == '__main__':
    test = con_at_steps({'dairy_lowdrain': 20,
                         'sbda_lowdrain': 10,
                         'dairy_highdrain': 20,
                         'sbda_highdrain': 5,
                         'lifestyle': 0},
                        True)
    print('done')
