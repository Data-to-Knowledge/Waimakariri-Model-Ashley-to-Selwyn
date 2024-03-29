# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/04/2018 10:08 AM

Non mt3d predictions -- superceeded
"""

from __future__ import division
from run_source_uncertainty import calc_all_ns, _np_describe
import pandas as pd
import env
import numpy as np
import os
import itertools
from glob import glob
from waimak_extended_boundry.model_run_tools.n_analysis_support.nitrate_at_key_receptors import \
    get_well_ids, get_str_ids
from warnings import warn

# break up the private and WDC wells into different folders for this process to avoid name overlaps
well_zones_to_modify = {
    # site name: modified n value

    # WDC wells
    # ## leave as is:
    # ###### fernside
    # ###### manderville
    # ###### oxford urban
    # ###### pegasus
    # ###### poyntzs rd
    # ###### waikuku

    # cust emmma calculation s median of 18% alpine
    'wdc_Cust': 6.5,

    # kaiapoi # from emma medain of 32% alpine
    'wdc_Kaiapoi': 6.03,

    # kairaki  #from emma medain from kaiakpoi of 32% alpine
    'wdc_Kairaki': 4.9,

    # ohoka # from EMMA median of 20%
    'wdc_Ohoka': 7.2,

    # rangiora #from emma medain from kaiakpoi of 32% alpine
    'wdc_Rangiora': 6.72,

    # west eyreton as west eyreton wells  # from emma 19% alpine
    'wdc_West Eyreton': 6.0,

    # private wells
    # ## leave as is:
    # 'Fernside',
    # 'Flaxton',
    # 'Horellville',
    # 'Mandeville',
    # 'North East Eyrewell_shallow',
    # 'North West Eyrewell_shallow',
    # 'Rangiora',
    # 'Swannanoa_shallow',
    # 'Waikuku',
    # 'Woodend - Tuahiwi',
    # 'West Eyreton_shallow',


    # Clarkville , scale similar to kaiapoi at island road waimakariri component 20% alpine river water
    'Clarkville': 7.6,

    # Cust, similar to cust WDC well emmma calculation s median of 18% alpine
    'Cust': 6.68,

    # 'Eyreton_deep', use our EMMA data 32% alpine river water
    'Eyreton_deep': 13.3,

    # 'Eyreton_shallow', as kaiapoi at island road component 20% alpine river water
    'Eyreton_shallow': 11,

    # 'North East Eyrewell_deep', # look for EMMA lump with NW eyrewell 28% alpine river water from emma
    'North East Eyrewell_deep': 6.9,

    # 'North West Eyrewell_deep', # look for EMMA lump with NE eyrewell 28% alpine river water from emma
    'North West Eyrewell_deep': 8.9,

    # 'Ohoka_deep', # as ohoka WDC well  from EMMA median of 20%
    'Ohoka_deep': 7.3,

    # 'Ohoka_shallow', # as ohoka stream # from emma we can expect a median of 12% waimak water
    'Ohoka_shallow': 6,

    # 'Springbank', #as cust (and cust wdc) emmma calculation s median of 18% alpine
    'Springbank': 6.5,

    # 'Summerhill', # look for emma 26% alpine river water
    'Summerhill': 9.5,

    # 'Swannanoa_deep', as WDC west eyreton # from emma 19% alpine
    'Swannanoa_deep': 8.9,

    # 'West Eyreton_deep', as swannanowa/WDC west eyreton # from emma 19% alpine
    'West Eyreton_deep': 6.55
}

stream_zones_to_modify = {
    # site name: modified n value

    # not modifing
    # cam_bramley_s: no change, we do not have EMMA data and it is unlikely to lead to any real changes (e.g. NOF bands)
    # cam at mashes as cam at bramleys
    # cust_skewbridge: no change, it is unlikely to lead to any real changes (e.g. NOF bands)
    # northbrook and southbrook at marshes: no change, it is unlikely to lead to any real changes (e.g. NOF bands)

    # the below modified from regressions created from n_vs_waimak_per values are in mg/l
    # courtaney at the kaiapoi
    # note the courtenay is possibly wrong zeb used SQ35169 which looks like a side trib
    #  when I think it should be SQ35170 which I think is the main courtenay

    'courtenay_kaiapoi_s': 6.33,  # from emma we expect a median of 8% waimak water

    # kaiapoi (silver stream) at harpers road
    'kaiapoi_harpers_s': 11.26,  # from emma we can expect a median of 23% waimak water

    # kaiapoi at island road
    'kaiapoi_island_s': 8.64,  # from emma we can expect a median of 20% waimak water

    # ohoka at island road
    'ohoka_island_s': 7.02  # from emma we can expect a median of 12% waimak water
}

stream_flows_gmp = {
    # key: percent flow under gmp from stocastic forward runs
    'cam_bramleys_s': .99,
    'cam_marshes_s': .99,
    'courtenay_kaiapoi_s': 1,
    'cust_skewbridge': 0.86,
    'kaiapoi_harpers_s': 1,  # link to silverstream at neeves
    'kaiapoi_island_s': 1,  # link to silverstream at neeves
    'northbrook_marshes_s': 1,
    'ohoka_island_s': 0.99,
    'southbrook_marshes_s': 0.99,
}


def apply_n_load_uncertainty(nvals, modifiers):
    nvals = np.atleast_1d(nvals)
    all_n = nvals[:, np.newaxis] * modifiers[np.newaxis, :]
    outdata = _np_describe(all_n)
    return outdata


def scale_gmp_for_flow_change(data, receptor_type):
    """
    scale the flows for gmp
    :param data: the averaged data for sites pd.Dataframe index site names, columns (see _np.describe)
    :param receptor_type: 'well' or 'stream'
    :return:
    """
    data = data.copy()
    if receptor_type == 'well':
        pass # decided instead to bracket with concentration and load... will run another mt3d run.
    elif receptor_type == 'stream':
        for site in data.index:
            try:
                modifier = stream_flows_gmp[site]
            except KeyError:
                modifier = 1
                warn('no key for {} in gmp flows!, setting modifier to 1'.format(site))

            data.loc[site] *= 1 / modifier
    else:
        raise ValueError('unexpected value {} for receptor type expected only "well" or "stream"'.format(receptor_type))
    return data


def output_actual_n_vals(outdir, mod_dir, gmp):
    """
    load the n modifiers, and applies them to the raw n data (calculated via Nitrate at key receptors) for both the
    stocastic set and AshOpt
    :param outdir: directory to save the output
    :param mod_dir: the direcotry that holds teh n_mods (generated by calc_all_ns)
    :param gmp: boolean, if True, adjust N values for GMP
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if os.path.exists(os.path.join(outdir, 'missing_sites.txt')):
        os.remove(os.path.join(outdir, 'missing_sites.txt'))
    # stocastic set
    # take raw data both N and modifier and put out values
    well_ids = get_well_ids()
    base_base_n_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results"
                               r"\ex_bd_va\n_results\n_results_at_points")

    # set up ashopt and stocastic set runs
    simsets = ['stocastic_set', 'ashopt']
    base_str_n_paths = [os.path.join(base_base_n_path, "raw_stocastic_set_str_data.csv"),
                        os.path.join(base_base_n_path, "AshOpt_stream_data.csv")]
    base_well_n_paths = [os.path.join(base_base_n_path, "raw_stocastic_set_well_data.csv"),
                         os.path.join(base_base_n_path, "AshOpt_grouped_well_data.csv")]

    for simset, base_well_n_path, base_str_n_path in zip(simsets, base_well_n_paths, base_str_n_paths):
        print(simset)
        if simset == 'ashopt':
            continue
            base_well_n = pd.read_csv(base_well_n_path, index_col=0, names=['n_con'])
        else:
            base_well_n = pd.read_csv(base_well_n_path, index_col=0).transpose()
        # well groups
        zone_sets = [set(well_ids.Zone[well_ids.Zone.notnull()]),
                     set(well_ids.Zone_1[well_ids.Zone_1.notnull()]),
                     set(well_ids.zone_2[well_ids.zone_2.notnull()])]
        outdata = {}
        for zone_set, key in zip(zone_sets, ['Zone', 'Zone_1', 'zone_2']):
            for zone in zone_set:
                try:
                    modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(zone)))
                except IOError as val:
                    with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                        f.write('missing: {}, exception: {}\n'.format(zone, val))
                        continue
                if key == 'Zone':  # handle some weird wdc overlap with private well names
                    use_zone = 'wdc_{}'.format(zone)
                else:
                    use_zone = zone

                if use_zone in well_zones_to_modify.keys():
                    nvals = well_zones_to_modify[use_zone]
                else:
                    idxs = well_ids.loc[well_ids[key] == zone].index.values
                    nvals = base_well_n.loc[idxs].mean().values
                outdata[use_zone] = apply_n_load_uncertainty(nvals, modifiers)
        outdata = pd.DataFrame(outdata).transpose()
        if not outdata.empty:
            if gmp:
                outdata = scale_gmp_for_flow_change(outdata, 'well')
            outdata.to_csv(os.path.join(outdir, 'full_{}_n_wells.csv'.format(simset)))

        # sfr_groups
        str_ids = get_str_ids() + ['waimakariri', 'conservative_interzone', 'highly_likely_interzone']
        base_str_n = pd.read_csv(base_str_n_path, index_col=0)
        if simset == 'ashopt':
            base_str_n = base_str_n.transpose()
        outdata = {}
        for str_id in str_ids:
            try:
                modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(str_id)))
            except IOError as val:
                with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                    f.write('missing: {}, exception: {}\n'.format(str_id, val))
                    continue
            if str_id in stream_zones_to_modify.keys():
                nvals = stream_zones_to_modify[str_id]
            elif str_id == 'waimakariri' or 'interzone' in str_id:
                with open(os.path.join(outdir,'READ_ME.txt'), 'a') as f:
                    f.write('there was no N data for {} set base value to 1 to look '
                            'at percentage changes'.format(str_id))
                nvals = np.array([1])
            else:
                nvals = base_str_n.loc[:, str_id].values
            outdata[str_id] = apply_n_load_uncertainty(nvals, modifiers)
        outdata = pd.DataFrame(outdata).transpose()
        if not outdata.empty:
            if gmp:
                scale_gmp_for_flow_change(outdata, 'stream')
            outdata.to_csv(os.path.join(outdir, 'full_{}_n_strs.csv'.format(simset)))


def run_all_nload_stuffs(base_outdir, szdirs, output_act_n=True):
    """
    wrapper to run all the nloads stocastics
    :param base_outdir: where to put the data
    :param szdirs: the directories for the source dirs
    :return:
    """
    if not os.path.exists(base_outdir):
        os.makedirs(base_outdir)
    szdirs = np.atleast_1d(szdirs)
    n_names = ['nload_cmp', 'nload_gmp', 'nconc_gmp', 'nconc_cmp']
    for sim_end in ['without_trans', 'with_trans']:  # with and without transition to CMP
        for sz_dir in szdirs:
            for n_name in n_names:
                if 'gmp' in n_name:
                    gmp = True
                else:
                    gmp = False
                print('starting N analysis for {} load, {} sims, and {} polygons'.format(n_name, sim_end,
                                                                                         os.path.basename(sz_dir)))
                outdir = os.path.join(base_outdir, sim_end, '{}_{}'.format(n_name, os.path.basename(sz_dir)))
                sims = pd.read_csv(env.gw_met_data(
                    "mh_modeling\stocastic_n_load_results\component_uncertainty_data_{}.csv".format(sim_end)),
                    index_col=0)
                if 'load' in n_name:
                    org_n_load_name = 'nload_cmp'
                elif 'con' in n_name:
                    org_n_load_name = 'nconc_cmp'
                else:
                    raise NotImplementedError
                # Calculate modifiers for the N load
                calc_all_ns(sims_org=sims, n_load_name=n_name, outdir=outdir, source_zone_dir=sz_dir,
                            org_n_load_name=org_n_load_name)

                if output_act_n: # hack to allow the waimakariri to run
                    output_actual_n_vals(outdir=outdir, mod_dir=outdir, gmp=gmp)
    if output_act_n:
        create_tabulated_results(base_outdir)


def create_tabulated_results(base_outdir):
    for trans, mset in itertools.product(['with_trans', 'without_trans'], ['stocastic_set']): # pulled out 'ashopt' for time
        gmp_path = glob(os.path.join(base_outdir, trans, '*load_gmp*', '*{}*.csv'.format(mset)))[0]
        cmp_path = glob(os.path.join(base_outdir, trans, '*load_cmp*', '*{}*.csv'.format(mset)))[0]
        gmp = pd.read_csv(gmp_path, index_col=0)
        cmp_n = pd.read_csv(cmp_path, index_col=0)
        outdata_load = pd.merge(gmp, cmp_n, right_index=True, left_index=True, suffixes=('_gmp_load', '_cmp_load'))
        gmp_path = glob(os.path.join(base_outdir, trans, '*conc_gmp*', '*{}*.csv'.format(mset)))[0]
        cmp_path = glob(os.path.join(base_outdir, trans, '*conc_cmp*', '*{}*.csv'.format(mset)))[0]
        gmp = pd.read_csv(gmp_path, index_col=0)
        cmp_n = pd.read_csv(cmp_path, index_col=0)
        outdata_con = pd.merge(gmp, cmp_n, right_index=True, left_index=True, suffixes=('_gmp_con', '_cmp_con'))
        outdata = pd.merge(outdata_load, outdata_con, right_index=True, left_index=True)
        outdata.to_csv(os.path.join(base_outdir, '{}_{}_overview.csv'.format(mset, trans)))


if __name__ == '__main__':
    run_all_nload_stuffs(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\stocastic_n_load_results\second_tranche_gmp_mod",
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche")
