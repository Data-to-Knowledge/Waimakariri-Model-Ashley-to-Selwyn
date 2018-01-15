# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:41:08 2017

@author: michaelek
"""
from collections import OrderedDict
from itertools import product

#######################################
### Dictionaries to relate the long names to the acronyms
### Need to be kept up-to-date compared to the SQL server
### If any new keys need to be added, add them to the end of the dictionaries!

feature_code_dict = OrderedDict((('River', 'river'), ('Aquifer', 'aq'), ('Atmosphere', 'atmos'), ('Soil', 'soil'), ('Lake', 'lake')))

mtype_code_dict = OrderedDict((('Water Level', 'wl'), ('Flow', 'flow'), ('Temperature', 'T'), ('Temperature Max', 'T_max'), ('Temperature Min', 'T_min'), ('Precipitation', 'precip'), ('Abstraction', 'abstr'), ('Net Radiation', 'R_n'), ('Shortwave Radiation', 'R_s'), ('Relative Humidity Min', 'RH_min'), ('Relative Humidity Max', 'RH_max'), ('Wind Speed', 'U_z'), ('Barometric Pressure', 'P'), ('ET Actual', 'ET'), ('ET Reference', 'ETo'), ('ET Potential', 'PET')))

msource_code_dict = OrderedDict((('Recorder', 'rec'), ('Manual Field', 'mfield'), ('Manual Lab', 'mlab')))

qual_state_code_dict = OrderedDict((('RAW', 'raw'), ('Quality Controlled', 'qc')))

## Resample functions
resample_mtype_fun = {'flow': 'mean', 'wl': 'mean', 'precip': 'sum', 'R_n': 'sum', 'R_s': 'sum', 'T_min': 'mean', 'T_max': 'mean', 'T': 'mean', 'RH_min': 'mean', 'RH_max': 'mean', 'U_z': 'mean', 'P': 'mean', 'abstr': 'sum', 'PET': 'sum'}

## Mtype codes with numeric values
#mtype_codes = {'river_wl_rec_qc': 1, 'river_flow_rec_qc': 2, 'river_temp_rec_qc': 3, 'aq_wl_rec_qc': 4, 'atmos_precip_rec_qc': 5, 'lake_wl_rec_qc': 6}

#######################################
### Reorganisation of the above dictionaries

key_fields = OrderedDict((('Feature', feature_code_dict), ('Mtype', mtype_code_dict), ('MeasurementSource', msource_code_dict), ('QualityState', qual_state_code_dict)))

code_fields = OrderedDict((i, key_fields[i].values()) for i in key_fields)
list_combo = list(product(*code_fields.values()))

feature_mtype = OrderedDict((('river', ('wl', 'flow', 'abstr', 'T')), ('aq', ('wl', 'abstr', 'T')), ('atmos', ('precip', 'R_n', 'R_s', 'T_min', 'T_max', 'RH_min', 'RH_max', 'U_z', 'P', 'T')), ('soil', ('T', 'PET')), ('lake', ('wl', 'T'))))

list_feature_mtype = [list(product([i], feature_mtype[i])) for i in feature_mtype]

fields_lst = [i for i in list_combo if i[1] in feature_mtype[i[0]]]

resample_fun = OrderedDict(('_'.join(i), resample_mtype_fun[i[1]]) for i in fields_lst)

legacy_mtypes = {'usage': 'sum'}

resample_fun.update(legacy_mtypes)

#all_mtypes = OrderedDict((i, '_'.join(fields_lst[i-1])) for i in range(1, len(fields_lst)))

all_mtypes = ['river_flow_cont_qc', 'river_flow_disc_qc', 'river_wl_cont_qc', 'aq_wl_cont_qc', 'aq_wl_disc_qc', 'atmos_precip_cont_qc', 'river_flow_cont_raw', 'atmos_precip_cont_raw', 'river_wl_cont_raw', 'aq_wl_cont_raw', 'usage']

######################################
### The main class

class hydro(object):
    """
    A class to handle environmental time series data where a site has a measurement type, a time series, and a location.
    """
    from core.classes.hydro.import_fun import add_geo_loc, missing_geo_loc_sites, add_data, _import_attr, add_geo_catch, _add_geo_data, rd_csv, rd_netcdf, _rd_hydro_mssql, combine, _rd_hydro_geo_mssql, _proc_hydro_sql
    from core.classes.hydro.indexing import sel_ts, sel_sites_by_poly, sel_ts_by_poly, sel, sel_by_poly, __getitem__, _comp_by_buffer, _comp_by_catch, sel_by_geo_attr
    from core.classes.hydro.misc import _check_mtypes_sites, _base_stats_fun, _mtype_check
    from core.classes.hydro.export_fun import to_csv, to_netcdf, to_shp
    from core.classes.hydro.ecan_import import get_geo_loc, _rd_hydstra, _rd_hydrotel, _rd_henry, get_data, get_site_geo_attr
    from copy import copy
    from core.classes.hydro.tools.sw import malf7d, flow_reg
    from core.classes.hydro.tools.general import resample, stats
    from core.classes.hydro.tools.plot import plot_hydrograph, plot_reg
    from core.classes.hydro.tools.gw import gwl_reg

#    @property
#    def _constructor(self):
#        return(hydro)
    ### General attributes

    ### Initial import and assignment function
    def __init__(self, data=None, time=None, sites=None, mtypes=None, values=None, dformat=None):
        if data is None:
            pass
        else:
            ## Read in data
            self.add_data(data=data, time=time, sites=sites, mtypes=mtypes, values=values, dformat=dformat)

    ### Call
#    def __call__(self, data=None, time=None, sites=None, mtypes=None, values=None, dformat=None):
#        self.add_data(data=data, time=time, sites=sites, mtypes=mtypes, values=values, dformat=dformat)

    ### What to return when the oject is called alone
    def __repr__(self):
        if hasattr(self, 'data'):
            if not hasattr(self, '_base_stats'):
                self._base_stats_fun()
            return(repr(self._base_stats))
        else:
            print("There's no data here. Add some in.")




