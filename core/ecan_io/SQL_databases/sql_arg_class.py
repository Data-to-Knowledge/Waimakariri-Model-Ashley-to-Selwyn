# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 11:10:04 2017

@author: MichaelEK
"""

#################################################################
### Helper functions


def _get_vars(dict1):
    v1 = [i for i in dict1 if '__' not in i]
    return(v1)

##################################################################
### Class dictionary containers


class consents(object):
    """
    Consents data.
    """
    crc_dates = {'server': 'SQL2012PROD03', 'database': 'DataWarehouse', 'table': 'F_ACC_Permit', 'col_names': ['B1_ALT_ID', 'fmDate', 'toDate', 'Longitude', 'Latitude', 'StatusType', 'B1_APPL_STATUS'], 'rename_cols': ['crc', 'from_date', 'to_date', 'wgs84_x', 'wgs84_y', 'status', 'status_details']}
    crc_relation = {'server': 'SQL2012PROD03', 'database': 'DataWarehouse', 'table': 'D_ACC_Relationships', 'col_names': ['ParentRecordNo', 'ChildRecordNo'], 'rename_cols': ['crc', 'crc_child']}
    crc_wap_allo = {'server': 'SQL2012PROD03', 'database': 'DataWarehouse', 'table': 'D_ACC_Act_Water_TakeWaterWAPAlloc', 'col_names': ['RecordNo', 'Activity', 'Allocation Block', 'ROW_INDEX', 'WAP', 'Max Rate for WAP (l/s)', 'Max Rate Pro Rata (l/s)', 'Max Vol Pro Rata (m3)', 'Consecutive Day Period', 'Include in SW Allocation?', 'First Stream Depletion Rate', 'From Month'], 'rename_cols': ['crc', 'take_type', 'allo_block', 'row_index', 'wap', 'max_rate_wap', 'max_rate', 'max_vol', 'return_period', 'in_sw_allo', 'sd', 'from_month']}
    crc_cav = {'server': 'SQL2012PROD03', 'database': 'DataWarehouse', 'table': 'D_ACC_Act_Water_TakeWaterConsent', 'col_names': ['RecordNo', 'Activity', 'Consented Annual Volume (m3/year)', 'Has a low flow restriction condition?'], 'rename_cols': ['crc', 'take_type', 'cav', 'min_flow']}
    crc_use_type = {'server': 'SQL2012PROD03', 'database': 'DataWarehouse', 'table': 'D_ACC_Act_Water_TakeWaterAllocData', 'col_names': ['RecordNo', 'Activity', 'Allocation Block', 'ROW_INDEX', 'Include in GW Allocation?', 'Water Use', 'Irrigation Area (ha)', 'Full Effective Annual Volume (m3/year)'], 'rename_cols': ['crc', 'take_type', 'allo_block', 'row_index', 'in_gw_allo', 'use_type', 'irr_area', 'feav']}


class water_use(object):
    """
    Water usage data.
    """
    wus_mon = {'server': 'SQL2012PROD03', 'database': 'WUS', 'table': 'vw_WUS_Fact_MonthlyUsageByUsageSite_MsAccess', 'col_names': ['UsageSite', 'Year', 'Month', 'Usage'], 'rename_cols': ['wap', 'year', 'month', 'usage']}
    wus_day = {'server': 'SQL2012PROD03', 'database': 'WUS', 'table': 'vw_WUS_Fact_DailyUsageByUsageSite', 'col_names': ['UsageSite', 'Day', 'Usage'], 'rename_cols': ['wap', 'date', 'usage']}


class GW(object):
    """
    Groundwater related tables.
    """
    sd = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'Well_StreamDepletion_Locations', 'col_names': ['Well_No', 'NZTMXwell', 'NZTMYwell', 'NZTMX', 'NZTMY', 'Distance', 'SD1_7', 'SD1_30', 'SD1_150', 'NZTMX2', 'NZTMY2', 'Distance2', 'SD2_7', 'SD2_30', 'SD2_150'], 'rename_cols': ['wap', 'wap_x', 'wap_y', 'sd1_x', 'sd1_y', 'dist1', 'sd1_7', 'sd1_30', 'sd1_150', 'sd2_x', 'sd2_y', 'dist2', 'sd2_7', 'sd2_30', 'sd2_150']}
    well_details = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'WELL_DETAILS', 'col_names': ['WELL_NO', 'WELL_TYPE', 'Well_Status', 'DEPTH', 'DIAMETER', 'NZTMX', 'NZTMY'], 'rename_cols': ['wap', 'well_type', 'well_status', 'depth', 'dia', 'NZTMX', 'NZTMY']}


class GIS(object):
    """
    Tables with geometry.
    """
    swaz_gis = {'server': 'SQL2012PROD05', 'database': 'GISPUBLIC', 'table': 'PLAN_NZTM_SURFACE_WATER_ALLOCATION_ZONES', 'col_names': ['ZONE_GROUP_NAME', 'ZONE_NAME'], 'rename_cols': ['swaz_grp', 'swaz'], 'geo_col': True}
    gwaz_gis = {'server': 'SQL2012PROD05', 'database': 'GISPUBLIC', 'table': 'PLAN_NZTM_GROUNDWATER_ALLOCATION_ZONES', 'col_names': ['Zone_Name'], 'rename_cols': ['gwaz'], 'geo_col': True}
    catch_gis = {'server': 'SQL2012PROD05', 'database': 'GISPUBLIC', 'table': 'SWATER_NZTM_CATCHMENTS', 'col_names': ['NIWACatchmentNumber', 'HydrologicalCatchmentName', 'CatchmentGroup', 'CatchmentGroupName', 'OldHydrologicalCatchmentNumber'], 'rename_cols': ['niwa_catch', 'catch_name', 'catch_grp', 'catch_grp_name', 'old_catch'], 'geo_col': True}
    catch_grp_gis = {'server': 'SQL2012PROD05', 'database': 'GISPUBLIC', 'table': 'SWATER_NZTM_CATCHMENT_GROUPS', 'col_names': ['CatchmentGroup', 'CatchmentGroupName'], 'rename_cols': ['catch_grp', 'catch_grp_name'], 'geo_col': True}
    cwms_gis = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'CWMS_NZTM_ZONES', 'col_names': ['ZONE_NAME'], 'rename_cols': ['cwms'], 'geo_col': True}
    crc_gis = {'server': 'SQL2012PROD05', 'database': 'GIS_Accela', 'table': 'vCOCOA_NZTM_ResourceConsents', 'col_names': ['ConsentNo', 'ConsentType', 'NZTMX', 'NZTMY'], 'rename_cols': ['crc', 'crc_type', 'NZTMX', 'NZTMY'], 'geo_col': True}
    rec_gis = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'MFE_NZTM_REC', 'col_names': ['NZREACH', 'NZFNODE', 'NZTNODE', 'ORDER_'], 'rename_cols': ['NZREACH', 'NZFNODE', 'NZTNODE', 'order'], 'geo_col': True}
    rec_catch_gis = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'MFE_NZTM_RECWATERSHEDCANTERBURY', 'col_names': ['NZREACH'], 'geo_col': True}


class lowflows(object):
    """
    Tables in the LowFlows database.
    """
    lowflow_restr_day = {'server': 'SQL2012PROD03', 'database': 'LowFlows', 'table': 'LowFlowSiteRestrictionDaily', 'col_names': ['SiteID', 'RestrictionDate', 'BandNo', 'BandAllocation'], 'rename_cols': ['lowflow_id', 'date', 'band', 'band_restr'], 'where_col': {'SnapshotType': ['Live']}}
    lowflow_restr_type = {'server': 'SQL2012PROD03', 'database': 'LowFlows', 'table': 'LowFlowSiteBandPeriodAllocation', 'col_names': ['SiteID', 'BandNo', 'PeriodNo', 'Allocation', 'Flow'], 'rename_cols': ['lowflow_id', 'band', 'mon', 'allo', 'min_flow']}
    lowflow_band_crc = {'server': 'SQL2012PROD03', 'database': 'LowFlows', 'table': 'vLowFlowConsents2', 'col_names': ['SiteID', 'BandNo', 'RecordNo', 'isCurrent'], 'rename_cols': ['lowflow_id', 'band', 'crc', 'active']}
    lowflow_gauge = {'server': 'SQL2012PROD03', 'database': 'LowFlows', 'table': 'LowFlowSite', 'col_names': ['SiteID', 'RefDBase', 'RefDBaseKey', 'isActive'], 'rename_cols': ['lowflow_id', 'DB', 'site', 'active']}


class SW(object):
    """
    Surface water tables.
    """
    hydro_sites = {'server': 'SQL2012PROD05', 'database': 'Bgauging', 'table': 'RSITES', 'col_names': ['SiteNumber', 'River', 'SiteName', 'RecorderSite', 'RainfallSite', 'FlowOrStage', 'Telemetered', 'MinimumFlowSite', 'NZTMX', 'NZTMY'], 'rename_cols': ['site', 'river', 'site_name', 'rec_site_bool', 'rain_site_bool', 'flow_stage', 'tel_bool', 'min_flow_bool', 'NZTMX', 'NZTMY']}



##################################################################
### Main class call

class sql_arg(object):
    """
    Class container to retreive dictionaries to pass to rd_sql.
    """
    from core.ecan_io.SQL_databases.sql_arg_class import _get_vars, consents, water_use, GW, SW, GIS, lowflows

    _code_dict = {'consents': consents, 'water_use': water_use, 'GW': GW, 'lowflows': lowflows, 'SW': SW, 'GIS': GIS}

    def __init__(self):
        classes = {i: _get_vars(self._code_dict[i].__dict__) for i in self._code_dict}
        codes = []
        [codes.extend(classes[i]) for i in classes]
        self.classes = classes
        self.codes = codes

    def get_dict(self, code=None, fields=None, minimum=False):
        """
        Main function to retreive dictionaries to pass to rd_sql.

        Arguments:
        code -- Special dictionary code to retreive the dictionary parameters for SQL (str).\n
        fields -- The dictionary keys that should be retreived (str or list). None returns all fields.\n
        minimum -- Should only the server, database, and table values be returned?
        """
        if isinstance(code, str):
            if code in self.codes:
                class1 = [i for i in self.classes if code in self.classes[i]][0]
                class2 = getattr(self, class1)
                dict1 = getattr(class2, code)
                if fields is not None:
                    if isinstance(fields, str):
                        fields = [fields]
                    if isinstance(fields, list):
                        dict1 = {i: dict1[i] for i in fields if i in dict1.keys()}
                if minimum:
                    fields = ['server', 'database', 'table']
                    dict1 = {i: dict1[i] for i in fields if i in dict1.keys()}

                return(dict1.copy())

    def __repr__(self):
        return(repr(self.classes))



































