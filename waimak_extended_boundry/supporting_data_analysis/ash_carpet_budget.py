from __future__ import division
import pandas as pd
from waimak_extended_boundry import smt
from waimak_extended_boundry.m_packages.drn_packages import _get_drn_spd
import flopy_mh as flopy


def ash_carpet_budget(path):
    drn_data = _get_drn_spd(1, 1)
    cbb = flopy.utils.CellBudgetFile(path)
    flux = cbb.get_data(kstpkper=cbb.get_kstpkper()[-1], text='drain', full3D=True)[0][0]
    ncarpetj = (smt.df_to_array(drn_data.loc[drn_data.group == 'ash_carpet'], 'j'))
    ncarpeti = (smt.df_to_array(drn_data.loc[drn_data.group == 'ash_carpet'], 'i'))

    outdata = {}
    for target_id, (imin, imax), (jmin, jmax) in zip(['ne_ash', 'nw_ash', 'se_ash', 'sw_ash'],
                                                     [(0, 48), (0, 48), (49, 300), (49, 300)],  # i's
                                                     [(295, 364), (0, 294), (295, 364), (0, 294)]):  # j's
        temp = (ncarpeti <= imax) & (ncarpeti >= imin) & (ncarpetj <= jmax) & (ncarpetj >= jmin)
        outdata[target_id] = flux[temp].sum()

    outdata = pd.DataFrame({'flux': outdata}) / 86400
    print(outdata)

if __name__ == '__main__':
    paths = [r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\VertUnstabA\AW20171022_i3_optver\i3\mf_aw_ex.cbc",
             r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\VertUnstabB\AW20171022_i4_optver\i4\mf_aw_ex.cbc"]
    for path in paths:
        print(path)
        ash_carpet_budget(path)