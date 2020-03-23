"""
 Author: Matt Hanson
 Created: 3/3/2020 11:43 AM
 """

# just a test # second test

import netCDF4 as nc
import numpy as np

if __name__ == '__main__':
    data = nc.Dataset(
        r"D:\ecan_data_delete\Numerical GW model\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    filters = [
        'run_mt3d',
        'emma_converge',
        'n_converge',
        'emma_no_wt',
        'emma_eq_wt',
        'emma_chch_wt',
        'emma_str_wt',
        'emma_ewf_wt',
    ]

    nsmc = np.array(data.variables['nsmc_num'])
    nsmc_nums = [5, 17, 18, 26, 37, 44, 72, 103, 117, 133, 142,
                 204, 233, 240, 258, 271, 278, 308, 314, 388, 391, 439,
                 441, 447, 488, 491, 552, 555, 604, 640, 666, 790, 794,
                 800, 808, 809, 871, 883, 932, 952, 955, 992, 1001, 1029,
                 1052, 1055, 1078, 1099, 1106, 1111, 1112, 1114, 1129, 1168, 1181,
                 1250, 1326, 1328, 1330, 1349, 1365, 1368, 1370, 1395, 1399, 1454,
                 1481, 1489, 1498, 1499, 1507, 1574, 1585, 1595, 1618, 1636, 1656,
                 1660, 1681, 1694, 1749, 1773, 1817, 1819, 1836, 1842, 1851, 1870,
                 1881, 1897, 1915, 1916, 1967, 2044, 2052, 2068, 2085, 2114, 2117,
                 2139, 2160, 2200, 2209, 2234, 2278, 2382, 2390, 2430, 2460, 2520,
                 2579, 2588, 2590, 2604, 2643, 2654, 2661, 2673, 2691, 2735, 2748,
                 2752, 2757, 2760, 2818, 2840, 2871, 2875, 2885, 2918, 2971, 2980,
                 2997, 3039, 3064, 3067, 3075, 3100, 3116, 3165, 3216, 3289, 3310,
                 3318, 3350, 3378, 3389, 3418, 3428, 3456, 3464, 3482, 3513, 3585,
                 3641, 3651, 3663, 3685, 3712, 3717, 3762, 3796, 3836, 3861, 3910]

    for f in filters:
        print(f)
        t = np.array(data.variables[f])
        t[t<0] = 0
        print((nsmc_nums==nsmc[t.astype(bool)]))
        print('\n\n')

        # turns out it was emma no wieght