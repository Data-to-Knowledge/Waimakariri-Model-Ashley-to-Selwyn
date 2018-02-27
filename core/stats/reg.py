# -*- coding: utf-8 -*-
"""
Functions for regressions.
"""


def lin_reg(x, y, log_x=False, log_y=False):
    """
    Function to estimate the linear regression of two data frames/series.

    Arguments:\n
    x, y -- Either data frames or series.\n
    log_axes -- Should the axes be log transformed for the regression?
    """
    from scipy import stats, log, exp, mean, sqrt
    from pandas import DataFrame, concat
    from numpy import nan

    col_names = ["X_loc", "Y_loc", "Slope", "Intercept", "R2", "MANE", "NMAE", "NRMSE", "RMSE", "p-value", "n-obs", "max_y"]
    x1 = DataFrame(x)
    y1 = DataFrame(y)
    x_names = x1.columns.values
    y_names = y1.columns.values
    x_len = len(x_names)
    y_len = len(y_names)

    output1 = [DataFrame(index=y_names, columns=col_names) for i in range(x_len)]

    for j in range(x_len):
        for i in range(y_len):
            both1 = concat([x1.iloc[:, j], y1.iloc[:, i]], axis=1).dropna()
            if both1.size > 0:
                if log_x and log_y:
                    both2 = both1.apply(log)
                    slope, inter, r_val, p_val, rmse = stats.linregress(both2.iloc[:, 0], both2.iloc[:, 1])
                    est_val = exp(slope*both2.iloc[:, 0] + inter)
                elif log_x:
                    x_log = both1.iloc[:, 0].apply(log)
                    slope, inter, r_val, p_val, rmse = stats.linregress(x_log, both1.iloc[:, 1])
                    est_val = slope*x_log + inter
                elif log_y:
                    y_log = both1.iloc[:, 1].apply(log)
                    slope, inter, r_val, p_val, rmse = stats.linregress(both1.iloc[:, 0], y_log)
                    est_val = exp(slope*both1.iloc[:, 0] + inter)
                else:
                    both2 = both1
                    slope, inter, r_val, p_val, rmse = stats.linregress(both2.iloc[:, 0], both2.iloc[:, 1])
                    est_val = slope*both2.iloc[:, 0] + inter

                mane = mean(abs(est_val - both1.iloc[:, 1])/both1.iloc[:, 1])
                nmae = mean(abs(est_val - both1.iloc[:, 1]))/mean(both1.iloc[:, 1])
                nrmse = sqrt(mean((est_val - both1.iloc[:, 1])**2))/mean(both1.iloc[:, 1])
                n_obs = len(both1.iloc[:, 1])
                max_y = max(both1.iloc[:, 1])

                output1[j].iloc[i, :] = [x_names[j], y_names[i], round(slope, 5),
                    round(inter, 5), round(r_val, 3), round(mane, 3),
                    round(nmae, 3), round(nrmse, 3), round(rmse, 3), round(p_val, 3), n_obs, max_y]
            else:
                 output1[j].iloc[i, :] = nan

    return(output1)
