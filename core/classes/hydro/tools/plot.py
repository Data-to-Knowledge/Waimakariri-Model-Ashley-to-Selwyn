# -*- coding: utf-8 -*-
"""
Functions to plot data.
"""

from core.ts.plot import hydrograph_plot, reg_plot

#############################################
#### Flow related


def plot_hydrograph(self, flow_sites=None, precip_sites=None, x_period='day', x_n_periods=1, time_format='%d-%m %H', x_rot=45, alpha=0.6, x_grid=False, start=None, end=None, export_path=None):

    if flow_sites is None:
        flow = self.sel_ts(mtypes='flow', pivot=True, start=start, end=end)
    elif isinstance(flow_sites, (int, str, list)):
        flow = self.sel_ts(mtypes='flow', sites=flow_sites, pivot=True, start=start, end=end)
#    else:
#        raise ValueError('If the hydro class object contains more than one site, then flow_sites must be passed as an int, str, or list.')

    if precip_sites is None:
        precip = None
    elif isinstance(precip_sites, (int, str, list)):
        precip = self.sel_ts(mtypes='precip', sites=precip_sites, pivot=True, start=start, end=end)

    plt1 = hydrograph_plot(flow=flow, precip=precip, x_period=x_period, x_n_periods=x_n_periods, time_format=time_format, x_rot=x_rot, alpha=alpha, x_grid=x_grid, export_path=export_path)

    return(plt1)


def plot_reg(self, x_mtype, x_site, y_mtype, y_site, freq='day', n_periods=1, fun='mean', min_ratio=0.75, digits=3, x_max=None, y_max=None, logs=False, export=False, export_path='reg.png'):
    from numpy import log

    x = self.sel_ts(mtypes=x_mtype, sites=x_site, pivot=True)
    y = self.sel_ts(mtypes=y_mtype, sites=y_site, pivot=True)

    plt1, reg1 = reg_plot(x, y, freq=freq, n_periods=n_periods, fun=fun, min_ratio=min_ratio, digits=digits, x_max=x_max, y_max=y_max, logs=logs, export=export, export_path=export_path)

    return(plt1, reg1)




#######################################
#### Documentation for the functions

plot_hydrograph.__doc__ = hydrograph_plot.__doc__
plot_reg.__doc__ = reg_plot.__doc__


