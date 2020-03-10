"""
Author: matth
Date Created: 11/05/2017 2:50 PM
"""

from __future__ import division
import flopy
import matplotlib.pyplot as plt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir

def main():
    model_path = '{}/well_flow_paths'.format(base_mod_dir)
    m = flopy.modflow.Modflow.load('{}/flow_paths.nam'.format(model_path))

    end_path = r"C:\Users\MattH\Desktop\Waimak_modeling\python_models\well_flow_paths\wdc_forward_wells_weak_s.mpend"
    path_path = r"C:\Users\MattH\Desktop\Waimak_modeling\python_models\well_flow_paths\wdc_forward_wells_weak_s.mppth"

    end_pts = flopy.utils.EndpointFile(end_path).get_alldata()
    path_lines = flopy.utils.PathlineFile(path_path).get_alldata()

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')

    ax.set_title('plot_array()')
    modelmap = flopy.plot.ModelMap(model=m,)
    modelmap.plot_bc("WEL")
    modelmap.plot_bc('STR')
    modelmap.plot_bc('DRN')
    quadmesh = modelmap.plot_ibound()
    linecollection = modelmap.plot_grid()
    modelmap.plot_endpoint(end_pts, direction='ending',
                               zorder=100)

    # plot the pathlines
    modelmap.plot_pathline(path_lines, layer='all', colors='red', travel_time=None)
    plt.show()
    print('done')


#load endpoint and

if __name__ == '__main__':
    main()


