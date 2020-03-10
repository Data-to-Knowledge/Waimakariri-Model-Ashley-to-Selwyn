"""
Author: matth
Date Created: 19/05/2017 9:01 AM
"""

from __future__ import division
from users.MH.Waimak_modeling.supporting_data_path import temp_file_dir
import users.MH.Waimak_modeling.model_tools as mt


def main():
    wells = mt.get_original_well_data()
    wells = wells[wells['flux'] < 0]
    m = mt.get_base_mf_ss()
    fig, ax, mmap = mt.plt_default_map(m,0,title='Abstraction wells')
    p = ax.scatter(wells['x'], wells['y'], c=wells['flux'],cmap='plasma_r', vmin=0, vmax=wells['flux'].quantile(0.05))
    fig.colorbar(p, ax=ax, extend='max')
    fig.savefig('{}/well_abstraction_map.png'.format(temp_file_dir))
    mt.plt.close(fig)



if __name__ == '__main__':
    main()