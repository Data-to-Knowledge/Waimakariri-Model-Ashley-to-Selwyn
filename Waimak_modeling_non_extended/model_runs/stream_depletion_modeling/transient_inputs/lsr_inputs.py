"""
Author: matth
Date Created: 22/05/2017 10:26 AM
"""

from __future__ import division
import pandas as pd
from users.MH.Waimak_modeling.supporting_data_path import sdp

def get_mean_monthly_lsr():
    """
    used a bucket model, zeb put together the average profile.
    :return: 
    """
    lsr_scale_path = "{}\inputs\LSR daily Thorpe & Scott Waimak zone average.xls".format(sdp)
    lsr_scale = pd.read_excel(lsr_scale_path,sheetname='for_model')
    lsr_scale = lsr_scale.set_index('Month')
    return lsr_scale['factor']

if __name__ == '__main__':
    lsr = get_mean_monthly_lsr()
    print 'done'