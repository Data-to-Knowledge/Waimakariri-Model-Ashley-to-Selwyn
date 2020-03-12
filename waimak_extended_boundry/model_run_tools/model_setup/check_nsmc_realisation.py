# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 21/12/2017 2:36 PM
"""

from __future__ import division
from realisation_id import get_model
import numpy as np
import pandas as pd

if __name__ == '__main__':
    m_org = get_model('NsmcBaseB')
    m_new = get_model('NsmcReal{:06d}'.format(-2), save_to_dir=True)

    # check well
    temp = np.isclose(pd.DataFrame(m_new.wel.stress_period_data.data[0]), pd.DataFrame(m_org.wel.stress_period_data.data[0])).all()
    print('are all well numbers the same: {}'.format(temp)) # some differences due to a rounding error in the LRZF

    # check rch
    temp = np.isclose(m_org.rch.rech.array, m_new.rch.rech.array).all()
    print('are all rch numbers the same: {}'.format(temp))

    # check kv
    temp = np.isclose(m_org.upw.vka.array, m_new.upw.vka.array).all()
    print('are all kv numbers the same: {}'.format(temp))

    # check kh
    temp = np.isclose(m_org.upw.hk.array, m_new.upw.hk.array).all()
    print('are all kh numbers the same: {}'.format(temp))

    # check sfr reach data
    temp = np.isclose(pd.DataFrame(m_org.sfr.reach_data), pd.DataFrame(m_new.sfr.reach_data)).all()
    print('are all sfr reach numbers the same: {}'.format(temp))

    # check sfr seg data
    temp = (m_org.sfr.segment_data[0] == m_new.sfr.segment_data[0]).all()
    print('are all sfr seg numbers the same: {}'.format(temp))

    # check drain data
    temp = np.isclose(pd.DataFrame(m_org.drn.stress_period_data.data[0]), pd.DataFrame(m_new.drn.stress_period_data.data[0])).all()
    print('are all drn numbers the same: {}'.format(temp))
