# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 8:37 AM
"""
from __future__ import division
from core import env
import os
import tempfile
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.percentage_reduction_maps import get_pa_reductions
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt


if __name__ == '__main__':
    'scenario     |     50%    |    95%\ncmp          |{:8.2f}    |{:8.2f}\ngmp          |{:8.2f}    |{:8.2f}\nchch8        |{:8.2f}    |{:8.2f}\ninterzone8   |{:8.2f}    |{:8.2f}\n'.format(1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1)

    'scenario     |     50%    |    95%\ncmp          |{:8.2f}    |{:8.2f}\ngmp          |{:8.2f}    |{:8.2f}\nchch8        |{:8.2f}    |{:8.2f}\ninterzone8   |{:8.2f}    |{:8.2f}\n'.format( float([cmp_50]), float([cmp_95]) , float([gmp_50]), float([gmp_95]), float([chch50]), float([chch95]), float([inter50]), float([inder95]))