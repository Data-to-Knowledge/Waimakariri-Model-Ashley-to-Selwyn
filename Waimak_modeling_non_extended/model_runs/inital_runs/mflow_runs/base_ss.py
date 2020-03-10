"""
Author: matth
Date Created: 18/04/2017 3:34 PM
"""

from __future__ import division
import users.MH.Waimak_modeling.model_tools as mt


test_dir = r"C:\Users\MattH\Desktop\test_dir"
m = mt.wraps.mflow.import_gns_model('test2', test_dir, safe_mode=False,mt3d_link=True)
m.write_input()
m.write_name_file()
print(m.run_model())





