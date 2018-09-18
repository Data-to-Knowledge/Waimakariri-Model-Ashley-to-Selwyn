from __future__ import division
from glob import glob
import shutil
import os
import pandas as pd

if __name__ == '__main__':
    base_dir = r"D:\Ashley-Waimakariri_model\final_stocastic_set"
    for m in base_dir:
        model_dir = os.path.join(base_dir, m)
        for path in os.listdir(model_dir):
            if 'NsmcReal' in path:
                continue
            os.remove(os.path.join(model_dir,path))
    print('done')