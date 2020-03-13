"""
 Author: Matt Hanson
 Created: 3/11/2020 11:05 AM
 """

from __future__ import division
from env.env_paths import  sci, transfers
import os
import socket

# todo document this when it is all said and done
sdp_required = r"D:\Waimakariri_model_input_data\required"

sdp_recommended = r"D:\Waimakariri_model_input_data\recommended"

temp_file_dir =r"C:\Users\Matt Hanson\Downloads\temp_waimak_files"


if not os.path.exists(temp_file_dir):
    os.makedirs(temp_file_dir)
