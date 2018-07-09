# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 5/07/2018 1:23 PM
"""

from __future__ import division
from core import env
import os
from fpdf import FPDF



base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\step_projections"
imagelist = [
    # streams
    "cam_bramleys_s.png",
    "courtenay_kaiapoi_s.png",
    "cust_skewbridge.png",
    "kaiapoi_harpers_s.png",
    "kaiapoi_island_s.png",

    # interzone
    "conservative_interzone.png",
    "highly_likely_interzone.png",

    # WDC wells
    "wdc_Cust.png",
    "wdc_Fernside.png",
    "wdc_Kaiapoi.png",
    "wdc_Kairaki.png",
    "wdc_Mandeville.png",
    "wdc_Ohoka.png",
    "wdc_Oxford Urban.png",
    "wdc_Pegasus.png",
    "wdc_Poyntzs Road.png",
    "wdc_Rangiora.png",
    "wdc_Waikuku.png",
    "wdc_West Eyreton.png",

    # private wells
    "Clarkville.png",
    "Cust.png",
    "Eyreton_deep.png",
    "Eyreton_shallow.png",
    "Fernside.png",
    "Flaxton.png",
    "Horellville.png",
    "Mandeville.png",
    "North East Eyrewell_deep.png",
    "North East Eyrewell_shallow.png",
    "North West Eyrewell_deep.png",
    "North West Eyrewell_shallow.png",
    "Ohoka_deep.png",
    "ohoka_island_s.png",
    "Ohoka_shallow.png",
    "Rangiora.png",
    "Springbank.png",
    "Summerhill.png",
    "Swannanoa_deep.png",
    "Swannanoa_shallow.png",
    "Waikuku.png",
    "West Eyreton_deep.png",
    "West Eyreton_shallow.png",
    "Woodend - Tuahiwi.png",

]

pdf = FPDF() # units defined to mm here
# imagelist is the list with all image filenames
imagelist = [os.path.join(base_dir,e) for e in imagelist]
for image in imagelist:
    pdf.add_page()
    pdf.image(image,w=180)
pdf.output(os.path.join(base_dir,"grouped_pdf.pdf"), "F")
