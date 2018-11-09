import os
import sys
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pandas as pd
import operator

CMAP_CSV_DIR = "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps"
NCOLORS=256

DATA_DICT = { "LtCO2" : 
                { "xco2" : {"data_field_name" : "xco2", "preprocessing" : False, "range": [380, 430], "cmap" : "viridis", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "xco2_relative" : {"data_field_name" : None, "preprocessing" : "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_trend_gl.txt", "range": [-6, 6], "cmap" : "RdBu_r", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "tcwv" : {"data_field_name" : "Retrieval/tcwv", "preprocessing" : False, "range": [0, 75], "cmap" : "Blues", "quality_info" : {}}, 
                },
             "LtSIF" : 
                { "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}, 
                  "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}, 
                  "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}
                }
            }


for key, v in DATA_DICT.items():
    for var in v.keys():

        cmap_name = os.path.join(CMAP_CSV_DIR, var + "_" + DATA_DICT[key][var]["cmap"] + "_" +  str(DATA_DICT[key][var]["range"][0]) + "to" + str(DATA_DICT[key][var]["range"][1]) + ".csv")
        data_step = (DATA_DICT[key][var]["range"][1] - DATA_DICT[key][var]["range"][0])/float(NCOLORS)

        data_crange_low = np.arange(DATA_DICT[key][var]["range"][0], DATA_DICT[key][var]["range"][1] + data_step, data_step, dtype=float)[:-1]
        data_crange_high = np.arange(DATA_DICT[key][var]["range"][0], DATA_DICT[key][var]["range"][1] + data_step, data_step, dtype=float)[1:]

        cmap_arr = np.array([[np.append(mpl.colors.to_rgba_array(cm.viridis(n)),[data_crange_low[n], data_crange_high[n]])] for n in range(NCOLORS)]).squeeze()
        cmap_df = pd.DataFrame(cmap_arr, columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"])
        
        cmap_df.to_csv(cmap_name, index=False)
