import os
import sys
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pandas as pd
import re
import math

CMAP_CSV_DIR = "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps"

DATA_DICT = { "xco2" : {"range": [380, 430], "cmap" : "cm.viridis", "binsize": 0.2}, 
              "xco2_relative" : {"range": [-6, 6], "cmap" : "cm.RdBu_r", "binsize": 0.05}, 
              "tcwv" : {"range": [0, 75], "cmap" : "cm.Blues", "binsize": 1/3.}, 
              "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}, 
              "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}, 
              "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}
            }


def truncate(n, d):
    return math.floor(n * 10 ** d) / 10 ** d

for var in DATA_DICT.keys():
    #print(var)
    cmap_name = os.path.join(CMAP_CSV_DIR, var + "_" + re.split("\.", DATA_DICT[var]["cmap"])[-1] + "_" +  str(DATA_DICT[var]["range"][0]) + "to" + str(DATA_DICT[var]["range"][1]) + ".csv")
    
    data_crange_low = [i*DATA_DICT[var]["binsize"] for i in np.arange(DATA_DICT[var]["range"][0] / DATA_DICT[var]["binsize"], DATA_DICT[var]["range"][1] / DATA_DICT[var]["binsize"] + DATA_DICT[var]["binsize"], 1)][:-1]
    data_crange_high = [i*DATA_DICT[var]["binsize"] for i in np.arange(DATA_DICT[var]["range"][0] / DATA_DICT[var]["binsize"], DATA_DICT[var]["range"][1] / DATA_DICT[var]["binsize"] + DATA_DICT[var]["binsize"], 1)][1:]
    
    ncolors = len(data_crange_low)

    cmap_arr = np.array([[np.append(np.round(255*mpl.colors.to_rgba_array(eval(DATA_DICT[var]["cmap"])(n))), [round(data_crange_low[n], 3), round(data_crange_high[n], 3)])] for n in range(ncolors)]).squeeze()
    cmap_df = pd.DataFrame(cmap_arr, columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"]).astype({"red": int, "green": int, "blue": int, "alpha": int})
    
    if var == "tcwv":
        cmap_df["data_lim_low"] = cmap_df["data_lim_low"].map(lambda x: truncate(x, 2))
        cmap_df["data_lim_high"] = cmap_df["data_lim_high"].map(lambda x: truncate(x, 2))

    cmap_df.to_csv(cmap_name, index=False)
