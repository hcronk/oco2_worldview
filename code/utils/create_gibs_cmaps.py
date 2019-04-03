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
              #"xco2_relative" : {"range": [-6, 6], "cmap" : "cm.RdBu_r", "binsize": 0.05}, 
              "xco2_relative" : {"range": [-8, 8], "cmap" : "cm.RdBu_r", "binsize": 0.065},
              #"xco2_relative" : {"range": [-10, 10], "cmap" : "cm.RdBu_r", "binsize": 0.08},
              "tcwv" : {"range": [0, 75], "cmap" : "cm.Blues", "binsize": 1/3.}, 
              "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}, 
              "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}, 
              "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [-1, 2], "cmap" : "cm.YlGn", "binsize": 0.015}
            }


def truncate(n, d):
    if np.isfinite(n):
        return math.floor(n * 10 ** d) / 10 ** d
    else:
        return n

for var in DATA_DICT.keys():
    #print(var)
    
    unpadded_cmap_name = os.path.join(CMAP_CSV_DIR, "unpadded", var + "_" + re.split("\.", DATA_DICT[var]["cmap"])[-1] + "_" +  str(DATA_DICT[var]["range"][0]) + "to" + str(DATA_DICT[var]["range"][1]) + ".csv")
    padded_cmap_name = os.path.join(CMAP_CSV_DIR, "padded", var + "_" + re.split("\.", DATA_DICT[var]["cmap"])[-1] + "_" +  str(DATA_DICT[var]["range"][0]) + "to" + str(DATA_DICT[var]["range"][1]) + ".csv")
        
    data_crange_low = [i*DATA_DICT[var]["binsize"] for i in np.arange(DATA_DICT[var]["range"][0] / DATA_DICT[var]["binsize"], (DATA_DICT[var]["range"][1] + DATA_DICT[var]["binsize"]) / DATA_DICT[var]["binsize"], 1)][:-1]
    data_crange_high = [i*DATA_DICT[var]["binsize"] for i in np.arange(DATA_DICT[var]["range"][0] / DATA_DICT[var]["binsize"], (DATA_DICT[var]["range"][1] + DATA_DICT[var]["binsize"]) / DATA_DICT[var]["binsize"], 1)][1:]
    
    ncolors = len(data_crange_low)
    
    cmap_df = pd.DataFrame([np.array([0, 0, 0, 0, np.nan, np.nan])], columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"])
    array_to_append = np.append(np.round(255*mpl.colors.to_rgba_array(eval(DATA_DICT[var]["cmap"])(0))), [np.NINF, data_crange_low[0]])
    cmap_df = cmap_df.append(pd.DataFrame([array_to_append], columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"]))
    
    for n in range(ncolors):
        #print(n)
        if round(data_crange_high[n], 3) > DATA_DICT[var]["range"][1]:
            array_to_append = np.append(np.round(255*mpl.colors.to_rgba_array(eval(DATA_DICT[var]["cmap"])(n+1))), [round(data_crange_low[n], 3), round(data_crange_high[n])])
            cmap_df = cmap_df.append(pd.DataFrame([array_to_append], columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"]))
            break
        array_to_append = np.append(np.round(255*mpl.colors.to_rgba_array(eval(DATA_DICT[var]["cmap"])(n+1))), [round(data_crange_low[n], 3), round(data_crange_high[n], 3)])
        cmap_df = cmap_df.append(pd.DataFrame([array_to_append], columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"]))
            
    array_to_append = np.append(np.round(255*mpl.colors.to_rgba_array(eval(DATA_DICT[var]["cmap"])(ncolors+2))), [round(data_crange_high[n]), np.inf])
    cmap_df = cmap_df.append(pd.DataFrame([array_to_append], columns = ["red", "green", "blue", "alpha", "data_lim_low", "data_lim_high"]))

    if var == "tcwv":
        cmap_df["data_lim_low"] = cmap_df["data_lim_low"].map(lambda x: truncate(x, 2))
        cmap_df["data_lim_high"] = cmap_df["data_lim_high"].map(lambda x: truncate(x, 2))
    
    #deal with known duplicates (due to converting colormap entries to bytescale integers)
    if var == "xco2":
        # [38,130,142,255,402.4,402.6] -> [38,131,142,255,402.4,402.6]
        cmap_df.loc[(cmap_df["data_lim_low"] == 402.4) & (cmap_df["data_lim_high"] == 402.6), ["green"]] = 131
        # [33,145,140,255,405.4,405.6] -> [33,145,141,255,405.4,405.6]
        cmap_df.loc[(cmap_df["data_lim_low"] == 405.4) & (cmap_df["data_lim_high"] == 405.6), ["blue"]] = 141
        # [32,146,140,255,405.6,405.8] -> [32,146,141,255,405.6,405.8]
        cmap_df.loc[(cmap_df["data_lim_low"] == 405.6) & (cmap_df["data_lim_high"] == 405.8), ["blue"]] = 141
        
    if var == "tcwv":
        # [238,245,252,255,3.66,4.0] -> [237,245,252,255,3.66,4.0]
        cmap_df.loc[(cmap_df["data_lim_low"] == 3.66) & (cmap_df["data_lim_high"] == 4.0), ["red"]] = 237

    if cmap_df.duplicated(subset=["red", "green", "blue", "alpha"]).any():
        print("Duplicate RGBAs still exist for " + var)
        print(cmap_df[cmap_df.duplicated(subset=["red", "green", "blue", "alpha"])])
        sys.exit()
    
    #unpadded colormap CSVs go to GIBS
    cmap_df = cmap_df.astype({"red": int, "green": int, "blue": int, "alpha": int})
    cmap_df.to_csv(unpadded_cmap_name, index=False)
    
    #pad to 256 colors for imagery generation colormaps
    cmap_df = cmap_df.append(cmap_df.iloc[[-1]*(256-ncolors-3)])
    cmap_df.to_csv(padded_cmap_name, index=False)
    
    
