import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

DATA_DICT = {
             "xco2" : {
		       "range": [380, 430], 
		       "cmap_file" : "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps/xco2_viridis_380to430.csv", 
		       "axis": [0.05, 0.8, 0.9, 0.08], 
		       "label": "XCO2 (ppm)"				  
		       }, 
     "xco2_relative" : {
		       "range": [-6, 6], 
		       "cmap_file" : "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps/xco2_relative_RdBu_r_-6to6.csv", 
		       "axis": [0.05, 0.6, 0.9, 0.08], 
		       "label": "Relative XCO2 (ppm)"
		      }, 
             "tcwv" : {
		       "range": [0, 75], 
		       "cmap_file" : "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps/tcwv_Blues_0to75.csv", 
		       "axis": [0.05, 0.4, 0.9, 0.08], 
		       "label": "TCWV (mm)"
		      }, 
              "sif" : {
		       "range": [-1, 2], 
		       "cmap_file" : "/home/hcronk/oco2_worldview/code/utils/gibs_cmaps/sif757_YlGn_-1to2.csv", 
		       "axis": [0.05, 0.2, 0.9, 0.08], 
		       "label": "SIF (Wm^(-2)sr^(-1)um^(-1))"
		      }
            }


def make_cmap(gibs_csv_file):

    cmap_df = pd.read_csv(gibs_csv_file)
    
    cmap_list = list(zip(cmap_df.red/255, cmap_df.green/255, cmap_df.blue/255))
    bounds_list = list(cmap_df.data_lim_low)
    bounds_list.append(cmap_df.data_lim_high.iloc[-1])
    
    cmap = mpl.colors.LinearSegmentedColormap.from_list("gibs_cmap", cmap_list, len(cmap_list))
    norm = mpl.colors.BoundaryNorm(bounds_list, cmap.N)

    return cmap, norm


# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(8, 4))

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
for var in DATA_DICT.keys():
    ax = fig.add_axes(DATA_DICT[var]["axis"])
    cmap, norm = make_cmap(DATA_DICT[var]["cmap_file"])
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal',
				    ticks=[DATA_DICT[var]["range"][0],DATA_DICT[var]["range"][1]])
    cb.set_label(DATA_DICT[var]["label"])
    #cb.ax.set_xticklabels([str(DATA_DICT[var]["range"][0]), str(DATA_DICT[var]["range"][1])])


fig.savefig("worldview_cbars.png")
