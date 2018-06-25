import os
import sys
code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
sys.path.append(code_dir)
from oco2_worldview_imagery import oco2_worldview_imagery
from glob import glob

#Global Variables
overwrite = False
lite_file_dirs = {"LtCO2": "/data6/OCO2/product/Lite/B8/LtCO2", 
                  "LtSIF": "/cloudsat/LtSIF"}
output_dir = 

data_dict = { "LtCO2" : 
                { "xco2" : {"data_field_name" : "xco2", "preprocessing" : False, "range": [380, 430], "cmap" : "jet", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "xco2_relative" : {"data_field_name" : None, "preprocessing" : "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_trend_gl.txt", "range": [-6, 6], "cmap" : "RdYlBu_r", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "tcwv" : {"data_field_name" : "Retrieval/tcwv", "preprocessing" : False, "range": [0, 75], "cmap" : "viridis", "quality_info" : {}}, 
                },
             "LtSIF" : 
                { "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}, 
                  "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}, 
                  "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}
                }
            }

resolution = "500m"

#These numbers came from the GIBS ICD
gibs_resolution_dict = {"2km" : 0.017578125, "1km" : 0.0087890625, "500m" : 0.00439453125, "250m" : 0.002197265625}

#South to North by 1/2 km bins, starting at -90 and ending at 89.
lat_bins = np.arange(-90, 90, gibs_resolution_dict[resolution], dtype=float)
#West to East by 1km bins
lon_bins = np.arange(-180, 180, gibs_resolution_dict[resolution], dtype=float)

#South to North, starting 1/2km North of the southern most bin line and ending 1/2 km North of the northern most bin line
lat_centers = np.arange(lat_bins[0] + gibs_resolution_dict[resolution] / 2., lat_bins[-1] + gibs_resolution_dict[resolution], gibs_resolution_dict[resolution], dtype=float)
#West to East, starting 1/2km East of the western most bin line and ending 1/2 km east of the easternmost bin line
lon_centers = np.arange(lon_bins[0] + gibs_resolution_dict[resolution] / 2., lon_bins[-1] + gibs_resolution_dict[resolution], gibs_resolution_dict[resolution], dtype=float)


#Default is to process 4 quadrants per day per variable                      
north_subset_grid_indices = np.where(lat_bins >= 0)[0]
east_subset_grid_indices =  np.where(lon_bins >= 0)[0]
south_subset_grid_indices = np.where(lat_bins <= 0)[0]
west_subset_grid_indices = np.where(lon_bins <= 0)[0]

north_subset_data_indices = np.where(lat_bins >= 0)[0]
east_subset_data_indices =  np.where(lon_bins >= 0)[0]
south_subset_data_indices = np.where(lat_bins < 0)[0]
west_subset_data_indices = np.where(lon_bins < 0)[0]


quadrant_dict = { "NE": { "grid_lat_indices" : north_subset_grid_indices, "grid_lon_indices" : east_subset_grid_indices, 
                          "data_lat_indices" : north_subset_data_indices, "data_lon_indices" : east_subset_data_indices,
                          "extent_box" : [0, 180, 0, 90]
                        },
                  "SE": { "grid_lat_indices" : south_subset_grid_indices, "grid_lon_indices" : east_subset_grid_indices, 
                          "data_lat_indices" : south_subset_data_indices, "data_lon_indices" : east_subset_data_indices,
                          "extent_box" : [0, 180, -90, 0]
                        },
                  "SW": { "grid_lat_indices" : south_subset_grid_indices, "grid_lon_indices" : west_subset_grid_indices, 
                          "data_lat_indices" : south_subset_data_indices, "data_lon_indices" : west_subset_data_indices,
                          "extent_box" : [-180, 0, -90, 0]
                        },
                  "NW": { "grid_lat_indices" : north_subset_grid_indices, "grid_lon_indices" : west_subset_grid_indices, 
                          "data_lat_indices" : north_subset_data_indices, "data_lon_indices" : west_subset_data_indices,
                          "extent_box" : [-180, 0, 0, 90]
                        }  
                 }


def build_config(oco2_file, variable, overwrite=overwrite, rgb=False, debug=False, verbose=False):
    
    config_dict = {}
    
    config_dict["rgb"] = rgb
    
    
    
    
if __name__ == "__main__":

    files_to_
    
    if custom_geo_box:
        custom_geo_box = [float(x.strip("[,]")) for x in custom_geo_box]
        if len(custom_geo_box) != 4:
            print("Custom geolocation box format: [lat_min, lat_max, lon_min, lon_max]")
            print("Exiting.")
            sys.exit()
        lon_ul = custom_geo_box[2]
        lon_lr = custom_geo_box[3]
        lat_lr = custom_geo_box[0]
        lat_ul = custom_geo_box[1]

    if custom_geo_box:
        lon_indices = np.where(np.logical_and(lon_centers >= lon_ul, lon_centers <= lon_lr))[0]
        lat_indices = np.where(np.logical_and(lat_centers >= lat_lr, lat_centers <= lat_ul))[0]
    else:
        #Plot quadrants
        layered_plot_names = []
        if rgb:
            print("RGB overlay is for testing and case study purposes only.")
            print("Overlaying the full geolocation quadrants on an RGB requires too much memory and exceeds the PIL pixel limit.")
            print("The quadrants will be plotted to 45 degrees N/S/E/W only for testing the stitching interface.")
            
            north_subset_grid_indices = np.where(np.logical_and(lat_bins >= 0, lat_bins <= 45))[0]
            east_subset_grid_indices = np.where(np.logical_and(lon_bins >= 0, lon_bins <= 45))[0]
            south_subset_grid_indices = np.where(np.logical_and(lat_bins <= 0, lat_bins >= -45))[0]
            west_subset_grid_indices = np.where(np.logical_and(lon_bins <= 0, lon_bins >= -45))[0]

            north_subset_data_indices = np.where(np.logical_and(lat_bins >= 0, lat_bins < 45))[0]
            east_subset_data_indices = np.where(np.logical_and(lon_bins >= 0, lon_bins < 45))[0]
            south_subset_data_indices = np.where(np.logical_and(lat_bins < 0, lat_bins >= -45))[0]
            west_subset_data_indices = np.where(np.logical_and(lon_bins < 0, lon_bins >= -45))[0]


            quadrant_dict = { "NE": { "grid_lat_indices" : north_subset_grid_indices, "grid_lon_indices" : east_subset_grid_indices, 
                                      "data_lat_indices" : north_subset_data_indices, "data_lon_indices" : east_subset_data_indices,
                                      "extent_box" : [0, 45, 0, 45]
                                    },
                              "SE": { "grid_lat_indices" : south_subset_grid_indices, "grid_lon_indices" : east_subset_grid_indices, 
                                      "data_lat_indices" : south_subset_data_indices, "data_lon_indices" : east_subset_data_indices,
                                      "extent_box" : [0, 45, -45, 0]
                                    },
                              "SW": { "grid_lat_indices" : south_subset_grid_indices, "grid_lon_indices" : west_subset_grid_indices, 
                                      "data_lat_indices" : south_subset_data_indices, "data_lon_indices" : west_subset_data_indices,
                                      "extent_box" : [-45, 0, -45, 0]
                                    },
                              "NW": { "grid_lat_indices" : north_subset_grid_indices, "grid_lon_indices" : west_subset_grid_indices, 
                                      "data_lat_indices" : north_subset_data_indices, "data_lon_indices" : west_subset_data_indices,
                                      "extent_box" : [-45, 0, 0, 45]
                                    }  
                             }
        
        else:
            #The geo grid needs overlap for plotting 
            #to complete the square grid box on the edges of the subset


    
    for p in out_in_product_map.keys():
	find_unprocessed_file(p)
    
