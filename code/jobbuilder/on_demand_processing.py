import os
import sys
import argparse
from glob import glob
import operator
import re
import numpy as np
from routine_processing import get_image_filename, check_processing_or_problem, build_config, get_intermediate_tif_filename
code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
sys.path.append(code_dir)
from oco2_worldview_imagery import stitch_quadrants

#Global Variables
DATA_DICT = { "LtCO2" : 
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
            
TILE_DICT = { "NE": {"extent_box" : [0, 180, 0, 90]
                    },
              "SE": {"extent_box" : [0, 180, -90, 0]
                    },
              "SW": {"extent_box" : [-180, 0, -90, 0]
                    },
              "NW": {"extent_box" : [-180, 0, 0, 90]
                    }
            }

LITE_FILE_REGEX = "oco2_(?P<product>[A-Za-z0-9]{5})_(?P<yymmdd>[0-9]{6})_(?P<version>B[0-9r]{,5})_[0-9]{12}s.nc4"

if __name__ == "__main__":
    
    code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
    
    parser = argparse.ArgumentParser(description="Flags", prefix_chars='@')
    parser.add_argument("@v", "@@verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("@f", "@@file", help="Full path to input file name for command line use", default=None)
    parser.add_argument("@l", "@@file_list", help="Full path to text file containing a list of files to process, one filename per line", default=None)
    parser.add_argument("@g", "@@geolocation", help="Lat/Lon extent box to plot. Fmt: [lat_min, lat_max, lon_min, lon_max]", nargs = '+', default=[])
    parser.add_argument("@r", "@@rgb", help="Overlay plots on MODIS RGB for testing and case study purposes", action="store_true")
    parser.add_argument("@o", "@@output_dir", help="Output directory for plots", default=os.path.join(code_dir, "plots"))
    parser.add_argument("@a", "@@vars", help="Variables to plot", nargs = '+', default=[])
    parser.add_argument("@w", "@@overwrite", help="Overwrite existing plots", action="store_true")
    parser.add_argument("@d", "@@debug", help="Just create job config file for debugging", action="store_true")
    parser.add_argument("@s", "@@stitch", help="Stitch together four quadrant plots", action="store_true")
    args = parser.parse_args()
    
    verbose = args.verbose
    lite_file = args.file
    file_list = args.file_list
    custom_geo_box = args.geolocation
    rgb = args.rgb
    out_plot_dir = args.output_dir
    user_defined_var_list = args.vars
    overwrite = args.overwrite
    debug = args.debug
    stitch = args.stitch
    
    if user_defined_var_list:
        user_defined_var_list = [x.strip("[,]") for x in user_defined_var_list]      
    
    if not os.path.exists(out_plot_dir):
        os.makedirs(out_plot_dir)
        
    if not lite_file and not file_list or lite_file and file_list:
        print("Please provide a single file to process with the @f flag OR a file containing a list of files to process with the @c flag.")
        print("Exiting.")
        sys.exit()
    
    if file_list:
        if glob(file_list):
            files = list(np.genfromtxt(file_list, dtype="str", delimiter="\n"))
        else:
            print(file_list + " DNE. Exiting.")
            sys.exit()
    
    if lite_file:
        if glob(lite_file):
            files = [lite_file]
        else:
            print(lite_file + " DNE. Exiting.")
            sys.exit()
    
    if custom_geo_box:
        custom_geo_box = [float(x.strip("[,]")) for x in custom_geo_box]
        if len(custom_geo_box) != 4 \
            or custom_geo_box[0] > custom_geo_box[1] \
            or custom_geo_box[2] > custom_geo_box[3] \
            or abs(custom_geo_box[0]) > 90. \
            or abs(custom_geo_box[1]) > 90. \
            or abs(custom_geo_box[2]) > 180. \
            or abs(custom_geo_box[3]) > 180.:
               print("Geolocation box format: [lat_min, lat_max, lon_min, lon_max] with Lat +/- 90 and Lon +/- 180")
               print("Exiting.")
               sys.exit()
        lon_ul = custom_geo_box[2]
        lon_lr = custom_geo_box[3]
        lat_lr = custom_geo_box[0]
        lat_ul = custom_geo_box[1]
        extent_box = [lon_ul, lon_lr, lat_lr, lat_ul]
    else:
        extent_box = []
        if stitch:
            print("Overlaying the full geolocation quadrants exceeds the PIL pixel limit.")
            print("The quadrants will be plotted to 45 degrees N/S/E/W only for testing the stitching interface.")
            #Update global variable
            TILE_DICT = { "NE": {"extent_box" : [0, 45, 0, 45]
                                },
                          "SE": {"extent_box" : [0, 45, -45, 0]
                                },
                          "SW": {"extent_box" : [-45, 0, -45, 0]
                                },
                          "NW": {"extent_box" : [-45, 0, 0, 45]
                                }
                        }

    if rgb:
        rgb = [os.path.join(code_dir, "GIBS_Aqua_MODIS_truecolor.xml")]
    else:
        rgb_list = False

    for lf in files:
        if verbose:
            print("Checking " + lf)  
        
        #lite_file_basename = os.path.basename(lf)
        #file_tokens = re.split("_", lite_file_basename)

        #product = file_tokens[1]
        #yymmdd = file_tokens[2]
        #version = file_tokens[3]
        
        lite_file_substring_dict = re.match(LITE_FILE_REGEX, os.path.basename(lf)).groupdict()

        date = "20" + lite_file_substring_dict["yymmdd"][:2] + "-" + lite_file_substring_dict["yymmdd"][2:4] + "-" + lite_file_substring_dict["yymmdd"][4:]

        plot_tags = lite_file_substring_dict["yymmdd"] + "_" + lite_file_substring_dict["version"] + ".png"
        
        if user_defined_var_list:
            var_list = user_defined_var_list
        else:
            var_list = DATA_DICT[lite_file_substring_dict["product"]].keys()  

        loop_list = list(var_list)
        for var in loop_list:
            if var not in DATA_DICT[lite_file_substring_dict["product"]].keys():
                if verbose:
                    print(var + " is not defined in the " + lite_file_substring_dict["product"] + " data dictionary. Please add it or check spelling.")
                var_list.remove(var)
        
        if not var_list:
            print("No variables to plot.")
            print("Exiting.")
            sys.exit()
        
        if extent_box:
            for var in var_list:
                if verbose:
                    print("Checking " + var)
                out_plot_name = get_image_filename(var, extent_box, plot_tags)
                layered_rgb_name = os.path.join(out_plot_dir, re.sub(var, var +"_onRGB", os.path.basename(out_plot_name)))
                if rgb:
                    intermediate_tif = [get_intermediate_tif_filename(extent_box, date)]
                    rgb_list = rgb + intermediate_tif
                if debug or overwrite or not glob(out_plot_name) or (rgb and not glob(layered_rgb_name)):
                    job_file = re.sub("png", "pkl", os.path.basename(out_plot_name))
                    processing_or_problem = check_processing_or_problem(job_file, verbose=verbose)
                    if not processing_or_problem:
                        if verbose:
                            print("Creating config file for " + var)
                        build_config(lf, lite_file_substring_dict["product"], var, extent_box, out_plot_name, job_file, rgb=rgb_list, debug=debug, verbose=verbose)
                else:
                    if verbose:
                        print(out_plot_name + " exists and will not be overwritten.")
                        print("To overwrite, change the value of 'overwrite' to True by setting the @w command line option")
                    continue
        else:
            for var in var_list:
                if verbose:
                    print("Checking " + var)
                if stitch:
                    plots_to_stitch = {}
                for t in TILE_DICT.keys():
                    if verbose:
                        print("Checking " + t + " tile")
                    out_plot_name = get_image_filename(var, TILE_DICT[t]["extent_box"], plot_tags)
                    layered_rgb_name = os.path.join(out_plot_dir, re.sub(var, var +"_onRGB", os.path.basename(out_plot_name)))
                    if stitch:
                        if rgb:
                            intermediate_tif = [get_intermediate_tif_filename(TILE_DICT[t]["extent_box"], date)]
                            rgb_list = rgb + intermediate_tif
                            plots_to_stitch[t] = layered_rgb_name
                        else:
                            plots_to_stitch[t] = out_plot_name
                    if debug or overwrite or not glob(out_plot_name) or (rgb and not glob(layered_rgb_name)):
                        job_file = re.sub("png", "pkl", os.path.basename(out_plot_name))
                        processing_or_problem = check_processing_or_problem(job_file, verbose=verbose)
                        if not processing_or_problem:
                            if verbose:
                                print("Creating config file for " + var)
                            build_config(lf, lite_file_substring_dict["product"], var, TILE_DICT[t]["extent_box"], out_plot_name, job_file, rgb=rgb_list, debug=debug, verbose=verbose)
                    else:
                        if verbose:
                            print(out_plot_name + " exists and will not be overwritten.")
                            print("To overwrite, change the value of 'overwrite' to True by setting the @w command line option")
                        continue
                if stitch:
                    tokens = re.split("_", plots_to_stitch["NE"])
                    absorb = [tokens.remove(token) for token in tokens if re.search("Lon", token)]
                    intermediate = ""
                    stitched_plot_name = "".join([intermediate + token + "_" if not re.search("Lat", token) else intermediate + "stitched_" for token in tokens])[:-1]
                    if debug or overwrite or not glob(stitched_plot_name):
                        success = stitch_quadrants(plots_to_stitch, os.path.join(out_plot_dir, stitched_plot_name))
            
