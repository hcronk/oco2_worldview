import warnings
warnings.filterwarnings("ignore")

import numpy as np
from functools import reduce
import argparse
import os
import operator
import datetime
import sys
import re
import h5py
import pandas as pd
from shapely.geometry import Polygon, Point
import pickle
from collections import namedtuple
from ftplib import FTP
from glob import glob
import xml.etree.ElementTree as ET
from jinja2 import Template
import matplotlib as mpl
from PIL import Image, ImagePalette
#NOTE: matplotlib.pyplot import is handled import_pyplot_appropriately function to support interactively selecting the backend based on usage

#Global Variables
CODE_DIR = os.path.dirname(os.path.realpath(__file__))
REFERENCE_XCO2_FILE_DICT = { "oco2" : 
                            os.path.join(CODE_DIR, "oco2_ref_co2_trend_gl.txt"),
                            "oco3" : 
                            os.path.join(CODE_DIR, "oco3_ref_co2_trend_gl.txt")
                            }
                            
REFERENCE_XCO2_FILE_HEADER = os.path.join(CODE_DIR, "ref_co2_trend_gl_header.txt")
INTERMEDIATE_REFERENCE_XCO2_FILE_DICT = { "oco2" : 
                                         os.path.join(CODE_DIR, "oco2_intermediate_reference_xco2.csv"),
                                         "oco3" : 
                                         os.path.join(CODE_DIR, "oco3_intermediate_reference_xco2.csv")
                                        }
                      

GEO_DICT = { "LtCO2" : {
                "lat" : "vertex_latitude",
                "lon" : "vertex_longitude",
                "sid" : "sounding_id",
                },
             "LtSIF" : {
                "lat" : "Latitude_Corners",
                "lon" : "Longitude_Corners",
                "sid" : "Metadata/SoundingId",
                }
            }          

LITE_FILE_REGEX = "(?P<satellite>[oco2|oco3]{4})_(?P<product>[A-Za-z0-9]{5})_(?P<yy>[0-9]{2})(?P<mm>[0-9]{2})(?P<dd>[0-9]{2})_(?P<version>[0-9A-Za-z]{5,8})_[0-9]{12}s.nc4"
FTP_REGEX = "ftp://(?P<ftp_host>([^/]*))(?P<ftp_cwd>/.*)/(?P<ftp_file>([^/]*)$)"

IMAGE_REGEX = "(?P<satellite>[oco2|oco3]{4})_(?P<var>[a-z0-9](.*))_(?P<latspan>[Lato\.-](.*))_(?P<lonspan>[Lton\.-](.*))_(?P<yymmdd>[0-9]{6})_(?P<version>[0-9A-Za-z]{5,8}).png"
METADATA_REGEX = "(?P<satellite>[oco2|oco3]{4})_(?P<var>[a-z0-9](.*))_(?P<latspan>[Lato\.-](.*))_(?P<lonspan>[Lton\.-](.*))_(?P<yymmdd>[0-9]{6})_(?P<version>[0-9A-Za-z]{5,8}).met"
ODL_METADATA_TEMPLATE = os.path.join(CODE_DIR, "template.met")
METADATA_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"
METADATA_DATE_FORMAT = "%Y-%m-%d"
METADATA_TIME_TRUNCATED_FORMAT = "%H:%M:%SZ"
METADATA_TIME_FORMAT = "%H:%M:%S.%fZ"
METADATA_DAY_FORMAT = "%Y%j"

RESOLUTION = "500m"
DPI = 10000
OCO2_MISSING = -999999

#These numbers came from the GIBS ICD
GIBS_RESOLUTION_DICT = {"2km" : 0.017578125, "1km" : 0.0087890625, "500m" : 0.00439453125, "250m" : 0.002197265625}

#South to North by 1/2 km
LAT_GRID_SOUTH = np.arange(-90, 90 + GIBS_RESOLUTION_DICT[RESOLUTION], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)[:-1]
LAT_GRID_NORTH = np.arange(-90, 90 + GIBS_RESOLUTION_DICT[RESOLUTION], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)[1:]

#West to East by 1/2km 
LON_GRID_WEST = np.arange(-180, 180 + GIBS_RESOLUTION_DICT[RESOLUTION], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)[:-1]
LON_GRID_EAST = np.arange(-180, 180 + GIBS_RESOLUTION_DICT[RESOLUTION], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)[1:]

#South to North, starting 1/2km North of the southern most bin line and ending 1/2 km North of the northern most bin line
LAT_CENTERS = np.arange(LAT_GRID_SOUTH[0] + GIBS_RESOLUTION_DICT[RESOLUTION] / 2., LAT_GRID_NORTH[-1], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)
#West to East, starting 1/2km East of the western most bin line and ending 1/2 km east of the easternmost bin line
LON_CENTERS = np.arange(LON_GRID_WEST[0] + GIBS_RESOLUTION_DICT[RESOLUTION] / 2., LON_GRID_EAST[-1], GIBS_RESOLUTION_DICT[RESOLUTION], dtype=float)

#Operational Functionality
def import_pyplot_appropriately(debug=False):
    """
    In matplotlib 2.x, the matplotlib.pyplot import chokes in certain scenarios if the 
    backend is not selected appropriately, so this function sets it to TKAgg if 
    an interactive plotting session is required for debugging, provided that 
    the code is running in an environment that supports Xwindows,
    and sets it to non-interactive agg otherwise.
    """    
    global plt
    
    if debug:
        if os.environ.get("DISPLAY") == None:
            print("No display found. Cannot produce interactive plots for debugging")
            sys.exit()
        else:
            mpl.use("TKAgg")
    else:
        mpl.use("agg")
    
    import matplotlib.pyplot as plt
    
    return True

def read_job_file(job_file):
    """
    Read job description configuration file created by the routine_processing.py or 
    on_demand_processing.py job builders.
    """
    if not job_file:
        print("Please supply the path and name a job description configuration file using -c")
        sys.exit()

    if not glob(job_file):
        print("The job description congfiguration file " + job_file + " DNE")
        sys.exit()

    with open(job_file, "rb") as jf:
        contents = pickle.load(jf)
    
    tuple_setup = namedtuple("job_info", sorted(contents))
    job_info = tuple_setup(**contents)
    return job_info

def get_hdf5_data(var, hdf5_file):
    """
    Extract given variable data from an HDF-5 file
    In operational path
    """
    f = h5py.File(hdf5_file, "r")
    try: 
        data = f[var][:]
    except:
        print("Problem retrieving " + var + " from " + hdf5_file)
        f.close()
        return
    f.close()
    return data

def ftp_pull(ftp_path, verbose=False):
    """
    Pulls data via FTP
    In operational path for pulling Reference XCO2 dataset from NOAA ESRL
    """
    global FTP_SUBSTRING_DICT
    
    FTP_SUBSTRING_DICT = re.match(FTP_REGEX, ftp_path).groupdict()
    
    count = 0
    tries = 5
    
    while count < tries:
        try:
            if verbose:
                print("FTP Connection attempt number " + str(count+1))
            cnx = FTP(FTP_SUBSTRING_DICT["ftp_host"])
            break
        except:
            print("Connection failed, trying again.")
        count += 1

    cnx.login()
    
    cnx.cwd(FTP_SUBSTRING_DICT["ftp_cwd"])
    
    with open(os.path.join(CODE_DIR, FTP_SUBSTRING_DICT["ftp_file"]), "wb") as f:
        cnx.retrbinary("RETR " + FTP_SUBSTRING_DICT["ftp_file"], f.write)
    
    if (glob(os.path.join(CODE_DIR, FTP_SUBSTRING_DICT["ftp_file"])) and 
      os.path.getsize(os.path.join(CODE_DIR, FTP_SUBSTRING_DICT["ftp_file"]))) == cnx.size(FTP_SUBSTRING_DICT["ftp_file"]):
        cnx.close()
        return True
    else:
        cnx.close()
        return False

def preprocess(var, lite_file, external_data_file=None, verbose=False):
    """
    Do any necessary processing to get data for the given variable
    In operational path
    """
    
    global REFERENCE_XCO2_TO_SAVE
    
    if var == "xco2_relative":
        success = ftp_pull(external_data_file, verbose=verbose)
        if success:
            df = pd.read_csv(os.path.join(CODE_DIR, FTP_SUBSTRING_DICT["ftp_file"]), comment="#", na_values=-99.99, header=None, \
                            names=['year', 'month', 'day', 'cycle', 'trend'], delim_whitespace=True)
            
            REFERENCE_XCO2_TO_SAVE = df.loc[(df["year"] == int("20" + LITE_FILE_SUBSTRING_DICT["yy"])) &
                    (df["month"] == int(LITE_FILE_SUBSTRING_DICT["mm"])) &
                    (df["day"] == int(LITE_FILE_SUBSTRING_DICT["dd"]))]
                            
            ref_xco2 = float(REFERENCE_XCO2_TO_SAVE["cycle"])
            oco2_xco2 = get_hdf5_data("xco2", lite_file)
            data = oco2_xco2 - ref_xco2
            del oco2_xco2
            del df
        else:
            print("Unable to retrieve " + external_data_file)
            sys.exit()  
    elif var == "sif_blended":
        sif757 = get_hdf5_data("Daily_SIF_757nm", lite_file)
        sif771 = get_hdf5_data("Daily_SIF_771nm", lite_file)
        data = 0.5 * (sif757 + 1.5 * sif771)
        del sif757
        del sif771
    else:
        print("No preprocessing required for " + var)
        return
    return data

def extent_box_to_indices(extent_box):
    """
    Translates the latitude and longitude from the provided extent box 
    into indices within the latitude and longitude arrays
    """

    lat_ul = extent_box[3]
    lon_ul = extent_box[0]
    lat_lr = extent_box[2]
    lon_lr = extent_box[1]

    lon_data_indices = np.where(np.logical_and(LON_CENTERS >= lon_ul, LON_CENTERS <= lon_lr))[0]
    lat_data_indices = np.where(np.logical_and(LAT_CENTERS >= lat_lr, LAT_CENTERS <= lat_ul))[0]
    
    return lat_data_indices, lon_data_indices       

def make_cnorm(gibs_csv_file):
    """
    Creates a matplotlib color normalization from the custom OCO-2 colormaps produced for GIBS.
    Returns the normalization and a list of red, green, and blue bytes
    """

    cmap_df = pd.read_csv(gibs_csv_file)

    rgb_list = cmap_df.red.tolist() + cmap_df.green.tolist() + cmap_df.blue.tolist()
    
    bounds_list = list(cmap_df.data_lim_low)
    n_colors = len(bounds_list)
    bounds_list.append(cmap_df.data_lim_high.iloc[-1])
    bounds_list = bounds_list[1:]
    
    norm = mpl.colors.BoundaryNorm(bounds_list, n_colors)  

    return norm, rgb_list
    
def color_idx_plot(grid, norm, rgb_list, out_plot_name, verbose=False):
    """
    Plot data to a paletted PNG
    """
        
    grid_norm = norm(grid)
    
    im = Image.fromarray(np.flipud(grid_norm.T))
    im = im.convert("P")
    im.putpalette(ImagePalette.ImagePalette(palette=rgb_list))
    
    if verbose:
        print("Saving plot to " + out_plot_name)
        
    im.save(out_plot_name, format="PNG", transparency=0)
    
    return True

def blank_plot(x, y, rgb_list, out_plot_name, verbose=False):
    """
    Create a blank PNG
    """
    
    im = Image.new("P", (x, y))
    im.putpalette(ImagePalette.ImagePalette(palette=rgb_list))
    
    if verbose:
        print("Saving plot to " + out_plot_name)
        
    im.save(out_plot_name, format="PNG", transparency=0)
    
    return True

def write_image_odl_metadata(start_ts, end_ts, out_plot_name, extent_box, lite_file):
    """
    Writes ODL formatted metadata. One per image
    """
    
    try:
        lat_ul = extent_box[3]
        lon_ul = extent_box[0]
        lat_lr = extent_box[2]
        lon_lr = extent_box[1]

        metadata_filename = re.sub("png", "met", out_plot_name)
        metadata_filename_dict = re.match(METADATA_REGEX, os.path.basename(metadata_filename)).groupdict()

        with open(ODL_METADATA_TEMPLATE, "r") as tf:
            template = Template(tf.read())

        render = template.render(image_id=os.path.splitext(os.path.basename(metadata_filename))[0],
                                 production_datetime=datetime.datetime.utcnow().strftime(METADATA_TIMESTAMP_FORMAT),
                                 lite_file=os.path.basename(lite_file), 
                                 latspan_lonspan=metadata_filename_dict["latspan"] + "_" + metadata_filename_dict["lonspan"],
                                 yyyyddd=start_ts.date().strftime(METADATA_DAY_FORMAT),
                                 data_start_timestamp_truncated=start_ts.replace(microsecond=0).strftime(METADATA_TIME_TRUNCATED_FORMAT),
                                 data_start_timestamp=start_ts.strftime(METADATA_TIME_FORMAT),
                                 data_end_date=end_ts.date().strftime(METADATA_DATE_FORMAT),
                                 data_start_date=start_ts.date().strftime(METADATA_DATE_FORMAT),
                                 data_end_timestamp=end_ts.strftime(METADATA_TIME_FORMAT),
                                 latmin=extent_box[2],
                                 lonmin=extent_box[0],
                                 latmax=extent_box[3],
                                 lonmax=extent_box[1])

        with open(metadata_filename, "w") as mf:
            mf.write(render)

        return True

    except:
        return False


def write_image_worldfile(x_indices, y_indices, out_plot_name):
    """
    Writes ESRI world file. Need one per image produced.
    """
    
    worldfile_name = re.sub("png", "pgw", out_plot_name)
    
    with open(worldfile_name, "w") as wf:
        #Pixel X size in map units/pixel
        wf.write(str(GIBS_RESOLUTION_DICT[RESOLUTION]) + "\n")
        #Rotation about the Y-axis
        wf.write(str(float(0)) + "\n")
        #Rotation about the X-axis
        wf.write(str(float(0)) + "\n")
        #Negative pixel Y size in map units/pixel
        wf.write(str(-1 * GIBS_RESOLUTION_DICT[RESOLUTION]) + "\n")
        #X coordinate of upper left pixel center
        wf.write(str(LON_CENTERS[x_indices[0]]) + "\n")
        #Y coordinate of upper left pixel center
        wf.write(str(LAT_CENTERS[y_indices[-1]]))
    
    return True

def get_lite_oco2_timestamps(sounding_id):
    """
    Produces start and end timestamp for the lite file from OCO-2 sounding IDs
    """

    start = (sounding_id[0])
    count=1
    while True:
        if start <= 0:
            start = (sounding_id[count])
            count+=1
        else:
            break

    end = (sounding_id[-1])
    count=2
    while True:
        if end <= 0:
            end = (sounding_id[-count])
            count+=1
        else:
            break

    start_yyyy = int(str(start)[0:4])
    start_mm = int(str(start)[4:6])
    start_dd = int(str(start)[6:8])
    start_hh = int(str(start)[8:10])
    start_mn = int(str(start)[10:12])
    start_ss = int(str(start)[12:14])

    start_ts = datetime.datetime(start_yyyy, start_mm, start_dd, start_hh, start_mn, start_ss)

    end_yyyy = int(str(end)[0:4])
    end_mm = int(str(end)[4:6])
    end_dd = int(str(end)[6:8])
    end_hh = int(str(end)[8:10])
    end_mn = int(str(end)[10:12])
    end_ss = int(str(end)[12:14])

    end_ts = datetime.datetime(end_yyyy, end_mm, end_dd, end_hh, end_mn, end_ss)

    #print(start_ts, end_ts)
    return start_ts, end_ts

def get_lite_sif_date_time_coverage(sif_lite_file):
    
    """
    Produces start and end timestamp for the B10+ SIF lite file from
    the date_time_coverage global attribute
    """ 
    
    f = h5py.File(sif_lite_file, "r")
    dt_attr_list = [x.strip("'[,]") for x in f.attrs["date_time_coverage"]]
    
    dt_regex = "(?P<yyyy>[0-9]{4})-(?P<mm>[0-9]{2})-(?P<dd>[0-9]{2})T(?P<hh>[0-9]{2}):(?P<mn>[0-9]{2}):(?P<ss>[0-9]{2})\.[0-9]{6}Z"
    
    start_dt_dict = re.match(dt_regex, dt_attr_list[0]).groupdict()
    end_dt_dict = re.match(dt_regex, dt_attr_list[1]).groupdict()
    
    start_ts = datetime.datetime(int(start_dt_dict["yyyy"]),
                                int(start_dt_dict["mm"]), int(start_dt_dict["dd"]),
                                int(start_dt_dict["hh"]), int(start_dt_dict["mn"]), 
                                int(start_dt_dict["ss"]))
    end_ts = datetime.datetime(int(end_dt_dict["yyyy"]),
                              int(end_dt_dict["mm"]), int(end_dt_dict["dd"]),
                              int(end_dt_dict["hh"]), int(end_dt_dict["mn"]), 
                              int(end_dt_dict["ss"]))
    #print(start_ts, end_ts)
    return start_ts, end_ts

def regrid_oco2(data, vertex_latitude, vertex_longitude, lat_centers_subset, lon_centers_subset, verbose=False, debug=False):
    
    """
    Put OCO-2 data on a regular grid
    In operational path
    """
    
    if not "matplotlib.pyplot" in sys.modules:
        success = import_pyplot_appropriately(debug)
    
    grid = np.full([len(lon_centers_subset), len(lat_centers_subset)], np.nan)
    
    #Create lat/lon corner pairs from vertices
    #Each element of this array is a 4x2 array of lat/lon points of the parallelogram corners (Order: LL, UL, UR, LR)
    poly = np.dstack([vertex_latitude, vertex_longitude])

    vlat_mins = vertex_latitude.min(axis=1)
    vlat_maxes = vertex_latitude.max(axis=1)
    vlon_mins = vertex_longitude.min(axis=1)
    vlon_maxes = vertex_longitude.max(axis=1)

    for n, vertices in enumerate(poly):
        #Create a shapely polygon from vertices (Point order: LL, UL, UR, LR, LL)
        pg = [Polygon((x, y) for x, y in vertices)][0]

        #Get the indexes of the center grid boxes where the lat/lon of the center is between the vertex min/max for this polygon
        lat_idx = np.where(np.logical_and(lat_centers_subset >= vlat_mins[n], lat_centers_subset <= vlat_maxes[n]))[0]
        lon_idx = np.where(np.logical_and(lon_centers_subset >= vlon_mins[n], lon_centers_subset <= vlon_maxes[n]))[0]

        #If there are no grid boxes inside this polygon, move on to the next polygon
        if len(lat_idx) == 0 or len(lon_idx) == 0:
            continue

        center_lat_subset = lat_centers_subset[lat_idx]
        center_lon_subset = lon_centers_subset[lon_idx]

        lat_m, lon_m = np.meshgrid(center_lat_subset, center_lon_subset)
            
        zip_it = zip(list(lat_m.flatten()), list(lon_m.flatten()))

        #For each grid center between the lat/lon bounds, create a shapely point and check if it's inside the shapely polygon

        for ll in zip_it:
            pt = Point(ll[0], ll[1])
            if pt.within(pg):
                x = np.where(ll[1] == lon_centers_subset)[0][0]
                y = np.where(ll[0] == lat_centers_subset)[0][0]
                if np.isnan(grid[x,y]):
                    grid[x,y] = data[n]
            else:
                if pg.exterior.distance(pt) <= 1e-3:
                    x = np.where(ll[1] == lon_centers_subset)[0][0]
                    y = np.where(ll[0] == lat_centers_subset)[0][0]
                    if np.isnan(grid[x,y]):
                        grid[x,y] = data[n]
                        
        if debug:
            #Plot polygon vertices and gridpoints to visualize/quality check
            dot_colors = []
            #Need a new iterator
            zip_it = zip(list(lat_m.flatten()), list(lon_m.flatten()))
            for ll in zip_it:
                pt = Point(ll[0], ll[1])
                if pt.within(pg):
                    dot_colors.append("cyan")
                else:
                    if pg.exterior.distance(pt) <= 1e-3:
                        dot_colors.append("cyan")
                    else:
                        dot_colors.append("black")
            fig = plt.figure(figsize=(10,8))
            ax = fig.add_subplot(111)
            plt.plot(np.append(vertices[:,1],vertices[0,1]), np.append(vertices[:,0], vertices[0,0]), "-o", c="red")
            for xy in zip(vertices[:,1], vertices[:,0]):
                ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize=8.5)
            plt.scatter(lon_m.flatten(), lat_m.flatten(), c=dot_colors, edgecolor='none')
            for xy in zip(np.round(lon_m.flatten(), 4), np.round(lat_m.flatten(), 4)):
                ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', rotation=-30, fontsize=8.5)
            plt.show()

    del poly
    
    if not np.any(np.isfinite(grid)):
        if verbose:
            print("There are no data points that fall within the given geolocation box")
        return False
    
    del lat_m
    del lon_m
    del zip_it
    
    return grid 
   
### End Operational Functionality ###

def oco2_worldview_imagery(job_file, update_db=False, verbose=False, debug=False):
    """
    Main code for generating gridded OCO-2 imagery for Worldview
    """
    
    global LITE_FILE_SUBSTRING_DICT
    
    if verbose:
        print("Processing " + job_file) 
    
    job_info = read_job_file(job_file)
    
    LITE_FILE_SUBSTRING_DICT = re.match(LITE_FILE_REGEX, os.path.basename(job_info.lite_file)).groupdict()

    date = "20" + LITE_FILE_SUBSTRING_DICT["yy"] + "-" + LITE_FILE_SUBSTRING_DICT["mm"] + "-" + LITE_FILE_SUBSTRING_DICT["dd"]

    if job_info.preprocessing:
        if type(job_info.preprocessing) == str:
            data = preprocess(job_info.var, job_info.lite_file, job_info.preprocessing, verbose=verbose)
        else:
             data = preprocess(job_info.var, job_info.lite_file, verbose=verbose)
    else:
        data = get_hdf5_data(job_info.data_field_name, job_info.lite_file)

    var_lat = get_hdf5_data(GEO_DICT[job_info.product]["lat"], job_info.lite_file)
    var_lon = get_hdf5_data(GEO_DICT[job_info.product]["lon"], job_info.lite_file)
    
    #Cut out the missing data and the data that crosses the date line
    
    #OCO-2 stores the lat and lon of vertices of soundings in an two [along-track X 4] arrays 
    #and for regridding, the geolocation of all vertices must exist. 
    #Therefore, if either lat or lon of one or more vertices are missing, that row should be masked for operations.
    #There was also a bug in some versions of the OCO-2 data that erroneously set the geolocation to 0.0
    #in places where it does not make sense, so those rows are masked too
    vertex_miss_mask = np.where(np.logical_not(np.any(var_lat == -999999, axis=1), np.any(var_lon == -999999, axis=1)))
    vertex_zero_mask = np.where(np.logical_not(np.any(var_lat == 0.0, axis=1), np.any(var_lon == 0.0, axis=1)))

    #Worldview segregates data by day at the dateline, 
    #so if the geolocation of any sounding vertices crosses the date line, 
    #that sounding it should not be included for the given day
    vertex_crossDL_mask = np.where(np.logical_not(np.any(var_lon <= -179.9, axis=1), np.any(var_lon >= 179.9, axis=1)))

    if job_info.quality_info:
        #"quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
        quality = get_hdf5_data(job_info.quality_info["quality_field_name"], job_info.lite_file)
        quality_mask = np.where(job_info.quality_info["qc_operator"](quality, job_info.quality_info["qc_val"]))
        total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask, quality_mask))
        del quality
        del quality_mask
    else:
        total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask))

    var_lat_gridding = np.squeeze(var_lat[total_mask, :])
    var_lon_gridding = np.squeeze(var_lon[total_mask, :])
    data = data[total_mask]

    lat_data_indices, lon_data_indices = extent_box_to_indices(job_info.extent_box)
    
    #Get the Lite File data located in each GIBS grid box
    grid = regrid_oco2(data, var_lat_gridding, var_lon_gridding, LAT_CENTERS[lat_data_indices], LON_CENTERS[lon_data_indices], verbose=verbose, debug=debug)
    
    del total_mask
    del var_lat_gridding
    del var_lon_gridding
    del data
    del var_lat
    del var_lon
    del vertex_miss_mask
    del vertex_zero_mask
    del vertex_crossDL_mask
    
    if verbose:
        print("Plotting")
        
    norm, rgb_list = make_cnorm(job_info.cmap_file)
        
    if isinstance(grid, bool):
        #Grid is empty (False returned from regrid_oco2) because there are no data points that fall within the given geolocation box
        success = blank_plot(len(LON_CENTERS[lon_data_indices]), len(LAT_CENTERS[lat_data_indices]), rgb_list, job_info.out_plot_name, verbose=verbose)
    else:    
        success = color_idx_plot(grid, norm, rgb_list, job_info.out_plot_name, verbose=verbose)
    
    del grid

    if not success:
        if verbose:
            print("Problem plotting!")
        return
    else:
        if job_info.product == "LtCO2":
            lite_sounding_id = get_hdf5_data(GEO_DICT[job_info.product]["sid"], job_info.lite_file)
            lite_start, lite_end = get_lite_oco2_timestamps(lite_sounding_id)
        elif job_info.product == "LtSIF":
            lite_start, lite_end = get_lite_sif_date_time_coverage(job_info.lite_file)
        success = write_image_odl_metadata(lite_start, lite_end, job_info.out_plot_name, job_info.extent_box, job_info.lite_file)
        if not success:
            if verbose:
                print("Problem writing metadata!")
            return
        success = write_image_worldfile(lon_data_indices, lat_data_indices, job_info.out_plot_name)
        if not success:
            if verbose:
                print("Problem writing world file!")
            return

    if job_info.rgb:
        if verbose:
            print("RGB Dealings")
        lat_ul = job_info.extent_box[3]
        lon_ul = job_info.extent_box[0]
        lat_lr = job_info.extent_box[2]
        lon_lr = job_info.extent_box[1]
        out_plot_dir = os.path.dirname(job_info.out_plot_name)
        just_plot_name = os.path.basename(job_info.out_plot_name)
        rgb_name = os.path.join(out_plot_dir, re.sub(job_info.var, "RGB", just_plot_name))
        
    if update_db and job_info.var == "xco2_relative":
        #official relative XCO2 imagery was produced; update the reference XCO2 file
        image_filename_dict = re.match(IMAGE_REGEX, os.path.basename(job_info.out_plot_name)).groupdict()
        latspan_lonspan = image_filename_dict["latspan"] + "_" + image_filename_dict["lonspan"]
        df = pd.read_csv(REFERENCE_XCO2_FILE_DICT[LITE_FILE_SUBSTRING_DICT["satellite"]], comment="#", na_values=-99.99, header=None, \
                         names=['year', 'month', 'day', 'cycle', 'trend', 'tile'], delim_whitespace=True)     
        REFERENCE_XCO2_TO_SAVE["tile"] = latspan_lonspan
        relevant_tile_vals = df["tile"][df.index[df["year"] == int("20" + LITE_FILE_SUBSTRING_DICT["yy"])] & df.index[df["month"] == int(LITE_FILE_SUBSTRING_DICT["mm"])] & df.index[df["day"] == int(LITE_FILE_SUBSTRING_DICT["dd"])]]        
        if latspan_lonspan in relevant_tile_vals.values:
            #Update the value of cycle and trend to the one that is available now
            df.loc[df.index[df["year"] == int("20" + LITE_FILE_SUBSTRING_DICT["yy"])] & 
                    df.index[df["month"] == int(LITE_FILE_SUBSTRING_DICT["mm"])] & 
                    df.index[df["day"] == int(LITE_FILE_SUBSTRING_DICT["dd"])] & 
                    df.index[df["tile"] == latspan_lonspan], "cycle"] = REFERENCE_XCO2_TO_SAVE["cycle"].values[0]
            df.loc[df.index[df["year"] == int("20" + LITE_FILE_SUBSTRING_DICT["yy"])] & 
                    df.index[df["month"] == int(LITE_FILE_SUBSTRING_DICT["mm"])] & 
                    df.index[df["day"] == int(LITE_FILE_SUBSTRING_DICT["dd"])] & 
                    df.index[df["tile"] == latspan_lonspan], "trend"] = REFERENCE_XCO2_TO_SAVE["trend"].values[0]
        elif latspan_lonspan not in relevant_tile_vals.values:
            # Add the cycle and trend for the given date and tile
            df = pd.concat([df, REFERENCE_XCO2_TO_SAVE]).sort_values(by=["year", "month", "day", "tile"])
        else:
            raise ValueError("Unexpected value in relevant_tile_vals:,", relevant_tile_vals)
        df = df.astype({"year": int, "month": int, "day": int})
        df = df[["year", "month", "day", "cycle", "trend", "tile"]]
        df.to_csv(INTERMEDIATE_REFERENCE_XCO2_FILE_DICT[LITE_FILE_SUBSTRING_DICT["satellite"]], header=False, index=False, sep="	")
        with open(REFERENCE_XCO2_FILE_DICT[LITE_FILE_SUBSTRING_DICT["satellite"]], "w") as rf:
            for fname in [REFERENCE_XCO2_FILE_HEADER, INTERMEDIATE_REFERENCE_XCO2_FILE_DICT[LITE_FILE_SUBSTRING_DICT["satellite"]]]:
                with open(fname) as infile:
                    rf.write(infile.read())
    
    return True
                    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="OCO-2 Worldview imagery generation job", prefix_chars='@')
    parser.add_argument("@c", "@@config_file", help="Path to the job configuration file", required=True)
    parser.add_argument("@v", "@@verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("@d", "@@debug", help="Plot polygon vertices and gridpoints to visualize/quality check", action="store_true")
    args = parser.parse_args()

    job_file = args.config_file
    verbose = args.verbose
    debug = args.debug
    
    success = import_pyplot_appropriately(debug=debug)
    
    success = oco2_worldview_imagery(job_file, verbose=verbose, debug=debug)
