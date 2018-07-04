import warnings
warnings.filterwarnings("ignore")

import numpy as np
from functools import reduce
import argparse
import os
import operator
import sys
import re
import h5py
import pandas as pd
from shapely.geometry import Polygon, Point, LineString
import matplotlib.pyplot as plt
#import json
import pickle
from ftplib import FTP
from glob import glob
import matplotlib.patches as mpatches
import xml.etree.ElementTree as ET
import matplotlib as mpl
from osgeo import gdal, osr
from PIL import Image
import cartopy
import cartopy.feature as cfeature
ccrs = cartopy.crs

#Global Variables
code_dir = os.path.dirname(os.path.realpath(__file__))

geo_dict = { "LtCO2" : {
                "lat" : "vertex_latitude",
                "lon" : "vertex_longitude"
                },
             "LtSIF" : {
                "lat" : "footprint_vertex_latitude",
                "lon" : "footprint_vertex_longitude"
                }
            }          

resolution = "500m"
dpi = 10000

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
    
    lite_file = contents["lite_file"]
    product = contents["product"]
    var = contents["var"]
    preprocessing = contents["preprocessing"]
    data_field_name = contents["data_field_name"]
    quality_info = contents["quality_info"]
    extent_box = contents["extent_box"]
    plot_name = contents["out_plot_name"]
    var_range = contents["range"]
    cmap = contents["cmap"]
    rgb = contents["rgb"]
    
    return lite_file, product, var, preprocessing, data_field_name, quality_info, extent_box, plot_name, var_range, cmap, rgb

def stitch_quadrants(quadrant_plot_names, result_plot):
    """
    Stitches four existing plots into one single plot
    """
    
    NE_plot = [n for n in quadrant_plot_names if re.search("NE", n)][0]
    SE_plot = [n for n in quadrant_plot_names if re.search("SE", n)][0]
    SW_plot = [n for n in quadrant_plot_names if re.search("SW", n)][0]
    NW_plot = [n for n in quadrant_plot_names if re.search("NW", n)][0]
    
    north_imgs = [Image.open(i) for i in [NW_plot, NE_plot]]
    min_shape = sorted([(np.sum(i.size), i.size) for i in north_imgs])[0][1]
    north = np.hstack((np.asarray(i.resize(min_shape)) for i in north_imgs))
    north_img = Image.fromarray(north)
    north_img.save("temp_north.png")
    
    south_imgs = [Image.open(i) for i in [SW_plot, SE_plot]]
    min_shape = sorted([(np.sum(i.size), i.size) for i in south_imgs])[0][1]
    south = np.hstack((np.asarray(i.resize(min_shape)) for i in south_imgs))
    south_img = Image.fromarray(south)
    south_img.save("temp_south.png")
    
    global_imgs = [Image.open(i) for i in ["temp_north.png", "temp_south.png"]]
    min_shape = sorted([(np.sum(i.size), i.size) for i in global_imgs])[0][1]
    global_stack = np.vstack((np.asarray(i.resize(min_shape)) for i in global_imgs))
    global_img = Image.fromarray(global_stack)
    global_img.save(result_plot)
    
    return True

def update_GIBS_xml(date, xml_file):
    """
    Puts the date of interest into the GIBS XML file.
    For research mode, not operations
    """
    
    tree = ET.parse(xml_file)
    root = tree.getroot()

    url = root[0][0].text

    old_date = re.split('/', url)[6]

    new_url = re.sub(old_date, date, url)

    root[0][0].text = new_url
    tree.write(xml_file)
    
    return True
    
def pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr, xml_file, tif_file, xsize=1200, ysize=1000):
    """
    Pulls the Aqua RGB imagery from WorldView using GIBS and puts it in specified tif file with associated metadata.
    For research mode, not operations
    """ 
    
    print("\nPulling RGB imagery from GIBS")
    gdal_path = os.popen("which gdal_translate").read().strip()
    cmd = gdal_path + " -of GTiff -outsize "+str(xsize)+" "+str(ysize)+" -projwin "+str(lon_ul)+" "+str(lat_ul)+" "+str(lon_lr)+" "+str(lat_lr)+" "+xml_file+" "+tif_file
    os.system(cmd)
    
    return True

def prep_RGB(rgb_name, extent, xpix, ypix, dpi):
    """
    Prepares the RGB geotiff for layering with the data and writes it to a png
    For research mode, not operations
    """
    
    fig = plt.figure(figsize=(xpix / dpi, ypix / dpi), dpi=dpi)

    img = plt.imread(code_dir+'/intermediate_RGB.tif')

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.outline_patch.set_visible(False)
    ax.imshow(img, origin='upper', transform=ccrs.PlateCarree(), extent=extent)
    
    fig.savefig(rgb_name, bbox_inches='tight', pad_inches=0, dpi=dpi)
    
    return True

def patch_plot(data, grid_lat, grid_lon, extent, data_limits, cmap, out_plot_name, xpix, ypix, dpi):
    """
    Plot data polygons on a lat/lon grid
    In operational path
    """
    #print "In the function"
    fig = plt.figure(figsize=(xpix / dpi, ypix / dpi), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.outline_patch.set_visible(False)

    xg, yg = np.nonzero(data)
    if not xg.size or not yg.size:
        return False
    
    valid_grid = data[xg,yg].astype(float)

    subset_lat_vertex = np.vstack([grid_lat[y], grid_lat[y], grid_lat[y + 1], grid_lat[y + 1]] for y in yg)
    subset_lon_vertex = np.vstack([grid_lon[x], grid_lon[x + 1], grid_lon[x + 1], grid_lon[x]] for x in xg)
    
    zip_it = np.dstack([subset_lon_vertex, subset_lat_vertex])

    patches = []

    for row in zip_it:
        polygon = mpatches.Polygon(row)
        patches.append(polygon)                 
    p = mpl.collections.PatchCollection(patches, cmap=cmap, edgecolor='none')
    p.set_array(valid_grid)
    p.set_clim(data_limits[0], data_limits[1])
    ax.add_collection(p)

    fig.savefig(out_plot_name, bbox_inches='tight', pad_inches=0, dpi=dpi)
    
    del xg
    del yg
    del valid_grid
    del subset_lat_vertex
    del subset_lon_vertex
    del zip_it
    del polygon
    del patches
    del p
    
    return True

def layer_rgb_and_data(rgb_name, data_plot_name, layered_plot_name):
    """
    Layers the data PNG on top of the RGB PNG
    For research mode, not operations
    """
    
    base = Image.open(rgb_name)
    top = Image.open(data_plot_name)
    pixel_dat = list(top.getdata())
    for i, p in enumerate(pixel_dat):
        if p[:3] == (255, 255, 255):
            pixel_dat[i] = (255, 255, 255, 0)
    top.putdata(pixel_dat)
    base_copy = base.copy()
    base_copy.paste(top, (0,0), top)
    base_copy.save(layered_plot_name)
    
    return True

def get_oco2_data(var, oco2_file):
    """
    Extract given variable data from the OCO-2 lite file (.h5 format)
    In operational path
    """
    f = h5py.File(oco2_file, "r")
    try: 
        data = f[var][:]
    except:
        print("Problem retrieving " + var + " from " + oco2_file)
        f.close()
        return
    f.close()
    return data

def ftp_pull(path, filename):
    """Pulls data via FTP
    In operational path for pulling Reference XCO2 dataset from NOAA ESRL
    """
    
    path_tokens = re.split("/", path)
    ftp_name = path_tokens[2]
    
    count = 0
    tries = 5
    
    while count < tries:
	try:
	    #print("Connection attempt number " + str(count+1) + " for " + product)
	    cnx = FTP(ftp_name)
	    break
	except:
	    print("Connection failed, trying again.")
	count += 1

    cnx.login()
    
    cwd = re.split(ftp_name, path)[-1]
    cnx.cwd(cwd)
    
    with open(filename, "wb") as f:
        cnx.retrbinary("RETR " + filename, f.write)
    
    if glob(filename) and os.path.getsize(filename) == cnx.size(filename):
        cnx.close()
        return True
    else:
        cnx.close()
        return False

def preprocess(var, oco2_file, external_data_file=None):
    """
    Do any necessary processing to get data for the given variable
    In operational path
    """
    
    if var == "xco2_relative":
        ftp_file = os.path.basename(external_data_file)
        ftp_path= os.path.dirname(external_data_file)
        success = ftp_pull(ftp_path, ftp_file)
        if success:
            df = pd.read_csv(ftp_file, comment="#", na_values=-99.99, header=None, \
                            names=['year', 'month', 'day', 'trend'], delim_whitespace=True)
            yyyy = int("20" + re.split("_", os.path.basename(oco2_file))[2][0:2])
            mm = int(re.split("_", os.path.basename(oco2_file))[2][2:4])
            dd = int(re.split("_", os.path.basename(oco2_file))[2][4:6])
            ref_xco2 = float(df["trend"][df.index[df["year"] == yyyy] & df.index[df["month"] == mm] & df.index[df["day"] == dd]])
            oco2_xco2 = get_oco2_data("xco2", oco2_file)
            data = oco2_xco2 - ref_xco2
            del oco2_xco2
            del df
        else:
            print("Unable to retrieve " + external_data_file)
            sys.exit()  
    elif var == "sif_blended":
        sif757 = get_oco2_data("SIF_757nm", oco2_file)
        sif771 = get_oco2_data("SIF_771nm", oco2_file)
        data = 0.5 * (sif757 + 1.5 * sif771)
        del sif757
        del sif771
    else:
        print("No preprocessing required for " + var)
        return
    
    return data

def extent_box_to_indices(extent_box):

    lat_ul = extent_box[3]
    lon_ul = extent_box[0]
    lat_lr = extent_box[2]
    lon_lr = extent_box[1]

    lon_data_indices = np.where(np.logical_and(lon_centers >= lon_ul, lon_centers <= lon_lr))[0]
    lat_data_indices = np.where(np.logical_and(lat_centers >= lat_lr, lat_centers <= lat_ul))[0]


    if lat_data_indices[0] == 0:
        low_lat_idx = 0
    else:
        low_lat_idx = min(lat_data_indices)

    if lat_data_indices[-1] == len(lat_centers) - 1:
        high_lat_idx = lat_data_indices[-1]
    else:
        high_lat_idx = max(lat_data_indices) + 1      

    lat_grid_indices = np.arange(low_lat_idx, high_lat_idx + 1, dtype=int)


    if lon_data_indices[0] == 0:
        low_lon_idx = 0
    else:
        low_lon_idx = min(lon_data_indices)

    if lon_data_indices[-1] == len(lon_centers) - 1:
        high_lon_idx = lon_data_indices[-1]
    else:
        high_lon_idx = max(lon_data_indices) + 1

    lon_grid_indices = np.arange(low_lon_idx, high_lon_idx + 1, dtype=int)
    
    return lat_data_indices, lon_data_indices, lat_grid_indices, lon_grid_indices        


def regrid_oco2(data, vertex_latitude, vertex_longitude, grid_lat_centers, grid_lon_centers, debug=False):
    """
    Put OCO-2 data on a regular grid
    In operational path
    """
    
    grid = np.empty([len(grid_lon_centers), len(grid_lat_centers)], dtype=np.object)
    
    #Create lat/lon corner pairs from vertices
    #Each element of this array is a 4x2 array of lat/lon points of the parallelogram corners (Order: LL, UL, UR, LR)
    poly = np.dstack([vertex_latitude, vertex_longitude])

    vlat_mins = vertex_latitude.min(axis=1)
    vlat_maxes = vertex_latitude.max(axis=1)
    vlon_mins = vertex_longitude.min(axis=1)
    vlon_maxes = vertex_longitude.max(axis=1)

    for n, vertices in enumerate(poly):
        #print n
        #print vertices
        #Create a shapely polygon from vertices (Point order: LL, UL, UR, LR, LL)
        pg = [Polygon((x, y) for x, y in vertices)][0]

        #Get the indexes of the center grid boxes where the lat/lon of the center is between the vertex min/max for this polygon
        lat_idx = np.where(np.logical_and(grid_lat_centers >= vlat_mins[n], grid_lat_centers <= vlat_maxes[n]))[0]
        lon_idx = np.where(np.logical_and(grid_lon_centers >= vlon_mins[n], grid_lon_centers <= vlon_maxes[n]))[0]

        #If there are no grid boxes inside this polygon, move on to the next polygon
        if len(lat_idx) == 0 or len(lon_idx) == 0:
            continue

        #Get the center lat/lon bounds of the grid boxes inside this polygon
        center_lat_subset = grid_lat_centers[lat_idx]
        center_lon_subset = grid_lon_centers[lon_idx]

        lat_m, lon_m = np.meshgrid(center_lat_subset, center_lon_subset)
        zip_it = zip(list(lat_m.flatten()), list(lon_m.flatten()))

        #For each grid center between the lat/lon bounds, create a shapely point and check if it's inside the shapely polygon

        for ll in zip_it:
            pt = Point(ll[0], ll[1])
            if pt.within(pg):
                x = np.where(ll[1] == lon_centers)[0][0]
                y = np.where(ll[0] == lat_centers)[0][0]
                if grid[x,y] is None:
                    grid[x,y] = [data[n]]
                    #grid[x,y] = [n]
                else:
                    grid[x,y].append(data[n])
                    #grid[x,y].append(n)
            else:
                if pg.exterior.distance(pt) <= 1e-3:
                    x = np.where(ll[1] == lon_centers)[0][0]
                    y = np.where(ll[0] == lat_centers)[0][0]
                    if grid[x,y] is None:
                        grid[x,y] = [data[n]]
                        #grid[x,y] = [n]
                    else:
                        grid[x,y].append(data[n])
                        #grid[x,y].append(n)
                        
#        #Plot polygon vertices and gridpoints to visualize/quality check
#        dot_colors = []
#        for ll in zip_it:
#            pt = Point(ll[0], ll[1])
#            if pt.within(pg):
#                dot_colors.append("cyan")
#            else:
#                if pg.exterior.distance(pt) <= 1e-3:
#                    dot_colors.append("cyan")
#                else:
#                    dot_colors.append("black")
#        fig = plt.figure(figsize=(10,8))
#        ax = fig.add_subplot(111)
#        #plt.scatter(vertices[:,1], vertices[:,0], c="red", edgecolor='none')
#        plt.plot(np.append(vertices[:,1],vertices[0,1]), np.append(vertices[:,0], vertices[0,0]), "-o", c="red")
#        for xy in zip(vertices[:,1], vertices[:,0]):
#            ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize=8.5)
#        plt.scatter(lon_m.flatten(), lat_m.flatten(), c=dot_colors, edgecolor='none')
#        for xy in zip(np.round(lon_m.flatten(), 4), np.round(lat_m.flatten(), 4)):
#            ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', rotation=-30, fontsize=8.5)
#        plt.show()

    del poly
    del lat_m
    del lon_m
    del zip_it
    
    return grid

def oco2_worldview_imagery(job_file, verbose=False, debug=False, stitch=False):
    
    if verbose:
        print("Processing " + job_file) 
    
    lite_file, product, var, preprocessing, data_field_name, quality_info, extent_box, plot_name, variable_plot_lims, cmap, rgb = read_job_file(job_file) 

    lite_file_basename = os.path.basename(lite_file)
    file_tokens = re.split("_", lite_file_basename)

    product = file_tokens[1]
    yymmdd = file_tokens[2]
    version = file_tokens[3]

    date = "20" + yymmdd[:2] + "-" + yymmdd[2:4] + "-" + yymmdd[4:]

    #regrid = True

    if preprocessing:
        if type(preprocessing) == str:
            data = preprocess(var, lite_file, preprocessing)
        else:
             data = preprocess(var, lite_file)
    else:
        data = get_oco2_data(data_field_name, lite_file)

    var_lat = get_oco2_data(geo_dict[product]["lat"], lite_file)
    var_lon = get_oco2_data(geo_dict[product]["lon"], lite_file)
    #Cut out the missing data and the data that crosses the date line
    vertex_miss_mask = np.where(np.logical_not(np.any(var_lat == -999999, axis=1), np.any(var_lon == -999999, axis=1)))
    vertex_zero_mask = np.where(np.logical_not(np.any(var_lat == 0.0, axis=1), np.any(var_lon == 0.0, axis=1)))
    vertex_crossDL_mask = np.where(np.logical_not(np.any(var_lon <= -179.9, axis=1), np.any(var_lon >= 179.9, axis=1)))

    total_gridding_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask))

    if quality_info:
        #"quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
        quality = get_oco2_data(quality_info["quality_field_name"], lite_file)
        quality_mask = np.where(quality_info["qc_operator"](quality, quality_info["qc_val"]))
        total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask, quality_mask))
        del quality
        del quality_mask
    else:
        total_mask = total_gridding_mask

    var_lat_gridding = np.squeeze(var_lat[total_mask, :])
    var_lon_gridding = np.squeeze(var_lon[total_mask, :])
    data = data[total_mask]

    #Get the Lite File indices in each GIBS grid box
    grid = regrid_oco2(data, var_lat_gridding, var_lon_gridding, lat_centers, lon_centers)
    
    del total_gridding_mask
    del total_mask
    del var_lat_gridding
    del var_lon_gridding
    del data
    del var_lat
    del var_lon
    del vertex_miss_mask
    del vertex_zero_mask
    del vertex_crossDL_mask

    x_action, y_action = np.nonzero(grid)

    for x, y in zip(x_action, y_action):
        if grid[x,y] is not None:
            grid[x,y] = np.mean(grid[x,y])

    del x_action
    del y_action
    
    lat_data_indices, lon_data_indices, lat_grid_indices, lon_grid_indices = extent_box_to_indices(extent_box)

    if verbose:
        print("Plotting")
    grid_subset = grid[int(lon_data_indices[0]) : int(lon_data_indices[-1] + 1), int(lat_data_indices[0]) : int(lat_data_indices[-1] + 1)]
    del grid
    lat_bin_subset = lat_bins[lat_grid_indices]
    lon_bin_subset = lon_bins[lon_grid_indices]
    success = patch_plot(grid_subset, lat_bin_subset, lon_bin_subset, extent_box, variable_plot_lims, cmap, plot_name, float(len(lon_data_indices)), float(len(lat_data_indices)), dpi)

    del grid_subset
    del lat_bin_subset
    del lon_bin_subset

    if not success:
        return

    if rgb:
        if verbose:
            print("RGB Dealings")
        lat_ul = extent_box[3]
        lon_ul = extent_box[0]
        lat_lr = extent_box[2]
        lon_lr = extent_box[1]
        out_plot_dir = os.path.dirname(plot_name)
        just_plot_name = os.path.basename(plot_name)
        rgb_name = os.path.join(out_plot_dir, re.sub(var, "RGB", just_plot_name))
        layered_rgb_name = os.path.join(out_plot_dir, re.sub(var, var +"_onRGB", just_plot_name))
        #rgb_name = os.path.join(out_plot_dir, "RGB_"+q+"_" + global_plot_name_tags)
        #layered_rgb_name = os.path.join(out_plot_dir, var + "_onRGB_"+q+"_" + global_plot_name_tags)

        success = update_GIBS_xml(date, rgb)
        success = pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr, rgb, code_dir+"/intermediate_RGB.tif")
        success = prep_RGB(rgb_name, extent_box, float(len(lon_data_indices)), float(len(lat_data_indices)), dpi)
        success = layer_rgb_and_data(rgb_name, plot_name, layered_rgb_name)

        #layered_plot_names.append(layered_rgb_name)
    #if rgb and stitch:            
    #    success = stitch_quadrants(layered_plot_names, os.path.join(out_plot_dir, var + "_onRGB_stitched_" + global_plot_name_tags))
    
    return True
                    
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description="OCO-2 Worldview imagery generation job", prefix_chars='-')
    parser.add_argument("-c", "--config_file", help="Path to the job configuration file", required=True)
    parser.add_argument("-v", "--verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("-d", "--debug", help="Plot polygon vertices and gridpoints to visualize/quality check", action="store_true")
    parser.add_argument("-s", "--stitch", help="Stitch a list of given files together", default=[])
    args = parser.parse_args()

    job_file = args.config_file
    verbose = args.verbose
    debug = args.debug
    stitch = args.stitch
    
    
    success = oco2_worldview_imagery(job_file, verbose=verbose)
