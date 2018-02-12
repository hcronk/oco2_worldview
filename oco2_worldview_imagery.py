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
import json
import urllib
from glob import glob
import matplotlib.patches as mpatches
import xml.etree.ElementTree as ET
import matplotlib as mpl
from osgeo import gdal, osr
from PIL import Image
import cartopy
import cartopy.feature as cfeature
ccrs = cartopy.crs

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
    Puts the date of interest into the GIBS XML file
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
    Pulls the Aqua RGB imagery from WorldView using GIBS and puts it in specified tif file with associated metadata
    """ 
    print("\nPulling RGB imagery from GIBS")
    gdal_path = os.popen("which gdal_translate").read().strip()
    cmd = gdal_path + " -of GTiff -outsize "+str(xsize)+" "+str(ysize)+" -projwin "+str(lon_ul)+" "+str(lat_ul)+" "+str(lon_lr)+" "+str(lat_lr)+" "+xml_file+" "+tif_file
    os.system(cmd)
    
    return True

def prep_RGB(rgb_name, extent, xpix, ypix, dpi):
    ### Pull in and prep RGB tif file ###
    ### Plot the RGB ###
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
    """
    #print "In the function"
    fig = plt.figure(figsize=(xpix / dpi, ypix / dpi), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.outline_patch.set_visible(False)

    xg, yg = np.nonzero(data)
    
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
    
    return True

def layer_rgb_and_data(rgb_name, data_plot_name, layered_plot_name):
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

def preprocess(var, oco2_file, external_data_file=None):
    """
    Do any necessary processing to get data for the given variable
    """
    if var == "xco2_relative":
        filename = os.path.basename(external_data_file)       
        urllib.urlretrieve(external_data_file, filename)
        df = pd.read_csv(filename, comment="#", na_values=-99.99, header=None, \
                        names=['year', 'month', 'decimal date', 'average', 'interpolated', 'seasonal corr trend', 'num days'], \
                        delim_whitespace=True)
        yyyy = int("20" + re.split("_", os.path.basename(oco2_file))[2][0:2])
        mm = int(re.split("_", os.path.basename(oco2_file))[2][2:4])
        ref_xco2 = float(df["average"][df.index[df["year"] == yyyy] & df.index[df["month"] == mm]])
        oco2_xco2 = get_oco2_data("xco2", oco2_file)
        data = oco2_xco2 - ref_xco2
        del oco2_xco2
        del df
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

def regrid_oco2(data, vertex_latitude, vertex_longitude, grid_lat_centers, grid_lon_centers):

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
    
    return grid

if __name__ == "__main__":

    code_dir = os.path.dirname(os.path.realpath(__file__))
    
    parser = argparse.ArgumentParser(description="Flags", prefix_chars='@')
    parser.add_argument("@v", "@@verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("@f", "@@file", help="Full path to input file name for command line use", default=None)
    parser.add_argument("@c", "@@config", help="Full path to text file containing a list of files to process, one filename per line", default=None)
    parser.add_argument("@g", "@@geolocation", help="Custom geolocation box for testing and case study purposes. Fmt: [lat_min, lat_max, lon_min, lon_max]", nargs = '+', default=[])
    parser.add_argument("@r", "@@rgb", help="Overlay plots on MODIS RGB for testing and case study purposes", action="store_true")
    args = parser.parse_args()
    verbose = args.verbose
    lite_file = args.file
    config_file = args.config
    custom_geo_box = args.geolocation
    rgb = args.rgb
    
    if not lite_file and not config_file or lite_file and config_file:
        print("Please provide a single file to process with the @f flag OR a file containing a list of files to process with the @c flag.")
        print("Exiting.")
        sys.exit()
    
    if config_file:
        if glob(config_file):
            files = list(np.genfromtxt(config_file, dtype="str", delimiter="\n"))
        else:
            print(config_file + " DNE. Exiting.")
            sys.exit()
    
    if lite_file:
        if glob(lite_file):
            files = [lite_file]
        else:
            print(lite_file + " DNE. Exiting.")
            sys.exit()
    
    if rgb:
        xml_file = os.path.join(code_dir, "GIBS_Aqua_MODIS_truecolor.xml")
    
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
    
    #Variables to be plotted, if not all of the ones available. Can be left as an empty list []
    user_defined_var_list = []
    #Output directory path
    out_plot_dir = "/home/hcronk/worldview/plots/for_feedback_1"
    #Overwrite existing plots in output directory, if applicable
    overwrite = True

    data_dict = { "LtCO2" : 
                    { "xco2" : {"data_field_name" : "xco2", "preprocessing" : False, "range": [380, 430], "cmap" : "jet", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                      "xco2_relative" : {"data_field_name" : None, "preprocessing" : "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt", "range": [-6, 6], "cmap" : "RdYlBu", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                      "tcwv" : {"data_field_name" : "Retrieval/tcwv", "preprocessing" : False, "range": [0, 75], "cmap" : "viridis", "quality_info" : {}}, 
                    },
                 "LtSIF" : 
                    { "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}, 
                      "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}, 
                      "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}
                    }
                }
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
    
    #grid_x_elem = int(360 / gibs_resolution_dict[resolution])
    #grid_y_elem = int(180 / gibs_resolution_dict[resolution])

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
    
    
    for lite_file in files:
        if verbose:
            print("Processing " + lite_file)  
        
        lite_file_basename = os.path.basename(lite_file)
        file_tokens = re.split("_", lite_file_basename)

        product = file_tokens[1]
        yymmdd = file_tokens[2]
        version = file_tokens[3]

        date = "20" + yymmdd[:2] + "-" + yymmdd[2:4] + "-" + yymmdd[4:]

        global_plot_name_tags = yymmdd + "_" + version + ".png"
        
        if user_defined_var_list:
            var_list = user_defined_var_list
        else:
            var_list = data_dict[product].keys()  

        if not overwrite:
            #double check there's something to do
            if custom_geo_box:
                for var in var_list:
                    if glob(os.path.join(out_plot_dir, var + "_Lat" + str(custom_geo_box[2]) + "to" + str(custom_geo_box[3]) + "_Lon" + str(custom_geo_box[0]) + "to" + str(custom_geo_box[1]) + "_" + global_plot_name_tags)):
                        var_list.remove(var)
            else:
                loop_list = list(var_list)
                for var in loop_list:
                    #print var
                    if glob(os.path.join(out_plot_dir, var + "_NE_" + global_plot_name_tags)) and glob(os.path.join(out_plot_dir, var + "_SE_" + global_plot_name_tags)) and glob(os.path.join(out_plot_dir, var + "_SW_" + global_plot_name_tags)) and glob(os.path.join(out_plot_dir, var + "_NW_" + global_plot_name_tags)):
                        var_list.remove(var)
            if not var_list:
                print("All plots exist. To overwrite, change the value of 'overwrite' to True")
                print("Exiting.")
                sys.exit()
        if verbose:
            print("Variables to be plotted: " + str(var_list))
            if overwrite:
                print("Any existing plots for these variables in " + out_plot_dir + " will be overwritten")
                print("This is the default behavior. To change, change the value of 'overwrite' to False")

        regrid = True
        #Parallelize variable processing?#
        for var in var_list:
            if verbose:
                print("Processing "+ var)
            if var not in data_dict[product].keys():
                print(var + " is not defined in the " + product + " data dictionary. Please add it or check spelling.")
                print("Exiting.")
                sys.exit()

            if data_dict[product][var]["preprocessing"]:
                if type(data_dict[product][var]["preprocessing"]) == str:
                    data = preprocess(var, lite_file, data_dict[product][var]["preprocessing"])
                else:
                     data = preprocess(var, lite_file)
            else:
                data = get_oco2_data(data_dict[product][var]["data_field_name"], lite_file)

            if var in data_dict[product] and regrid:
                var_lat = get_oco2_data(geo_dict[product]["lat"], lite_file)
                var_lon = get_oco2_data(geo_dict[product]["lon"], lite_file)
                #Cut out the missing data and the data that crosses the date line
                vertex_miss_mask = np.where(np.logical_not(np.any(var_lat == -999999, axis=1), np.any(var_lon == -999999, axis=1)))
                vertex_zero_mask = np.where(np.logical_not(np.any(var_lat == 0.0, axis=1), np.any(var_lon == 0.0, axis=1)))
                vertex_crossDL_mask = np.where(np.logical_not(np.any(var_lon <= -179.9, axis=1), np.any(var_lon >= 179.9, axis=1)))

                total_gridding_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask))

                if data_dict[product][var]["quality_info"]:
                    #"quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                    quality = get_oco2_data(data_dict[product][var]["quality_info"]["quality_field_name"], lite_file)
                    quality_mask = np.where(data_dict[product][var]["quality_info"]["qc_operator"](quality, data_dict[product][var]["quality_info"]["qc_val"]))
                    total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask, quality_mask))
                    del quality
                    del quality_mask
                else:
                    total_mask = total_gridding_mask

    #            var_lat_gridding = np.squeeze(var_lat[total_gridding_mask, :])
    #            var_lon_gridding = np.squeeze(var_lon[total_gridding_mask, :])

                var_lat_gridding = np.squeeze(var_lat[total_mask, :])
                var_lon_gridding = np.squeeze(var_lon[total_mask, :])
                data = data[total_mask]

                #print lat_bins[0], lat_bins[-1]
                #print lon_bins[0], lon_bins[-1]
                #print np.min(var_lat_gridding), np.max(var_lat_gridding)
                #print np.min(var_lon_gridding), np.max(var_lon_gridding)
                #sys.exit()

                #Get the LtCO2 indices in each GIBS grid box
                grid = regrid_oco2(data, var_lat_gridding, var_lon_gridding, lat_centers, lon_centers)
                #regrid = False
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

    #        if data_dict[product][var]["preprocessing"]:
    #            data = preprocessing(var, lite_file)
    #        else:
    #            data = get_oco2_data(data_dict[product][var]["data_field_name"], lite_file)

    #        if data_dict[product][var]["quality_info"]:
    #            #"quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
    #            quality = get_oco2_data(data_dict[product][var]["quality_info"]["quality_field_name"], lite_file)
    #            quality_mask = np.where(data_dict[product][var]["quality_info"]["qc_operator"](quality, data_dict[product][var]["quality_info"]["qc_val"]))
    #            total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask, quality_mask))
    #        else:
    #            total_mask = total_gridding_mask

            #data_grid = np.empty_like(grid)
            x_action, y_action = np.nonzero(grid)

            for x, y in zip(x_action, y_action):
                if grid[x,y] is not None:
                    grid[x,y] = np.mean(grid[x,y])
            
            del x_action
            del y_action

    #        for x, y in zip(x_action, y_action):
    #            #print x,y
    #            #print grid[x,y]
    #            #print data[grid[x,y]]
    #            #print type(data[grid[x,y]])
    #            #print np.float(np.mean([data[g] for g in grid[x,y]]))
    #            if grid[x,y] is not None and np.any(grid[x,y] in total_mask):
    #                #print "Valid"
    #                #print data[grid[x,y]]
    #                #print [(g, g in total_mask) for g in grid[x,y]]
    #                data_grid[x,y] = np.mean([data[g] for g in grid[x,y] if g in total_mask])
    #                #print data_grid[x,y]

            variable_plot_lims = data_dict[product][var]["range"]
            cmap = data_dict[product][var]["cmap"]
            #global_plot_name = os.path.join(var + "_" + global_plot_name_tags)


            if custom_geo_box:
                if verbose:
                    print("Plotting")
                plot_name = os.path.join(out_plot_dir, var + "_Lat" + str(custom_geo_box[2]) + "to" + str(custom_geo_box[3]) + "_Lon" + str(custom_geo_box[0]) + "to" + str(custom_geo_box[1]) + "_" + global_plot_name_tags)
                grid_subset = grid[int(lon_indices[0]) : int(lon_indices[-1] + 1), int(lat_indices[0]) : int(lat_indices[-1] + 1)]
                del grid
                lat_bin_subset = lat_bins[int(lat_indices[0]) : ]
                lon_bin_subset = lon_bins[int(lon_indices[0]) : ] 
                success = patch_plot(grid_subset, lat_bin_subset, lon_bin_subset, [lon_ul, lon_lr, lat_lr, lat_ul], variable_plot_lims, cmap, plot_name, float(len(lon_indices)), float(len(lat_indices)), dpi)
                
                del grid_subset
                del lat_bin_subset
                del lon_bin_subset
                
                if rgb:
                    if verbose:
                        print("RGB Dealings")
                    rgb_name = os.path.join(out_plot_dir, "RGB_Lat" + str(custom_geo_box[2]) + "to" + str(custom_geo_box[3]) + "_Lon" + str(custom_geo_box[0]) + "to" + str(custom_geo_box[1]) + "_" + global_plot_name_tags)
                    layered_rgb_name = os.path.join(out_plot_dir, var + "_onRGB_Lat" + str(custom_geo_box[2]) + "to" + str(custom_geo_box[3]) + "_Lon" + str(custom_geo_box[0]) + "to" + str(custom_geo_box[1]) + "_" + global_plot_name_tags)

                    success = update_GIBS_xml(date, xml_file)
                    success = pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr, xml_file, code_dir+"/intermediate_RGB.tif")
                    success = prep_RGB(rgb_name, [lon_ul, lon_lr, lat_lr, lat_ul], float(len(lon_indices)), float(len(lat_indices)), dpi)
                    success = layer_rgb_and_data(rgb_name, plot_name, layered_rgb_name)

            else:
                #Plot quadrants
                #Parallelize?#
                for q in quadrant_dict.keys():

                    plot_name = os.path.join(out_plot_dir, var + "_"+q+"_" + global_plot_name_tags)

                    if verbose:
                        print("Working on the " + q + " plotting quadrant")
                        if not glob(plot_name):
                            print("Creating " + plot_name)
                        if glob(plot_name) and overwrite:
                            print("Overwriting " + plot_name)
                        if glob(plot_name) and not overwrite:
                            print(plot_name + " exists and overwrite is not set. Moving on.")

                    if not glob(plot_name) or glob(plot_name) and overwrite:
                        grid_subset = grid[int(quadrant_dict[q]["data_lon_indices"][0]) : int(quadrant_dict[q]["data_lon_indices"][-1] + 1), int(quadrant_dict[q]["data_lat_indices"][0]) : int(quadrant_dict[q]["data_lat_indices"][-1] + 1)]
                        lat_bin_subset = lat_bins[quadrant_dict[q]["grid_lat_indices"]]
                        lon_bin_subset = lon_bins[quadrant_dict[q]["grid_lon_indices"]]
                        success = patch_plot(grid_subset, lat_bin_subset, lon_bin_subset, quadrant_dict[q]["extent_box"], variable_plot_lims, cmap, plot_name, float(len(quadrant_dict[q]["data_lon_indices"])), float(len(quadrant_dict[q]["data_lat_indices"])), dpi)
                        
                        del grid_subset
                        del lat_bin_subset
                        del lon_bin_subset
                        
                        if rgb:
                            lat_ul = quadrant_dict[q]["extent_box"][3]
                            lon_ul = quadrant_dict[q]["extent_box"][0]
                            lat_lr = quadrant_dict[q]["extent_box"][2]
                            lon_lr = quadrant_dict[q]["extent_box"][1]
                            rgb_name = os.path.join(out_plot_dir, "RGB_"+q+"_" + global_plot_name_tags)
                            layered_rgb_name = os.path.join(out_plot_dir, var + "_onRGB_"+q+"_" + global_plot_name_tags)

                            success = update_GIBS_xml(date, xml_file)
                            success = pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr, xml_file, code_dir+"/intermediate_RGB.tif")
                            success = prep_RGB(rgb_name, quadrant_dict[q]["extent_box"], float(len(quadrant_dict[q]["data_lon_indices"])), float(len(quadrant_dict[q]["data_lat_indices"])), dpi)
                            success = layer_rgb_and_data(rgb_name, plot_name, layered_rgb_name)

                            layered_plot_names.append(layered_rgb_name)
                del grid
                if rgb:            
                    success = stitch_quadrants(layered_plot_names, os.path.join(out_plot_dir, var + "_onRGB_stitched_" + global_plot_name_tags))
