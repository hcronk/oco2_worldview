#from oco2_modis_vistool.OCO2FileOps import LiteCO2File, LiteSIFFile

import numpy as np
from functools import reduce
import argparse
import os
import operator
import sys
import re
import h5py
from shapely.geometry import Polygon, Point, LineString
import matplotlib.pyplot as plt
import json
from glob import glob
import matplotlib.patches as mpatches
import matplotlib as mpl
from PIL import Image
import cartopy
import cartopy.feature as cfeature
ccrs = cartopy.crs

def stitch_east_west(west_plot, east_plot, global_plot):
    """
    Stitches two existing plots (west_plot, east_plot) into one single plot
    """
    if not glob(west_plot):
        print(west_plot + " DNE. Check path. Exiting")
        sys.exit()
    if not glob(east_plot):
        print(east_plot + " DNE. Check path. Exiting")
        sys.exit()
        
    imgs = map(Image.open, [west_plot, east_plot])
    w, h = zip(*(i.size for i in imgs))
    total_w = sum(w)
    max_h = max(h)
    
    stitched = Image.new('RGB', (total_w, max_h))
    
    offset = 0
    for i in imgs:
        stitched.paste(i, (offset, 0))
        offset += i.size[0]
    
    stitched.save(global_plot)
    
    return

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
        oco2_xco2 = get_oco2_data("xco2", oco2_file)
        ref_xco2 = 400.
        data = oco2_xco2 - ref_xco2
    elif var == "sif_blended":
        sif757 = get_oco2_data("SIF_757nm", oco2_file)
        sif771 = get_oco2_data("SIF_771nm", oco2_file)
        data = 0.5 * (sif757 + 1.5 * sif771)
    else:
        print("No preprocessing required for " + var)
        return
    
    return data

def regrid_oco2(vertex_latitude, vertex_longitude, grid_lat_centers, grid_lon_centers):
    
    grid = np.empty([len(grid_lon_centers), len(grid_lat_centers)], dtype=np.object)
    
    #Create lat/lon corner pairs from vertices
    #Each element of this array is a 4x2 array of lat/lon points of the parallelogram corners (Order: LL, UL, UR, LR)
    poly = np.dstack([vertex_latitude, vertex_longitude])

    vlat_mins = vertex_latitude.min(axis=1)
    vlat_maxes = vertex_latitude.max(axis=1)
    vlon_mins = vertex_longitude.min(axis=1)
    vlon_maxes = vertex_longitude.max(axis=1)

    for n, vertices in enumerate(poly):
        print n
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
                    grid[x,y] = [n]
                else:
                    grid[x,y].append(n)
            else:
                if pg.exterior.distance(pt) <= 1e-3:
                    x = np.where(ll[1] == lon_centers)[0][0]
                    y = np.where(ll[0] == lat_centers)[0][0]
                    if grid[x,y] is None:
                        grid[x,y] = [n]
                    else:
                        grid[x,y].append(n)
    return grid

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Flags", prefix_chars='-')
    parser.add_argument("-v", "--verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("-f", "--file", help="Full path to input file name for command line use", default=None)
    args = parser.parse_args()
    verbose = args.verbose
    lite_file = args.file
    
    # Local paths
    xco2_lite_file_dir = "/data6/OCO2/product/Lite/B8/LtCO2/"
    sif_lite_file_dir = "/cloudsat/LtSIF/"
    out_plot_dir = "/home/hcronk/worldview/plots"

    #Variables to be plotted, if not all of the ones available. Can be left as an empty list []
    user_defined_var_list = ["xco2"]
    #Overwrite existing plots, if applicable
    overwrite = True

    ### TESTING ###
    ##Test Case 1: Act America Pennsylvania
    #oco2_file = "oco2_LtCO2_160727_B8100r_171007102957s.nc4"
    ##Test Case 2: Somalia, GL_2014-12-30_02638_033
    oco2_file = "oco2_LtCO2_141230_B8100r_171012035906s.nc4"
    ##Test Case 3: Chicago SIF
    #oco2_file = "oco2_LtSIF_150820_B8100r_171011083216s.nc4"
    #This will be what is given in loop/command line
    lite_file= os.path.join(xco2_lite_file_dir, oco2_file)
    #lite_file= os.path.join(sif_lite_file_dir, oco2_file)
    ######

    lite_file_basename = os.path.basename(lite_file)
    file_tokens = re.split("_", lite_file_basename)
    
    product = file_tokens[1]
    yymmdd = file_tokens[2]
    version = file_tokens[3]

    date = "20" + yymmdd[:2] + "-" + yymmdd[2:4] + "-" + yymmdd[4:]

    global_plot_name_tags = yymmdd + "_" + version + ".png"

    data_dict = { "LtCO2" : 
                    { "xco2" : {"data_field_name" : "xco2", "preprocessing" : False, "range": [395, 408], "cmap" : "jet", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                      "xco2_relative" : {"data_field_name" : None, "preprocessing" : True, "range": [-6, 6], "cmap" : "seismic", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                      "tcwv" : {"data_field_name" : "Retrieval/tcwv", "range": [0, 75], "cmap" : "viridis", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                    },
                 "LtSIF" : 
                    { "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}, 
                      "sif771" : {"data_field_name" : "SIF_771nm", "range": [-6, 6], "cmap" : "jet", "quality_info" : {}}, 
                      "sif_blended" : {"data_field_name" : None, "range": [0, 2], "cmap" : "jet", "quality_info" : {}}
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

    if user_defined_var_list:
        var_list = user_defined_var_list
    else:
        var_list = data_dict[product].keys()

    if verbose:
        print("Variables to be plotted: " + str(var_list))
        if overwrite:
            print("Any existing plots for these variables in " + out_plot_dir + " will be overwritten")
            print("This is the default behavior. To change, change the value of 'overwrite' to False")
            

    resolution = "500m"
    dpi = 10000
    regrid = True

    #These numbers came from the GIBS ICD
    gibs_resolution_dict = {"2km" : 0.017578125, "1km" : 0.0087890625, "500m" : 0.00439453125, "250m" : 0.002197265625}

    #South to North by 1km bins, starting at -90 and ending at 89.
    lat_bins = np.arange(-90, 90, gibs_resolution_dict[resolution], dtype=float)
    #West to East by 1km bins
    lon_bins = np.arange(-180, 180, gibs_resolution_dict[resolution], dtype=float)

    #South to North, starting 1/2km North of the southern most bin line and ending 1/1 km North of the northern most bin line
    lat_centers = np.arange(lat_bins[0] + gibs_resolution_dict[resolution] / 2., lat_bins[-1] + gibs_resolution_dict[resolution], gibs_resolution_dict[resolution], dtype=float)
    #West to East, starting 1/2km East of the western most bin line and ending 1/2 km east of the easternmost bin line
    lon_centers = np.arange(lon_bins[0] + gibs_resolution_dict[resolution] / 2., lon_bins[-1] + gibs_resolution_dict[resolution], gibs_resolution_dict[resolution], dtype=float)

    for var in var_list:
        if verbose:
            print("Processing "+ var)
        if var not in data_dict[product].keys():
            print(var + " is not defined in the " + product + " data dictionary. Please add it or check spelling.")
            print("Exiting.")
            sys.exit()
        if var in data_dict[product] and regrid:
            var_lat = get_oco2_data(geo_dict[product]["lat"], lite_file)
            var_lon = get_oco2_data(geo_dict[product]["lon"], lite_file)
            #Cut out the missing data and the data that crosses the date line
            vertex_miss_mask = np.where(np.logical_not(np.any(var_lat == -999999, axis=1), np.any(var_lon == -999999, axis=1)))
            vertex_zero_mask = np.where(np.logical_not(np.any(var_lat == 0.0, axis=1), np.any(var_lon == 0.0, axis=1)))
            vertex_crossDL_mask = np.where(np.logical_not(np.any(var_lon <= -179.9, axis=1), np.any(var_lon >= 179.9, axis=1)))

            total_mask = reduce(np.intersect1d, (vertex_miss_mask, vertex_zero_mask, vertex_crossDL_mask))

            var_lat = np.squeeze(var_lat[total_mask, :])
            var_lon = np.squeeze(var_lon[total_mask, :])

            #Get the LtCO2 indices in each GIBS grid box
            grid = regrid_oco2(var_lat, var_lon, lat_centers, lon_centers)
            regrid = False
            print "success"
    sys.exit()
    good_quality_xco2_grid_scheme
    global_plot_name = var + "_" + global_plot_name_tags
    if data_dict[product][var]["preprocessing"]:
        data = preprocessing(var, lite_file)
    else:
        data = preprocessing(var, lite_file)

    variable_plot_lims = data_dict[product][var]["range"]


    sys.exit()
    #Read in data
    lf = LiteCO2File(oco2_file_path)
    lf.open_file()
    latitude = lf.get_lat()
    longitude = lf.get_lon()
    vertex_latitude = lf.get_vertex_lat()
    vertex_longitude = lf.get_vertex_lon()
    data = lf.get_xco2()
    qf = lf.get_qf()
    lf.close_file()

    #Cut out the bad quality data
    quality_mask = np.where(qf == 0)
    vertex_miss_mask = np.intersect1d(quality_mask, np.where(np.logical_not(np.any(vertex_latitude == -999999, axis=1), np.any(vertex_longitude == -999999, axis=1))))
    missing_pixels = list(np.intersect1d(quality_mask, np.where(np.logical_and(np.any(vertex_latitude == -999999, axis=1), np.any(vertex_longitude == -999999, axis=1)))[0]))
    vertex_zero_mask = np.intersect1d(quality_mask, np.where(np.logical_not(np.any(vertex_latitude == 0.0, axis=1), np.any(vertex_longitude == 0.0, axis=1))))
    zero_pixels = list(np.intersect1d(quality_mask, np.where(np.logical_and(np.any(vertex_latitude == 0.0, axis=1), np.any(vertex_longitude == 0.0, axis=1)))[0]))

    total_mask = np.intersect1d(vertex_miss_mask, vertex_zero_mask)

    if len(missing_pixels) != 0: 
        with open("good_quality_missing_vertex_geolocation.json", "r") as qc_file:
            try: 
                records = json.load(qc_file)
            except:
                records = {}
        if oco2_file_path not in records.keys():
            records[oco2_file_path] = missing_pixels
            with open("good_quality_missing_vertex_geolocation.json", "w") as qc_file:
                json.dump(records, qc_file, indent=4)       

    if len(zero_pixels) != 0: 
        with open("good_quality_zero_vertex_geolocation.json", "r") as qc_file:
            try: 
                records = json.load(qc_file)
            except:
                records = {}
        if oco2_file_path not in records.keys():
            records[oco2_file_path] = zero_pixels
            with open("good_quality_zero_vertex_geolocation.json", "w") as qc_file:
                json.dump(records, qc_file, indent=4)       

    latitude = latitude[total_mask]
    longitude = longitude[total_mask]
    vertex_latitude = np.squeeze(vertex_latitude[total_mask, :])
    vertex_longitude = np.squeeze(vertex_longitude[total_mask, :])
    data = data[total_mask]




    #    #Plot polygon vertices and gridpoints to visualize/quality check
    #    fig = plt.figure(figsize=(10,8))
    #    ax = fig.add_subplot(111)
    #    #plt.scatter(vertices[:,1], vertices[:,0], c="red", edgecolor='none')
    #    plt.plot(np.append(vertices[:,1],vertices[0,1]), np.append(vertices[:,0], vertices[0,0]), "-o", c="red")
    #    for xy in zip(vertices[:,1], vertices[:,0]):
    #        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize=8.5)
    #    plt.scatter(lon_m.flatten(), lat_m.flatten(), c="blue", edgecolor='none')
    #    for xy in zip(np.round(lon_m.flatten(), 4), np.round(lat_m.flatten(), 4)):
    #        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', rotation=-30, fontsize=8.5)
    #    plt.show()

    x_action, y_action = np.nonzero(grid)

    for x, y in zip(x_action, y_action):
        if grid[x,y] is not None:
            grid[x,y] = np.mean(grid[x,y])

    # Plot the Eastern Hemisphere
    #fig = plt.figure(figsize=(0.5 * grid_x_elem / dpi, 0.5 * grid_y_elem / dpi), dpi=dpi)
    ##ax = plt.subplot(111, projection=ccrs.PlateCarree())
    #ax = plt.axes(projection=ccrs.PlateCarree())
    ##ax.set_global()
    #ax.outline_patch.set_visible(False)
    ##ax = plt.subplot(111)
    ##ax.set_extent([
    #ax.coastlines(resolution='10m', color='black', linewidth=1)
    #ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    ##fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection': ccrs.PlateCarree()})

    xg, yg = np.nonzero(grid)
    valid_grid = grid[xg,yg].astype(float)

    #print valid_grid.shape

    #subset_lat_vertex = np.vstack([lat_bins[y], lat_bins[y], lat_bins[y + 1], lat_bins[y + 1]] for y in yg)
    #subset_lon_vertex = np.vstack([lon_bins[x], lon_bins[x + 1], lon_bins[x + 1], lon_bins[x]] for x in xg)

    east_subset_indices =  np.where(lon_centers >= 0)
    west_subset_indices = np.where(lon_centers < 0)

    #subset_lat_centers = lat_centers[y_subset_indices]
    #subset_lon_centers = lon_centers[x_subset_indices]
    east_lon_bins = lon_bins[east_subset_indices]
    west_lon_bins = lon_bins[west_subset_indices]

    east_grid_subset = grid[int(east_subset_indices[0][0]) : int(east_subset_indices[0][-1] + 1), :]
    west_grid_subset = grid[int(west_subset_indices[0][0]) : int(west_subset_indices[0][-1] + 1), :]

    ###East####

    fig = plt.figure(figsize=(0.5 * grid_x_elem / dpi, 0.5 * grid_y_elem / dpi), dpi=dpi)
    #ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_global()
    ax.set_extent([0, 180, -90, 90])
    ax.outline_patch.set_visible(False)
    #ax = plt.subplot(111)
    #ax.coastlines(resolution='10m', color='black', linewidth=1)
    #ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    #fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection': ccrs.PlateCarree()})

    xg, yg = np.nonzero(east_grid_subset)
    valid_grid = east_grid_subset[xg,yg].astype(float)

    subset_lat_vertex = np.vstack([lat_bins[y], lat_bins[y], lat_bins[y + 1], lat_bins[y + 1]] for y in yg)
    subset_lon_vertex = np.vstack([lon_bins[x], lon_bins[x + 1], lon_bins[x + 1], lon_bins[x]] for x in east_subset_indices[0][xg])


    zip_it = np.dstack([subset_lon_vertex, subset_lat_vertex])

    patches = []

    for row in zip_it:
        polygon = mpatches.Polygon(row)
        patches.append(polygon)                 
    p = mpl.collections.PatchCollection(patches, cmap='jet', edgecolor='none')
    p.set_array(valid_grid)
    p.set_clim(variable_plot_lims[0], variable_plot_lims[1])
    ax.add_collection(p)

    #plt.axis('off')
    #fig.axes.get_xaxis().set_visible(False)
    #fig.axes.get_yaxis().set_visible(False)
    #ax.set_axis_off()
    #plt.show()
    fig.savefig(global_plot_east, bbox_inches='tight', pad_inches=0, dpi=dpi)

    ###West###

    fig = plt.figure(figsize=(0.5 * grid_x_elem / dpi, 0.5 * grid_y_elem / dpi), dpi=dpi)
    #ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_global()
    ax.set_extent([-180, 0, -90, 90])
    ax.outline_patch.set_visible(False)
    #ax = plt.subplot(111)
    #ax.set_extent([
    #ax.coastlines(resolution='10m', color='black', linewidth=1)
    #ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    #fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection': ccrs.PlateCarree()})

    xg, yg = np.nonzero(west_grid_subset)
    valid_grid = west_grid_subset[xg,yg].astype(float)

    subset_lat_vertex = np.vstack([lat_bins[y], lat_bins[y], lat_bins[y + 1], lat_bins[y + 1]] for y in yg)
    subset_lon_vertex = np.vstack([lon_bins[x], lon_bins[x + 1], lon_bins[x + 1], lon_bins[x]] for x in west_subset_indices[0][xg])


    zip_it = np.dstack([subset_lon_vertex, subset_lat_vertex])

    patches = []

    for row in zip_it:
        polygon = mpatches.Polygon(row)
        patches.append(polygon)                 
    p = mpl.collections.PatchCollection(patches, cmap='jet', edgecolor='none')
    p.set_array(valid_grid)
    p.set_clim(variable_plot_lims[0], variable_plot_lims[1])
    ax.add_collection(p)

    #plt.axis('off')
    #fig.axes.get_xaxis().set_visible(False)
    #fig.axes.get_yaxis().set_visible(False)
    #ax.set_axis_off()
    #plt.show()
    fig.savefig(global_plot_west, bbox_inches='tight', pad_inches=0, dpi=dpi)

    if glob(global_plot_east) and glob(global_plot_west) and not glob(global_plot):
        stitch_east_west(global_plot_west, global_plot_east, global_plot)

main()
