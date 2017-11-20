from oco2_modis_vistool.oco2_modis_vistool import do_modis_overlay_plot
from oco2_modis_vistool.OCO2FileOps import LiteCO2File

import numpy as np
import os
import sys
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

def stitch_east_west(east_plot, west_plot, global_plot):
    imgs = map(Image.open, [east_plot, west_plot])
    w, h = zip(*(i.size for i in imgs))
    total_w = sum(w)
    max_h = max(h)
    
    stitched = Image.new('RGB', (total_w, max_h))
    
    offset = 0
    for i in imgs:
        stitched.paste(i, (offset, 0))
        offset += i.size[0]
    
    stitched.save(global_plot)

resolution = "500m"
dpi = 10000
oco2_file_dir = "/home/codell/OCO2_results/b70/lite_B7Br"
#oco2_file_dir = "/home/codell/OCO2_results/b70/lite_test_20170724/"
out_data_dir = "/home/hcronk/worldview/data"
out_plot_dir = "/home/hcronk/worldview/plots"

##Detroit
#date = "2015-03-24"
#geo_upper_left = [42.8,-83.4]
#geo_lower_right = [41.2,-82.5]
#variable_plot_lims = [395, 408]
#
#plot_name = "Detroit_" + resolution + ".png"
#out_data_name = "Detroit_" + resolution + ".h5"
#
#global_plot_east_name = "Detroit_global_east_" + resolution + ".png"
#global_plot_west_name = "Detroit_global_west_" + resolution + ".png"
#global_plot_name = "Detroit_global_" + resolution + ".png"
#
#baseline_plot_name = "Detroit_baseline.png"
#baseline_data_name = "Detroit_baseline.h5"
#
#oco2_file = "oco2_LtCO2_150324_B7305Br_160712115928s.nc4"

#---#

#Act America Pennsylvania
date = "2016-07-27"
geo_upper_left = [43.5, -78.9]
geo_lower_right = [39.6, -77.5]
variable_plot_lims = [395, 408]

plot_name = "AAPenn_" + resolution + ".png"
out_data_name = "AAPenn_" + resolution + ".h5"

global_plot_name = "AAPenn_global_" + resolution + ".png"
global_plot_east_name = "AAPenn_global_east_" + resolution + ".png"
global_plot_west_name = "AAPenn_global_west_" + resolution + ".png"

baseline_plot_name = "AAPenn_baseline.png"
baseline_data_name = "AAPenn_baseline.h5"

oco2_file = "oco2_LtCO2_160727_B7305Br_160923172049s.nc4"

#---#

##Tommy Case Study: GL_2016-09-07_11629_023
#date = "2016-09-07"
#geo_upper_left = [-0.02819999999999823,-103.757698059082]
#geo_lower_right = [-2.028199999999998,-103.3143615722656]
#variable_plot_lims = [399., 403.]
#
#plot_name = "GL_2016-09-07_11629_023_" + resolution + ".png"
#out_data_name = "GL_2016-09-07_11629_023_" + resolution + ".h5"
#
#global_plot_name = "GL_2016-09-07_11629_global_" + resolution + ".png"
#
#baseline_plot_name = "GL_2016-09-07_11629_023_baseline.png"
#baseline_data_name = "GL_2016-09-07_11629_023_baseline.h5"
#
#oco2_file = "oco2_LtCO2_20160907_B7302rb_r02_COoffline.nc4"

#---#

##Tommy Case Study: GL_2016-09-07_11629_024
#date = "2016-09-07"
#geo_upper_left = [1.971800000000002,-104.2035446166992]
#geo_lower_right = [-0.02819999999999823,-103.7544708251953]
#variable_plot_lims = [399., 403.]
#
#plot_name = "GL_2016-09-07_11629_024_" + resolution + ".png"
#out_data_name = "GL_2016-09-07_11629_024_" + resolution + ".h5"
#
#global_plot_name = "GL_2016-09-07_11629_global_" + resolution + ".png"
#
#baseline_plot_name = "GL_2016-09-07_11629_024_baseline.png"
#baseline_data_name = "GL_2016-09-07_11629_024_baseline.h5"
#
#oco2_file = "oco2_LtCO2_20160907_B7302rb_r02_COoffline.nc4"

#---#

##Somalia: GL_2014-12-30_02638_033
#date = "2014-12-30"
#geo_upper_left = [0.6, 42.45]
#geo_lower_right = [-1.0, 42.9]
#variable_plot_lims = [395, 402]
#
#plot_name = "GL_2014-12-30_02638_033_" + resolution + ".png"
#out_data_name = "GL_2014-12-30_02638_033_" + resolution + ".h5"
#
#global_plot_name = "GL_2014-12-30_02638_global_" + resolution + ".png"
#
#baseline_plot_name = "GL_2014-12-30_02638_033_baseline.png"
#baseline_data_name = "GL_2014-12-30_02638_033_baseline.h5"
#
#oco2_file = "oco2_LtCO2_141230_B7305Br_160712143805s.nc4"

#----#

out_plot = os.path.join(out_plot_dir, plot_name)
out_data = os.path.join(out_data_dir, out_data_name)
global_plot_east = os.path.join(out_plot_dir, global_plot_east_name)
global_plot_west = os.path.join(out_plot_dir, global_plot_west_name)
global_plot = os.path.join(out_plot_dir, global_plot_name)

baseline_plot = os.path.join(out_plot_dir, baseline_plot_name)
baseline_data = os.path.join(out_data_dir, baseline_data_name)

if glob(global_plot_east) and glob(global_plot_west) and not glob(global_plot):
    print "Just the global plot"
    stitch_east_west(global_plot_east, global_plot_west, global_plot)

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

#print lat_bins[0], lat_centers[0]
#print lat_bins[-1], lat_centers[-1]
#print lon_bins[0], lon_centers[0]
#print lon_bins[-1], lon_centers[-1]

grid_x_elem = int(360 / gibs_resolution_dict[resolution])
grid_y_elem = int(180 / gibs_resolution_dict[resolution])

grid = np.empty([grid_x_elem, grid_y_elem], dtype=np.object)

oco2_file_path = os.path.join(oco2_file_dir, oco2_file)

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

if not glob(baseline_plot):
    do_modis_overlay_plot(geo_upper_left,
                          geo_lower_right, 
		          date, vertex_latitude, vertex_longitude, data, var_lims = variable_plot_lims,
                          out_plot = baseline_plot, out_data = baseline_data)

#Create lat/lon corner pairs from vertices
#Each element of this array is a 4x2 array of lat/lon points of the parallelogram corners (Order: LL, UL, UR, LR)
poly = np.dstack([vertex_latitude, vertex_longitude])

vlat_mins = vertex_latitude.min(axis=1)
vlat_maxes = vertex_latitude.max(axis=1)
vlon_mins = vertex_longitude.min(axis=1)
vlon_maxes = vertex_longitude.max(axis=1)

for n, vertices in enumerate(poly):
    #print n
    #Create a shapely polygon from vertices (Point order: LL, UL, UR, LR, LL)
    pg = [Polygon((x, y) for x, y in vertices)][0]
    
    #Get the indexes of the center grid boxes where the lat/lon of the center is between the vertex min/max for this polygon
    lat_idx = np.where(np.logical_and(lat_centers >= vlat_mins[n], lat_centers <= vlat_maxes[n]))[0]
    lon_idx = np.where(np.logical_and(lon_centers >= vlon_mins[n], lon_centers <= vlon_maxes[n]))[0]
    
    #If there are no grid boxes inside this polygon, move on to the next polygon
    if len(lat_idx) == 0 or len(lon_idx) == 0:
        continue
    
    #Get the center lat/lon bounds of the grid boxes inside this polygon
    center_lat_subset = lat_centers[lat_idx]
    center_lon_subset = lon_centers[lon_idx]
    
    lat_m, lon_m = np.meshgrid(center_lat_subset, center_lon_subset)
    zip_it = zip(list(lat_m.flatten()), list(lon_m.flatten()))
    
    #For each grid center between the lat/lon bounds, create a shapely point and check if it's inside the shapely polygon
    
    for ll in zip_it:
        pt = Point(ll[0], ll[1])
        #print pg
        #print pt
        if pt.within(pg):
            x = np.where(ll[1] == lon_centers)[0][0]
            y = np.where(ll[0] == lat_centers)[0][0]
            if grid[x,y] is None:
                grid[x,y] = [data[n]]
            else:
                grid[x,y].append(data[n])
        else:
            if pg.exterior.distance(pt) <= 1e-3:
                x = np.where(ll[1] == lon_centers)[0][0]
                y = np.where(ll[0] == lat_centers)[0][0]
                if grid[x,y] is None:
                    grid[x,y] = [data[n]]
                else:
                    grid[x,y].append(data[n])                  
    
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
    stitch_east_west(global_plot_east, global_plot_west, global_plot)

## Plot the subset
#lat_ul = geo_upper_left[0]
#lon_ul = geo_upper_left[1]
#lat_lr = geo_lower_right[0]
#lon_lr = geo_lower_right[1]
#
#x_subset_indices = np.where(np.logical_and(lon_centers >= lon_ul, lon_centers <= lon_lr))
#y_subset_indices = np.where(np.logical_and(lat_centers >= lat_lr, lat_centers <= lat_ul))
#
#subset_lat_centers = lat_centers[y_subset_indices]
#subset_lon_centers = lon_centers[x_subset_indices]
#subset_lat_bins = lat_bins[y_subset_indices]
#subset_lon_bins = lon_bins[x_subset_indices]
#grid_subset = grid[int(x_subset_indices[0][0]) : int(x_subset_indices[0][-1] + 1), int(y_subset_indices[0][0]) : int(y_subset_indices[0][-1] + 1)]
#
#xg, yg = np.nonzero(grid_subset)
#valid_grid = grid_subset[xg,yg].astype(float)
#
#subset_lat_vertex = np.vstack([lat_bins[y], lat_bins[y], lat_bins[y + 1], lat_bins[y + 1]] for y in y_subset_indices[0][yg])
#subset_lon_vertex = np.vstack([lon_bins[x], lon_bins[x + 1], lon_bins[x + 1], lon_bins[x]] for x in x_subset_indices[0][xg])
#
#do_modis_overlay_plot(geo_upper_left,
#                      geo_lower_right, 
#		      date, subset_lat_vertex, subset_lon_vertex, valid_grid, var_lims = variable_plot_lims,
#                      out_plot = out_plot, out_data = out_data)
