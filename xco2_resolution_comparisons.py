from oco2_modis_vistool.oco2_modis_vistool import do_modis_overlay_plot
from oco2_modis_vistool.OCO2FileOps import LiteCO2File

import numpy as np
import os
import sys
from shapely.geometry import Polygon, Point, MultiPoint
import matplotlib.pyplot as plt


resolution = "1km"
date = "2015-03-24"
geo_upper_left = [42.8,-83.4]
geo_lower_right = [41.2,-82.5]

out_plot_dir = "/home/hcronk/worldview/plots"
plot_name = "Detroit_vistool_" + resolution + ".png"
out_data_dir = "/home/hcronk/worldview/data"
out_data_name = "Detroit_vistool_" + resolution + ".h5"
out_plot = os.path.join(out_plot_dir, plot_name)
out_data = os.path.join(out_data_dir, out_data_name)

oco2_file_dir = "/home/codell/OCO2_results/b70/lite_B7Br"
oco2_file = "oco2_LtCO2_150324_B7305Br_160712115928s.nc4"

#These numbers came from the GIBS ICD
gibs_resolution_dict = {"2km" : 0.017578125, "1km" : 0.0087890625, "500m" : 0.00439453125, "250m" : 0.002197265625}

#lat_bins = np.arange(90, -90, -gibs_resolution_dict[resolution], dtype=float)
#lon_bins = np.arange(-180, 180, gibs_resolution_dict[resolution], dtype=float)
#
#lat_centers = np.arange(lat_bins[0] - gibs_resolution_dict[resolution] / 2., lat_bins[-1] - gibs_resolution_dict[resolution], -gibs_resolution_dict[resolution], dtype=float)
#lon_centers = np.arange(lon_bins[0] + gibs_resolution_dict[resolution] / 2., lon_bins[-1] + gibs_resolution_dict[resolution], gibs_resolution_dict[resolution], dtype=float)

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

print os.path.join(oco2_file_dir, oco2_file)

#Read in data
lf = LiteCO2File(os.path.join(oco2_file_dir, oco2_file))
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
latitude = latitude[quality_mask]
longitude = longitude[quality_mask]
vertex_latitude = np.squeeze(vertex_latitude[quality_mask, :])
vertex_longitude = np.squeeze(vertex_longitude[quality_mask, :])
data = data[quality_mask]

latlon_subset_mask = np.logical_and(
            np.logical_and(latitude <= 42.52, latitude >= 42.48), 
            np.logical_and(longitude <= -83.0, longitude >= -83.3))

latitude = latitude[latlon_subset_mask]
longitude = longitude[latlon_subset_mask]
vertex_latitude = np.squeeze(vertex_latitude[latlon_subset_mask, :])
vertex_longitude = np.squeeze(vertex_longitude[latlon_subset_mask, :])
data = data[latlon_subset_mask]

#do_modis_overlay_plot(geo_upper_left,
#                      geo_lower_right, 
#		      date, vertex_latitude, vertex_longitude, data, var_lims = [395, 408],
#                      out_plot = "data_check.png", out_data = "data_check.h5")                     
#sys.exit()

#Create lat/lon corner pairs from vertices
#Each element of this array is a 4x2 array of lat/lon points of the parallelogram corners (Order: LL, UL, UR, LR)
poly = np.dstack([vertex_latitude, vertex_longitude])

vlat_mins = vertex_latitude.min(axis=1)
vlat_maxes = vertex_latitude.max(axis=1)
vlon_mins = vertex_longitude.min(axis=1)
vlon_maxes = vertex_longitude.max(axis=1)

data_check = []
lat_check = []
lon_check = []
trigger = False

fig = plt.figure()
ax = fig.add_subplot(111)

for n, vertices in enumerate(poly):
    print n
    #Create a shapely polygon from vertices (Point order: LL, UL, UR, LR, LL)
    pg = [Polygon((x,y) for x, y in vertices)][0]
    
    #Get the indexes of the center grid boxes where the lat/lon of the center is between the vertex min/max for this polygon
    lat_idx = np.where(np.logical_and(lat_centers >= vlat_mins[n], lat_centers <= vlat_maxes[n]))[0]
    lon_idx = np.where(np.logical_and(lon_centers >= vlon_mins[n], lon_centers <= vlon_maxes[n]))[0]

#    lat_idx = np.where(np.logical_and(lat_centers >= vlat_mins[n] - 0.01, lat_centers <= vlat_maxes[n] + 0.01))[0]
#    lon_idx = np.where(np.logical_and(lon_centers >= vlon_mins[n] - 0.01, lon_centers <= vlon_maxes[n] + 0.01))[0]
    
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
        if pt.within(pg):
            print pt, "is within", pg
            x = np.where(ll[1] == lon_centers)[0][0]
            y = np.where(ll[0] == lat_centers)[0][0]
            if grid[x,y] is None:
                grid[x,y] = [data[n]]
            else:
                grid[x,y].append(data[n])

    if n in [3,4,10]: 
        plt.plot(np.append(vertices[:,1],vertices[0,1]), np.append(vertices[:,0], vertices[0,0]), "-o", c="red")
        for xy in zip(vertices[:,1], vertices[:,0]):
            ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize=8.5)
        plt.scatter(lon_m.flatten(), lat_m.flatten(), c="blue", edgecolor='none')
        for xy in zip(np.round(lon_m.flatten(), 4), np.round(lat_m.flatten(), 4)):
            if np.all(np.round(xy, 4) == [-83.2368, 42.4995]):
                ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', rotation=-30, fontsize=8.5)
plt.show()
    
#    #Plot polygon vertices and gridpoints to visualize/quality check
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    #plt.scatter(vertices[:,1], vertices[:,0], c="red", edgecolor='none')
#    plt.plot(np.append(vertices[:,1],vertices[0,1]), np.append(vertices[:,0], vertices[0,0]), "-o", c="red")
#    for xy in zip(vertices[:,1], vertices[:,0]):
#        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize=8.5)
#    plt.scatter(lon_m.flatten(), lat_m.flatten(), c="blue", edgecolor='none')
#    for xy in zip(np.round(lon_m.flatten(), 4), np.round(lat_m.flatten(), 4)):
#        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', rotation=-30, fontsize=8.5)
#    plt.show()
#    print
#    #sys.exit()

x_action, y_action = np.nonzero(grid)

for x, y in zip(x_action, y_action):
    if grid[x,y] is not None:
        grid[x,y] = np.mean(grid[x,y])
        

# Plot the subset
lat_ul = geo_upper_left[0]
lon_ul = geo_upper_left[1]
lat_lr = geo_lower_right[0]
lon_lr = geo_lower_right[1]

x_subset_indices = np.where(np.logical_and(lon_centers >= lon_ul, lon_centers <= lon_lr))
y_subset_indices = np.where(np.logical_and(lat_centers >= lat_lr, lat_centers <= lat_ul))

subset_lat_centers = lat_centers[y_subset_indices]
subset_lon_centers = lon_centers[x_subset_indices]
subset_lat_bins = lat_bins[y_subset_indices]
subset_lon_bins = lon_bins[x_subset_indices]
grid_subset = grid[int(x_subset_indices[0][0]) : int(x_subset_indices[0][-1] + 1), int(y_subset_indices[0][0]) : int(y_subset_indices[0][-1] + 1)]

xg, yg = np.nonzero(grid_subset)
valid_grid = grid_subset[xg,yg].astype(float)

subset_lat_vertex = np.vstack([lat_bins[y], lat_bins[y], lat_bins[y + 1], lat_bins[y + 1]] for y in y_subset_indices[0][yg])
subset_lon_vertex = np.vstack([lon_bins[x], lon_bins[x + 1], lon_bins[x + 1], lon_bins[x]] for x in x_subset_indices[0][xg])

#do_modis_overlay_plot(geo_upper_left,
#                      geo_lower_right, 
#		      date, subset_lat_vertex, subset_lon_vertex, valid_grid, var_lims = [395, 408],
#                      out_plot = out_plot, out_data = out_data)

do_modis_overlay_plot(geo_upper_left,
                      geo_lower_right, 
		      date, subset_lat_vertex, subset_lon_vertex, valid_grid, var_lims = [395, 408],
                      out_plot = "data_check.png", out_data = "data_check.h5")

