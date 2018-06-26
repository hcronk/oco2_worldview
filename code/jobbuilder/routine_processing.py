import os
import sys
code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
sys.path.append(code_dir)
from oco2_worldview_imagery import oco2_worldview_imagery
from glob import glob
import operator
import re

#Global Variables
overwrite = False
lite_file_dirs = {"LtCO2": "/data6/OCO2/product/Lite/B8/LtCO2", 
                  "LtSIF": "/cloudsat/LtSIF"}
out_plot_dir = "/home/hcronk/worldview/plots/jobbuilder_testing"
lockfile_dir = "/home/hcronk/worldview/processing_status"
try_threshold = 3 #(number of times to try to process before moving to issues for analysis)
try_wait = 3600 #(number of seconds to wait before trying to reprocess a failed job)

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
            
tile_dict = { "NE": {"extent_box" : [0, 180, 0, 90]
                    },
              "SE": {"extent_box" : [0, 180, -90, 0]
                    },
              "SW": {"extent_box" : [-180, 0, -90, 0]
                    },
              "NW": {"extent_box" : [-180, 0, 0, 90]
                    }
            }


def find_unprocessed_file(lite_product, verbose=False):
    for f in glob(os.path.join(lite_file_dirs[lite_product], "*")):
        if verbose:
            print(f)
        lite_file_basename = os.path.basename(f)
        file_tokens = re.split("_", lite_file_basename)
        yymmdd = file_tokens[2]
        version = file_tokens[3]
        
        plot_tags = yymmdd + "_" + version + ".png"
        for v in data_dict[lite_product].keys():
            if verbose:
                print(v)
            for t in tile_dict.keys():
                if verbose:
                    print(t)
                out_plot_name = get_image_filename(v, tile_dict[t]["extent_box"], plot_tags)
                if not glob(out_plot_name):
                    job_file = re.sub("png", "json", os.path.basename(out_plot_name))
                    processing_or_problem = check_processing_or_problem(job_file)
                    if not processing_or_problem:
                        build_config(f, lite_product, v, out_plot_name, job_file)
            

def build_config(oco2_file, lite_product, var, out_plot_name, job_file, overwrite=overwrite, rgb=False, debug=False, verbose=False):
    
    global lockfile
    
    config_dict = {}
    
    config_dict["lite_file"] = oco2_file
    config_dict["var"] = data_dict[lite_product][var]
    config_dict["rgb"] = rgb
    
    with open(job_file, "w") as config_file:
        json.dump(config_dict, config_file, indent=4, sort_keys=True)
    
    if debug:
        print(job_file + " created for subsequent debugging")
	os.remove(lockfile)
	sys.exit()
	
    run_job(job_file, verbose=verbose)


    
    
def get_image_filename(var, extent_box, plot_name_tags):
    """
    Build the filename of the output image
    """
    
    return os.path.join(out_plot_dir, var + "_Lat" + str(extent_box[2]) + "to" + str(extent_box[3]) + "_Lon" + str(extent_box[0]) + "to" + str(extent_box[1]) + "_" + plot_name_tags)

def check_processing_or_problem(job_file):
    """
    Check if the job is already processing or if there is a problem with it.
    problem == the job has been run and faild more times than the try threshold defined globally
    """
    
    global lockfile
    global issue_file
    
    #Check for / create lockfile
    basename = re.sub("json", "proc", os.path.basename(job_file))
    lockfile = os.path.join(lockfile_dir, "processing", basename)
    issue_file = os.path.join(lockfile_dir, "problem", basename)

    if glob(issue_file):
	return True

    if glob(lockfile):
	with open(lockfile, "r") as lf:
	    tries = lf.readlines()
	
	latest_try = tries[-1].rstrip("\n")
	latest_try_dt = datetime.datetime.strptime(re.split("\.", latest_try)[0], "%Y-%m-%d %H:%M:%S")
	today = datetime.datetime.now()
	delta_t = (today - latest_try_dt).total_seconds()
	if delta_t <= try_wait:
	    #Give it some time before trying to process again
	    return True
	
	if len(tries) > try_threshold:
	    shutil.copy(lockfile, issue_file)
	    if glob(issue_file) and os.stat(lockfile).st_size == os.stat(issue_file).st_size:
		os.remove(lockfile)
	    return True
	else:
	    with open(lockfile, "a") as lf:
		lf.write(str(datetime.datetime.now()) + "\n")
	    return False
    else:
	#Create lockfile
	with open(lockfile, "w") as lf:
	    lf.write(str(datetime.datetime.now()) + "\n")
    
        return False

    
if __name__ == "__main__":

    for p in data_dict.keys():
	find_unprocessed_file(p)
    
