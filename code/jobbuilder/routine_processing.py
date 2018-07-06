import os
import sys
code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
sys.path.append(code_dir)
from oco2_worldview_imagery import oco2_worldview_imagery
from glob import glob
import operator
import re
import datetime
#import json
import pickle
from multiprocessing import Process

#Global Variables
lite_file_dirs = {"LtCO2": "/data6/OCO2/product/Lite/B8/LtCO2", 
                  "LtSIF": "/cloudsat/LtSIF"}
out_plot_dir = "/home/hcronk/worldview/plots/jobbuilder_testing"
lockfile_dir = "/home/hcronk/worldview/processing_status"
try_threshold = 3 #(number of times to try to process before moving to issues for analysis)
try_wait = 3600 #(number of seconds to wait before trying to reprocess a failed job)
overwrite = False

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
                out_plot_name = get_image_filename(out_plot_dir, v, tile_dict[t]["extent_box"], plot_tags)
                if not glob(out_plot_name) or overwrite:
                    #job_file = re.sub("png", "json", os.path.basename(out_plot_name))
                    job_file = re.sub("png", "pkl", os.path.basename(out_plot_name))
                    processing_or_problem = check_processing_or_problem(job_file)
                    if not processing_or_problem:
                        build_config(f, lite_product, v, tile_dict[t]["extent_box"], out_plot_name, job_file)
            

def build_config(oco2_file, lite_product, var, extent_box, out_plot_name, job_file, rgb=False, debug=False, verbose=False):
    
    global lockfile
    
    config_dict = data_dict[lite_product][var]
    
    config_dict["lite_file"] = oco2_file
    config_dict["product"] = lite_product
    config_dict["var"] = var
    config_dict["extent_box"] = extent_box
    config_dict["out_plot_name"] = out_plot_name
    config_dict["rgb"] = rgb
    
#    with open(job_file, "w") as config_file:
#        json.dump(config_dict, config_file, indent=4, sort_keys=True)
    
    with open(job_file, "wb") as config_file:
        pickle.dump(config_dict, config_file)
    
    if debug:
        print(job_file + " created for subsequent debugging")
	os.remove(lockfile)
	sys.exit()

    #sys.exit()
    #start process here, give arguments, after runs and returns, daemon=true to make sure children die if parent does _process = Process(number of cores), return true/false
    #not multithreading
    process = Process(target=run_job, args=(job_file, verbose))
    process.daemon = True
    process.start()
    process.join()
    if process.is_alive():
        process.terminate()
    #run_job(job_file, verbose=verbose)    
    
def get_image_filename(out_plot_dir, var, extent_box, plot_name_tags):
    """
    Build the filename of the output image
    """
    
    return os.path.join(out_plot_dir, var + "_Lat" + str(extent_box[2]) + "to" + str(extent_box[3]) + "_Lon" + str(extent_box[0]) + "to" + str(extent_box[1]) + "_" + plot_name_tags)

def check_processing_or_problem(job_file, verbose=False):
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
	if verbose:
            print("There is a problem with the job " + job_file)
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
            if verbose:
                print(job_file + " is already processing")
	    return True
	
	if len(tries) > try_threshold:
	    shutil.copy(lockfile, issue_file)
	    if glob(issue_file) and os.stat(lockfile).st_size == os.stat(issue_file).st_size:
		os.remove(lockfile)
	    if verbose:
                print("There is a problem with the job " + job_file)
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

def run_job(job_file, verbose=False):
    
    global lockfile
    global issue_file

    success = oco2_worldview_imagery(job_file, verbose=verbose)
    
    if success:
        with open(job_file, "rb") as jf:
            contents = pickle.load(jf)
    
        plot_name = contents["out_plot_name"]
        rgb = contents["rgb"]
        var = contents["var"]
        
        job_worked = check_job_worked(plot_name, var, rgb)
    
    if success and job_worked:
        os.remove(job_file)
	os.remove(lockfile)
	if glob(issue_file):
	    os.remove(issue_file)
        if rgb:
            out_plot_dir = os.path.dirname(plot_name)
            just_plot_name = os.path.basename(plot_name)
            rgb_name = os.path.join(out_plot_dir, re.sub(var, "RGB", just_plot_name))
            os.remove(rgb_name)        

def check_job_worked(plot_name, var, rgb=False):
 
    #Check the file exists
    if rgb:
        just_plot_name = os.path.basename(plot_name)
        layered_rgb_name = os.path.join(out_plot_dir, re.sub(var, var +"_onRGB", just_plot_name))
        if glob(plot_name) and glob(layered_rgb_name):
            return True
        else:
            return False
    else:
        if glob(plot_name):
            return True
        else:
            return False

    
if __name__ == "__main__":

    for p in data_dict.keys():
        find_unprocessed_file(p)
    
