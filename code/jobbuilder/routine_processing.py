import os
import sys
CODE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
sys.path.append(CODE_DIR)
from oco2_worldview_imagery import oco2_worldview_imagery
from glob import glob
import operator
import re
import datetime
import shutil
import pickle
import errno
from multiprocessing import Process
import sqlite3

#Global Variables
LITE_FILE_DIRS = {"LtCO2": "/data/oco2/scf/product/Lite/B9003r/r02", 
                  "LtSIF": "/data/oco2/scf/product/Lite/B8100r/r02"}
OUT_PLOT_DIR = "/home/hcronk/oco2_worldview/test_plots"
LOCKFILE_DIR = "/home/hcronk/oco2_worldview/processing_status"
TRY_THRESHOLD = 3 #(number of times to try to process before moving to issues for analysis)
TRY_WAIT = 3600 #(number of seconds to wait before trying to reprocess a failed job)
OVERWRITE = False
CONN = sqlite3.connect(os.path.join(CODE_DIR, "oco2_worldview_imagery.db"))
CUR = CONN.cursor()

DATA_DICT = { "LtCO2" : 
                { "xco2" : {"data_field_name" : "xco2", "preprocessing" : False, "range": [380, 430], "cmap" : "viridis", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "xco2_relative" : {"data_field_name" : None, "preprocessing" : "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_trend_gl.txt", "range": [-6, 6], "cmap" : "RdBu_r", "quality_info" : {"quality_field_name" : "xco2_quality_flag", "qc_val" :  0, "qc_operator" : operator.eq }}, 
                  "tcwv" : {"data_field_name" : "Retrieval/tcwv", "preprocessing" : False, "range": [0, 75], "cmap" : "Blues", "quality_info" : {}}, 
                },
             "LtSIF" : 
                { "sif757" : {"data_field_name" : "SIF_757nm", "preprocessing" : False, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}, 
                  "sif771" : {"data_field_name" : "SIF_771nm", "preprocessing" : False, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}, 
                  "sif_blended" : {"data_field_name" : None, "preprocessing" : True, "range": [0, 2], "cmap" : "YlGn", "quality_info" : {}}
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


def find_unprocessed_file(lite_product, verbose=False):
    """
    Crawl Lite File directories and check if all expected output imagery exists and
    initiate a job for any missing imagery
    """
    
    for root, subdirs, files in os.walk(LITE_FILE_DIRS[lite_product]):
        subdirs[:] = [d for d in subdirs if not d[0] == "."]
        files = [f for f in files if not f[0] == "."]
        for just_filename in files:
            f = os.path.join(root, just_filename)
            if verbose:
                print(f)
            sys.exit()
            lite_file_substring_dict = re.match(LITE_FILE_REGEX, just_filename).groupdict()
        
            plot_tags = lite_file_substring_dict["yymmdd"] + "_" + lite_file_substring_dict["version"] + ".png"

            for v in DATA_DICT[lite_product].keys():
                if verbose:
                    print(v)
                for t in TILE_DICT.keys():
                    if verbose:
                        print(t)
                    out_plot_name = get_image_filename(OUT_PLOT_DIR, v, TILE_DICT[t]["extent_box"], plot_tags)
                    db_entry = CUR.execute("SELECT filename FROM created_imagery WHERE filename=?", (os.path.basename(out_plot_name),)).fetchall()
                    if not db_entry or OVERWRITE:
                        #job_file = re.sub("png", "json", os.path.basename(out_plot_name))
                        job_file = re.sub("png", "pkl", os.path.basename(out_plot_name))
                        processing_or_problem = check_processing_or_problem(job_file)
                        if not processing_or_problem:
                            build_config(f, lite_product, v, TILE_DICT[t]["extent_box"], out_plot_name, job_file)
            

def build_config(oco2_file, lite_product, var, extent_box, out_plot_name, job_file, rgb=False, debug=False, verbose=False):
    """
    Create a pickle (.pkl) formatted job configuration file with 
    the necessary parameters to create OCO-2 Worldview imagery
    """  
    config_dict = DATA_DICT[lite_product][var]
    
    config_dict["lite_file"] = oco2_file
    config_dict["product"] = lite_product
    config_dict["var"] = var
    config_dict["extent_box"] = extent_box
    config_dict["out_plot_name"] = out_plot_name
    config_dict["rgb"] = rgb
    
    with open(job_file, "wb") as config_file:
        pickle.dump(config_dict, config_file)
    
    if debug:
        print(job_file + " created for subsequent debugging")
        os.remove(LOCKFILE)
        sys.exit()

    #Run as a multiprocessing Process even though just one process runs now to allow for multiprocessing
    #and also to ensure memory is released between jobs and to make sure jobs are killed if the parent process
    #is killed
    process = Process(target=run_job, args=(job_file, verbose))
    process.daemon = True
    process.start()
    process.join()
    if process.is_alive():
        process.terminate()   
    
def get_image_filename(image_dir, var, extent_box, plot_name_tags):
    """
    Build the filename of the output image
    """
    
    return os.path.join(image_dir, var + "_Lat" + str(extent_box[2]) + "to" + str(extent_box[3]) + "_Lon" + str(extent_box[0]) + "to" + str(extent_box[1]) + "_" + plot_name_tags)

def get_GIBS_xml_filename(date):
    """
    Build the filename of the GIBS XML configuration file
    """
    
    return os.path.join(CODE_DIR, "GIBS_Aqua_MODIS_truecolor_" + date + ".xml")

def get_intermediate_tif_filename(tif_dir, extent_box, date):
    """
    Build the filename of the intermediate TIFF imagery filename
    """
    
    return os.path.join(tif_dir, "intermediate_RGB_Lat" + str(extent_box[2]) + "to" + str(extent_box[3]) + "_Lon" + str(extent_box[0]) + "to" + str(extent_box[1]) + "_" + date + ".tif")

def check_processing_or_problem(job_file, verbose=False):
    """
    Check if the job is already processing or if there is a problem with it.
    problem == the job has been run and failed more times than the try threshold defined globally
    """
    
    global LOCKFILE
    global ISSUE_FILE
    
    #Check for / create lockfile
    basename = re.sub("json", "proc", os.path.basename(job_file))
    LOCKFILE = os.path.join(LOCKFILE_DIR, "processing", basename)
    ISSUE_FILE = os.path.join(LOCKFILE_DIR, "problem", basename)

    if glob(ISSUE_FILE):
        if verbose:
            print("There is a problem with the job " + job_file)
        return True

    if glob(LOCKFILE):
        with open(LOCKFILE, "r") as lf:
            tries = lf.readlines()
        
        latest_try = tries[-1].rstrip("\n")
        latest_try_dt = datetime.datetime.strptime(re.split("\.", latest_try)[0], "%Y-%m-%d %H:%M:%S")
        today = datetime.datetime.now()
        delta_t = (today - latest_try_dt).total_seconds()
        if delta_t <= TRY_WAIT:
            #Give it some time before trying to process again
            if verbose:
                print(job_file + " is already processing")
            return True
        
        if len(tries) > TRY_THRESHOLD:
            shutil.copy(LOCKFILE, ISSUE_FILE)
            if glob(ISSUE_FILE) and os.stat(LOCKFILE).st_size == os.stat(ISSUE_FILE).st_size:
                os.remove(LOCKFILE)
            if verbose:
                print("There is a problem with the job " + job_file)
            return True
        else:
            with open(LOCKFILE, "a") as lf:
                lf.write(str(datetime.datetime.now()) + "\n")
            return False
    else:
        #Create lockfile
        with open(LOCKFILE, "w") as lf:
            lf.write(str(datetime.datetime.now()) + "\n")
    
        return False

def run_job(job_file, verbose=False):
    """
    Create imagery using specifications in job file
    and if the imagery is produced, get rid of the lockfile
    and all intermediate files.
    """
    
    success = oco2_worldview_imagery(job_file, verbose=verbose)
    
    with open(job_file, "rb") as jf:
        contents = pickle.load(jf)
    
    if success:        
        job_worked = check_job_worked(contents["out_plot_name"], contents["var"], contents["rgb"])
    
    if success and job_worked:
        if verbose:
            print("Success!")
        db_entry = CUR.execute("SELECT * FROM created_imagery WHERE filename=?", (os.path.basename(contents["out_plot_name"]),)).fetchall()
        if not db_entry:
            if verbose:
                print("Updating database")
            #Update database
            lite_file_substring_dict = re.match(LITE_FILE_REGEX, os.path.basename(contents["lite_file"])).groupdict()
            CUR.execute("INSERT INTO created_imagery (filename, var, date, input_product, input_file) VALUES (?, ?, ?, ?, ?)", 
                        ((os.path.basename(contents["out_plot_name"])), contents["var"], lite_file_substring_dict["yymmdd"], contents["product"], os.path.basename(contents["lite_file"])))
            CONN.commit()
        else:
            if verbose:
                print("This file is already in the database")
                print(db_entry)
        
        os.remove(job_file)
        os.remove(LOCKFILE)
        if glob(ISSUE_FILE):
            os.remove(ISSUE_FILE)
    else:
        #if the job was unsuccessful, get rid of the output plot if it exists in any form
        #but leave the job file for debugging
        if verbose:
            print("There was a problem with the job")
            print("The job file has been left for debugging: " + job_file)            
        silent_remove(contents["out_plot_name"])
    
    #Remove intermediate files no matter what
    if contents["rgb"]:
        just_plot_name = os.path.basename(contents["out_plot_name"])
        just_plot_dir = os.path.dirname(contents["out_plot_name"])
        rgb_name = os.path.join(just_plot_dir, re.sub(contents["var"], "RGB", just_plot_name))
        silent_remove(rgb_name)
        silent_remove(rgb["xml"])
        silent_remove(rgb["intermediate_tif"])   

def check_job_worked(plot_name, var, rgb=False):
    """
    Check if the expected output files exist
    """ 
    metadata_name = re.sub("png", "met", plot_name)
    worldfile_name = re.sub("png", "pgw", plot_name)
    
    if rgb:
        if (glob(plot_name) and os.path.getsize(plot_name) > 0 and 
            glob(rgb["layered_rgb_name"]) and os.path.getsize(rgb["layered_rgb_name"]) > 0):
            return True
        else:
            return False
    else:
        if (glob(plot_name) and os.path.getsize(plot_name) > 0 and
        glob(metadata_name) and os.path.getsize(metadata_name) > 0 and 
        glob(worldfile_name) and os.path.getsize(worldfile_name) > 0):
            return True
        else:
            return False

def silent_remove(filename):
    """
    Remove a file if it exists, and keep quiet if it doesn't
    But if there's some other problem with the operation, yell
    """
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: 
            raise
    
if __name__ == "__main__":

    if not os.path.exists(OUT_PLOT_DIR):
        os.makedirs(OUT_PLOT_DIR)
    
    for p in DATA_DICT.keys():
        find_unprocessed_file(p)
    
    if CONN:
        #print("Closing database connection from routine_processing.py")
        CONN.close()
