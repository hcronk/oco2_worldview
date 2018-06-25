import os
import sys
import argparse
from glob import glob

if __name__ == "__main__":
    
    code_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
    
    parser = argparse.ArgumentParser(description="Flags", prefix_chars='@')
    parser.add_argument("@v", "@@verbose", help="Prints some basic information during code execution", action="store_true")
    parser.add_argument("@f", "@@file", help="Full path to input file name for command line use", default=None)
    parser.add_argument("@c", "@@config", help="Full path to text file containing a list of files to process, one filename per line", default=None)
    parser.add_argument("@g", "@@geolocation", help="Custom geolocation box for testing and case study purposes. Fmt: [lat_min, lat_max, lon_min, lon_max]", nargs = '+', default=[])
    parser.add_argument("@r", "@@rgb", help="Overlay plots on MODIS RGB for testing and case study purposes", default=os.path.join(code_dir, "GIBS_Aqua_MODIS_truecolor.xml"))
    parser.add_argument("@o", "@@output_dir", help="Output directory for plots", default=os.path.join(code_dir, "plots"))
    parser.add_argument("@a", "@@vars", help="Variables to plot", nargs = '+', default=[])
    args = parser.parse_args()
    verbose = args.verbose
    lite_file = args.file
    config_file = args.config
    custom_geo_box = args.geolocation
    rgb = args.rgb
    out_plot_dir = args.output_dir
    user_defined_var_list = args.vars
    
    if user_defined_var_list:
        user_defined_var_list = [x.strip("[,]") for x in user_defined_var_list]
    
    if not os.path.exists(out_plot_dir):
        os.makedirs(out_plot_dir)
        
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
    
    if custom_geo_box:
        custom_geo_box = [float(x.strip("[,]")) for x in custom_geo_box]
        if len(custom_geo_box) != 4:
            print("Custom geolocation box format: [lat_min, lat_max, lon_min, lon_max]")
            print("Exiting.")
            sys.exit()
