import os
import sys
import subprocess
from datetime import datetime
from datetime import timedelta
import numpy as np
import math

# Helper function to map arbitrary number between 0 and 2*pi
# using modular arithmetic
def mod_circle(angle):
    m_pi = math.pi
    
    return angle-((2*m_pi)*math.floor(angle/(2*m_pi)));

# Function to compute the projection from lat/lon grid
# to cartesian grid using CPP projection.
# Output grid points rounded to nearest meter
#def compute_projection(lat1, lat2, lon1, lon2, ref_lat, ref_lon):
def compute_projection(input_options):
   
    # Define constants
    Rearth = 6378206.4
    m_pi = math.pi

    # Get values from dictionary
    lat1, lat2 = float(input_options['lat1']), float(input_options['lat2'])
    lon1, lon2 = float(input_options['lon1']), float(input_options['lon2'])
    ref_lat    = float(input_options['reference_lat'])
    ref_lon    = float(input_options['reference_lon'])

    # Set reference values in radians
    phi0 = ref_lat*(m_pi/180.)
    lam0 = mod_circle(ref_lon*(m_pi/180.))
    
    # Location of split determines negative or positive x-values
    lam_split = mod_circle(lam0+m_pi)

    # Convert lat/lon from degrees to radians
    phi = lat1*(m_pi/180.)              # Latitude of first grid point
    lam = mod_circle(lon1*(m_pi/180.))  # Longitude of first grid point

    # Determine first point (E->W, N->S)
    x1 = Rearth * math.cos(phi0) * (mod_circle(lam-lam_split) - m_pi)
    y1 = Rearth * (phi - phi0)
   
    # Determine second point
    x2 = x1 + (lon2-lon1) * (m_pi/180.) * Rearth * math.cos(phi0);
    y2 = y1 - (lat1-lat2) * (m_pi/180.) * Rearth;

    # Round values to meters
    x1 = round(x1)
    x2 = round(x2)
    y1 = round(y1)
    y2 = round(y2)
    
    return x1, x2, y1, y2

# Parse the inventory file for the time in hours of each entry.
# Also record the number of entries for wind velocity pairs and
# determine which wind velocity is given first (u_first=0 or 1)
def parse_inventory(inventory_file, ref_time):

    # Try to open inventory file
    try:
        f = open(inventory_file, 'r')
    except:
        print("Could not open inventory file.")
        exit
   
    # Create array for storing entry times
    entry_times = []

    # Loop over the lines of the inventory file
    for count, line in enumerate(f, start=1):
        
        # Split line for distinct info
        split_line = line.split(":")
        
        # Check which direction is first (u- and v- winds should alternate)
        if count == 1:
            if split_line[3] == "10U":
                u_first = 1
            elif split_line[3] == "10V":
                u_first = 0
            else:
                print("Warning. Parse inventory function assumes variables are \
10m winds in u- and v-directions. \
\nIn order to use new variables, parts of this utility must be updated. \
\nExiting without performing interpolation.")
                exit()
        
        # Skip every other line since we're only interested in dates
        # IMPORTANT!: u- and v- winds assumed to have matching dates
        # and alternate in the grib file.
        if count%2 == 0:
            continue

        # Read in the date and the hours of the record
        date = split_line[2][2:]
        hour = int(split_line[8][3:])
        
        # Convert to from YY to YYYY
        # We assume date is pre-2050 here
        if int(date[0]) > 5:
            date = "19"+date
        else:
            date = "20"+date
            
        # Create datetime object for date of entry with hours
        datetime_object = datetime.strptime(
                date, '%Y%m%d%H')
        hours_added = timedelta(hours = hour)
        datetime_object = datetime_object + hours_added

        # Find difference between entry and reference dates
        tdelta = datetime_object - ref_time

        # Convert to total hours
        hours_from_ref = int(tdelta.seconds/3600) + int(tdelta.days*24)

        # Append to list
        entry_times.append(hours_from_ref)

    # Record number of times at which u/v wind vels recorded
    num_times = count/2
    
    # Close inventory file
    f.close()

    return entry_times, num_times, u_first

# Read the input file
def read_input_file(filename):

    # Dictionary for storing the input options from input file
    input_options = {}

    # Open the file
    try:
        f = open(filename, "r")
    except:
        print("Error. Could not open ",filename," Exiting.")
        exit()

    # List of necessary input options; we will check against
    # this list to ensure user gave all necessary options
    var_names_dict = {'wgrib_exe','interpolator', \
            'grib_file','adh_root', \
            'lat1','lat2','lon1','lon2','ni','nj', \
            'reference_time','reference_lat','reference_lon', \
            'binary_format','variables_string'}

    # Read lines of file
    options = f.readlines()

    # Loop over lines and extract variables (ignoring '#' comments)
    for line in options:
        if line.startswith('#'):
            continue
        try:
            key, value = line.split(":",1)
            input_options[key.strip()] = value.strip()
        except:
            continue
   
    # Close input file
    f.close()

    # Check that we have all necessary input options
    missing_input_variable=False
    for var in var_names_dict:
        if var not in input_options:
            print("Missing ", var, " in input file.")
            missing_input_variable=True

    # If missing an option quit
    if missing_input_variable:
        print("Exiting without creating interpolated winds file.")
        exit()

    return input_options


### Main function ###
def main():

    if len(sys.argv) < 2:
        print("Error. Need input filename with extension.")
        exit()

    # Set input filename from command line
    input_filename = sys.argv[1]

    # Read the input options
    input_options = read_input_file(input_filename)

    ## Set the reference time for AdH
    # Ref time is in YYYY, MM, DD, HH
    # We assume date is pre-2050 here
    reference_time = datetime.strptime(input_options['reference_time'],"%Y-%m-%d-%H")

    # Inventory, binary, and time files given the same name
    # as the GRIB file but with different extension
    grib_file = input_options['grib_file']
    inventory_file = grib_file+".inventory"
    binary_file    = grib_file+".bin"
    time_file      = grib_file+".time"

    # Add .grib extension to GRIB file
    grib_file      = grib_file+".grib"

    # Create inventory file using wgrib
    # Only select records of 10m u- and v- wind velocity
    os.system(input_options['wgrib_exe']  + " " + grib_file
            + " | grep -E " + input_options['variables_string']
            + " > " + inventory_file)

    # Use inventory file to create binary file with
    # values of selected parameters
    os.system("cat " + inventory_file + " | "
            + input_options['wgrib_exe'] + " -i " + grib_file
            + " -o " + binary_file
            + " 1> /dev/null")

    # Read the inventory file for the times of u/v wind vel records.
    # Also get the number of pairs of records and determine which
    # direction appears first.
    entry_times, num_times, u_first = parse_inventory(inventory_file, reference_time)

    # Create file with the times in hours past reference date
    f = open(time_file, "w")
    for time in entry_times:
        f.write(str(time) + "\n")
    f.close()

    # Compute the projection from lat-lon to cartesian
    x1, x2, y1, y2 = compute_projection(input_options)

    cartesian_grid = str(x1) + " " + str(x2) + " " \
            + str(y1) + " " + str(y2) + " " \
            + input_options['ni'] + " " + input_options['nj']
    
    # Call interpolation routine to interpolate
    # ERA5 values from binary_file to the AdH grid
    os.system(input_options['interpolator'] + " "
            + input_options['adh_root'] + " "
            + time_file                 + " "
            + binary_file               + " "
            + str(u_first)              + " "
            + str(num_times)            + " "
            + cartesian_grid            + " "
            + input_options['binary_format'])


if __name__ == "__main__":
    main()
