#!/bin/python

import argparse
import os
import glob
import math
import time
import trim
from datetime import datetime
import numpy as np
from natsort import natsorted


parser = argparse.ArgumentParser()
parser.add_argument(u"-nprocs", type=str, required=True,
                    help=u"Give the number of processors")
parser.add_argument(u"-radius", type=float,default=1.5, 
                    help=u"Give the radius manually")
parser.add_argument(u"-basename", "-bn", type=str, required=True, 
                    help=u"Give the job basename (e.g. bem for'bem.1.inp' ")
args = parser.parse_args()
nprocs = args.nprocs
basename = args.basename
# Store current working dir as variable

present_dir = os.getcwd()
#
#basename = trim.get_basename()

# Read the input file <input.txt>
trim_iter, trim_meth, dCd_dcoll, dCq_dcoll, avg_over, \
trim_meth, extra_force, cm_tail, gridgen_path, tclfile, tail_line, rtype = trim.readInput("input.txt")

# Get number of rotors based on the inputs
num_rotors = trim.get_num_rotors(rtype)

# Write the date to output file
output = open("trim_output.log", "a")
print >> output, "\n"
print >> output, datetime.now()
output.close()

# Get the force and moment history files
fomoco_file = trim.get_output_file('fomoco', basename)
if rtype == 'bem':
    rdisk_file = trim.get_output_file('rdisk_trim', basename)
elif rtype == 'blades':
    hist_files = natsorted(glob.glob(u"*.hist.txt"))

# Read the OVERFLOW input file so I can dimensionalize forces
vinf, vref, alpha, beta, omega, a, tinf, rho = trim.readOverflowInputFile(basename + '.1.inp', args.radius)

# Read the fomoco file to determine how many components I have
fomoco_obj = trim.number_components(fomoco_file)

# Read mixsur to get aref, lref, and cg of each component
mixsur = 'mixsur.inp'
lref, aref, cg_location, comp_refs = trim.read_mixsur(mixsur, fomoco_obj)

# The number of trim iterations are set in the input file
for i in range(trim_iter):

    # Case #1 or #2 (2 is just restart from 1)
    case = trim.casenum(fomoco_file)

    # Comment this line out if I don't want to run OVERFLOW
    print('overrunmpi -np ' + nprocs + ' ' + basename + ' ' + str(case))
    #os.system('overrunmpi -np ' + str(nprocs) + ' ' + str(basename) + ' ' + str(case))
    
    # Get the FOMOCO forces for all components
    fomoco_avg, minimum, maximum, comp_names = trim.fomoco_wind(
            fomoco_file, avg_over, fomoco_obj, omega, args.radius, lref, aref, rho)

    # Determine which fomoco component is the 'TARGET'
    fomoco_target, lref_target, aref_target = trim.get_target(fomoco_avg, comp_names, lref, aref)
    
    # Get the non-dimensionalized extra force from input file
    nondim_extra = trim.nondim_force(aref_target, lref_target, rho, vref, vinf, extra_force, num_rotors)
    
    # Assign target values to simpler variables
    CD_target = fomoco_target[0] + nondim_extra
    CM_target = fomoco_target[4]    
    
    # Output current iteration to log file
    output = open("trim_output.log", "a")
    print >> output, "\n########################\n     ITERATION: #" + str(i + 1)
    print >> output, "########################\n"
    output.close()

    if rtype == 'bem':
    # I'm getting average rdisk forces and returning them in order from r1 - r6
        rdisk_avg = trim.readRdisk(rdisk_file, num_rotors, avg_over)
        
    elif rtype == 'blades':
    # Rotors are automatically sorted based on the natsorted hist_files command above
#        rdisk_avg = trim.rotor_fomoco(fomoco_avg) # Remember these are in the wind frame, so I'm not sure I can use that...
        comp_range = np.arange(-4, -1)
        rdisk_avg = fomoco_avg[comp_range, :]
   
    if abs(CD_target) > 0.0010 or abs(CM_target) > 0.0010:
        dr1_coll, dr2_coll, dr3_coll, dr4_coll, dr5_coll, dr6_coll, r5_cq, \
        r6_cq = trim.trimRotors(trim_meth, fomoco_target, rdisk_avg, dCq_dcoll, dCd_dcoll, nondim_extra, num_rotors, rtype)

        # Write the rdisk input files based on the delta
        coll_old, coll_new, coll_delta = trim.writeRotorInp(num_rotors, trim_meth, 
                                                            dr1_coll, dr2_coll, dr3_coll, \
                                                            dr4_coll, dr5_coll, dr6_coll, \
                                                            fomoco_target, rdisk_avg, r5_cq, \
                                                            r6_cq, nondim_extra, rtype)

        if abs(CD_target) < 0.00001:
            dtail = trim.trimTail(cm_tail, fomoco_target, trim_meth)    
            # Make the tail grid by running the bash script
            told, dtail, tnew = trim.makeTail(num_rotors, dtail, gridgen_path, tclfile, tail_line)
        else:
            dtail = 0
            told = trim.get_old_tail(tclfile, tail_line)
            tnew = dtail + told

        
        trim.write_plot_out(i, CD_target, CM_target, r5_cq, coll_old, coll_delta, coll_new, told, dtail, tnew)
   
    else:
        output = open("trim_output.log","a")
        print >> output, "\n *** Reached the trim tolerance. Exiting the loop. *** \n"
        output.close()
        break


