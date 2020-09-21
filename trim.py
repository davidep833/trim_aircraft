#!/usr/bin/python

# Import the usual modules - some may not be needed
import glob
import os
import f90nml
import math
import numpy as np
from scipy.constants import convert_temperature
import natsort
import string
import re
import unicodedata

# GAS CONSTANT IN UNITS J/(KG-K)
gas_constant = 287.05

# RATIO OF SPECIFIC HEATS
gamma = 1.4

# ASSUMPTIONS
# 1: "TARGET" IS THE NAME OF THE COMPONENT TO TRIM
# 2: NFOMO IS SET TO 1 SO THAT THE OUTPUT MATCHES THE RDISK_TRIM.OUT SAVE FREQUENCY
# 3: EXTRA FORCE/DRAG IS IN M^2
# 4: IF TRIM ITERATIONS IN "INPUT.TXT" IS != 0 THEN RUN THROUGH 3 ITERATIONS


def number_components(fomoco_file):
    # READ THE FOMOCO FILE AND GET COEFFICIENTS FOR TRIMMING
    # DECLARE ARRAYS SO THAT I CAN POPULATE THEM WITH FOMOCO DATA
    raw_data = []
    count = 0

    #############################
    # READ OVERFLOW FOMOCO FILE #
    #############################

    fomoco_raw = open(fomoco_file, u"r")
    for lines in fomoco_raw:
        junk = lines.split()
        for columns in junk:
            junk = columns.split()
            junk = u''.join(unicode(x) for x in junk)  # CONVERT LISTS TO STRINGS
            raw_data.append(junk)
    length = len(raw_data)
    fomoco_raw.close()

    ###############################################
    # ASSIGN DATA IN FOMOCO.OUT FILE TO VARIABLES #
    ###############################################

    for n in xrange(0, length, 39):
        Comp = raw_data[n + 2]
        CompNum = u''.join(unicode(x) for x in Comp)
        CompNum = abs(int(CompNum))
        count = count + 1

    if count != abs(CompNum):
        fomoco_obj = abs(CompNum)

    return fomoco_obj


# READ THE OVERFLOW INPUT FILE AND RETURN FLOINP AND NQT
def readOverflowInputFile(basename, radius):
    parser = f90nml.Parser()
    overflow_input = parser.read(basename)
    alpha = overflow_input[u'floinp'][u'alpha'] # DEGREES
    beta = overflow_input[u'floinp'][u'beta']  # DEGREES
    fsmach = overflow_input[u'floinp'][u'fsmach']
    rey = overflow_input[u'floinp'][u'rey']   # RE/(GRID-UNIT)

    if u'tinf' in overflow_input[u'floinp']:
        tinf_tmp = overflow_input[u'floinp'][u'tinf']  # RANKINE
    else:
        tinf_tmp = 518.7
    if u'refmach' in overflow_input[u'floinp']:
        refmach = overflow_input[u'floinp'][u'refmach']
    else:
        refmach = fsmach

# CALCULATE OMEGA AND VINF BASED ON NON-DIMENSIONAL INPUTS
    tinf = convert_temperature(tinf_tmp, u'Rankine', u'Kelvin')
    a = math.sqrt(gamma*gas_constant*tinf)  # M/S
    omega = refmach*a/radius                         # RAD/S
    vinf = fsmach*a                            # M/S
    vref = refmach*a
    nu = (1E-03/rey)*(refmach*a)

# CALCULATE THE DENSITY USING SUTHERLAND'S FORMULA FOR VISCOSITY
    mu0 = 1.716E-05 # lb-s/ft^2
    T0 = 273.15    # Rankine
    mu = mu0 * (tinf/T0)**1.5 * ((T0 + 110.4)/(tinf + 110.4))
    rho = mu/nu
    print u"rho = " + str(rho)

    return vinf, vref, alpha, beta, omega, a, tinf, rho


def readInput(input_filename):
    #######################
    # READ THE INPUT FILE #
    #######################
    inputs = {}
    # Make the input file a dictionary to read from
    with open(input_filename, "r") as file:
        comment_less = filter(None, (line.split('#')[0].strip() for line in file.readlines()))
        inputs = {item.split('=')[0].strip(): item.split('=')[1].strip() for item in comment_less}
        
    avg_over = int(inputs['avg_over'])
    trim_meth = int(inputs['trim_method'])
    dCd_dcoll = float(inputs['Cdcoll'])
    dCq_dcoll = float(inputs['Cqcoll'])

    cm_tail = float(inputs['Cmtail'])
    trim_iter = int(inputs['trim_iters'])  # define number of trim iterations
    extra_force = float(inputs['extra_force'])  # climb thrust/extra cooling,excresence drag     
    gridgen_path = inputs['gridgen_path']
    tclfile = inputs['gridgen_file']
    tail_line = int(inputs['tail_line'])
    rtype = inputs['rotor_type'].lower() # Make lowercase for the if statement

    return trim_iter, trim_meth, dCd_dcoll, dCq_dcoll, avg_over, \
           trim_meth, extra_force, cm_tail, gridgen_path, tclfile, tail_line, rtype
  
def get_num_rotors(rtype):
    if rtype == 'bem':
        rotor_inp = natsort.natsorted(glob.glob('rotor.r*.inp'))
    elif rtype == 'blades':
        rotor_inp = natsort.natsorted(glob.glob(u"*.hist.txt"))
    num_rotors = len(rotor_inp)
    return num_rotors
            
            
def nondim_force(aref_target, lref_target, rho, vref, vinf, extra_force, num_rotors):
    drag = 0.5 * rho * vinf**2 * extra_force
    
    if num_rotors == 3:
        aref = aref_target*2
    else:
        aref = aref_target
    nondim_force = drag  / (0.5 * rho * aref * vref**2)
    return nondim_force

# GET FOMOCO COEFFICIENTS IN BODY FRAME
def fomoco_wind(fomoco_file, avg_over, fomoco_obj, omega, radius, lref, aref, rho):
    # READ THE FOMOCO FILE AND GET COEFFICIENTS FOR TRIMMING
    # DECLARE ARRAYS SO THAT I CAN POPULATE THEM WITH FOMOCO DATA
    raw_data = []
    fomoco_data = []
    comp_names = []
    name, Iter, Comp, Area, Area_X, Area_Y, Area_Z = [], [], [], [], [], [], []
    CxPres, CyPres, CzPres, CxVis, CyVis, CzVis = [], [], [], [], [], []
    CxMom, CyMom, CzMom, CLPres, CDPres, CSPres = [], [], [], [], [], []
    CLvis, CDvis, CSvis, CLmom, CDmom, CSmom = [], [], [], [], [], []
    CmRoll, CmPitch, CmYaw, MassF, Time = [], [], [], [], []
    CmRoll, CmPitch, CmYaw, MassF, Time = [], [], [], [], []
    CMXp, CMYp, CMZp, CMXv, CMYv, CMZv = [], [], [], [], [], []
    CMXm, CMYm, CMZm = [], [], []
    count = 0

#############################
# READ OVERFLOW FOMOCO FILE #
#############################

    fomoco_raw = open(fomoco_file, u"r")
    for lines in fomoco_raw:
        junk = lines.split()
        for columns in junk:
            junk = columns.split()
            junk = u''.join(unicode(x) for x in junk)  # CONVERT LISTS TO STRINGS
            raw_data.append(junk)
    length = len(raw_data)
    fomoco_raw.close()

###############################################
# ASSIGN DATA IN FOMOCO.OUT FILE TO VARIABLES #
###############################################

    for n in xrange(0, length, 39):
        name = unicode(raw_data[n])
        Comp = raw_data[n + 2]
        CxPres = raw_data[n + 7]
        CyPres = raw_data[n + 8]
        CzPres = raw_data[n + 9]
        CxVis = raw_data[n + 10]
        CyVis = raw_data[n + 11]
        CzVis = raw_data[n + 12]
        CxMom = raw_data[n + 13]
        CyMom = raw_data[n + 14]
        CzMom = raw_data[n + 15]
        CLPres = raw_data[n + 16]
        CDPres = raw_data[n + 17]
        CSPres = raw_data[n + 18]
        CLvis = raw_data[n + 19]
        CDvis = raw_data[n + 20]
        CSvis = raw_data[n + 21]
        CLmom = raw_data[n + 22]
        CDmom = raw_data[n + 23]
        CSmom = raw_data[n + 24]
        CmRoll = raw_data[n + 25]
        CmPitch = raw_data[n + 26]
        CmYaw = raw_data[n + 27]
        CompNum = u''.join(unicode(x) for x in Comp)
        CmRoll = u''.join(unicode(x) for x in CmRoll)
        CmPitch = u''.join(unicode(x) for x in CmPitch)
        CmYaw = u''.join(unicode(x) for x in CmYaw)
        Fx_wind = float(CDPres) + float(CDvis) + float(CDmom)
        Fy_wind = float(CSPres) + float(CSvis) + float(CSmom)
        Fz_wind = float(CLPres) + float(CLvis) + float(CLmom)
        CompNum = int(CompNum)
        count = count + 1
  
        fomoco_data.append([Fx_wind, Fy_wind, Fz_wind, float(CmRoll), 
                            float(CmPitch), float(CmYaw)])
        if count <= fomoco_obj:
            comp_names.append(name)

    fomoco_data = np.array(fomoco_data)
    count = count - avg_over * fomoco_obj - 1
    size = (avg_over, fomoco_obj, 6)
    fomoco_data2 = np.zeros(size)

    # CREATE 3-D MATRIX WITH
    for k in xrange(avg_over):
        for j in xrange(fomoco_obj):
            count = count + 1
            for i in xrange(6):
                fomoco_data2[k, j, i] = fomoco_data[count, i]

    fomoco_avg = 0.0
    size = (fomoco_obj, 6)
    fomoco_avg = np.zeros(size)
    fomoco_min = np.zeros(size)
    fomoco_max = np.zeros(size)

    # Calculate the average force and moment coeffients for all the components.
    for i in xrange(fomoco_obj):
        for j in xrange(6):
            fomoco_avg[i, j] = np.mean(fomoco_data2[:, i, j])
            fomoco_min[i, j] = np.min(fomoco_data2[:, i, j])
            fomoco_max[i, j] = np.max(fomoco_data2[:, i, j])

    return fomoco_avg, fomoco_min, fomoco_max, comp_names


def read_mixsur(mixsur, fomoco_obj):
    # Declare lists
    data, cg_location, lref, aref, xmc, ymc, zmc = [], [], [], [], [], [], []
    tmp_refs, tmp2_refs, tmp3_refs, comp_refs = [], [], [], []

    abc = string.ascii_lowercase
    ABC = string.ascii_uppercase
    alphabet = list(abc + ABC)
    
    count = 0
    with open(mixsur, mode='r') as mixsurfile:
        junk = mixsurfile.readlines()
        for line in junk:
            x = re.split(", |\s |,", line)
        # This part is great for getting aref, lref, and moment centers, but not the component moment centers
            if len(x) >= 1:
                data.append(x)
        data2 = []
        for x in data:
            y, tmp = [], []
            for y in x:
                if y != '':
                    tmp.append(y)
            data2.append(tmp)
        
        count = 0
        
        # determine where the reference components line starts
        for line in data2:
#            if len(line) >= 8:
            if line[0][0] in alphabet:
                end_line = count
                break
            count += 1
        # I want one list with all of the mixsur file lines in it
        tmp_refs.extend(data2)
        data = data2
        nref = int(data[1][0])
        for i in xrange(nref):
            xmc.append(float(data[i+2][2]))
            ymc.append(float(data[i+2][3]))
            zmc.append(float(data[i+2][4]))
            lref.append(float(data[i+2][0]))
            aref.append(float(data[i+2][1]))
    cg_location = (xmc, ymc, zmc)
    cg_location = np.array(cg_location)
    cg_location = np.transpose(cg_location)
    lref = np.array(lref)
    aref = np.array(aref)
    cg_location_tmp = cg_location / 1E+03
    lref_tmp = lref / 1E+03 # convert from mm to m
    aref_tmp = aref / 1E+06 # convert from mm^2 to m^2
    
    # I need to convert all the data to str or int instead of unicode to properly read the rest of the file
    for row in tmp_refs:
        z = []
        for column in row:
            txt = column
            z.append(txt)           
        tmp2_refs.append(z)

    # Read through the components and get the reference conditions for each integration
    # by taking the NREF after the string.
    # e.g. After "TARGET" get the reference component.
    tmp3_refs = tmp2_refs[(end_line-1):]
    for i in range(len(tmp3_refs)):
        if tmp3_refs[i][0][0] in alphabet:
            comp_refs.append(int(tmp3_refs[i+1][1]))
                
    
    cg_location, aref, lref = [], [], []
    for i in range(fomoco_obj):
    # convert to python indexing with 0 being the first fomoco object
        tmp_ref = int(comp_refs[i]) - 1
        cg_location.append(cg_location_tmp[tmp_ref].tolist())
        aref.append(aref_tmp[tmp_ref].tolist())
        lref.append(lref_tmp[tmp_ref].tolist())
    return lref, aref, cg_location, comp_refs

# I'm getting average rdisk forces and returning them in order from r1 - r6
def readRdisk(rdisk_file, num_rotors, avg_over):
    ############################
    # READ RDISK_TRIM.OUT FILE #
    ############################

    # READ RDISK_TRIM.OUT FILE. RDISK SHOULD WRITE OUT EVERY TIMESTEP (SO SHOULD FOMOCO)

    rdisk_filename = open(rdisk_file, "r")
    rdisk_data_raw = []
    count = 0
    for lines in rdisk_filename:
        junk = lines.split()
        count = count + 1
        for columns in junk:
            junk = columns.split()
            junk = ''.join(str(x) for x in junk)  # CONVERT LISTS TO STRINGS
            rdisk_data_raw.append(junk)

    rdisk_filename.close()
    length = len(rdisk_data_raw)
    size = (num_rotors, 11)
    rdisk_avg = np.zeros(size)
    size = (int(count / num_rotors), num_rotors, 11)
    rdisk_data = np.empty(size)
    rdisk_data_tmp = []

    for n in range(0, length, 13):
        if n > length - 13:
            break
        rotor = rdisk_data_raw[n + 2]
        istep = rdisk_data_raw[n + 3]
        nit = rdisk_data_raw[n + 4]
        omega = rdisk_data_raw[n + 5]
        a0 = rdisk_data_raw[n + 6]
        a1 = rdisk_data_raw[n + 7]
        b1 = rdisk_data_raw[n + 8]
        CT = rdisk_data_raw[n + 9]
        CQ = rdisk_data_raw[n + 10]
        CMS = rdisk_data_raw[n + 11]
        CMR = rdisk_data_raw[n + 12]
        rdisk_data_tmp.append(
            [int(rotor), int(istep), int(nit), float(omega), float(a0), float(a1), float(b1), float(CT), float(CQ),
             float(CMS), float(CMR)])

    # MAKE THE DATA A NUMPY ARRAY SO I CAN MANIPULATE IT EASIER

    rdisk_data_tmp = np.array(rdisk_data_tmp)
    count = len(rdisk_data_tmp)
    num = 0

    for i in range(int(count / num_rotors)):
        for j in range(num_rotors):
            rdisk_data[i, j, :] = rdisk_data_tmp[num, :]
            num = num + 1


    # Just writing out the rotor numbers as positive numbers
    for i in range(num_rotors):
        rdisk_avg[i, 0] = abs(rdisk_data[0, i, 0])
    
    # Calculate average rotor thrust and torque
    for i in range(int(count / num_rotors) - avg_over, int(count) // num_rotors):
        for j in range(num_rotors):
            rdisk_avg[j, 1:] = rdisk_avg[j, 1:] + rdisk_data[i, j, 1:] / avg_over
            
    # Sort the array in ascending order from R1 to R6
    average = rdisk_avg[np.argsort(rdisk_avg[:, 0])]

    return average


###################
# TRIM THE ROTORS #
###################

def trimRotors(trim_meth, fomoco_target, rdisk_avg, dCq_dcoll, dCd_dcoll, extra_force, num_rotors, rtype):
    output = open("trim_output.log", "a")        
    ''' 
    Here I'm cycling through each variable in the rdisk_avg to find out 
    which rotor it is. I don't know if the rotors will be in order from 1-6 or listed 
    as 1, 2, 5, 3, 4, and 6 in the rdisk_trim.out file.
    '''
    if num_rotors != 6:
        r3_avg = np.zeros_like(rdisk_avg[0, :])
        r4_avg = np.zeros_like(rdisk_avg[0, :])
        r6_avg = np.zeros_like(rdisk_avg[0, :])
    
    '''
    Use the average angular rate variable "omega" from "rdisk_trim.out" 
    file to determine the rotational direction of the tail rotor. 
    Return those values in a list in ascending rotor order.
    '''
    rdir = []
    if rtype == 'bem':
        for val in rdisk_avg[:,3]:
            if val > 0:
                rdir.append(1)
            else: 
                rdir.append(-1)
                
    elif rtype == 'blades':
        # rotor direction is set for R1, R2, and R5
        rdir = [-1, 1, -1]
    
    if trim_meth == 1:
        print >> output, "\nTrim Method #1 - Trim R1-4 to drag/extra force and R5/6 to zero torque."  # TRIM R1 & R2 TO DRAG/EXTRA FORCE, TRIM R5 TO 0 TORQUE
        
        # Get the information to trim the tail rotors to torque
        if num_rotors == 3:
            r5_avg = rdisk_avg[2, :]
            r6_avg = np.zeros_like(r5_avg)
            r5_dir = rdir[2]
            r6_dir = 1 # dummy number that's not used in trim from 1/2 aircraft
            num_front = 2
        elif num_rotors == 6:
            r5_avg = rdisk_avg[4, :]
            r6_avg = rdisk_avg[5, :]
            r5_dir = rdir[4]
            r6_dir = rdir[5]
            num_front = 4
        
        if rtype == 'bem':
            col = 8
        elif rtype == 'blades':
            # Hist outputs are in the inertial frame, so we want Mx = 0, not Mz = 0
            col = 3
        
        r5_cq = r5_avg[col]
        r6_cq = r6_avg[col]
        
        dr5_coll = (r5_cq * r5_dir / dCq_dcoll)  # NEED DIRECTION TO KNOW IF ABOUT +/- TORQUE
        dr6_coll = (r6_cq * r6_dir / dCq_dcoll)  # NEED DIRECTION TO KNOW IF ABOUT +/- TORQUE
        
        dpitches = [dr5_coll, dr6_coll]
        # apply limiter, so the collective doesn't change too quickly
        n = 1.0 # This is the limiter in degrees
        for i in range(len(dpitches)):
            if dpitches[i] != 0.0:
                sign = dpitches[i] / abs(dpitches[i])
                dpitches[i] = min(n, abs(dpitches[i])) * sign
            else:
                dpitches[i] = 0.0
        dr5_coll, dr6_coll = dpitches[:]
        
        dr1_coll = ((fomoco_target[0] + extra_force) / num_front - dCd_dcoll * (dr5_coll + dr6_coll) / num_front) / dCd_dcoll
        dr2_coll = dr1_coll
        dr3_coll = dr1_coll
        dr4_coll = dr1_coll
        
        dpitches = [dr1_coll, dr2_coll, dr3_coll, dr4_coll]
        # apply limiter, so the collective doesn't change too quickly
        n = 1.0 # This is the limiter in degrees
        for i in range(len(dpitches)):
            if dpitches[i] != 0.0:
                sign = dpitches[i] / abs(dpitches[i])
                dpitches[i] = min(n, abs(dpitches[i])) * sign
            else:
                dpitches[i] = 0.0
                
        dr1_coll, dr2_coll, dr3_coll, dr4_coll = dpitches[:]
    
    elif trim_meth == 2:
        print >> output, "\nTrim Method #2 - Trimming all rotors to same pitch"  # TRIM ALL ROTORS TO SAME COLLECTIVE PITCH
        dr1_coll = ( (fomoco_target[0] + extra_force) / num_rotors) / dCd_dcoll
        dr2_coll = dr1_coll
        dr3_coll = dr1_coll
        dr4_coll = dr1_coll
        dr5_coll = dr1_coll
        dr6_coll = dr1_coll
    
    elif trim_meth == 3:
        # Trim the a/c in hover
        print >> output, "\nTrim Method #3 - Trimming to hover thrust and pitching moment" 
        
        '''        
        Trimmming in hover to the aircraft weight and pitching moment. From
        Jeremy: "... trim hover to a value and PM to zero by varying the mean 
            collective and the delta (r5-r2) ..."
        The target lift would be sum of the forces of weight and rotor thrust
            - fomoco_target[2] = lift of rotors + a/c
            - dCd_dcoll is the change in thrust per change in degree collective
        '''            
        # Get information from rdisk_avg to trim to hover thrust and pitching moment
        r2_avg = rdisk_avg[1, :]
        r5_avg = rdisk_avg[4, :]
        
        # This is the thrust target we want where extra force is input in Newtons.
        target_force = fomoco_target[2] + extra_force
        
        # The current pitching moment of the a/c, cmy
        cm_pitch = fomoco_target[4]
        mean_coll = 0
       
        # Take the average of the collectives
        for i in rdisk_avg[:,4]:
            mean_coll = mean_coll + i/num_rotors
        
        # Newton's method to trim thrust based on mean collective pitch
        dcoll = target_force / num_rotors / dCd_dcoll
        dr_old = r5_avg[4] - r2_avg[4]
        dcmy_dr = 0.01
        dr_new = dr_old - cm_pitch/dcmy_dr

        # Now update the deltas based on mean collective and cmy
        dr1_coll = dcoll
        dr2_coll = r2_avg[4] + abs(dr_new)/2 + dcoll
        dr5_coll = r5_avg[4] - abs(dr_new)/2 + dcoll
        dr4_coll = dr1_coll
        dr3_coll = dr2_coll
        dr6_coll = dr5_coll
            
        
        
    elif trim_meth == 4:
        # Trim all props to zero torque
        print >> output, "\nTrim Method #4 - Trimming All props to zero torque"

        # Get information from rdisk_avg to trim to hover thrust and pitching moment
        r1_avg = rdisk_avg[0, :]
        r2_avg = rdisk_avg[1, :]
        r5_avg = rdisk_avg[2, :]

 # using s4 dir's other dir values not in this function

        dr1_coll = (r1_avg[8] * -1 / dCq_dcoll)  # NEED DIRECTION TO KNOW IF ABOUT +/- TORQUE
        dr2_coll = (r2_avg[8]   / dCq_dcoll)  # NEED DIRECTION TO KNOW IF ABOUT +/- TORQUE
        dr5_coll = (r5_avg[8] * -1 / dCq_dcoll)  # NEED DIRECTION TO KNOW IF ABOUT +/- TORQUE


        dr4_coll = dr1_coll
        dr3_coll = dr2_coll
        dr6_coll = dr5_coll
        r5_cq = r5_avg[8]
        r6_cq = r6_avg[8]

        
    output.close()

    return dr1_coll, dr2_coll, dr3_coll, dr4_coll, dr5_coll, dr6_coll, r5_cq, r6_cq


#############################
# TRIM THE TAIL RUDDERVATOR #
#############################

def trimTail(cm_tail, fomoco_target, trim_meth):
    if trim_meth == 1 or trim_meth == 2 or trim_meth == 4:    
        delta = fomoco_target[4] / cm_tail
        n = 1 # This is the limiter
        sign = delta/abs(delta) # We need the sign of the delta
        dtail = min(n, abs(delta)) * sign
    elif trim_meth == 3:
    # Don't use the tail to trim in hover
        dtail = 0.0
    return dtail

def get_target(fomoco_avg, comp_names, lref, aref):
    target = 'TARGET' 
    index = 0
    for i in comp_names:
        if i == target:
            num_target = index
        index = index + 1
    fomoco_target = fomoco_avg[num_target]
    aref_target = aref[num_target]
    lref_target = lref[num_target]
    return fomoco_target, lref_target, aref_target
        

def get_old_rdisk(num_rotors):
    parser = f90nml.Parser()
    
    # Create list of rotor input files. Read old collectives from these files.
    rotor_files = natsort.natsorted(glob.glob("rotor.r*.inp"))
    rotor_names = []
    for file in rotor_files:
        junk = file.replace('rotor.','')
        junk = junk.replace('.inp','')
        rotor_names.append(junk)
        
    coll_old = []
    # Read all the old rdisk files to get the old collectives so I can upate to the new collectives
    for i in range(num_rotors):
        rdisk_input = parser.read(rotor_files[i])
        coll_old.append(rdisk_input['rdisk_trim']['A0'])
    
    return coll_old, rotor_files, rotor_names

def readMotiontxt(motion):
    data = []
    with open(motion, u'r') as inputfile:
        for lines in inputfile:
            x = lines.split()
            data.append(x)
    nspan = int(data[1][0])
    npsi = int(data[1][1])
    data2 = np.zeros((nspan, npsi))
    count = nspan+3
    data1 = data[count:-1]
    count = 0
    for i in xrange(npsi-1):
        for n in xrange(nspan):
            data2[n, i] = data1[(i+1)*3+count][3]
            count = count + 1
    data2[:, npsi-1] = data2[:, 1]
    pitch = np.mean(data2)
    return pitch

def get_old_blades(num_rotors):
    coll_old, rotor_names = [], []
    motion_files = glob.glob(u"motion.r*.txt")
    motion_files.remove("motion.ref.txt")
    for file in motion_files:
        if file != u'motion.ref.txt' and file != u'motion.txt_history':
            coll_old.append(readMotiontxt(file))
            junk = file.replace('motion.','')
            junk = junk.replace('.txt','')
            rotor_names.append(junk)

    return coll_old, rotor_names

def writeRotorInp(num_rotors, trim_meth, dr1_coll, dr2_coll, dr3_coll, dr4_coll, \
                  dr5_coll, dr6_coll, fomoco_target, rdisk_avg, r5_cq, r6_cq, nondim_extra, rtype):
    #############################################################
    # NOW THE TRIMMING_WTAIL.PY SCRIPT IS AMMENDED TO THIS CODE #
    #############################################################
    parser = f90nml.Parser()
    
    if rtype == 'bem':
        # Get rdisk old collectives
        coll_old, rotor_files, rotor_names = get_old_rdisk(num_rotors)
    elif rtype == 'blades':
        coll_old, rotor_names = get_old_blades(num_rotors)
        
    
    # Make a list of the delta collectives for full & half aircraft
    if num_rotors == 3:
        coll_delta = list([dr1_coll, dr2_coll, dr5_coll])
        
        # Open the output file so that I can write to the file.
        output = open("trim_output.log", "a")
        if trim_meth == 1 or trim_meth == 2:
            print >> output, "fomoco_avg(0) [Total Drag] = " + str(fomoco_target[0] + nondim_extra)
            print >> output, "fomoco_avg(4) [Total PM]   = " + str(fomoco_target[4])
            print >> output, "r5_avg [Torque]            = " + str(r5_cq)
            print >> output, "\n"
        elif trim_meth == 3:
            print >> output, "fomoco_target(0) [Total A/C Lift] = " + str(fomoco_target[0])
            print >> output, "fomoco_target(4) [Total A/C PM] = " + str(fomoco_target[4])
            print >> output, "\n"
            
        
        print >> output, 'Old collectives'
        for i in range(num_rotors):
            print >> output, str(rotor_names[i]) + ": " + str(coll_old[i])
        print >> output, "\n"
    
        print >> output, 'Delta collectives'
        for i in range(num_rotors):
            print >> output, str(rotor_names[i]) + ": " + str(coll_delta[i])
        print >> output, "\n"
       
        coll_new = []
        for i in range(num_rotors):
            coll_new.append(coll_old[i] + coll_delta[i])
    
        print >> output, 'New collectives'
        for i in range(num_rotors):
             print >> output, str(rotor_names[i]) + ": " + str(coll_new[i])
        print >> output, "\n"     
          
        if rtype == 'bem':
            # Write new rotor files using f90nml
            os.system('rm *.inp2')
            for i in range(len(rotor_files)):
                rdisk_input = parser.read(rotor_files[i])
                rdisk_input['rdisk_trim']['a0'] = coll_new[i]
                name = rotor_files[i].replace('.inp','.inp2')
                rdisk_input.write(name)
                os.system('cp ' + name + ' ' + str(rotor_files[i]))
                   
            # Copy temporary files for the next iteration
            os.system('cp rotor.r1.inp2 rotor.r1.inp')
            os.system('cp rotor.r2.inp2 rotor.r2.inp')
            os.system('cp rotor.r5.inp2 rotor.r5.inp')
        elif rtype == 'blades':
            os.system("sed 's/coll/{p}/' motion.ref.txt > motion.r1.tmp".format(p=coll_new[0]))
            os.system("sed 's/coll/{p}/' motion.ref.txt > motion.r2.tmp".format(p=coll_new[1]))
            os.system("sed 's/coll/{p}/' motion.ref.txt > motion.r5.tmp".format(p=coll_new[2]))
            os.system("mv motion.r1.tmp motion.r1.txt")
            os.system("mv motion.r2.tmp motion.r2.txt")
            os.system("mv motion.r5.tmp motion.r5.txt")
        output.close()
        
    ###########################################################################
    ################# Full aircraft functionality here ########################
    ###########################################################################
    
    elif num_rotors == 6:
        coll_delta = list([dr1_coll, dr2_coll, dr3_coll, dr4_coll, dr5_coll, dr6_coll])

    
        output = open("trim_output.log", "a")
        r5_avg = rdisk_avg[4,8]
        r6_avg = rdisk_avg[5,8]
        
        if trim_meth == 1 or trim_meth == 2:
            print >> output, "CD_TARGET [Total Drag]   = " + str(fomoco_target[0] + nondim_extra)
            print >> output, "CM_TARGET [Total PM]     = " + str(fomoco_target[4])
            print >> output, "r5_avg [Torque]          = " + str(r5_cq)
            print >> output, "r6_avg [Torque]          = " + str(r6_cq)
            print >> output, "\n"
        elif trim_meth == 3:
            print >> output, "fomoco_target(0) [Total A/C Lift] = " + str(fomoco_target[0])
            print >> output, "fomoco_target(4) [Total A/C PM]   = " + str(fomoco_target[4])
            print >> output, "\n"
    
        print >> output, 'Old collectives'
        for i in range(num_rotors):
            print >> output, str(rotor_names[i]) + ": " + str(coll_old[i])
        print >> output, "\n"
    
        print >> output, 'Delta collectives'
        for i in range(num_rotors):
            print >> output, str(rotor_names[i]) + ": " + str(coll_delta[i])
        print >> output, "\n"
    
        coll_new = []
        for i in range(num_rotors):
            coll_new.append(coll_old[i] + coll_delta[i])
    
        print >> output, 'New collectives'
        for i in range(num_rotors):
             print >> output, str(rotor_names[i]) + ": " + str(coll_new[i])
    
        # Write new rotor files using f90nml 
        os.system('rm *.inp2')
        for i in range(len(rotor_files)):
            rdisk_input = parser.read(rotor_files[i])
            rdisk_input['rdisk_trim']['a0'] = coll_new[i]
            name = rotor_files[i].replace('.inp','.inp2')
            rdisk_input.write(name)
            os.system('cp ' + name + ' ' + str(rotor_files[i]))
              
        # COPY TEMPORARY FILES FOR NEXT INTERATION
        os.system('cp rotor.r1.inp2 rotor.r1.inp')
        os.system('cp rotor.r2.inp2 rotor.r2.inp')
        os.system('cp rotor.r3.inp2 rotor.r3.inp')
        os.system('cp rotor.r4.inp2 rotor.r4.inp')
        os.system('cp rotor.r5.inp2 rotor.r5.inp')
        os.system('cp rotor.r6.inp2 rotor.r6.inp')
    
        output.close()

    return coll_old, coll_new, coll_delta

def get_old_tail(tclfile, tail_line):
    # Read the old tail deflection angle
    tail_file = open(tclfile, "r")
    junk = tail_file.readlines()
    tmp = []
    for lines in junk:
        x = lines.split(' ')
        tmp.append(x)
    
    tail_last = float(tmp[tail_line-1][-1])
    return tail_last

def make_tcl_script (new_tail, tclfile, symmetric):
    # Write the new grid generation file
    refname = tclfile + '_ref'
    os.system('cp {ref} {tcl}'.format(ref=refname, tcl=tclfile))
    os.system("sed -i s/TAIL_ANGLE/{tail}/ {tcl}".format(tail=new_tail, tcl=tclfile))
    # Is the grid half a/c or full? 1 = full ac, 0 = half ac
    os.system("sed -i s/SYMMETRIC_REF/{sym}/ {tcl}".format(sym=symmetric, tcl=tclfile))
    return 

def makeTail(num_rotors, dtail, gridgen_path, tclfile, tail_line):
    present_dir = os.getcwd()
    
    # Am I making a symmetric grid or not?
    if num_rotors == 3:
        symmetric = 0
    elif num_rotors == 6:
        symmetric = 1
    
    tail_last = get_old_tail(tclfile, tail_line)

    # This is the same for both sides of the ruddervator
    os.chdir(present_dir)
    
    # Calculate the new tail deflection angle with Newton iteration
    new_tail = dtail + tail_last
    
    # Write results to file for debugging
    output = open("trim_output.log", "a")
    print >> output, "Tail Old  : " + str(tail_last)
    print >> output, "Tail Delta: " + str(dtail)
    print >> output, "Tail New  : " + str(new_tail)
    output.close()
    
    # Make the new TCL script by changing reference values
    make_tcl_script(new_tail, tclfile, symmetric)
    
    # Copy that file into the grid generation directory
    os.system("cp {tcl} {dir}/.".format(tcl=tclfile, dir=gridgen_path))
    
    # Go into directory with all grids and make the new ruddervator grid
    os.chdir(gridgen_path)
    os.system('./{tcl}'.format(tcl=tclfile))
    os.chdir(present_dir)
    
    # Export the gridgen path as a BASH variable
    os.putenv("GRIDGEN", gridgen_path)
    
    # This script with replace the new tail in the grid.in 
    os.system('./replace_tail_grid {halffull}'.format(halffull=symmetric))
            
    return tail_last, dtail, new_tail
    
def get_output_file(outname, basename):
    dout = outname + '.out'
    files = glob.glob("*{out}*".format(out=outname))
    if (basename + '.{out}'.format(out=outname)) in files:
        filename = basename + '.' + outname
    elif dout in files:
        filename = dout
    return filename

def casenum(fomoco_file):
    rst = 0
    # Now I'm just determining if this is the beginning of the loop or a restart
    for root, dirs, files in os.walk("."):
        if 'fomoco.out' in files:
            rst = 1
        elif fomoco_file in files:
            rst = 1
    
    # If not a restart, start with case=1 
    if rst == 0:
        case = 1
    elif rst == 1:
        case = 2
    
    return case
    
def get_basename():
    files = glob.glob("*.*.out")
    files.sort(key=os.path.getmtime)
    basefile = files[-1]
    basename = basefile.split('.')[0]
    return basename


def read_hist(hist_files, avg_over):
    size = (len(hist_files), 6)
    rotor_avg = np.zeros(size)
    for f in range(len(hist_files)):
        data_tmp = np.loadtxt(hist_files[f])
        rotor = data_tmp[-avg_over:, 1:]
        for i in range(6):
            rotor_avg[f, i] = np.mean(rotor[:, i])
    return rotor_avg


def change_nondim_thrust(rotor_avg, radius, aref_target, lref_target, rho, vref, rtype):     
    # Convert to the aircraft non-dimensional variables
    Aref_helo = np.pi * radius**2
    Aref_targ = aref_target*2
    
    fdim_helo = rho * Aref_helo * vref**2
    mdim_helo = rho * Aref_helo * vref**2 * radius
    fdim_targ = 0.5 * rho * Aref_targ * vref**2
    mdim_targ = 0.5 * rho * Aref_targ * vref**2 * lref_target
    
    fdim_ratio = fdim_targ/fdim_helo
    mdim_ratio = mdim_targ/mdim_helo
    
    if rtype == 'bem':
        col_thrust = 7
    else:
        col_thrust = 3
    col_thrust -= 1 # change to python indexing
    rotor_avg[:, col_thrust] = rotor_avg[:, col_thrust] * fdim_ratio
#    rotor_avg[:, 3:] = rotor_avg[:, 3:] * mdim_ratio
    
    newdata = rotor_avg
    return newdata


def rotor_fomoco(fomoco_avg):
    rotor_range = range(5, 8)
    rdisk_forces = fomoco_avg[rotor_range, :]
    return rdisk_forces
    
def write_plot_out (iteration, drag, cmy, r5_cq, coll_old, coll_delta, coll_new, tail_old, tail_delta, tail_new):
    iteration += 1
    fileout = 'TrimData.txt'
    if len(coll_old) == 3:
        names = [1, 2, 5]
    else:
        names = range(1,7)
    for i in range(len(coll_old)):
        output = [iteration, drag, cmy, r5_cq, names[i], coll_old[i], coll_delta[i], coll_new[i], tail_old, tail_delta, tail_new]
        with open(fileout, 'a') as file:
            for val in output:
                file.write(str(val))
                file.write(',')
            file.write('\n')
            file.close()
    return
    
    
    
