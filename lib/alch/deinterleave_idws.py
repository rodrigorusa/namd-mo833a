#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Input: fepout file from an interleaved double-wide sampling run from 0 to 1
# (with the last window from 1 to the previous lambda value)
# Output:
# - <base>_fwd.fepout (forward deltaE only)
# - <base>_bwd.fepout (backward deltaE only)
# - <base>_bwd_r.fepout (backward data, with order of windows reversed to simulate a
#   reverse sequential transformation (0 -> 1 -> 0)

# Notes on further analysis:
# Missing samples of deltaE must be interpolated to recover the same number of samples in each window
# This is necessary for ParseFEP up to version 2.1
# Hopefully later versions will fix this

# Note:
# The first window samples forward only, the last window backward only
#
# Jérôme Hénin <henin@ibpc.fr> (2018)
# with contributions from Tom Joseph (UPenn)

import sys
import os.path
import argparse
import re
import math
import fileinput

try:
    from scipy.interpolate import InterpolatedUnivariateSpline
    interp_avail = True
except ImportError:
    interp_avail = False

def interpolated_window(steps_window, deltaE_window, lambda1, lambda2, lines=[]):
    equil = True
    acc_exp = 0
    acc_N = 0
    buf = ''
    # Create interpolation function for deltaE data
    func_dE = InterpolatedUnivariateSpline(steps_window, deltaE_window)
    for t in range(len(steps)):
        s = steps[t]
        if equil and s >= equilSteps:
            buf += ('#{} STEPS OF EQUILIBRATION AT LAMBDA {} COMPLETED\n#STARTING COLLECTION OF ENSEMBLE AVERAGE\n').format(equilSteps, lambda1)
            equil = False
            acc_exp = 0
            acc_N = 0
        dE = float(func_dE(s))
        acc_exp += math.exp(-beta * dE)
        if acc_exp == 0: # catch floating-point underflow due to very high dE (clash)
            acc_exp = 1e-16
        acc_N += 1
        buf += 'FepEnergy: {:6d}  {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(s, 0, 0, 0, 0, dE, 0., T[t], -RT*math.log(acc_exp / acc_N))
    return buf, -RT*math.log(acc_exp / acc_N)

# Note: this doesn't subsample anymore, as I've given up the subsampling trick, which is frustratingly complicated to get right
# and wouldn't be necessary if ParseFEP was slightly less ****
# which hopefully it will be soon
def subsampled_window(steps_window, deltaE_window, lambda1, lambda2, lines):
    equil = True
    acc_exp = 0
    acc_N = 0
    buf = ''
    first = steps_window[0]

    stride = 1
    # Disabling the trick blow, as this does not reliably produce even-sampled windows
    # and then parseFEP has a fit
    # find first gap in data
    #for i in range(len(steps_window)-1):
    #    diff = steps_window[i+1] - steps_window[i]
    #    if diff != alchOutFreq:
    #        if diff % alchOutFreq != 0:
    #            sys.exit('Error: gap between data points is ' + str(diff) + ' instead of a multiple of ' + str(alchOutFreq))
    #        else:
    #            stride = diff // alchOutFreq - 1
    #        break

    for t in range(0, len(steps_window), stride):
        s = steps_window[t]
        if equil and s >= equilSteps:
            buf += ('#{} STEPS OF EQUILIBRATION AT LAMBDA {} COMPLETED\n#STARTING COLLECTION OF ENSEMBLE AVERAGE\n').format(equilSteps, lambda1)
            equil = False
            acc_exp = 0
            acc_N = 0
        dE = deltaE_window[t]
        acc_exp += math.exp(-beta * dE)
        if acc_exp == 0: # catch floating-point underflow due to very high dE (clash)
            acc_exp = 1e-16
        acc_N += 1
        buf += 'FepEnergy:' + lines[t]
        #buf += 'FepEnergy: {:6d}  {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(s, 0, 0, 0, 0, dE, 0., T[t], -RT*math.log(acc_exp / acc_N))
    return buf, -RT*math.log(acc_exp / acc_N)


def process_window():
    # Create interpolated function from data and write new data for the whole window
 
    if not last_window and len(steps_f) > 0:
        print("Forward data: " + str(lambda1) + " to " + str(lambda2) + " (" + str(len(steps_f)) + " points)")
        outf.write("#NEW FEP WINDOW: LAMBDA SET TO {} LAMBDA2 {}\n".format(lambda1, lambda2))
        buf, dG = resample(steps_f, deltaE_f, lambda1, lambda2, lines_f)
        outf.write(buf)
        if line.startswith('#Free'):
            outf.write(line)
        else:
            outf.write('#Free energy change for lambda window [ {} {} ] is {} ; net change until now is 0.\n'.format(lambda1, lambda2, dG))

    if len(steps_b) > 0: # only if we do have backward data in this window
        # create new window in the backward output
        # for every window but the first one, which is forward-only
        print("Backward data: " + str(lambda1) + " to " + str(lambdaIDWS) + " (" + str(len(steps_b)) + " points)")
        outb_buffer.append("#NEW FEP WINDOW: LAMBDA SET TO {} LAMBDA2 {}\n".format(lambda1, lambdaIDWS))
        buf, dG = resample(steps_b, deltaE_b, lambda1, lambdaIDWS, lines_b)
        outb_buffer[-1] += buf
        outb_buffer[-1] += '#Free energy change for lambda window [ {} {} ] is {} ; net change until now is 0.\n'.format(lambda1, lambdaIDWS, dG)


parser = argparse.ArgumentParser(description='Process fepout data from Interleaved Double-Wide Sampling')
parser.add_argument('filenames', nargs='+', help='fepout file name(s) with lambda from 0 to 1')
parser.add_argument('-T', '--temperature', type=float, default=300., help='temperature for FEP estimate, in K')
parser.add_argument('-r', '--reverse', type=bool, default=True, help='reverse order of backward-sampling windows')
parser.add_argument('-i', '--interpolate', default=True, action='store_true', help='interpolate data rather than subsampling it')
args = parser.parse_args()

# Store the output in the same directory as the first file
basename = os.path.splitext(args.filenames[0])[0]
outf = open(basename + '_fwd.fepout', 'w')

interpolating = False
if args.interpolate:
    if interp_avail:
        resample = interpolated_window
        interpolating = True
    else:
        print('Interpolation requires the scipy.interpolate package. Defaulting to subsampling.')
        resample = subsampled_window
else:
    resample = subsampled_window

# Buffer to hold the contents of the backwards fepout file, by window. Will be written in reverse.
outb_buffer = []

first_window = True
last_window = False

RT = 1.98720650096e-3 * args.temperature #in kcal/mol
beta = 1. / RT

# Global step counter, because we don't depend on the fepout files to have consistent step numbers
step_counter = -1
alchOutFreq = -1
equilSteps = 0

deltaE_f = []
steps_f = []
lines_f = []
deltaE_b = []
steps_b = []
lines_b = []
steps = []
T = []

equil_re = re.compile('#(\d*) STEPS OF EQUILIBRATION')

for line in fileinput.input(args.filenames):
    fields = re.split(' +', line.strip())
    if line.startswith('#NEW FEP WINDOW:'):
        # Process previous window, if any
        process_window()
        # data will contain all dE data for the window
        # make two copies, forward and backward
        deltaE_f.clear()
        steps_f.clear()
        lines_f.clear()
        deltaE_b.clear()
        steps_b.clear()
        lines_b.clear()
        # Common data: step number, temperature
        steps.clear()
        T.clear()
        step_counter = -1
        first_step = -1
        alchOutFreq = -1
        equilSteps = 0

        # Extract lambda values
        lambda1 = float(fields[6])
        lambda2 = float(fields[8])
        if len(fields) == 11:
            lambdaIDWS = float(fields[10])
            IDWS_on = True
        else:
            lambdaIDWS = -1.
            IDWS_on = False

        # Check that we start from an end-point
        if first_window:
            if lambda1 == 1:
                start_from_one = True
            elif lambda1 == 0:
                start_from_one = False
            else:
                sys.exit('Lambda should start at zero or one, found ' + str(lambda1))

        # Skip the last window for the "forward" output, as it is backward from 1
        if ((not start_from_one and lambda1 == 1) or (start_from_one and lambda1 == 0)):
            last_window = True

        if (last_window):
            # special case, last window is backward
            lambdaIDWS = lambda2
            lambda2 = -1

        # Done processing the header of the first window
        first_window = False

    match = equil_re.match(line)
    if match != None:
        equilSteps = int(match.group(1))

    # Collect all timestep numbers (both fwd and bwd) for interpolation
    if line.startswith('FepE'):
        l = len(fields)
        if l == 8:
            i_temp = 7
        elif l == 10:
            i_temp = 8
        else:
            print("Wrong number of fields (expected 8 or 10) in line:\n" + line)

        # Detect first step and value of alchOutFreq
        if step_counter == -1:
            first_step = int(fields[1])
            step_counter = first_step
        elif step_counter == first_step:
            alchOutFreq = int(fields[1]) - first_step
            step_counter += alchOutFreq

        steps.append(step_counter)
        T.append(float(fields[i_temp]))

        if line.startswith('FepEnergy:') and not last_window:
            steps_f.append(step_counter)
            deltaE_f.append(float(fields[6]))
            if not interpolating:
                lines_f.append(line[10:])

        if line.startswith('FepE_back:') or (last_window and line[:10] == 'FepEnergy:'):
            steps_b.append(step_counter)
            deltaE_b.append(float(fields[6]))
            if not interpolating:
                lines_b.append(line[10:])

        if alchOutFreq > 0:
            step_counter += alchOutFreq

# Process last window
process_window()


outf.close()

outb = open(basename + '_bwd.fepout', 'w')
if args.reverse:
    # Write backward windows in reverse order
    for string in reversed(outb_buffer):
        outb.write(string)
else:
    # Write output in same order as run
    for string in outb_buffer:
        outb.write(string)
outb.close()

