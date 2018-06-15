#!/usr/bin/python3

# Input: fepout file from an interleaved double-wide sampling run from 0 to 1
# (with the last window from 1 to the previous lambda value)
# Output:
# - <base>_fwd.fepout (forward deltaE only)
# - <base>_bwd.fepout (backward deltaE only)
# - <base>_bwd_r.fepout (backward data, with order of windows reversed to simulate a
#   reverse sequential transformation (0 -> 1 -> 0)
# Process fwd.fepout and bwd_r.fepout with ParseFEP for SOS and BAR estimators
# Missing samples of deltaE are interpolated
#
# Note:
# The first window samples forward only, the last window backward only
#
# Jérôme Hénin <henin@ibpc.fr> (2018)

import sys
import os.path
import argparse
import re
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline


def interpolated_window(steps_window, deltaE_window, lambda1, lambda2):
  equil = True
  acc_exp = 0
  acc_N = 0
  buf = ''
  # Create interpolation function for deltaE data
  func_dE = InterpolatedUnivariateSpline(steps_window, deltaE_window)
  for t in range(len(steps)):
    s = steps[t]
    if equil and s >= equilSteps:
      buf += '#' + str(s) + ' STEPS OF EQUILIBRATION AT LAMBDA ' + str(0) + ' COMPLETED\n#STARTING COLLECTION OF ENSEMBLE AVERAGE\n'
      equil = False
      acc_exp = 0
      acc_N = 0
    dE = float(func_dE(s))
    acc_exp += math.exp(-beta * dE)
    acc_N += 1
    buf += 'FepEnergy: {:6d}  {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(s, 0, 0, 0, 0, dE, 0., T[t], -RT*math.log(acc_exp / acc_N))
  return buf, -RT*math.log(acc_exp / acc_N)


parser = argparse.ArgumentParser(description='Process fepout data from Interleaved Double-Wide Sampling')
parser.add_argument('filename', help='fepout file name with lambda from 0 to 1')
parser.add_argument('-T', type=float, default=300., help='temperature for FEP estimate')
parser.add_argument('-r', type=bool, default=True, help='reverse order of backward-sampling windows')
args = parser.parse_args()

basename = os.path.splitext(args.filename)[0]
outf = open(basename + '_fwd.fepout', 'w')

# Only interpolate deltaE data
quick = True

# Buffer to hold the contents of the backwards fepout file, by window. Will be written in reverse.
outb_buffer = []

first_window = True
last_window = False

l_l1 = []
l_l2 = []
l_l3 = []

RT = 1.98720650096e-3 * args.T # in kcal/mol
beta = 1. / RT

with open(args.filename, 'r') as fepout:
  for line in fepout:
    fields = re.split(' +', line.strip())
    if line[:16] == '#NEW FEP WINDOW:':

      # data will contain all dE data for the window
      # make two copies, forward and backward
      deltaE_f = []
      steps_f = []
      deltaE_b = []
      steps_b = []
      # Common data: step number, temperature
      steps = []
      T = []

      lambda1 = float(fields[6])
      lambda2 = float(fields[8])
      if len(fields) == 11:
        lambdaIDWS = float(fields[10])
      else:
        lambdaIDWS = -1.
      print("Lambda window: " + str(lambda1) + " to " + str(lambda2) + (" IDWS: " + str(lambdaIDWS) if lambdaIDWS >= 0 else ""))

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
      if (not last_window):
        outf.write("#NEW FEP WINDOW: LAMBDA SET TO {} LAMBDA2 {}\n".format(lambda1, lambda2))
      else:
        # special case, last window is backward
        lambdaIDWS = lambda2
        lambda2 = -1

      if not first_window:
        # create new window in the backward output
        # for every window but the first one, which is forward-only
        outb_buffer.append("#NEW FEP WINDOW: LAMBDA SET TO {} LAMBDA2 {}\n".format(lambda1, lambdaIDWS))

      # Done processing the header of the first window
      first_window = False

    if line.strip() == '#STARTING COLLECTION OF ENSEMBLE AVERAGE':
      equilSteps = steps[-1] + 1 # actually a lower bound of equilSteps if fepOutFreq > 1

    if line[:10] == 'FepEnergy:' and not last_window:
      steps_f.append(int(fields[1]))
      deltaE_f.append(float(fields[6]))

    if line[:10] == 'FepE_back:' or (last_window and line[:10] == 'FepEnergy:'):
      steps_b.append(int(fields[1]))
      deltaE_b.append(float(fields[6]))

    # Collect all timestep numbers (both fwd and bwd) for interpolation
    if line[:4] == 'FepE':
      steps.append(int(fields[1]))
      T.append(float(fields[8]))

    if line[:19] == '#Free energy change':
      # End of window
      # Create interpolated function from data and write new data for the whole window

      if not last_window:
        buf, dG = interpolated_window(steps_f, deltaE_f, lambda1, lambda2)
        outf.write(buf)
        outf.write(line)

      if len(steps_b) > 0: # only if we do have backward data in this window
        buf, dG = interpolated_window(steps_b, deltaE_b, lambda1, lambdaIDWS)
        outb_buffer[-1] += buf
        outb_buffer[-1] += '#Free energy change for lambda window [ {} {} ] is {} ; net change until now is 0.\n'.format(lambda1, lambda2, dG)


outf.close()

if args.r:
  # Write backward windows in reverse order
  outb_reversed = open(basename + '_bwd.fepout', 'w')
  for string in reversed(outb_buffer):
    outb_reversed.write(string)
  outb_reversed.close()
else:
  # Write output in same order as run
  outb = open(basename + '_bwd.fepout', 'w')
  for string in outb_buffer:
    outb.write(string)
  outb.close()

