#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# simple.py
# Purpose: mimic operation of couple/simple/simple.cpp via Python
# Syntax:  simple.py in.liggghts
#          in.liggghts = LIGGGHTS input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print "Syntax: simple.py in.liggghts"
  sys.exit()

infile = sys.argv[1]

me = 0

# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from liggghts import liggghts
lmp = liggghts()

# run infile one line at a time

lines = open(infile,'r').readlines()
for line in lines: lmp.command(line)

# run 10 more steps
# get coords from LIGGGHTS
# change coords of 1st atom
# put coords back into LIGGGHTS
# run a single step with changed coords

lmp.command("run 10")
x = lmp.gather_atoms("x",1,3)
epsilon = 0.1
x[0] += epsilon
lmp.scatter_atoms("x",1,3,x)
lmp.command("run 1");

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
