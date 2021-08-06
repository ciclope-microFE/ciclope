#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Process txt file and produce force-displacement plot.
More info: https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts
"""

import pylab
import numpy
import argparse

parser = argparse.ArgumentParser(description='Process txt file and produce force-displacement plot.')
parser.add_argument('-i', '--filein', type=str, help='input .TXT filename produced from dat2txt.py')
parser.add_argument('-o', '--fileout', type=str, default=None, help='output plot filename')

args = parser.parse_args()

# load file
print('Loading file: {}'.format(args.filein))

# de = numpy.genfromtxt("total internal energy_EDRAHT.txt")
dm = numpy.genfromtxt(args.filein)
pylab.plot(dm[:,0],dm[:,3],'r')
pylab.grid(True)
pylab.xlim([0,1])
pylab.xlabel("t")
pylab.ylabel("F_z")
pylab.legend([args.filein],loc=0)
pylab.savefig(args.fileout)