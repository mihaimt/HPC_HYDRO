#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os

import argparse



parser = argparse.ArgumentParser(description='Convert the output of a sim (vts files) to a movie via pngs')
parser.add_argument('-i', '--input', type=str, required=True,
                    metavar='vts_files/',
                    help='a relative path to the output folder of the simulation')
parser.add_argument('-o', '--output', type=str,required=False, default='movie.mpg',
                    metavar='movie.mpg',
                    help='the relative path to the output movie')
parser.add_argument('-n', '--nsteps', type=int, required=False, default=-1,
                    help='the number of timesteps to actually use [default = -1 = all]')
parser.add_argument('-d', '--delay', type=int, required=False, default=100,
                    help='the delay between frames in ms')


args = parser.parse_args()
print os.path.abspath(args.input)

try:
    args.input = os.path.abspath(args.input)
    if not os.path.isdir(args.input):
        raise KeyError
except:
    parser.error("--input is not a valid path")

