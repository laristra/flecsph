#!/usr/bin/env python
 
# Author : Alexander Kaltenborn
# Date : 06/04/2019
#

# package loading
import os, sys, h5py, argparse
import numpy as np

my_description = """
Extracts single step from an h5part multistep output file:
  large.h5part   ==>   large.h5part<STEP#####>.h5part
"""
my_usage = ("   %(prog)s <infile> [-o|--outfile <filename>]" +
            " [-s|--step <val> | -l|--last] [-h|--help]")

parser = argparse.ArgumentParser(description=my_description, 
         usage=my_usage, formatter_class=argparse.RawTextHelpFormatter,
         epilog="example: ./%(prog)s -f sodtube_evolution.h5part -s 10")
parser.add_argument('infile', type=str,
      help='input file in the FleCSPH hdf5 format')
parser.add_argument("-o", "--outfile", action="store", type=str, default="",
      help="output file in the FleCSPH hdf5 format", dest="outfile")

group = parser.add_mutually_exclusive_group()
group.add_argument("-s", "--step", action="store", type=int, default=0, 
      help="step number (default: 0)", dest="step")
group.add_argument("-l", "--last", action="store_true", default=False, 
      dest="last", help="selects the last step of the file")
args = parser.parse_args()

# read the input file
try:
  h5file = h5py.File(args.infile,'r')
except:
  sys.exit ("ERROR: cannot read input file %s" % args.infile)

# by default, read the first timestep
if(args.last == True):
  iterhold = len(list(h5file.keys())) - 1
else:
  iterhold = args.step
key_step = ("Step#%d" % iterhold)
if not key_step in h5file:
  sys.exit ("ERROR: cannot find step '%s' in your file" % key_step)

# check if the output file already exists; 
# prompt the user if they want to overwrite it
outfile_name = args.outfile
if (outfile_name == ""):
  namesplit = args.infile.split(".h5part")
  outfile_name = "%s_%05d.h5part" % (namesplit[0], iterhold)

if os.path.isfile(outfile_name):
  ans = None
  while ans not in ("Y", "y", "N", "n"):
    ans = input ("File %s exists. Overwrite? [y/N] " % outfile_name)
    if ans == "n" or ans == "N":
      exit(0)

# open output h5part file
outfile = h5py.File(outfile_name,'w')
out_step = 0

grp = outfile.create_group("/Step#"+str(out_step))
dset   = h5file[key_step+"/P"]
grp.create_dataset("P",data=dset)
dsetX  = h5file[key_step+"/x"]
grp.create_dataset("x",data=dsetX)
dsetY  = h5file[key_step+"/y"]
grp.create_dataset("y",data=dsetY)
dsetZ  = h5file[key_step+"/z"]
grp.create_dataset("z",data=dsetZ)
dset   = h5file[key_step+"/vx"]
grp.create_dataset("vx",data=dset)
dset   = h5file[key_step+"/vy"]
grp.create_dataset("vy",data=dset)
dset   = h5file[key_step+"/vz"]
grp.create_dataset("vz",data=dset)
dset   = h5file[key_step+"/ax"]
grp.create_dataset("ax",data=dset)
dset   = h5file[key_step+"/ay"]
grp.create_dataset("ay",data=dset)
dset   = h5file[key_step+"/az"]
grp.create_dataset("az",data=dset)
dsetD  = h5file[key_step+"/rho"]
grp.create_dataset("rho",data=dsetD)
dsetM  = h5file[key_step+"/m"]
grp.create_dataset("m",data=dsetM)
dsetH  = h5file[key_step+"/h"]
grp.create_dataset("h",data=dsetH)
dsetU  = h5file[key_step+"/u"]
grp.create_dataset("u",data=dsetU)
dsetT  = h5file[key_step+"/type"]
grp.create_dataset("type",data=dsetT)
dsetID = h5file[key_step+"/id"]
grp.create_dataset("id",data=dsetID)

outfile.close()
h5file.close()

# report the name of the output file
print("Extracted step %05d to output file: %s" % (iterhold, outfile_name))
print("bye bye")
