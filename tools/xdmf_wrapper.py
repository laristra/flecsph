#!/usr/bin/env python

# Copyright (c) 2019 Triad National Security, LLC
# All rights reserved

import re, os, sys, h5py, argparse
import numpy as np
from math import *

# ===========================================================================
#
# ## FUNCTION TO OUTPUT TIMESTEP TO XDMF FILE
#    Params:
#    - in_h5file: HDF5 input file
#    - out_xmdfile: output XMDF file
#    - Npart: number of particles
#    - key_step: timestep
#
def write_xdmf_timestep (in_h5file, out_xdmfile, in_fname, key_step, ndim):
  # -- Dictionary of variables  -----
  vps = {'x'    : {'type':'Float','size':8},
         'y'    : {'type':'Float','size':8},
         'z'    : {'type':'Float','size':8},
         'vx'   : {'type':'Float','size':8},
         'vy'   : {'type':'Float','size':8},
         'vz'   : {'type':'Float','size':8},
         'ax'   : {'type':'Float','size':8},
         'ay'   : {'type':'Float','size':8},
         'az'   : {'type':'Float','size':8},
         'rho'  : {'type':'Float','size':8},
         'P'    : {'type':'Float','size':8},
         'u'    : {'type':'Float','size':8},
         'dt'   : {'type':'Float','size':8},
         'h'    : {'type':'Float','size':8},
         'm'    : {'type':'Float','size':8},
         'id'   : {'type':'Integer','size':8},
         'key'  : {'type':'Integer','size':8},
         'rank' : {'type':'Integer','size':8},
         'type' : {'type':'Integer','size':4}}
  # open the <Grid> section in the file
  out_xdmfile.write("""
    <Grid>""")
  dset = in_h5file[key_step]
  Npart = dset['x'].len()

  # -- 1D --------- ---------------------
  if ndim == 1:
    out_xdmfile.write("""
      <Topology TopologyType="Polyvertex" NumberOfElements="%d"/>
      <Geometry Type="X_Y">
         <DataItem Name="x" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/x
         </DataItem>
         <DataItem Name="y" Format="XML" 
          Dimensions="%d" NumberType="%s" Precision="%d">
                            """ % (Npart,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step,
       Npart,vps['x']['type'],vps['x']['size']))
    for i in range(Npart): 
      out_xdmfile.write("0.0 ")
      if i % 20 == 0: out_xdmfile.write("\n")
    out_xdmfile.write("""
         </DataItem>
      </Geometry>""")
    for var in ['x','vx','ax']:
      if var in dset.keys(): out_xdmfile.write("""
      <Attribute Name="%s">
        <DataItem Format="HDF" Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/%s
        </DataItem>
      </Attribute>""" % (var,
       Npart,vps[var]['type'],vps[var]['size'],
       in_fname,key_step,var)) 

  # -- 2D --------- ---------------------
  if ndim == 2:
    out_xdmfile.write("""
      <Topology TopologyType="Polyvertex" NumberOfElements="%d"/>
      <Geometry Type="X_Y">
         <DataItem Name="x" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/x
         </DataItem>
         <DataItem Name="y" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/y
         </DataItem>
      </Geometry>""" % (Npart,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step))

    for var in ['x','y','vx','vy','ax','ay']:
      if var in dset.keys(): out_xdmfile.write("""
      <Attribute Name="%s">
        <DataItem Format="HDF" Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/%s
        </DataItem>
      </Attribute>""" % (var,
       Npart,vps[var]['type'],vps[var]['size'],
       in_fname,key_step,var))

  # -- 3D --------- ---------------------
  if ndim == 3:
    out_xdmfile.write("""
      <Topology TopologyType="Polyvertex" NumberOfElements="%d"/>
      <Geometry Type="X_Y_Z">
         <DataItem Name="x" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/x
         </DataItem>
         <DataItem Name="y" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/y
         </DataItem>
         <DataItem Name="z" Format="HDF" 
          Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/z
         </DataItem>
      </Geometry>""" % (Npart,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step,
       Npart,vps['x']['type'],vps['x']['size'],
       in_fname,key_step))
    for var in ['x','y','z','vx','vy','vz','ax','ay','az']:
      if var in dset.keys(): out_xdmfile.write("""
      <Attribute Name="%s">
        <DataItem Format="HDF" Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/%s
        </DataItem>
      </Attribute>""" % (var,
       Npart,vps[var]['type'],vps[var]['size'],
       in_fname,key_step,var))


  # -- Other scalar variables   --------
  for var in ['rho','h','m','P','u','dt','id','key','rank','type']:
    if var in dset.keys(): out_xdmfile.write("""
      <Attribute Name="%s">
        <DataItem Format="HDF" Dimensions="%d" NumberType="%s" Precision="%d">
            %s:/%s/%s
        </DataItem>
      </Attribute>""" % (var,
       Npart,vps[var]['type'],vps[var]['size'],
       in_fname,key_step,var))

  # -- Output timestamp    --------------
  try:
    dset = in_h5file[key_step]
    tm = dset.attrs['time'][0]
  except:
    tm = 0.0
  out_xdmfile.write("""
      <Time Value="%e" />""" % tm)

  # -- Footer ---------------------------
  xdmfile.write("""
  </Grid>""")


# 
# ## FUNCTION TO SORT STEP LABELS
#    because alphabetical sort won't work, for example
#    'Step2 > Step100' # incorrect!
#
def stepLabelSort (stepLabel):
  # Step#23 -> 23
  return int(stepLabel[5:])

########################
my_desc = """
Creates an XDMF wrapper for one or multiple H5Part files"""
my_usag = """
   %(prog)s inprefix[.h5part] [-h] [-o <fname.xdmf>]
"""
parser = argparse.ArgumentParser ( description=my_desc, usage=my_usag,
    epilog="example: ./%(prog)s sodtube_evolution.h5part -s 10 -v2 u")

parser.add_argument('inprefix', type=str,
    help='input file prefix (or full path) to the FleCSPH hdf5 output')
parser.add_argument('-o', type=str,
    help='output file prefix (by default, same as input)')
args = parser.parse_args()

# input prefix: strip file extension and 5-digit suffix, if any:
infile_prefix = os.path.basename(args.inprefix)
dirname = os.path.dirname(args.inprefix)
if dirname == '': dirname = '.'
else:    os.chdir(dirname)
if infile_prefix[-7:] == '.h5part': infile_prefix = infile_prefix[:-7]
m = re.search('^(.*)_\d\d*$',infile_prefix)
ifname = infile_prefix + ".h5part"

# -- Check if file(s) exist   -----------------------------------------------
#
# m != None?
# |
# +-- yes: user wants multiple files - check for them; 
# |        +-- found: proceed with multiple files
# |        +-- none exists: complain and exit
# |
# +-- no:  check for single file
#          +-- found: proceed with single file
#          +-- not found: check for multiple files
#              +-- not found: complain and exit
#              +-- found: proceed with multiple files
# 
multiple_files_mode = False
if m != None: 
  if not os.path.isfile(ifname):
    sys.exit("ERROR: cannot find input file: " + ifname)
  multiple_files_mode = True
  infile_prefix = m.groups()[0]

else:
  if not os.path.isfile(ifname):
    # check for multiple files
    for f in os.listdir('.'):
      if re.search(infile_prefix+'_\d\d*.h5part',f) != None:
        multiple_files_mode = True
        break
    if not multiple_files_mode:
      sys.exit("ERROR: cannot find input file(s): " + infile_prefix)


# -- Compile and sort list of h5part snapshots   ----------------------------
mflist = []
if multiple_files_mode:
  for f in os.listdir('.'):
    if re.search(infile_prefix+'_\d\d*.h5part',f) != None:
      mflist.append(f)
  mflist.sort()
  ifname = mflist[0]

# -- Open the file and determine the dimension  -----------------------------
try:
  h5file = h5py.File(ifname)
except:
  sys.exit ("ERROR: cannot read input file %s" % ifname)

# determine the dimension
try:
  ndim = h5file.attrs['ndim'][0]
except:
  try:
    ndim = h5file.attrs['dimension'][0]
  except:
    sys.exit ("ERROR: cannot determine file dimensionality")


# -- Check if output file exists; if so, delete it    -----------------------
xdmfname = infile_prefix + ".xdmf"
if os.path.isfile(xdmfname):
  os.remove(xdmfname)


# -- Write header -----------------------------------------------------------
xdmfile = open(xdmfname,'w')
xdmfile.write("""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain><Grid GridType="Collection" CollectionType="Temporal" >""")

# -- Output snapshots in either multiple-files, or single-file mode  --------
if multiple_files_mode:
  for ifname in mflist:
    h5file.close()
    try:
      h5file = h5py.File(ifname,'r')
    except:
      sys.exit ("ERROR: cannot read input file %s" % ifname)
    steps = h5file.keys()
    write_xdmf_timestep (h5file, xdmfile, ifname, steps[0], ndim)

else:
  steps = h5file.keys()
  steps.sort(key=stepLabelSort)
  for key_step in steps:
    write_xdmf_timestep (h5file, xdmfile, ifname, key_step, ndim)

# -- Write footer, close files  ---------------------------------------------
xdmfile.write("""
 </Grid></Domain>
</Xdmf>
""")
xdmfile.close()
h5file.close()


