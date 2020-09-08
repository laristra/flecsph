#!/usr/bin/env python3
"""
Takes a(/two) single star h5part file(/s) and produces binary system at a specified orbital separation:
"""

#  IMPORTANT if there are asymmetries in initial files, orient this way:
#            M1 ->            COM             <- M2
#           ***                                  ***
#        **********                          **********
#       *************          x           *************
#        **********                          **********
#           ***                                  ***
#

# Store data in h5part format
import h5py, sys, argparse, os
import numpy as np

# IMPORTANT CONSTANTS
G_newt  = 6.67430e-8 #cm3/g/s2
M_solar = 1.988435e33 #g
R_solar = 6.957e10 #cm
c_light = 2.99792458e10 #cm/s

ln = np.log

my_description = """
Takes a(/two) single star h5part file(/s) and produces binary system at a specified orbital separation:
  star.h5part (left-justified) [+ star2.h5part (right-justified)]  ==>   binary.h5part"""
my_usage = """
    %(prog)s [-f|--file <filename(s)>] [-i|--identical] [-pm|--pointmass <val>] [-a|--orbsep <val>] [-h|--help]"""

parser = argparse.ArgumentParser(description=my_description, usage=my_usage, epilog="EXAMPLE: $python %(prog)s -f wd.h5part -i -a 1.e9")

parser.add_argument("-f",   "--file",          action="store",      type=argparse.FileType('r'),                 help="input file(s) in FleCSPH hdf5 format",         nargs="+",  dest="infile")
parser.add_argument("-i",   "--identical",     action="store_true",                              default=False,  help="set this if you want an equal mass binary",                dest="ident")
parser.add_argument("-pm",  "--pointmass",     action="store",      type=float,                  default=-1.,    help="total mass of the point star (default: -1.)",              dest="pm")
parser.add_argument("-a",   "--orbsep",        action="store",      type=float,                  default=-1.,    help="orbital separation in cm.",                                dest="orbsep")
parser.add_argument("-pmd", "--pointmassrho",  action="store",      type=float,                  default=4.e+14, help="density of the point mass (default: 4e14)",                dest="pmd")
parser.add_argument("-pmh", "--pointmassh",    action="store",      type=float,                  default=1.,     help="smoothing len of the point mass (default: 1)",             dest="pmh")
parser.add_argument("-pmu", "--pointmassu",    action="store",      type=float,                  default=0.,     help="int. energy of the point mass (default: 0)",               dest="pmu")
parser.add_argument("-pmt", "--pointmasstype", action="store",      type=float,                  default=2.,     help="type of the point mass particle (default: 2)",             dest="pmt")
parser.add_argument("-pmp", "--pointmasspres", action="store",      type=float,                  default=2.e+28, help="pressure of the point mass (default: 2e28)",               dest="pmp")
parser.add_argument("-dir", "--direction_2",   action="store",      type=int,                    default=1,      help="direction of the 2nd star (+/-1)*x_pos (default:+1)",      dest="dir2")
#parser.add_argument("-t", "--type", action="store", type=str, help="type of the star 'wd' or 'ns' ('WD' or 'NS')", dest="type")
args = parser.parse_args()

#def get_WDradius(WD_Mass):
#  Rwd = 1.e9*(WD_Mass/0.7/M_solar)**(-1./3.)*(1.-(WD_Mass/1.44/M_solar)**(4./3.))**(1./2.)*(2./2.)**(-5./3.)
#  return Rwd

def main():
  # read the input file
  try:
    h5_in = h5py.File(args.infile[0].name,'r')
  except:
    sys.exit ("ERROR: cannot read first input file %s" % args.infile)
  if(len(args.infile) == 2):
    try:
      h5_in2 = h5py.File(args.infile[1].name,'r')
    except:
      sys.exit ("ERROR: cannot read second input file %s" % args.infile)

  try:
    os.remove("binary.h5part")
  except:
    pass

  if (args.orbsep == -1.0):
    sys.exit ("ERROR: no orbital separation supplied in arguments")
  # Open the output file HDF5
  output_name = "binary.h5part"
  out = h5py.File(output_name,'w')

  # Open the first input file HDF5
  dataset = ("Step#%d" % (len(list(h5_in.keys()))-1))
  dsetP  = h5_in[dataset+"/P"]
  dsetX  = h5_in[dataset+"/x"]
  dsetY  = h5_in[dataset+"/y"]
  dsetZ  = h5_in[dataset+"/z"]
  dsetVX = h5_in[dataset+"/vx"]
  dsetVY = h5_in[dataset+"/vy"]
  dsetVZ = h5_in[dataset+"/vz"]
  dsetAX = h5_in[dataset+"/ax"]
  dsetAY = h5_in[dataset+"/ay"]
  dsetAZ = h5_in[dataset+"/az"]
  dsetD  = h5_in[dataset+"/rho"]
  dsetM  = h5_in[dataset+"/m"]
  dsetH  = h5_in[dataset+"/h"]
  dsetU  = h5_in[dataset+"/u"]

  size   = len(dsetP[()])
  try:
    dsett  = h5_in[dataset+"/type"]
  except:
    dsett  = np.zeros((size))
  print("calculating total mass star 1")
  M_star = sum(dsetM[()])
  print("mass of star 1 in solar masses: {0:3.2f}".format(M_star/M_solar))
  #R_star = get_WDradius(M_star)

  # Open the second input file HDF5
  if(len(args.infile) == 2):
    orient2 = args.dir2
    dataset2 = ("Step#%d" % (len(list(h5_in2.keys()))-1))
    dsetP2  = h5_in2[dataset2+"/P"]
    dsetX2  = h5_in2[dataset2+"/x"]
    dsetY2  = h5_in2[dataset2+"/y"]
    dsetZ2  = h5_in2[dataset2+"/z"]
    dsetVX2 = h5_in2[dataset2+"/vx"]
    dsetVY2 = h5_in2[dataset2+"/vy"]
    dsetVZ2 = h5_in2[dataset2+"/vz"]
    dsetAX2 = h5_in2[dataset2+"/ax"]
    dsetAY2 = h5_in2[dataset2+"/ay"]
    dsetAZ2 = h5_in2[dataset2+"/az"]
    dsetD2  = h5_in2[dataset2+"/rho"]
    dsetM2  = h5_in2[dataset2+"/m"]
    dsetH2  = h5_in2[dataset2+"/h"]
    dsetU2  = h5_in2[dataset2+"/u"]
    dsett2  = h5_in[dataset+"/type"]
    size2   = len(dsetP2[()])
    try:
      dsett2  = h5_in[dataset+"/type"]
    except:
      dsett2  = np.zeros((size))
    print("calculating total mass star 2")
    M_star2 = sum(dsetM2[()])
    print("mass of star 2 in solar masses: {0:3.2f}".format(M_star2/M_solar))
    #R_star2 = get_WDradius(M_star2)

    newsize = size + size2
    Mtot = M_star + M_star2
    sep = args.orbsep
    x_offset  = sep*M_star2/Mtot
    x_offset2 = sep*M_star/Mtot

    aran = np.arange
    arra = np.array
    zero = np.zeros
    part_id = aran(newsize)
    tempX   = zero((newsize))
    tempY   = zero((newsize))
    tempR   = zero((newsize))
    tempP   = zero((newsize))
    tempM   = zero((newsize))
    tempD   = zero((newsize))
    tempO   = zero((newsize))
    tempVY  = zero((newsize))
    tempVX  = zero((newsize))
    temp    = zero((newsize))
    omega = np.sqrt(G_newt*(Mtot)/(sep**3.0))
    state = np.hstack((np.full(dsetX.shape,1,dtype=int),np.full(dsetX2.shape,2,dtype=int)))

    print("calculating X and Y coordinates")
    tempX[:size] = dsetX[()] - x_offset
    tempX[size:] = orient2*dsetX2[()] + x_offset2
    tempY[:size] = dsetY[()]
    tempY[size:] = dsetY2[()]
    print("calculating Vx and Vy velocities")
    tempR = np.sqrt(tempX*tempX +tempY*tempY)
    tempO = np.arctan2(tempY,tempX)
    tempVX[:size] = dsetVX[()] - omega*tempR[:size]*np.sin(tempO[:size])
    tempVX[size:] = orient2*dsetVX2[()] - omega*tempR[size:]*np.sin(tempO[size:])
    tempVY[:size] = dsetVY[()] + omega*tempR[:size]*np.cos(tempO[:size])
    tempVY[size:] = dsetVY2[()] + omega*tempR[size:]*np.cos(tempO[size:])
    grp = out.create_group("/Step#0")
    print("done")
    print("setting P and x")
    temp[:size] = dsetP[()]
    temp[size:] = dsetP2[()]
    grp.create_dataset("P",data=temp)
    grp.create_dataset("x",data=tempX)
    print("setting y")
    grp.create_dataset("y",data=tempY)
    print("setting z")
    temp[:size] = dsetZ[()]
    temp[size:] = dsetZ2[()]
    grp.create_dataset("z",data=temp)
    print("setting vx")
    grp.create_dataset("vx",data=tempVX)
    print("setting vy")
    grp.create_dataset("vy",data=tempVY)
    print("setting vz")
    temp[:size] = dsetVZ[()]
    temp[size:] = dsetVZ2[()]
    grp.create_dataset("vz",data=temp)
    print("setting ax")
    temp[:size] = 0.0
    temp[size:] = 0.0
    grp.create_dataset("ax",data=temp)
    print("setting ay")
    temp[:size] = 0.0
    temp[size:] = 0.0
    grp.create_dataset("ay",data=temp)
    print("setting az")
    temp[:size] = 0.0
    temp[size:] = 0.0
    grp.create_dataset("az",data=temp)
    print("setting rho")
    temp[:size] = dsetD[()]
    temp[size:] = dsetD2[()]
    grp.create_dataset("rho",data=temp)
    print("setting m")
    temp[:size] = dsetM[()]
    temp[size:] = dsetM2[()]
    grp.create_dataset("m",data=temp)
    print("setting h")
    temp[:size] = dsetH[()]
    temp[size:] = dsetH2[()]
    grp.create_dataset("h",data=temp)
    print("setting u")
    tempP[:size] = dsetP[()]
    tempP[size:] = dsetP2[()]
    tempM[:size] = dsetM[()]
    tempM[size:] = dsetM2[()]
    tempD[:size] = dsetD[()]
    tempD[size:] = dsetD2[()]
    temp = 3./2.*tempP*tempM/tempD
    grp.create_dataset("u",data=temp)
    print("setting type")
    temp[:size] = 0
    temp[size:] = 0
    grp.create_dataset("type",data=temp)
    print("setting id")
    grp.create_dataset("id",data=part_id)
    print("setting state")
    grp.create_dataset("state",data=state)
  else:
    if(args.ident):
      newsize = 2*size
      print("using identical mass binary")
      M_star2 = M_star
      print("mass of star 2 in solar masses: {0:3.2f}".format(M_star2/M_solar))
    else:
      newsize = size+1
      M_star2 = args.pm
      print("using point mass secondary")
      print("mass of star 2 in solar masses: {0:3.2f}".format(M_star2/M_solar))

    Mtot = M_star + M_star2
    sep = args.orbsep

    x_offset  = sep*M_star2/Mtot
    x_offset2 = sep*M_star/Mtot
    aran = np.arange
    arra = np.array
    zero = np.zeros
    part_id = aran(newsize)
    tempX   = zero((newsize))
    tempY   = zero((newsize))
    tempVX  = zero((newsize))
    tempVY  = zero((newsize))
    tempR   = zero((newsize))
    tempP   = zero((newsize))
    tempM   = zero((newsize))
    tempD   = zero((newsize))
    tempO   = zero((newsize))
    temp    = zero((newsize))
    omega = np.sqrt(G_newt*(Mtot)/(sep**3.0))
    state = np.hstack((np.full(dsetX.shape,1,dtype=int),np.full(dsetX.shape,2,dtype=int)))
    
    if(args.ident):
      print("calculating X and Y coordinates")
      tempX[:size] = dsetX[()] - x_offset
      tempX[size:] = -dsetX[()] + x_offset2
      tempY[:size] = dsetY[()]
      tempY[size:] = dsetY[()]
      print("calculating Vx and Vy velocities")
      tempR = np.sqrt(tempX*tempX +tempY*tempY)
      tempO = np.arctan2(tempY,tempX)
      tempVX[:size] = dsetVX[()] - omega*tempR[:size]*np.sin(tempO[:size])
      tempVX[size:] = -dsetVX[()] - omega*tempR[size:]*np.sin(tempO[size:])
      tempVY[:size] = dsetVY[()] + omega*tempR[:size]*np.cos(tempO[:size])
      tempVY[size:] = dsetVY[()] + omega*tempR[size:]*np.cos(tempO[size:])
      grp = out.create_group("/Step#0")
      print("done")
      print("setting P and x")
      temp[:size] = dsetP[()]
      temp[size:] = dsetP[()]
      grp.create_dataset("P",data=temp)
      grp.create_dataset("x",data=tempX)
      print("setting y")
      grp.create_dataset("y",data=tempY)
      print("setting z")
      temp[:size] = dsetZ[()]
      temp[size:] = dsetZ[()]
      grp.create_dataset("z",data=temp)
      print("setting vx")
      grp.create_dataset("vx",data=tempVX)
      print("setting vy")
      grp.create_dataset("vy",data=tempVY)
      print("setting vz")
      temp[:size] = dsetVZ[()]
      temp[size:] = dsetVZ[()]
      grp.create_dataset("vz",data=temp)
      print("setting ax")
      temp[:size] = dsetAX[()]
      temp[size:] = -dsetAX[()]
      grp.create_dataset("ax",data=temp)
      print("setting ay")
      temp[:size] = dsetAY[()]
      temp[size:] = dsetAY[()]
      grp.create_dataset("ay",data=temp)
      print("setting az")
      temp[:size] = dsetAZ[()]
      temp[size:] = dsetAZ[()]
      grp.create_dataset("az",data=temp)
      print("setting rho")
      temp[:size] = dsetD[()]
      temp[size:] = dsetD[()]
      grp.create_dataset("rho",data=temp)
      print("setting m")
      temp[:size] = dsetM[()]
      temp[size:] = dsetM[()]
      grp.create_dataset("m",data=temp)
      print("setting h")
      temp[:size] = dsetH[()]
      temp[size:] = dsetH[()]
      grp.create_dataset("h",data=temp)
      print("setting u")
      tempP[:size] = dsetP[()]
      tempP[size:] = dsetP2[()]
      tempM[:size] = dsetM[()]
      tempM[size:] = dsetM2[()]
      tempD[:size] = dsetD[()]
      tempD[size:] = dsetD2[()]
      temp = 3./2.*tempP*tempM/tempD
      grp.create_dataset("u",data=temp)
      print("setting type")
      temp[:size] = 0
      temp[size:] = 0
      grp.create_dataset("type",data=temp)
      print("setting id")
      grp.create_dataset("id",data=part_id)
      print("setting state")
      grp.create_dataset("state",data=state)
    else:
      state = np.hstack((np.full(dsetX.shape,1,dtype=int),np.full(dsetX.shape,3,dtype=int)))
      print("calculating X and Y coordinates")
      tempX[:size] = dsetX[()] - x_offset
      tempX[size:] = x_offset2
      tempY[:size] = dsetY[()]
      tempY[size:] = 0.
      print("calculating Vx and Vy velocities")
      tempR = np.sqrt(tempX*tempX +tempY*tempY)
      tempO = np.arctan2(tempY,tempX)
      tempVX[:size] = dsetVX[()] - omega*tempR[:size]*np.sin(tempO[:size])
      tempVX[size:] = 0.
      tempVY[:size] = dsetVY[()] + omega*tempR[:size]*np.cos(tempO[:size])
      tempVY[size:] = omega*tempX[size:]
      grp = out.create_group("/Step#0")
      print("done")
      print("setting P and x")
      temp[:size] = dsetP[()]
      temp[size:] = args.pmp
      grp.create_dataset("P",data=temp)
      grp.create_dataset("x",data=tempX)
      print("setting y")
      grp.create_dataset("y",data=tempY)
      print("setting z")
      temp[:size] = dsetZ[()]
      temp[size:] = 0.0
      grp.create_dataset("z",data=temp)
      print("setting vx")
      grp.create_dataset("vx",data=tempVX)
      print("setting vy")
      grp.create_dataset("vy",data=tempVY)
      print("setting vz")
      temp[:size] = dsetVZ[()]
      temp[size:] = 0.0
      grp.create_dataset("vz",data=temp)
      print("setting ax")
      temp[:size] = dsetAX[()]
      temp[size:] = 0.0
      grp.create_dataset("ax",data=temp)
      print("setting ay")
      temp[:size] = dsetAY[()]
      temp[size:] = 0.0
      grp.create_dataset("ay",data=temp)
      print("setting az")
      temp[:size] = dsetAZ[()]
      temp[size:] = 0.0
      grp.create_dataset("az",data=temp)
      print("setting rho")
      temp[:size] = dsetD[()]
      temp[size:] = args.pmd
      grp.create_dataset("rho",data=temp)
      print("setting m")
      temp[:size] = dsetM[()]
      temp[size:] = M_star2
      grp.create_dataset("m",data=temp)
      print("setting h")
      temp[:size] = dsetH[()]
      temp[size:] = args.pmh
      grp.create_dataset("h",data=temp)
      print("setting u")
      tempP[:size] = dsetP[()]
      tempM[:size] = dsetM[()]
      tempD[:size] = dsetD[()]
      temp[:size] = 3./2.*tempP[:size]*tempM[:size]/tempD[:size]
      temp[size:] = args.pmu
      grp.create_dataset("u",data=temp)
      print("setting type")
      temp[:size] = 0
      temp[size:] = 0 #args.pmt
      grp.create_dataset("type",data=temp)
      print("setting id")
      grp.create_dataset("id",data=part_id)
      print("setting state")
      grp.create_dataset("state",data=state)

  print("Done creating hdf5 file")
  out.close()
  h5_in.close()
  if(len(args.infile) == 2):
    h5_in2.close()

if __name__ == "__main__":
    main()
