#!/usr/bin/env python
"""
Takes a(/two) single star h5part file and produces binary system at RLOF:
"""

# Store data in h5part format
import h5py, sys, argparse, os
import numpy as np

# IMPORTANT CONSTANTS
G_newt  = 6.67259e-8 #cm3/g/s2
M_solar = 1.99e33 #g
R_solar = 6.96e10 #cm
c_light = 2.99792458e10 #cm/s

ln = np.log

my_description = """
Takes a(/two) single star h5part file(s) and produces binary system at RLOF or given orbital separation:
  star.h5part (+ star2.h5part)  ==>   binary.h5part"""
my_usage = """
    %(prog)s [-f|--file <filename(s)>] [-i|--identical] [-m|--pointmass <val>] [-mrf|--MRfile <filename>] [-a|--orbsep <val>] [-h|--help]"""

parser = argparse.ArgumentParser(description=my_description, usage=my_usage, epilog="EXAMPLE: $python %(prog)s -f wd.h5part -i -a 1.e9")

parser.add_argument("-f",   "--file",      action="store",      type=argparse.FileType('r'),                help="input file(s) in FleCSPH hdf5 format",          nargs="+", dest="infile")
parser.add_argument("-i",   "--identical", action="store_true",                              default=False, help="set this if you want an equal mass binary",                dest="ident")
parser.add_argument("-pm",  "--pointmass", action="store",      type=float,                  default=-1.,   help="total mass of the point star (default: -1.)",              dest="pm")
parser.add_argument("-mrf", "--MRfile",    action="store",      type=argparse.FileType('r'),                help="Mass Radius file for supplied star type",                  dest="MRfile")
parser.add_argument("-a",   "--orbsep",    action="store",      type=float,                  default=-1.,   help="orbital separation in cm. default RLOF orbsep",            dest="orbsep")
#parser.add_argument("-t", "--type", action="store", type=str, help="type of the star 'wd' or 'ns' ('WD' or 'NS')", dest="type")
args = parser.parse_args()

def get_WDradius(WD_Mass):
  Rwd = 1.e9*(WD_Mass/0.7/M_solar)**(-1./3.)*(1.-(WD_Mass/1.44/M_solar)**(4./3.))**(1./2.)*(2./2.)**(-5./3.)
  return Rwd

def main():
  # read the input file
  try:
    h5_in = h5py.File(args.infile[0].name)
  except:
    sys.exit ("ERROR: cannot read first input file %s" % args.infile)
  if(len(args.infile) == 2):
    try:
      h5_in2 = h5py.File(args.infile[1].name)
    except:
      sys.exit ("ERROR: cannot read second input file %s" % args.infile)

  try:
    os.remove("binary.h5part")
  except:
    pass

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
  print("calculating total mass star 1")
  M_star = sum(dsetM[()])
  print("mass of star 1 in solar masses: ", M_star/M_solar)
  R_star = get_WDradius(M_star)

  # Open the second input file HDF5
  if(len(args.infile) == 2):
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
    size2   = len(dsetP2[()])
    print("calculating total mass star 2")
    M_star2 = sum(dsetM2[()])
    print("mass of star 2 in solar masses: ", M_star/M_solar)
    R_star2 = get_WDradius(M_star2)


    newsize = size + size2
    Mtot = M_star + M_star2
    if (args.orbsep == -1.0):
      if (M_star <= M_star2):
        q = M_star/M_star2
        sep = R_star*(0.6*q**(2./3.)+ln(1.+q**(1./3.)))/(0.49*q**(2./3.))
      else:
        q = M_star2/M_star
        sep = R_star2*(0.6*q**(2./3.)+ln(1.+q**(1./3.)))/(0.49*q**(2./3.))
    else:
      sep = args.orbsep
    x_offset  = sep*M_star2/Mtot
    x_offset2 = sep*M_star/Mtot

    aran = np.arange
    arra = np.array
    zero = np.zeros
    part_id = aran(newsize)
    tempX   = zero((newsize))
    tempVY  = zero((newsize))
    temp    = zero((newsize))
    omega = np.sqrt(G_newt*(M_star + M_star2)/(sep**3.0))

    print("calculate tempX and VY")
    tempX[:size] = x_offset + dsetX[()]
    tempX[size:] = -x_offset2 + dsetX2[()]
    tempVY[:size] = dsetVY[()] + omega*tempX[:size]
    tempVY[size:] = dsetVY2[()] + omega*tempX[size:]
    grp = out.create_group("/Step#0")
    print("done")
    print("set P and x")
    temp[:size] = dsetP[()]
    temp[size:] = dsetP2[()]
    grp.create_dataset("P",data=temp)
    grp.create_dataset("x",data=tempX)
    print("set y")
    temp[:size] = dsetY[()]
    temp[size:] = dsetY2[()]
    grp.create_dataset("y",data=temp)
    print("set z")
    temp[:size] = dsetZ[()]
    temp[size:] = dsetZ2[()]
    grp.create_dataset("z",data=temp)
    print("set vx")
    temp[:size] = dsetVX[()]
    temp[size:] = dsetVX2[()]
    grp.create_dataset("vx",data=temp)
    print("set vy")
    grp.create_dataset("vy",data=tempVY)
    print("set vz")
    temp[:size] = dsetVZ[()]
    temp[size:] = dsetVZ2[()]
    grp.create_dataset("vz",data=temp)
    print("set ax")
    temp[:size] = dsetAX[()]
    temp[size:] = dsetAX2[()]
    grp.create_dataset("ax",data=temp)
    print("set ay")
    temp[:size] = dsetAY[()]
    temp[size:] = dsetAY2[()]
    grp.create_dataset("ay",data=temp)
    print("set az")
    temp[:size] = dsetAZ[()]
    temp[size:] = dsetAZ2[()]
    grp.create_dataset("az",data=temp)
    print("set rho")
    temp[:size] = dsetD[()]
    temp[size:] = dsetD2[()]
    grp.create_dataset("rho",data=temp)
    print("set m")
    temp[:size] = dsetM[()]
    temp[size:] = dsetM2[()]
    grp.create_dataset("m",data=temp)
    print("set h")
    temp[:size] = dsetH[()]
    temp[size:] = dsetH2[()]
    grp.create_dataset("h",data=temp)
    print("set u")
    temp[:size] = dsetU[()]
    temp[size:] = dsetU2[()]
    grp.create_dataset("u",data=temp)
    print("set id")
    grp.create_dataset("id",data=part_id)
  else:
    if(args.ident):
      newsize = 2*size
      print("using identical mass binary")
      M_star2 = M_star
      print("mass of star 2 in solar masses: ", M_star/M_solar)
    else:
      newsize = size+1
      M_star2 = args.pm

    Mtot = M_star + M_star2
    q = M_star/M_star2
    sep = R_star*(0.6*q**(2./3.)+ln(1.+q**(1./3.)))/(0.49*q**(2./3.))

    x_offset  = sep*M_star2/Mtot
    x_offset2 = sep*M_star/Mtot
    aran = np.arange
    arra = np.array
    zero = np.zeros
    part_id = aran(newsize)
    tempX   = zero((newsize))
    tempVY  = zero((newsize))
    temp    = zero((newsize))
    omega = np.sqrt(G_newt*(M_star + M_star2)/(sep**3.0))
    if(args.ident):
      print("calculate tempX and VY")
      tempX[:size] = x_offset + dsetX[()]
      tempX[size:] = -x_offset2 + dsetX[()]
      tempVY[:size] = dsetVY[()] + omega*tempX[:size]
      tempVY[size:] = dsetVY[()] + omega*tempX[size:]
      grp = out.create_group("/Step#0")
      print("done")
      print("set P and x")
      temp[:size] = dsetP[()]
      temp[size:] = dsetP[()]
      grp.create_dataset("P",data=temp)
      grp.create_dataset("x",data=tempX)
      print("set y")
      temp[:size] = dsetY[()]
      temp[size:] = dsetY[()]
      grp.create_dataset("y",data=temp)
      print("set z")
      temp[:size] = dsetZ[()]
      temp[size:] = dsetZ[()]
      grp.create_dataset("z",data=temp)
      print("set vx")
      temp[:size] = dsetVX[()]
      temp[size:] = dsetVX[()]
      grp.create_dataset("vx",data=temp)
      print("set vy")
      grp.create_dataset("vy",data=tempVY)
      print("set vz")
      temp[:size] = dsetVZ[()]
      temp[size:] = dsetVZ[()]
      grp.create_dataset("vz",data=temp)
      print("set ax")
      temp[:size] = dsetAX[()]
      temp[size:] = dsetAX[()]
      grp.create_dataset("ax",data=temp)
      print("set ay")
      temp[:size] = dsetAY[()]
      temp[size:] = dsetAY[()]
      grp.create_dataset("ay",data=temp)
      print("set az")
      temp[:size] = dsetAZ[()]
      temp[size:] = dsetAZ[()]
      grp.create_dataset("az",data=temp)
      print("set rho")
      temp[:size] = dsetD[()]
      temp[size:] = dsetD[()]
      grp.create_dataset("rho",data=temp)
      print("set m")
      temp[:size] = dsetM[()]
      temp[size:] = dsetM[()]
      grp.create_dataset("m",data=temp)
      print("set h")
      temp[:size] = dsetH[()]
      temp[size:] = dsetH[()]
      grp.create_dataset("h",data=temp)
      print("set u")
      temp[:size] = dsetU[()]
      temp[size:] = dsetU[()]
      grp.create_dataset("u",data=temp)
      print("set id")
      grp.create_dataset("id",data=part_id)
    else:
      print("calculate tempX and VY")
      tempX[:size] = x_offset + dsetX[()]
      tempX[size:] = -x_offset2
      tempVY[:size] = dsetVY[()] + omega*tempX[:size]
      tempVY[size:] = omega*tempX[size:]
      grp = out.create_group("/Step#0")
      print("done")
      print("set P and x")
      temp[:size] = dsetP[()]
      temp[size:] = 0.0
      grp.create_dataset("P",data=temp)
      grp.create_dataset("x",data=tempX)
      print("set y")
      temp[:size] = dsetY[()]
      temp[size:] = 0.0
      grp.create_dataset("y",data=temp)
      print("set z")
      temp[:size] = dsetZ[()]
      temp[size:] = 0.0
      grp.create_dataset("z",data=temp)
      print("set vx")
      temp[:size] = dsetVX[()]
      temp[size:] = 0.0
      grp.create_dataset("vx",data=temp)
      print("set vy")
      grp.create_dataset("vy",data=tempVY)
      print("set vz")
      temp[:size] = dsetVZ[()]
      temp[size:] = 0.0
      grp.create_dataset("vz",data=temp)
      print("set ax")
      temp[:size] = dsetAX[()]
      temp[size:] = 0.0
      grp.create_dataset("ax",data=temp)
      print("set ay")
      temp[:size] = dsetAY[()]
      temp[size:] = 0.0
      grp.create_dataset("ay",data=temp)
      print("set az")
      temp[:size] = dsetAZ[()]
      temp[size:] = 0.0
      grp.create_dataset("az",data=temp)
      print("set rho")
      temp[:size] = dsetD[()]
      temp[size:] = 0.0
      grp.create_dataset("rho",data=temp)
      print("set m")
      temp[:size] = dsetM[()]
      temp[size:] = M_star2
      grp.create_dataset("m",data=temp)
      print("set h")
      temp[:size] = dsetH[()]
      temp[size:] = 0.0
      grp.create_dataset("h",data=temp)
      print("set u")
      temp[:size] = dsetU[()]
      temp[size:] = 0.0
      grp.create_dataset("u",data=temp)
      print("set type")
      temp[:size] = 0
      temp[size:] = 2
      grp.create_dataset("type",data=temp)
      print("set id")
      grp.create_dataset("id",data=part_id)

  print("Done creating hdf5")
  out.close()
  h5_in.close()
  if(len(args.infile) == 2):
    h5_in2.close()

if __name__ == "__main__":
    main()
