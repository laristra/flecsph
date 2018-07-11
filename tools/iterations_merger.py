"""
Step mergers for h5part files 
This merger takes in input several h5part and merge them regarding the 
iterations:  f_0001.h5part f_0002.h5part ==> merged.h5part
THIS MERGER DOES NOT MERGE FILES FROM DIFFERENT PROCESSES
"""

# Store data in h5part format
import h5py
import sys 
import argparse

parser = argparse.ArgumentParser()      
parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
args = parser.parse_args()

# Open the output file HDF5
output_name = "merged.h5part"
out = h5py.File(output_name,'w')
out_step = 0

for f in args.file:
	in_step = 0
	done = False
	print("Merging: "+f.name)
	
	h5_in = h5py.File(f.name,'r')
	
	while(not(done)):
		dataset = "/Step#"+str(in_step)
		if dataset not in h5_in.keys():
			done = True
			continue
		grp = out.create_group("/Step#"+str(out_step))
		dset = h5_in[dataset+"/P"]
		grp.create_dataset("P",data=dset)
		dsetX = h5_in[dataset+"/x"]
		grp.create_dataset("x",data=dsetX)
		dsetY = h5_in[dataset+"/y"]
		grp.create_dataset("y",data=dsetY)
		dsetZ = h5_in[dataset+"/z"]
		grp.create_dataset("z",data=dsetZ)
		dset = h5_in[dataset+"/vx"]
		grp.create_dataset("vx",data=dset)
		dset = h5_in[dataset+"/vy"]
		grp.create_dataset("vy",data=dset)
		dset = h5_in[dataset+"/vz"]
		grp.create_dataset("vz",data=dset)
		dset = h5_in[dataset+"/ax"]
		grp.create_dataset("ax",data=dset)
		dset = h5_in[dataset+"/ay"]
		grp.create_dataset("ay",data=dset)
		dset = h5_in[dataset+"/az"]
		grp.create_dataset("az",data=dset)
		dsetD = h5_in[dataset+"/density"]
		grp.create_dataset("density",data=dsetD)
		dsetM = h5_in[dataset+"/mass"]
		grp.create_dataset("mass",data=dsetM)
		dsetH = h5_in[dataset+"/h"]
		grp.create_dataset("h",data=dsetH)
		dsetU = h5_in[dataset+"/u"]
		grp.create_dataset("u",data=dsetU)
		dsetID = h5_in[dataset+"/id"]
		grp.create_dataset("id",data=dsetID)

		out_step = out_step + 1
		in_step = in_step + 1
	print("Done file: "+f.name+" at step: "+str(in_step)+"/"+str(out_step))
	h5_in.close()

print("Merging done: "+str(out_step)+" steps total")
out.close()
