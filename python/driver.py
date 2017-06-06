import h5py
import numpy as np
import random

from hdf5ParticleIO import HDF5ParticleIO

# Create file
testIO = HDF5ParticleIO("test.h5part")

# Create dataset attribute
testIO.writeDatasetAttribute("numParticles", 80)
testIO.writeDatasetAttribute("gravitational constant", 6.67300E-11)
testIO.writeDatasetAttribute("numDims", 3)
testIO.writeDatasetAttribute("use_fixed_timestep", 0)
testIO.writeDatasetAttribute("sim_name", "random_data")






for ts in range(0,4):

	testIO.setTimeStep(ts)

	testIO.writeTimestepAttribute("sim_time", ts/10.0)
	testIO.writeTimestepAttribute("step_time", ts)

	data_x  = np.random.random(80)
	data_y  = np.random.random(80)
	data_z  = np.random.random(80)

	data_pressure = []
	for j in range(80):
		data_pressure.append( ts )

	testIO.writeVariable("x", data_x)
	testIO.writeVariable("y", data_y)
	testIO.writeVariable("z", data_z)
	testIO.writeVariable("pressure", data_pressure)

	data_pressure = []








#sudo -H pip --proxy http://proxyout.lanl.gov:8080 install h5py