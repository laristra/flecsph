#include "hdf5ParticleIO.h"

#include <boost/python.hpp>

using namespace boost::python;
using namespace Flecsi_Sim_IO;


BOOST_PYTHON_MODULE(hdf5PythonWrapper)
{
	class_<HDF5ParticleIO>("HDF5ParticleIO")
		.def("closeFile", &HDF5ParticleIO::closeFile)
		.def("setFilename", &HDF5ParticleIO::setFilename)
		.def("getFilename", &HDF5ParticleIO::getFilename)
		.def("createDataset", &HDF5ParticleIO::getFilename)
	;
}
