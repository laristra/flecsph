import h5py
import numpy

class HDF5ParticleIO:
	outputFileName = ""
	currentTimestep = -1

	def __init__(self, _outputFileName):
		self.createDataset(_outputFileName)

	#def __del__(self):
	#	self.closeFile()


	def createDataset(self, _outputFileName):
		self.outputFileName = _outputFileName
		self.dataFile = h5py.File(self.outputFileName,'w')

	def closeFile(self):
		self.dataFile.close()
	

	def setTimeStep(self, _ts):
		self.currentTimestep = self.currentTimestep + 1
		self.currentGroup = self.dataFile.create_group("Step#" + str(self.currentTimestep) )


	def writeDatasetAttribute(self, _name, _value):
		if type(_value) is str:
			_value = numpy.string_(_value)
		self.dataFile.attrs.create(_name, _value)

	def writeTimestepAttribute(self, _name, _value):
		self.currentGroup.attrs.create(_name, _value)

	def writeVariable(self, _varName, _data):
		self.currentGroup.create_dataset(_varName, data=_data)




