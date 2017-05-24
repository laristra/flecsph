#ifndef _SIM_IO_H_
#define _SIM_IO_H_

#include <string>
#include <vector>

namespace Flecsi_Sim_IO
{


//
// Enums
enum VariableType
{ 
	point               = 0,
	grid_cellCentered   = 1,
	grid_vertexCentered = 2,
	grid_faceCentered   = 3
};


enum OutputType
{
	structuredGrid 	 = 0,
	unStructuredGrid = 1
};


enum TreeType
{
	octree = 0,
	kdtree = 1
};


enum Endianness
{
	little  = 0,
	big = 1
};



struct Variable
{
	std::string name;				// Variable name
	VariableType varType;			// point or ...
	std::string dataType;			// datatype: int, char, float, ...
	int numElements;
	void *data;						// the actual data

	Variable(){ data = NULL; }
	~Variable()
	{
		numElements = 0;

		if (data != NULL) return;
		
		if ( dataType == "float")			delete [] (float *) data;
		else if ( dataType == "double")		delete [] (double *) data;
		else if ( dataType == "int")		delete [] (int *) data;
		else if ( dataType == "int8_t")		delete [] (int8_t *) data;
		else if ( dataType == "int16_t")	delete [] (int16_t *) data;
		else if ( dataType == "int32_t")	delete [] (int32_t *) data;
		else if ( dataType == "int64_t")	delete [] (int64_t *) data;
		else if ( dataType == "uint8_t")	delete [] (uint8_t *) data;
		else if ( dataType == "uint16_t")	delete [] (uint16_t *) data;
		else if ( dataType == "uint32_t")	delete [] (uint32_t *) data;
		else if ( dataType == "uint64_t")	delete [] (uint64_t *) data;
		else { }
	}

	void createVariable(std::string _name, VariableType _varType, std::string _dataType, int _numElements, void *_data)
	{
		name = _name;
		varType = _varType;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}
};



class SimIO
{
  public:  
	std::string outputFileName;
	Endianness endian;
	OutputType datasetType;
	int numVars;
	int numDims;					// 1 or 2 or 3 ...
	
  public:
  	std::vector<Variable> vars;
  	std::vector<int> gridDims;		// Size of each dimension
	std::vector<double> extents;	// (min, max) pair for each dimension

  	SimIO(){ endian=little; }
	SimIO(std::string _outputFileName):outputFileName(_outputFileName){ endian=little; }
	~SimIO(){};

	void createDataset(OutputType _datasetType, int _numVars, int _numDims)
			{ datasetType=_datasetType;  numVars=_numVars;  numDims=_numDims; }
	virtual int writeGridData() = 0;
	virtual int writePointData() = 0;
	virtual int writePointData(int ts) = 0;

	void setFilename(std::string _name){ outputFileName=_name; }
	void setEndianness(Endianness _en){ endian=_en; }
	void setDatasetType(OutputType _datasetType){ datasetType=_datasetType; }
	void setNumVars(int _numVars){ numVars=_numVars; }
	void setNumDims(int _numDims){ numDims=_numDims; }
};


} // end Flecsi_Sim_IO namespace

#endif
