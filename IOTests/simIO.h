#ifndef _SIM_IO_H_
#define _SIM_IO_H_

#include <string>
#include <vector>

//#include "field.h"
//#include "octree.h"
//#include "cellSet.h"

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
	void *data;						// the actual data

	Variable(){ data = NULL; }
	~Variable(){
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
};



class SimIO
{
  public:  
	std::string outputFileName;
	Endianness endian;

	OutputType datasetType;
	int numDims;					// 1 or 2 or 3 ...

	std::vector<int> gridDims;		// Size of each dimension
	std::vector<double> extents;	// (min, max) pair for each dimension

	std::vector<Variable> vars;

  public:
  	SimIO(){ endian = little; }
	SimIO(std::string _outputFileName):outputFileName(_outputFileName){ endian = little; }
	~SimIO(){};

	void setFilename(std::string _name){ outputFileName = _name; }
	void setEndianness(Endianness _en){ endian = _en; }


	virtual int createDataset() = 0;

	virtual int writeGridData() = 0;
	virtual int writePointData() = 0;
};


} // end Flecsi_Sim_IO namespace

#endif
