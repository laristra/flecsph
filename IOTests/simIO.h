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


enum AttributeType
{ 
	global   = 0,
	timestep = 1
};


enum OutputType
{
	structuredGrid 	 = 0,
	unStructuredGrid = 1
};


enum Operation
{
	READING = 0,
	WRITING = 1
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



struct Attribute
{
	std::string name;				// Variable name
	AttributeType attributeType;	// global or timestep
	std::string dataType;			// datatype: int, char, float, ...
	int numElements;
	void *data;						// the actual data	

	Attribute(){ data = NULL; }

	template <class T>
	Attribute(std::string _name, AttributeType _attributeType, std::string _dataType, T _data)
	{
		name = _name;
		attributeType = _attributeType;
		dataType = _dataType;
		numElements = 1;

		if ( dataType == "float")       { data = new float[numElements];	((float *)  data)[0] = _data; }
		else if ( dataType == "double") { data = new double[numElements];	((double *) data)[0] = _data; }
		else if ( dataType == "int32_t"){ data = new int32_t[numElements];	((int32_t *)data)[0] = _data; }
		else if ( dataType == "int64_t"){ data = new int64_t[numElements];	((int64_t *)data)[0] = _data; }
		else if ( dataType == "string") { data = new char[numElements];		((char *)   data)[0] = _data; }
		else { }
	}

	~Attribute(){ numElements = 0; }


	template <class T>
	void createAttribute(std::string _name, AttributeType _attributeType, std::string _dataType, T _data)
	{
		name = _name;
		attributeType = _attributeType;
		dataType = _dataType;
		numElements = 1;

		if ( dataType == "float")       { data = new float[numElements];	((float *)  data)[0] = _data; }
		else if ( dataType == "double") { data = new double[numElements];	((double *) data)[0] = _data; }
		else if ( dataType == "int32_t"){ data = new int32_t[numElements];	((int32_t *)data)[0] = _data; }
		else if ( dataType == "int64_t"){ data = new int64_t[numElements];	((int64_t *)data)[0] = _data; }
		else if ( dataType == "string") { data = new char[numElements];		((char *)   data)[0] = _data; }
		else { }
	}

	void createAttributeArray(std::string _name, AttributeType _attributeType, std::string _dataType, int _numElements, void *_data)
	{
		name = _name;
		attributeType = _attributeType;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}


	template <class T>
	T getAttributeValue()
	{
		if ( dataType == "float")       { return ((float *)  data)[0]; }
		else if ( dataType == "double") { return ((double *) data)[0]; }
		else if ( dataType == "int32_t"){ return ((int32_t *)data)[0]; }
		else if ( dataType == "int64_t"){ return ((int64_t *)data)[0]; }
		else if ( dataType == "string") { return ((char *)   data)[0]; }
		else { return -1; }
	}


	// void getAttributeArray(void *_data)
	// {
	// 	if ( dataType == "float")       { return ((float *)  data)[0]; }
	// 	else if ( dataType == "double") { return ((double *) data)[0]; }
	// 	else if ( dataType == "int32_t"){ return ((int32_t *)data)[0]; }
	// 	else if ( dataType == "int64_t"){ return ((int64_t *)data)[0]; }
	// 	else if ( dataType == "string") { return ((char *)   data)[0]; }
	// 	else { return -1; }
	// }
};

					

struct Variable
{
	std::string name;				// Variable name
	VariableType varType;			// point or ...
	std::string dataType;			// datatype: int, char, float, ...
	int numElements;
	void *data;						// the actual data
				

	Variable(){ data = NULL; }
	Variable(std::string _name, VariableType _varType, std::string _dataType, int _numElements, void *_data)
	{ 
		name = _name;
		varType = _varType;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}

	~Variable()
	{
		numElements = 0;
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


struct Timestep
{
	std::vector<Attribute> attributes;
	std::vector<Variable> vars;
};


class SimIO
{
  public:  
	std::string outputFileName;
	Endianness endian;
	OutputType datasetType;
	int numTimesteps;
	int numDatasetAttributes;
	int numVars;
	int numDims;					// 1 or 2 or 3 ...
	
  public:
  	std::vector<int> gridDims;		// Size of each dimension
	std::vector<double> extents;	// (min, max) pair for each dimension

  	std::vector<Variable> vars;
  	std::vector<Attribute> timestepAttributes;

  	std::vector<Attribute> datasetAttributes;
  	std::vector<Timestep> timesteps;
  	

  	SimIO(){ init(); }
	SimIO(std::string _outputFileName):outputFileName(_outputFileName){ init(); }
	~SimIO(){};

	void init();

	void setFilename(std::string _name){ outputFileName=_name; }
	void setEndianness(Endianness _en){ endian=_en; }
	void setDatasetType(OutputType _datasetType){ datasetType=_datasetType; }
	void setNumVars(int _numVars){ numVars=_numVars; }
	void setNumDims(int _numDims){ numDims=_numDims; }

	std::string getFilename(){ return outputFileName; }

	int getNumTimesteps(){ return numTimesteps; }
	int getNumDatasetAttributes(){ return numDatasetAttributes; }

	
	void addTimeStepAttribute(Attribute _a){ timestepAttributes.push_back(_a); }
	void addVariable(Variable _v){  vars.push_back(_v); }
	
	void clearTimestepAttributes(){ timestepAttributes.clear(); }
	void clearTimestepVariables() { vars.clear();  }
	void clearTimestepData()	  { vars.clear(); timestepAttributes.clear(); }
};


inline void SimIO::init()
{ 
	numDims = 0;
	endian = little; 
	numTimesteps = 0;
	numDatasetAttributes = 0;
	
}



} // end Flecsi_Sim_IO namespace

#endif
