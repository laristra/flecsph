#ifndef _SIM_IO_H_
#define _SIM_IO_H_

#include <string>
#include <vector>
#include <iostream>

namespace Flecsi_Sim_IO
{


//
// Enums
enum Operation
{
	READING = 0,
	WRITING = 1
};



//
// Structures
struct Attribute
{
	std::string name;				// Variable name
	std::string dataType;			// datatype: int, char, float, ...
	int numElements;
	void *data;						// the actual data	

	Attribute(){ data = NULL; }

	template <class T>
	Attribute(std::string _name, std::string _dataType, T _data)
	{
		createAttribute(_name, _dataType, _data);
	}


	Attribute(std::string _name, std::string _dataType, int _numElements, void *_data)
	{
		createAttributeArray(_name, _dataType, _numElements, _data);
	}

	~Attribute(){ numElements = 0; }


	template <class T>
	void createAttribute(std::string _name, std::string _dataType, T _data)
	{
		name = _name;
		dataType = _dataType;
		numElements = 1;

		if ( dataType == "float")       { data = new float[numElements];	((float *)  data)[0] = _data; }
		else if ( dataType == "double") { data = new double[numElements];	((double *) data)[0] = _data; }
		else if ( dataType == "int32_t"){ data = new int32_t[numElements];	((int32_t *)data)[0] = _data; }
		else if ( dataType == "int64_t"){ data = new int64_t[numElements];	((int64_t *)data)[0] = _data; }
		else if ( dataType == "string") { data = new char[numElements];		((char *)   data)[0] = _data; }
		else { }
	}


	void createAttributeArray(std::string _name, std::string _dataType, int _numElements, void *_data)
	{
		name = _name;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}


	template <class T>
	T getAttributeValue()
	{
		if ( dataType != "string")
			return static_cast<T *>(data)[0];
	}


	template <class T>
	void getAttributeArray(T _data[])
	{
		std::cout << numElements << std::endl;
		for (int i=0; i<numElements; i++)
		{
			std::cout << static_cast<T *>(data)[i] << std::endl;
			_data[i] = static_cast<T *>(data)[i];
		}
	}


	void deleteData()
	{

	}
};
					

struct Variable
{
	std::string name;				// Variable name
	std::string dataType;			// datatype: int, char, float, ...
	int numElements;
	void *data;						// the actual data
				

	Variable(){ data = NULL; }
	Variable(std::string _name, std::string _dataType, int _numElements, void *_data)
	{ 
		name = _name;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}

	~Variable()
	{
		numElements = 0;
	}

	void createVariable(std::string _name, std::string _dataType, int _numElements, void *_data)
	{
		name = _name;
		dataType = _dataType;
		numElements = _numElements;
		data = _data;
	}

	void deleteData()
	{

	}
};


struct Timestep
{
	int timeStep;

	int numAttributes;
	std::vector<Attribute> attributes;

	int numVariables;
	std::vector<Variable> vars;

	Timestep(){ timeStep=0; numAttributes=0; numVariables=0;}
	Timestep(int _ts){ timeStep = _ts; }
	void addAttribute(Attribute _a){ attributes.push_back(_a); numAttributes = attributes.size(); }
	void addVariable(Variable _v){ vars.push_back(_v); numVariables = vars.size(); }
};



class SimIO
{
  public:
	//std::string outputFileName;
	std::string fileName;
	int numTimesteps;
	int numDatasetAttributes;

  public:
  	std::vector<Attribute> datasetAttributes;
  	std::vector<Timestep> timesteps;

  	std::vector<Variable> vars;					// To remove
  	std::vector<Attribute> timestepAttributes;	// To remove

  		

  	SimIO(){ init(); }
	SimIO(std::string _fileName):fileName(_fileName){ init(); }
	~SimIO(){};

	void init();

	void addDatasetAttribute(Attribute _a){  datasetAttributes.push_back(_a); }
	void addTimestep(Timestep _t){ timesteps.push_back(_t); }


	void setFilename(std::string _name){ fileName=_name; }

	std::string getFilename(){ return fileName; }
	int getNumTimesteps(){ return numTimesteps; }
	int getNumDatasetAttributes(){ return numDatasetAttributes; }



	// To remove
	void addTimeStepAttribute(Attribute _a){ timestepAttributes.push_back(_a); }
	void addVariable(Variable _v){  vars.push_back(_v); }
};


inline void SimIO::init()
{ 
	numTimesteps = 0;
	numDatasetAttributes = 0;
	
	datasetAttributes.clear();
	timesteps.clear();
}


} // end Flecsi_Sim_IO namespace

#endif
