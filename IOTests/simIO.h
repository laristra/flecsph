#ifndef _SIM_IO_H_
#define _SIM_IO_H_

#include <string>
#include <vector>


namespace Flecsi_Sim_IO
{

enum VariableType
{ 
	point               = 0,
	grid_cellCentered   = 1,
	grid_vertexCentered = 2,
	grid_faceCentered   = 3
};

enum TreeType
{
	octree  = 0,
	k-dtree = 1
};


enum Endianness
{
	little  = 0,
	big = 1
};


struct Variable
{
	std::string name;				// Variable name
	int numDims;					// 1 or 2 or 3 ...
	std::vector<int> gridDims;		// 
	std::vector<double> extents;	// (min, max) pair for each dimension
	VariableType varType;			// point or ...
	std::string dataType;			// datatype: int, char, float, ...
	void* data;						// the actual data
};

struct timeStep
{
	std::vector<Variable> vars;
	int timestep;
	double timeStamp;
};



struct dataOrganization
{	
	bool available;
	TreeType type;
	void *data;
};



class SimIO
{
  public:  
	std::string outputFileName;
	Endianness endian;
	std::vector<timeStep> timeVariables;
	dataOrganization layout;				// storing the organization as just one of the variables

  public:
	SimIO(std::string _outputFileName){ outputFileName = _outputFileName; }
	~SimIO(){}{ }

	virtual int createDataset() = 0;
	virtual int writeData(int _timeStep) = 0;
};


} // end Flecsi_Sim_IO namespace

#endif
