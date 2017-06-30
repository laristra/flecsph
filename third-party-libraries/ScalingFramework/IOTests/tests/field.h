#ifndef _FIELD_H_
#define _FIELD_H_

#include "simIO.h"

class Field
{
  public:
	std::string name;		// name of field; e.g pressure, temperature
	std::string dataType;	// datatype: int, char, float, ...
	FieldType association;	// cell-centered, face-centered, ...
	size_t numElements;		// # elements

	void *data;				// the actual data

  public:
  	Field(){};
  	Field(std::string _name):name(_name){}
  	~Field(){};

  	void setName(std::string _name){ name=_name; }
  	void setDataype(std::string _dataType){ dataType=_dataType; }
  	void setAssociation(FieldType _association){ association = _association; }
  	void setNumElements(size_t _numElements){ numElements=_numElements; }
};

#endif