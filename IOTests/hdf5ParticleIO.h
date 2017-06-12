#ifndef _H5_PARTICLE_IO_H
#define _H5_PARTICLE_IO_H

#include "simIO.h"
#include "hdf5.h"
#include "H5hut.h"

#include <iostream>
#include <mpi.h>


namespace Flecsi_Sim_IO
{


class HDF5ParticleIO: public SimIO
{
    h5_file_t *dataFile;
    int myRank;

  public:
    HDF5ParticleIO():SimIO(){ dataFile=NULL; }
    HDF5ParticleIO(std::string _outputFileName, MPI_Comm _comm):SimIO(_outputFileName){ createDataset(_outputFileName, _comm); }
    ~HDF5ParticleIO(){ closeFile(); }

    int createDataset(std::string _outputFileName, MPI_Comm _comm);
    int openFile(MPI_Comm _comm);
    void closeFile();

    // Reading
    int readFile(std::string inputFileName);

    int getNumTimestepAttributes(int ts);
    int getNumDatasetAttributes();
    int getNumVariables(int ts);

    void displayStepAttributes(int ts);
    void displayFileAttributes();

    template <class T>
    T readDatasetAttribute(std::string _varName);
    void readDatasetAttributeArray(std::string _name, std::string _dataType, void *_data, int numElements=1);

    template <class T>
    T readTimestepAttribute(std::string _varName);
    void readTimestepAttributeArray(std::string _name, std::string _dataType, void *_data, int numElements=1);



    // Writing 
    void setTimeStep(int ts){ H5SetStep(dataFile, ts); }

    template <class T>
    int writeDatasetAttribute(std::string _name, std::string _dataType, T _data);
    int writeDatasetAttributeArray(std::string _name, std::string _dataType, void *_data, int numElements=1);

    template <class T>
    int writeTimestepAttribute(std::string _name, std::string _dataType, T _data);
    int writeTimestepAttributeArray(std::string _name, std::string _dataType, void * _data, int numElements=1);

    int writeTimestepAttributes();
    int writeVariables();
};



int HDF5ParticleIO::readFile(std::string inputFileName, int &numTs, int &numAttrs)
{
    dataFile = H5PartOpenFile(inputFileName.c_str(),H5_O_RDONLY);
    if (!dataFile)
        return -1;

    numTimesteps = H5GetNumSteps(dataFile);
    numDatasetAttributes = H5GetNumFileAttribs(dataFile);

    return 0;
}


int HDF5ParticleIO::getNumTimestepAttributes(int ts)
{
    H5PartSetStep(dataFile, ts);
    return H5GetNumStepAttribs(ts);
}

int HDF5ParticleIO::getNumDatasetAttributes()
{

    return H5GetNumFileAttribs(dataFile);
}

int HDF5ParticleIO::getNumVariables(int ts)
{
    H5PartSetStep(dataFile, ts);
    return H5PartGetNumDatasets(ts);
}

























































































void HDF5ParticleIO::readDatasetAttribute(std::string _name, std::string _dataType, void *_data)
{
    if ( _dataType == "int32_t" )
    {
        H5ReadFileAttribInt32(dataFile, _name.c_str(), (int32_t *)_data);
    }
    else if ( _dataType == "int64_t" )
    {
        H5ReadFileAttribInt64(dataFile, _name.c_str(), (int64_t *)_data);
    }
    else if ( _dataType == "float" )
    {
        H5ReadFileAttribFloat32(dataFile, _name.c_str(), (float *)_data);
    }
    else if ( _dataType == "double" )
    {
        H5ReadFileAttribFloat64(dataFile, _name.c_str(), (double *)_data);
    }
    else if ( _dataType == "string" )
    {
        H5ReadFileAttribString(dataFile, _name.c_str(), (char *)_data);
    }
    else
        return -1;
}



void HDF5ParticleIO::getStepAttributes(int ts)
{
    H5PartSetStep(dataFile, ts);

    const int MAX_LEN = 256;
    h5_int64_t type;
    h5_int64_t numElem;
    int numTimestepAttributes = H5PartGetNumStepAttribs(ts);
    for (int i=0; i<numTimestepAttributes; i++)
    {
        h5_int64_t type;
        char fileAttrib[MAX_LEN];

        H5GetStepAttribInfo(dataFile, i, fileAttrib, MAX_LEN, &type, &numElem);
        std::cout << fileAttrib << " " << numElem << " ";
        if (type == H5_INT64_T)
        {
            std::cout << "H5_INT64_T" << std::endl;
        }
        elseif (type == H5_INT32_T)
        {
            std::cout << "H5_INT32_T" << std::endl;
        }
        elseif (type == H5_FLOAT64_T)
        {
            std::cout << "H5_INT64_T" << std::endl;
        }
        elseif (type == H5_FLOAT32_T)
        {
            std::cout << "H5_FLOAT32_T" << std::endl;
        }
        elseif (type == H5_STRING_T)
        {
            std::cout << "H5_STRING_T" << std::endl;
        }
        else
        {
            std::cout << "unknown!!!" << std::endl;
        }
    }
}

void HDF5ParticleIO::getFileAttributes()
{
    const int MAX_LEN = 256;
    h5_int64_t type;
    h5_int64_t numElem;
    int numDatasetAttributes = H5GetNumFileAttribs(dataFile);
    for (int i=0; i<numDatasetAttributes; i++)
    {
        h5_int64_t type;
        char fileAttrib[MAX_LEN];

        H5GetFileAttribInfo(dataFile, i, fileAttrib, MAX_LEN, &type, &numElem);
        std::cout << fileAttrib << " " << numElem << " ";
        if (type == H5_INT64_T)
        {
            std::cout << "H5_INT64_T" << std::endl;
        }
        elseif (type == H5_INT32_T)
        {
            std::cout << "H5_INT32_T" << std::endl;
        }
        elseif (type == H5_FLOAT64_T)
        {
            std::cout << "H5_INT64_T" << std::endl;
        }
        elseif (type == H5_FLOAT32_T)
        {
            std::cout << "H5_FLOAT32_T" << std::endl;
        }
        elseif (type == H5_STRING_T)
        {
            std::cout << "H5_STRING_T" << std::endl;
        }
        else
        {
            std::cout << "unknown!!!" << std::endl;
        }
    }
}





int HDF5ParticleIO::createDataset(std::string _outputFileName, MPI_Comm _comm)
{ 
    MPI_Comm_rank(_comm, &myRank);

    outputFileName = _outputFileName;
    dataFile = H5OpenFile(outputFileName.c_str(), H5_O_RDWR, _comm);
    if (!dataFile)
        return -1;

    return 0;
}


int HDF5ParticleIO::openFile(MPI_Comm  _comm)
{
    dataFile = H5OpenFile(outputFileName.c_str(), H5_O_RDWR, _comm);
    if (!dataFile)
        return -1;

    return 0;
}

void HDF5ParticleIO::closeFile()
{ 
    if (dataFile!=NULL)
    {
        H5CloseFile(dataFile);
        dataFile=NULL;
    }
}

template <class T>
inline int HDF5ParticleIO::writeDatasetAttribute(std::string _name, std::string _dataType, T _data)
{
    if ( _dataType == "int32_t" )
    {
        int32_t data = (int32_t)_data;
        H5WriteFileAttribInt32(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "int64_t" )
    {
        int64_t data = (int64_t)_data;
        H5WriteFileAttribInt64(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "float" )
    {
        float data = (float)_data;
        H5WriteFileAttribFloat32(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "double" )
    {
        double data = (double)_data;
        H5WriteFileAttribFloat64(dataFile, _name.c_str(), &data, 1);
    }
    else
        return -1;

    return 0;
}

inline int HDF5ParticleIO::writeDatasetAttributeArray(std::string _name, std::string _dataType, void *_data, int numElements)
{
    if ( _dataType == "int32_t" )
    {
        H5WriteFileAttribInt32(dataFile, _name.c_str(), (int32_t *)_data, numElements);
    }
    else if ( _dataType == "int64_t" )
    {
        H5WriteFileAttribInt64(dataFile, _name.c_str(), (int64_t *)_data, numElements);
    }
    else if ( _dataType == "float" )
    {
        H5WriteFileAttribFloat32(dataFile, _name.c_str(), (float *)_data, numElements);
    }
    else if ( _dataType == "double" )
    {
        H5WriteFileAttribFloat64(dataFile, _name.c_str(), (double *)_data, numElements);
    }
    else if ( _dataType == "string" )
    {
        H5WriteFileAttribString(dataFile, _name.c_str(), (char *)_data);
    }
    else
        return -1;
}



template <class T>
inline int HDF5ParticleIO::writeTimestepAttribute(std::string _name, std::string _dataType, T _data)
{
    if ( _dataType == "int32_t" )
    {
        int32_t data = (int32_t)_data;
        H5WriteStepAttribInt32(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "int64_t" )
    {
        int64_t data = (int64_t)_data;
        H5WriteStepAttribInt64(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "float" )
    {
        float data = (float)_data;
        H5WriteStepAttribFloat32(dataFile, _name.c_str(), &data, 1);
    }
    else if ( _dataType == "double" )
    {
        double data = (double)_data;
        H5WriteStepAttribFloat64(dataFile, _name.c_str(), &data, 1);
    }
    else 
        return -1;
}

inline int HDF5ParticleIO::writeTimestepAttributeArray(std::string _name, std::string _dataType, void * _data, int numElements)
{
    if ( _dataType == "int32_t" )
    {
        H5WriteStepAttribInt32(dataFile, _name.c_str(), (int32_t *)_data, numElements);
    }
    else if ( _dataType == "int64_t" )
    {
        H5WriteStepAttribInt64(dataFile, _name.c_str(), (int64_t *)_data, numElements);
    }
    else if ( _dataType == "float" )
    {
        H5WriteStepAttribFloat32(dataFile, _name.c_str(), (float *)_data, numElements);
    }
    else if ( _dataType == "double" )
    {
        H5WriteStepAttribFloat64(dataFile, _name.c_str(), (double *)_data, numElements);
    }
    else if ( _dataType == "string" )
    {
        H5WriteStepAttribString(dataFile, _name.c_str(), (char *)_data);
    }
    else 
        return -1;
}



inline int HDF5ParticleIO::writeTimestepAttributes()
{
    for (size_t i=0; i<timestepAttributes.size(); i++)
    {
        if ( timestepAttributes[i].dataType == "int32_t" )
            H5WriteStepAttribInt32(dataFile, (timestepAttributes[i].name).c_str(), (int32_t *)timestepAttributes[i].data, timestepAttributes[i].numElements);
        else if ( timestepAttributes[i].dataType == "int64_t" )
            H5WriteStepAttribInt64(dataFile, (timestepAttributes[i].name).c_str(), (int64_t *)timestepAttributes[i].data, timestepAttributes[i].numElements);
        else if ( timestepAttributes[i].dataType == "float" )
            H5WriteStepAttribFloat32(dataFile, (timestepAttributes[i].name).c_str(), (float *)timestepAttributes[i].data, timestepAttributes[i].numElements);
        else if ( timestepAttributes[i].dataType == "double" )
            H5WriteStepAttribFloat64(dataFile, (timestepAttributes[i].name).c_str(), (double *)timestepAttributes[i].data, timestepAttributes[i].numElements);
        else
            return -1;
    }

    timestepAttributes.clear();
    return 0;
}


inline int HDF5ParticleIO::writeVariables()
{
    for (size_t i=0; i<vars.size(); i++)
    {
        H5PartSetNumParticles(dataFile, vars[i].numElements);
        
        if ( vars[i].dataType == "float" )
            H5PartWriteDataFloat32(dataFile, (vars[i].name).c_str(), (float *)vars[i].data);
        else if ( vars[i].dataType == "double" )
            H5PartWriteDataFloat64(dataFile, (vars[i].name).c_str(), (double *)vars[i].data);
        else if ( vars[i].dataType == "int32_t" )
            H5PartWriteDataInt32(dataFile, (vars[i].name).c_str(), (int32_t *)vars[i].data);
        else if ( vars[i].dataType == "int64_t" )
            H5PartWriteDataInt64(dataFile, (vars[i].name).c_str(), (int64_t *)vars[i].data);
        else
            return -1;
    }

    vars.clear();
    return 0;
}


} // end Flecsi_Sim_IO namespace

#endif