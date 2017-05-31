#ifndef _CELL_SET_H_
#define _CELL_SET_H_

class CellSet
{
	enum CellSetType{ point=0, structuredGrid=1}
};


class StructuredGrid: public CellSet
{
  public:
	std::vector<int> dims;          // e.g. 8 x 8 x 4
	std::vector<float> extents;     // minX, maxX   minY, maxY,   minZ, maxZ
};


class UnstructuredGrid: public CellSet
{
  public:
	int numPoints;
};


class Point: public CellSet
{
  public:
	int numPoints;
};

#endif