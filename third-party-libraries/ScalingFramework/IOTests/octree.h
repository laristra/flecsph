#ifndef _OCTREE_H_
#define _OCTREE_H_

#include <vector>
#include <list>
#include <iostream>

#include "simIO.h"

//
// Octree Encoding:
//
//	morton order (z order)
//  
//  .....
//   . 
//    .
//     .
//  .....
//  bottom left  ->  bottom right -> top left  -> top right
//
//
//  010 011   bottom
//  000 001
//
//  110 111   top
//  100 101
//
//
//  Levels:
//	0:		  X				
//	1:	X   X   X   X		  010 (3 bits)
//  2:        ...         000 010 (6 bits)
 
// Directions:
// 0:left  1:right  2:bottom  3:top  4:front  5:back



class Octree
{
	int maxLevel;
	std::list<std::string> node;

  public:
  	static const int BIT_STRING_SIZE = 3;
  	enum direction{LEFT=0, RIGHT=1, BOTTOM=2, TOP=3, FRONT=4, BACK=5};

  	Octree(){  node.clear(); maxLevel=0; }
  	~Octree(){ node.clear(); maxLevel=0; }

  	void init();
  	
  	void getSpatialExtent(std::string _node, float array[]);
  	int findNeighborNode(std::string _node, direction _dir, std::string resultNode);

  	int buildTree(int _numLevels);
  	int splitRootNode();
  	int splitNode(std::string _node="");
  	void splitNode(std::list<std::string>::iterator it, std::string _node);


  	int getMaxLevel() { return maxLevel; }
  	std::string getNodeAt(size_t pos);
  	int getNodeLevel(std::string _node){ return _node.size()/BIT_STRING_SIZE; }

  	std::string print();
};



void Octree::init()
{
	maxLevel = 0;
}


int Octree::buildTree(int _numLevels)
{
	int numNodes = 8;
	splitNode();	// creates the first 8 nodes for level 1


	for (int _lvl=1; _lvl<_numLevels; _lvl++)
	{

		auto it=node.begin();
		while (it != node.end())
		{
			std::string nodeName = *it;
			splitNode(it, nodeName);

			// Skip 8 nodes ahead to account for creating 8 children of a node
			for (int i=0; i<8; i++)
				++it;
		}
	}

	return numNodes;
}


int Octree::splitRootNode()
{
	auto it = node.begin();
	node.insert(it, "000" ); 	node.insert(it, "001" );
	node.insert(it, "010" ); 	node.insert(it, "011" );
	 
	node.insert(it, "100" ); 	node.insert(it, "101" ); 
	node.insert(it, "110" ); 	node.insert(it, "111" );
	
	maxLevel = 1;
}


void Octree::splitNode(std::list<std::string>::iterator it, std::string _node)
{
	// Ease it from the octree
	auto itErase = it;	it++;
	node.erase(itErase);

	// Insert children
	node.insert(it, "000" + _node ); 	node.insert(it, "001" + _node ); 
	node.insert(it, "010" + _node ); 	node.insert(it, "011" + _node ); 
	
	node.insert(it, "100" + _node ); 	node.insert(it, "101" + _node );
	node.insert(it, "110" + _node ); 	node.insert(it, "111" + _node );
	
	int nodeLevel = (_node.size())/BIT_STRING_SIZE + 1;
	maxLevel = std::max(nodeLevel, maxLevel);
}


int Octree::splitNode(std::string _node)
{
	// If first 8 nodes
	if (_node == "" && node.size() == 0)
	{
		splitRootNode();
		return 0;
	}


	// Navigate to the node
	auto it = node.begin();
	while (*it != _node)
		it++;	

	// Leave if node is not found
	if (it == node.end())
		return -1;

	// Split the node found
	splitNode(it, _node);


	return 0;
}


void Octree::getSpatialExtent(std::string _node, float array[])
{
	int nodeLevel = (_node.size())/BIT_STRING_SIZE;

	array[0] = 0;	array[1] = 1.0;		// x: min, max
	array[2] = 0;	array[3] = 1.0;		// y: min, max
	array[4] = 0;	array[5] = 1.0;		// z: min, max

	std::string _nodeCopy = _node;
	int _level = 1;
	while (_level < nodeLevel+1)
	{
		// z: front or back
		if (_node[0] == '0')	
			array[5] = (array[5] - array[4])/2.0 + array[4];
		else
			array[4] = (array[5] - array[4])/2.0 + array[4];

		// y: top or bottom
		if (_node[1] == '0')
			array[3] = (array[3] - array[2])/2.0 + array[2];
		else
			array[2] = (array[3] - array[2])/2.0 + array[2];

		//x: left or right
		if (_node[2] == '0')
			array[1] = (array[1] - array[0])/2.0 + array[0];
		else
			array[0] = (array[1] - array[0])/2.0 + array[0];


		_nodeCopy = _node.substr(BIT_STRING_SIZE);
		_node = _nodeCopy;

		_level++;
	}
}


std::string Octree::getNodeAt(size_t pos)
{
	size_t count = 0;
	for (auto it : node)
	{
		if (pos == count)
			return it;

		count++;
	}
}


std::string Octree::print()
{
	std::string outputString = "Num nodes: " + std::to_string(node.size()) + ", max level: " + std::to_string(maxLevel) + "\n";

	for (auto it=node.begin(); it!=node.end(); it++)
		outputString += *it + "  ";

	outputString += "\n\n";

	return outputString;
}

// 10 11
// 00 01


//  010 011   bottom
//  000 001
//
//  110 111   top
//  100 101

int Octree::findNeighborNode(std::string _node, direction _dir, std::string resultNode)
{
	int strLength = _node.size();

	if (_dir ==  LEFT)
		if (_node[strLength-1] == '1')
		{
			resultNode[strLength-1] = '0';
			return 0;
		}

	if (_dir ==  RIGHT)
		if (_node[strLength-1] == '0')
		{
			resultNode[strLength-1] = '1';
			return 0;
		}

	
	if (_dir ==  TOP)
		if (_node[strLength-2] == '0')
		{
			resultNode[strLength-2] = '1';
			return 0;
		}

	if (_dir ==  BOTTOM)
		if (_node[strLength-2] == '1')
		{
			resultNode[strLength-2] = '0';
			return 0;
		}


	if (_dir ==  TOP)
		if (_node[strLength-2] == '0')
		{
			resultNode[strLength-2] = '1';
			return 0;
		}

	if (_dir ==  BOTTOM)
		if (_node[strLength-2] == '1')
		{
			resultNode[strLength-2] = '0';
			return 0;
		}
}

#endif
