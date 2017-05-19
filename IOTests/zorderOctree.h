#ifndef _Z_ORDER_OCTREE_H_
#define _Z_ORDER_OCTREE_H_

#include <vector>
#include <list>
#include <iostream>


//
// Oncteee Encoding:
//
//  001 011   bottom
//  000 010
//
//  101 111   top
//  100 110
//
//
//  Levels:
//	0:		  X				
//	1:	X   X   X   X		  000 (3 bits)
//  2:        ...         000 000 (6 bits)


const int BIT_STRING_SIZE = 3;

std::string toBinary(int num)
{
    std::string binStr = "";

    while (num>0)
	{
	    binStr = std::to_string(num % 2) + binStr;
	    num /= 2;
	}

	return binStr;
}


int toDecimal(std::string binary)
{
	int num = 0;
	int pow = 1;

	for (int i=binary.size()-1; i>=0; i--)
	{
		num += ( (binary[i] == '0') ? 0 : 1 ) * pow;
		pow *= 2;
	}

	return num;
}


class ZorderOctree
{
	int maxLevel;
	std::list<int> code;
	std::list<int> level;


  public:
  	ZorderOctree(){ initEmpty(); };
  	~ZorderOctree(){ initEmpty(); };

  	void initEmpty();
  	void splitRoot();
  	void splitNode(int _node);

  	void printOctree();
};


void ZorderOctree::initEmpty()
{
	maxLevel = 0;
	code.clear();
	level.clear();
}


void ZorderOctree::splitRoot()
{
	initEmpty();
	maxLevel = 1;
	for (int i=0; i<8; i++)
	{
		code.push_front(i);
		level.push_front(1);
	}
}


void ZorderOctree::splitNode(int _node)
{
	auto itNode = code.begin();
	auto itLevel = level.begin();

	while (*itNode != _node)
	{
		itNode++;	
		itLevel++;
	}
	

	//
	// Make the string size correct
	std::string binString = toBinary(_node);
	int stringSize 		= binString.size();
	int bitStringSize 	= nodeLevel * BIT_STRING_SIZE;
	while (stringSize < bitStringSize)
	{
		binString = '0' + binString;
		stringSize++;
	}


	//
	// Do the split
	std::vector< std::string> octants;

	octants.push_back( "000" + binString );
	octants.push_back( "001" + binString );
	octants.push_back( "010" + binString );
	octants.push_back( "011" + binString );
	octants.push_back( "100" + binString );
	octants.push_back( "101" + binString );
	octants.push_back( "110" + binString );
	octants.push_back( "111" + binString );
	

	// Insert in the correct location
	code.remove(_node);
	code.insert(itNode, toDecimal(octants[0]) );	itNode++;
	code.insert(itNode, toDecimal(octants[1]) );	itNode++;
	code.insert(itNode, toDecimal(octants[2]) );	itNode++;
	code.insert(itNode, toDecimal(octants[3]) );	itNode++;
	code.insert(itNode, toDecimal(octants[4]) );	itNode++;
	code.insert(itNode, toDecimal(octants[5]) );	itNode++;
	code.insert(itNode, toDecimal(octants[6]) );	itNode++;
	code.insert(itNode, toDecimal(octants[7]) );

	level.remove(nodeLevel);
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );	itLevel++;
	level.insert(itLevel, nodeLevel+1 );

	maxLevel = std::max(nodeLevel+1, maxLevel);
}


void ZorderOctree::printOctree()
{
	std::cout << "Num nodes: " << code.size() << std::endl;

	for (int i=0; i<maxLevel; i++)
	{
		for (int j=0; j<level; j++)
			if (level[j] == i)
				std::cout << code[j];
			else
				std::cout << "  ";
		
		std::cout << "\n";
	}
}


#endif
