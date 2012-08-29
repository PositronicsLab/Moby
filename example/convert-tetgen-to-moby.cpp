#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    std::cerr << "syntax: convert-tetgen-to-moby <node file> <element file> <output file>" << std::endl;
    return -1; 
  }

  // open the node file
  std::ifstream in_node(argv[1]);
  if (in_node.fail())
  {
    std::cerr << " -- unable to open " << argv[1] << " for reading!" << std::endl;
    return -1;
  }

  // open the element file
  std::ifstream in_ele(argv[2]);
  if (in_ele.fail())
  {
    std::cerr << " -- unable to open " << argv[2] << " for reading!" << std::endl;
    return -1;
  }

  // open the tetra file for writing
  std::ofstream out(argv[3]);
  if (out.fail())
  {
    std::cerr << " -- unable to open " << argv[3] << " for writing!" << std::endl;
    return -1;
  }
 
  // read in the number of nodes
  unsigned nnodes, ncols;
  in_node >> nnodes;
  in_node >> ncols;
  if (ncols != 3)
  {
    std::cerr << " -- format of " << argv[1] << " is incorrect; terminating!" << std::endl;
    return -1;
  }

  // read in and discard garbage
  unsigned garbage;
  in_node >> garbage;
  if (garbage != 0)
  {
    std::cerr << " -- format of " << argv[1] << " is incorrect; terminating!" << std::endl;
    return -1;
  }
  in_node >> garbage;
  if (garbage != 0)
  {
    std::cerr << " -- format of " << argv[1] << " is incorrect; terminating!" << std::endl;
    return -1;
  }

  // read in vertices, simultaneously writing them out
  for (unsigned i=0; i< nnodes; i++)
  {
    // read and discard the vertex index
    in_node >> garbage;
    double number;
    in_node >> number;
    out << "v " << number << " ";
    in_node >> number;
    out << number << " ";
    in_node >> number;
    out << number << std::endl;
  }

  // close the node file
  in_node.close();

  // read the number of elements
  unsigned nele;
  in_ele >> nele;
  in_ele >> ncols;
  if (ncols != 4)
  {
    std::cerr << " -- format of " << argv[2] << " is incorrect; terminating!" << std::endl;
    return -1;
  }

  // read in and discard garbage
  in_ele >> garbage;
  if (garbage != 0)
  {
    std::cerr << " -- format of " << argv[2] << " is incorrect; terminating!" << std::endl;
    return -1;
  }

  // read in tetrahedra, simultaneously writing them out
  for (unsigned i=0; i< nele; i++)
  {
    // read and discard the tetrahedron index
    in_ele >> garbage;
    unsigned index;
    in_ele >> index;
    out << "t " << index << " ";
    in_ele >> index;
    out << index << " ";
    in_ele >> index;
    out << index << " ";
    in_ele >> index;
    out << index << std::endl;
  }

  // close files
  in_ele.close();
  out.close();
}

