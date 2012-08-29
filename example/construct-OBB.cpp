#include <sstream>
#include <iostream>
#include <fstream>
#include <Moby/XMLReader.h>
#include <Moby/XMLTree.h>
#include <Moby/OBB.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/VExactTOIFinder.h>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::cerr;
using std::endl;
using std::map;

int main(int argc, char* argv[])
{
  // check that syntax ok
  if (argc < 2)
  {
    cerr << "syntax: construct-obb <xml file>" << endl;
    return -1;
  }

  // load the simulation
  map<std::string, BasePtr> read_map = XMLReader::read(std::string(argv[argc-1]));

  // get the VExactTOIFinder
  shared_ptr<VExactTOIFinder> coldet;
  for (map<std::string, BasePtr>::const_iterator i = read_map.begin(); i != read_map.end(); i++)
    if ((coldet = dynamic_pointer_cast<VExactTOIFinder>(i->second)))
      break;
  if (!coldet)
  {
    cerr << "no VExactTOIFinder found!" << endl;
    return -1;
  }
  
  // get the OBBs
  const map<CollisionGeometryPtr, OBBPtr>& obbs = coldet->get_OBBs();

  // write the OBBs
  for (map<CollisionGeometryPtr, OBBPtr>::const_iterator i = obbs.begin(); i != obbs.end(); i++)
  {
    // create the filename
    std::ostringstream fname;
    fname << i->first->id << ".obb.xml";

    // open the filename for writing
    std::ofstream out(fname.str().c_str());

    // save the OBB to an XML tree
    XMLTreePtr tree = i->second->save_to_xml_tree();
    out << *tree << endl;

    // close the file
    out.close();
  }
}

