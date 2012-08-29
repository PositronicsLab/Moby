#include <errno.h>
#include "ObjectFormImpl.h"
#include <qthread.h>
#include <qmutex.h>
#include <qlabel.h>
#include <qaction.h>
#include <qmessagebox.h>
#include <q3listbox.h>
#include <q3filedialog.h>
#include <qmainwindow.h>
#include <Moby/DynamicBody.h>
#include <Moby/EventDrivenSimulator.h>
#include <Moby/XMLWriter.h>

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::list;
using std::map;
using std::set;
using std::make_pair;
using std::pair;
using std::stack;
using namespace Moby;

extern Vector2 g_size, g_loc;
extern Vector3 g_up, g_target, g_position;

ObjectFormImpl::ObjectFormImpl(SimulatorPtr sim, QMutex* iv_mutex, void (*read_simulator)(const std::string& filename), std::list<std::list<Triangle> >* colliding_tris, QWidget* parent) : QMainWindow(parent)
{
  // setup the UI
  ui.setupUi(this);

  // init the kinematic form
  _kf = NULL;

  // save the mutex
  _mutex = iv_mutex;

  // save the function pointer
  this->read_simulator = read_simulator;

  // save pointer to colliding triangle list
  this->colliding_tris = colliding_tris;

  // set the simulator
  set_simulator(sim);

  // connect the listbox selection to the slot
  connect(ui._lb, SIGNAL(selectionChanged()), this, SLOT(body_selected()));
  connect( (QObject*) ui.actionOpen, SIGNAL( activated() ), this, SLOT( fileOpen() ) );
  connect( (QObject*) ui.actionSave, SIGNAL( activated() ), this, SLOT( fileSave() ) );
  connect( (QObject*) ui.actionSave_As, SIGNAL( activated() ), this, SLOT( fileSaveAs() ) );
  connect( (QObject*) ui.actionQuit, SIGNAL( activated() ), this, SLOT( fileExit() ) );
  connect((QObject*) ui.actionDraw_Colliding_Tris, SIGNAL(activated()), this, SLOT(toggleCollidingTris()));

  // sort the items in the list box
  ui._lb->sort();
}

ObjectFormImpl::~ObjectFormImpl()
{
  if (_kf)
    delete _kf;
}

/// sets up the XML attributes for the driver node
XMLTreePtr ObjectFormImpl::create_driver_node()
{
  // create the top level node
  XMLTreePtr node(new XMLTree("Driver"));

  // create the camera node
  XMLTreePtr camera_node(new XMLTree("camera"));
  camera_node->attribs.insert(XMLAttrib("position", g_position));
  camera_node->attribs.insert(XMLAttrib("target", g_target));
  camera_node->attribs.insert(XMLAttrib("up", g_up));

  // create the window node
  XMLTreePtr window_node(new XMLTree("window"));
  window_node->attribs.insert(XMLAttrib("location", g_loc));
  window_node->attribs.insert(XMLAttrib("size", g_size));

  // setup the tree
  node->add_child(camera_node);
  node->add_child(window_node);
  return node;
}

/// Serializes the given object (and all of its dependencies) to XML
void ObjectFormImpl::serialize_to_xml(const std::string& fname, BaseConstPtr object)
{
  // get the filename
  std::string filename = fname;

  // get the current path
  size_t BUFSIZE = 8192;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "XMLWriter::serialize_to_xml() - unable to allocate sufficient memory!" << std::endl;
      return;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != std::string::npos)
  {
    // get the new working path
    std::string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,std::string::npos);
  }

  // create a new XMLTree
  XMLTreePtr topnode(new XMLTree("XML"));

  // create a node for driver
  XMLTreePtr driver_node = create_driver_node();

  // create a node for Moby
  XMLTreePtr node(new XMLTree("Moby"));
  topnode->add_child(driver_node);
  topnode->add_child(node);

  // setup a list of shared objects
  std::list<BaseConstPtr> shared_objects;

  // add the desired object to the list of shared objects
  shared_objects.push_back(object);

  // init a set of serialized objects
  std::set<BaseConstPtr> serialized;
  
  // develop the XML tree until there is nothing more to serialize
  while (!shared_objects.empty())
  {
    // get the object off of the front of the queue
    BaseConstPtr obj = shared_objects.front();
    assert(obj);
    shared_objects.pop_front();

    // if this object has already been serialized, skip it
    if (serialized.find(obj) != serialized.end())
      continue;

    // create a new node for this object under the parent
    XMLTreePtr new_node(new XMLTree(""));
    node->add_child(new_node);

    // serialize to this new node
    obj->save_to_xml(new_node, shared_objects);

    // indicate that the node has been serialized
    serialized.insert(obj);
  }

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the tree to the file
  out << *topnode << std::endl;

  // close the file
  out.close();

  // change back to the initial working directory
  chdir(cwd.get());
}

/// Toggles drawing colliding triangles
void ObjectFormImpl::drawCollidingTris()
{
  // determine whether we activate or deactivate the feature
  if (ui.actionDraw_Colliding_Tris->isOn())
  {
    // make sure that this is a EventDrivenSimulator
    shared_ptr<EventDrivenSimulator> cs = dynamic_pointer_cast<EventDrivenSimulator>(_sim); 

    // if it's not a EventDrivenSimulator, disable the toggle
    if (!cs)
    {
      ui.actionDraw_Colliding_Tris->setOn(false);
      return;
    }

    // clear the list of colliding triangles
    colliding_tris->clear();

    // setup a map of triangle pointers to triangles
    typedef pair<const IndexedTriArray*, unsigned> UniqueTri;
    map<UniqueTri, Triangle> tris;

    // we want a spanning edges of all triangles that are colliding 
    map<UniqueTri, set<UniqueTri> > edges;

    // get colliding triangles, if any
    BOOST_FOREACH(shared_ptr<CollisionDetection> coldet, cs->collision_detectors)
    {
      _mutex->lock();
      coldet->mode = CollisionDetection::eAllContacts;
      coldet->is_collision();
      list<CollidingTriPair> ctris = coldet->colliding_tris;
      BOOST_FOREACH(const CollidingTriPair& cpair, ctris)
      {
        // get the geometries
        CollisionGeometryPtr g1 = cpair.geom1;
        CollisionGeometryPtr g2 = cpair.geom2;
        if (!g1 || !g2)
          continue;

        // get the two triangles
        Triangle t1 = cpair.mesh1->get_triangle(cpair.tri1);
        Triangle t2 = cpair.mesh2->get_triangle(cpair.tri2);

        // transform the triangles
        t1 = Triangle::transform(t1, g1->get_transform());
        t2 = Triangle::transform(t2, g2->get_transform());

        // save the triangles
        tris[make_pair(cpair.mesh1, cpair.tri1)] = t1;
        tris[make_pair(cpair.mesh2, cpair.tri2)] = t2;

        // add connections between the triangles
        edges[make_pair(cpair.mesh1, cpair.tri1)].insert(make_pair(cpair.mesh2, cpair.tri2));
        edges[make_pair(cpair.mesh2, cpair.tri2)].insert(make_pair(cpair.mesh1, cpair.tri1));
      }
      _mutex->unlock();
    }

    // setup triangles that should be in each group
    set<pair<const IndexedTriArray*, unsigned> > processed;
    while (!edges.empty())
    {
      // get first node from graph
      stack<UniqueTri> s;
      s.push(edges.begin()->first);

      // add a new list
      colliding_tris->push_back(list<Triangle>());

      // process all nodes connected to this node using DFS
      while (!s.empty())
      {
        // get the node
        UniqueTri node = s.top();
        s.pop();

        // add all nodes connected to this one for processing
        if (edges.find(node) == edges.end())
          continue; // node has already been processed
        else
          BOOST_FOREACH(const UniqueTri& ut, edges[node])
          {
            // we don't want to process node again
            if (edges.find(ut) != edges.end())
              edges[ut].erase(node);
            s.push(ut);
          }

        // remove the edges emanating from this node
        edges.erase(node);

        // push the colliding triangle to the list
        colliding_tris->back().push_back(tris[node]);
      }
    }
  }
  else
  {
    // don't draw any colliding triangles
   _mutex->lock();
   colliding_tris->clear();
   _mutex->unlock();
  }
}

/// Quits
void ObjectFormImpl::fileExit()
{
  if (QMessageBox::question(this, tr("Confirmation"), tr("Do you really want to exit?"), QMessageBox::Ok, QMessageBox::Cancel, QMessageBox::NoButton) == QMessageBox::Ok)
    exit(0); 
}

/// Opens a simulator file
void ObjectFormImpl::fileOpen()
{
  // open the file dialog
  QString s = Q3FileDialog::getOpenFileName(".", "XML files (*.xml)", this, 
"Open File", "Choose a file");

  // make sure a file was selected
  if (s == "")
    return;

  // attempt to open the file
  std::cout << "-- opening " << s.ascii() << std::endl;
  (*read_simulator)(std::string(s.ascii()));

  // clear the save file name
  _save_fname = "";
}

/// Saves the file to a new filename
void ObjectFormImpl::fileSaveAs()
{
  // make sure that we have a simulator open
  if (!_sim)
  {
    std::cerr << "-- no XML file opened!  nothing to save to..." << std::endl;
    return;
  }

  // save the filename used to save
  QString s = Q3FileDialog::getSaveFileName(".", "XML files (*.xml)", this, "Save file", "Choose a filename");

  // make sure a file was specified
  if (s == "")
    return;

  // see whether the file already exists and user wants to overwrite it
  if (QFile::exists(s) && QMessageBox::question(this, tr("Overwrite file?"), tr("The file already exists.  Overwrite it?"), QMessageBox::Ok, QMessageBox::Cancel, QMessageBox::NoButton) == QMessageBox::Cancel)
    return;

  // save this filename as the last used
  _save_fname = std::string(s.ascii());

  // save
  std::cout << "-- saved simulator to " << _save_fname << std::endl;
  serialize_to_xml(_save_fname, _sim);
}

/// Saves the file to the last used filename for saving (if any)
void ObjectFormImpl::fileSave()
{
  // make sure that we have a simulator open
  if (!_sim)
  {
    std::cerr << "-- no XML file opened!  nothing to save to..." << std::endl;
    return;
  }

  // if there is no saved file, use save as
  if (_save_fname == "")
    fileSaveAs();
  else
  {
    serialize_to_xml(_save_fname, _sim);
    std::cout << "-- saved simulator to " << _save_fname << std::endl;
  }
}

/// Called when the colliding tris option is toggled
void ObjectFormImpl::toggleCollidingTris()
{
  drawCollidingTris();
}

/// Sets the simulator (if any)
void ObjectFormImpl::set_simulator(shared_ptr<Simulator> sim)
{
  // store the pointer to the simulator
  _sim = sim;

  // clear all elements from the list box
  ui._lb->clear();

  // get the bodies in the simulator (if any)
  if (sim)
  {
    const std::list<DynamicBodyPtr>& bodies = sim->get_dynamic_bodies();
    BOOST_FOREACH(DynamicBodyPtr body, bodies)
      ui._lb->insertItem(QString(body->id.c_str()));
  }

  // sort the items in the list box
  ui._lb->sort();
}

void ObjectFormImpl::body_selected()
{
  // kill any launched forms...
  if (_kf)
  {
    _kf->hide();
    delete _kf;
    _kf = NULL;
  }

  // verify that a body is selected
  Q3ListBoxItem* lbi = ui._lb->selectedItem();
  if (!lbi)
    return;

  // see whether the simulator is a contact simulator
  shared_ptr<EventDrivenSimulator> csim = dynamic_pointer_cast<EventDrivenSimulator>(_sim);
  list<shared_ptr<CollisionDetection> > coldets;
  if (csim)
    coldets = csim->collision_detectors;

  // get the body id()
  const std::string& id(lbi->text().ascii());

  // look for the body in the simulator
  const std::list<DynamicBodyPtr>& bodies = _sim->get_dynamic_bodies();
  BOOST_FOREACH(DynamicBodyPtr body, bodies)
    if (body->id == id)
    {
      _kf = new KinematicFormImpl(body, coldets, this, _mutex);
      _kf->show();
      return;
    }

  // should not be here!
  assert(false);
}

