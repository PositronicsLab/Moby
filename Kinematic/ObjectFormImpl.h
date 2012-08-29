#ifndef _OBJECT_FORM_IMPL_H_
#define _OBJECT_FORM_IMPL_H_

#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>
#include <Moby/Triangle.h>
#include "KinematicFormImpl.h"
#include "ObjectForm.h"
//Added by qt3to4:
#include <Q3BoxLayout>
#include <q3mainwindow.h>

class Q3BoxLayout;
class Q3ListBox;
class QMutex;
class ObjectFormImpl : public QMainWindow 
{
  friend class KinematicFormImpl;

    Q_OBJECT

  public:
    ObjectFormImpl(boost::shared_ptr<Moby::Simulator> s, QMutex* iv_mutex, void (*read_simulator)(const std::string&), std::list<std::list<Moby::Triangle> >* colliding_tris, QWidget* parent = 0);
    ~ObjectFormImpl();
    void set_simulator(boost::shared_ptr<Moby::Simulator> s);

  public slots:
    void body_selected();
    virtual void fileExit();
    virtual void fileOpen();
    virtual void fileSave();
    virtual void fileSaveAs();
    virtual void toggleCollidingTris();

  private:
    Moby::XMLTreePtr create_driver_node();
    void serialize_to_xml(const std::string& fname, Moby::BaseConstPtr object);
    void drawCollidingTris();
    void (*read_simulator)(const std::string& filename);
    std::list<std::list<Moby::Triangle> >* colliding_tris;
    KinematicFormImpl* _kf;
    QMutex* _mutex;
    boost::shared_ptr<Moby::Simulator> _sim;
    std::string _save_fname;
    Ui::MainWindow ui; 
};

#endif

