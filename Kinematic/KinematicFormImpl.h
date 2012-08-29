#ifndef _KINEMATIC_TRUE_H
#define _KINEMATIC_TRUE_H

#include <string>
#include <map>
#include <Moby/DynamicBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/CollisionDetection.h>
#include <Moby/RigidBody.h>
#include <boost/shared_ptr.hpp>
#include "KinematicForm.h"

class QMutex;
class QSlider;
class QWidget;
class QComboBox;
class QStatusBar;
class QwtDial;
class QwtCounter;
class ObjectFormImpl;
class Q3Table;

class KinematicFormImpl : public QDialog
{
  // Needed because of Qt signal / slot declarations
  Q_OBJECT

  public:
    KinematicFormImpl(boost::shared_ptr<Moby::DynamicBody> body, const std::list<boost::shared_ptr<Moby::CollisionDetection> > collision_detectors, ObjectFormImpl* oform, QMutex* iv_mutex, QWidget* parent = 0);
    ~KinematicFormImpl();

  public slots:
    void DOF_changed(double);
    void DOF_selected(const QString&);

  private:

    ObjectFormImpl* _oform;
    QMutex* _mutex;
    boost::shared_ptr<Moby::RCArticulatedBody> _abody;
    boost::shared_ptr<Moby::RigidBody> _rbody;
    std::list<boost::shared_ptr<Moby::CollisionDetection> > _coldet;
    QwtDial* _dof_dial;
    QwtCounter* _dof_counter;
    QComboBox* _combobox;
    Q3Table* _table;
    bool _update_DOFs;
    Ui::KinematicForm ui;
    double _alpha, _beta, _gamma;
};

#endif

