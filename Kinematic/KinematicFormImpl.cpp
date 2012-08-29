#include "KinematicFormImpl.h"
#include "ObjectFormImpl.h"
#include <iostream>
#include <cmath>
#include <Moby/Joint.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/CollisionDetection.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/Quat.h>
#include <qmutex.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qslider.h>
#include <qdial.h>
#include <q3table.h>
#include <qcombobox.h>
#include <qwt/qwt_counter.h>
#include <qwt/qwt_slider.h>
#include <qwt/qwt_dial.h>
#include <qwt/qwt_dial_needle.h>

using namespace Moby;
using boost::shared_ptr;

// converts euler angles to a rotation matrix
Matrix3 ea_to_matrix(Real alpha, Real beta, Real gamma)
{
  const unsigned X = 0, Y = 1, Z = 2;

  Real ch = std::cos(alpha);
  Real sh = std::sin(alpha);
  Real ca = std::cos(beta);
  Real sa = std::sin(beta);
  Real cb = std::cos(gamma);
  Real sb = std::sin(gamma);

  Matrix3 R;
  R(X,X) = ch*ca;   R(X,Y) = -ch*sa*cb + sh*sb;   R(X,Z) = ch*sa*sb + sh*cb;
  R(Y,X) = sa;      R(Y,Y) = ca*cb;               R(Y,Z) = -ca*sb;
  R(Z,X) = -sh*ca;  R(Z,Y) = sh*sa*cb + ch*sb;    R(Z,Z) = -sh*sa*sb + ch*cb; 

  return R;
}

// converts a rotation matrix to euler angles
void matrix_to_ea(const Matrix3& R, Real& alpha, Real& beta, Real& gamma)
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (R(Y,X) > (Real) 1.0 - NEAR_ZERO)
  {
    alpha = std::atan2(R(X,Z), R(Z,Z));
    beta = M_PI_2;
    gamma = (Real) 0.0;
  }
  else if (R(Y,X) < (Real) 1.0 + NEAR_ZERO)
  {
    alpha = std::atan2(R(X,Z),R(Z,Z));
    beta = -M_PI_2;
    gamma = (Real) 0.0;
  }
  else
  {
    alpha = std::atan2(-R(Z,X),R(X,X));
    beta = std::atan2(-R(Y,Z),R(Y,Y));
    gamma = std::asin(R(Y,X));
  }
}

KinematicFormImpl::KinematicFormImpl(DynamicBodyPtr body, const std::list<shared_ptr<CollisionDetection> > collision_detectors, ObjectFormImpl* oform, QMutex* iv_mutex, QWidget* parent) : QDialog(parent) 
{
  // setup the ui
  ui.setupUi(this);

  // store the object form
  _oform = oform;

  // store the mutex
  _mutex = iv_mutex;

  // copy list of collision detectors 
  _coldet = collision_detectors;

  // indicate that we want to update DOFs
  _update_DOFs = true;

  // determine whether the body is articulated or rigid
  _abody = boost::dynamic_pointer_cast<RCArticulatedBody>(body);
  if (!_abody)
    _rbody = boost::dynamic_pointer_cast<RigidBody>(body);
  assert(_abody || _rbody);

  // create a combo box to hold everything -- make it read-only
  _combobox = new QComboBox(false, this);
  _combobox->setGeometry(QRect(10,10,200,35));

  // create a dial for changing DOFs
  _dof_dial = new QwtDial(this);
  _dof_dial->setGeometry(QRect(300,10,250,250));
  QwtDialSimpleNeedle* needle = new QwtDialSimpleNeedle(QwtDialSimpleNeedle::Arrow);
  _dof_dial->setNeedle(needle);
  _dof_dial->setLineWidth(10);
  _dof_dial->setScale(12,36,0.0);

  // create a counter for precise control 
  _dof_counter = new QwtCounter(this);
  _dof_counter->setGeometry(QRect(10,130,200,35));
  _dof_counter->setNumButtons(3);
  _dof_counter->setIncSteps(QwtCounter::Button1, 1);
  _dof_counter->setIncSteps(QwtCounter::Button2, 10);
  _dof_counter->setIncSteps(QwtCounter::Button3, 100);

  // create a table for displaying pairs in collision
  _table = new Q3Table(this);
  _table->setGeometry(QRect(10,280,width()-10,height()-300));
  _table->setVScrollBarMode(Q3ScrollView::Auto);
  _table->verticalHeader()->hide();
  _table->setNumCols(2);
  _table->horizontalHeader()->setLabel(0, QString("Body1"));
  _table->horizontalHeader()->setLabel(1, QString("Body2"));
  _table->setColumnWidth(0,(width()-20)/2 - 20);
  _table->setColumnWidth(1,(width()-20)/2 - 20);

  // add any joints 
  QStringList jointids;
  if (_abody)
  {
    const std::vector<JointPtr>& joints = _abody->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
    {
      // if the joint has exactly one DOF, add it unchanged
      if (joints[i]->num_dof() == 1)
        jointids.append(QString(joints[i]->id.c_str()));
      // ignore fixed joints
      else if (joints[i]->num_dof() > 1)
      {
        // add a string for each DOF
        assert(joints[i]->num_dof() < 10);
        for (unsigned j=0; j< joints[i]->num_dof(); j++)
        {
          char buffer[5];
          sprintf(buffer," [%u]", j);
          QString s = QString(joints[i]->id.c_str()) + QString(buffer);
          jointids.append(s);
        }
      }
    }
  }

  // sort the list
  jointids.sort();

  // add the items to the combo box
  _combobox->insertStringList(jointids);

  // add free controls
  if (!_abody || _abody->is_floating_base())
  {
    _combobox->insertItem(QString("z-rotation"), 0);
    _combobox->insertItem(QString("y-rotation"), 0);
    _combobox->insertItem(QString("x-rotation"), 0);
    _combobox->insertItem(QString("z-translation"), 0);
    _combobox->insertItem(QString("y-translation"), 0);
    _combobox->insertItem(QString("x-translation"), 0);
  }

  // connect the dial and the combobox to this
  QObject::connect(_combobox, SIGNAL(activated(const QString&)), this, SLOT(DOF_selected(const QString&)));
  QObject::connect(_dof_counter, SIGNAL(valueChanged(double)), this, SLOT(DOF_changed(double)));
  QObject::connect(_dof_dial, SIGNAL(valueChanged(double)), this, SLOT(DOF_changed(double)));

  // setup the Euler angles for free movement..
  Matrix3 R;
  RigidBodyPtr rbody = _rbody;
  if (!rbody)
    rbody = RigidBodyPtr(_abody->get_base_link());
  rbody->get_transform().get_rotation(&R);
  matrix_to_ea(R, _alpha, _beta, _gamma);

  // select the first item in the list box
  _combobox->setCurrentItem(0);
  DOF_selected(_combobox->currentText());
}

KinematicFormImpl::~KinematicFormImpl()
{
  delete _dof_dial;
  delete _dof_counter;
  delete _combobox;
}

/// Called when a DOF in the list box is selected
void KinematicFormImpl::DOF_selected(const QString& dof)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the name of the DOF 
  std::string DOF(dof.ascii());

  // DOF index is zero by default, but change it if present between []
  unsigned dof_idx = 0;
  if (DOF.find('[') != std::string::npos)
  {
    dof_idx = DOF[DOF.find('[')+1] - '0';
    DOF = DOF.substr(0,DOF.find('[')-1);
  }

  // indicate we don't want to update any DOFs
  _update_DOFs = false;

  // look for one of the joints, if it is an articulated body
  if (_abody)
  {
    const std::vector<JointPtr>& joints = _abody->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
      if (joints[i]->id == DOF)
      {  
        // get the lower and upper limits
        Real lolimit = joints[i]->lolimit[dof_idx];
        Real hilimit = joints[i]->hilimit[dof_idx];

        // set a range on the lower and upper limits
        lolimit = std::max((Real) -1000.0, lolimit);
        hilimit = std::min((Real) 1000.0, hilimit);

        // set the limits on the dial
        _dof_dial->setRange(lolimit,hilimit,.01,10);
        _dof_counter->setRange(lolimit,hilimit,.001);

        // get the current value of the dial
        Real q = joints[i]->q[dof_idx];

        // update the value on the dial (counter will be updated automatically) 
        _dof_dial->setValue(q);

        // indicate we now want to update any DOFs
        _update_DOFs = true;

        // repainting is necessary for dials...
        _dof_dial->repaint();  
      }
  }

  // check for a free-body control
  RigidBodyPtr rbody = (_rbody) ? _rbody : RigidBodyPtr(_abody->get_base_link());
  if (DOF == "x-translation" || DOF == "y-translation" || DOF == "z-translation")
  {
    _dof_dial->setRange(-10.0,10.0,0.001,10);
    _dof_counter->setRange(-10.0,10.0,0.001);
    if (DOF[0] == 'x')
      _dof_dial->setValue(rbody->get_position()[X]);    
    else if (DOF[0] == 'y')
      _dof_dial->setValue(rbody->get_position()[Y]);    
    else
      _dof_dial->setValue(rbody->get_position()[Z]);    
    _dof_dial->repaint();  
  }
  else if (DOF == "x-rotation" || DOF == "y-rotation" || DOF == "z-rotation")
  {
    _dof_dial->setRange(-M_PI,M_PI,0.1,10);
    _dof_counter->setRange(-M_PI,M_PI,0.01);
    if (DOF[0] == 'x')
      _dof_dial->setValue(_alpha);    
    else if (DOF[0] == 'y')
      _dof_dial->setValue(_beta);    
    else
      _dof_dial->setValue(_gamma);
    _dof_dial->repaint();  
  }

  // update the DOF value
  DOF_changed(_dof_dial->value());

  // indicate we now want to update any DOFs
  _update_DOFs = true;
}

/// Method for handling when one of the DOF sliders has been changed
void KinematicFormImpl::DOF_changed(double value)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // change the value on the other control (if hasn't already been done)
  if (_dof_counter->value() != value)
  {
    _dof_counter->setValue(value);
    return;
  }
  else if (_dof_dial->value() != value)
  {
    _dof_dial->setValue(value);
    return;
  }

  // the following is only for updating DOFs
  if (_update_DOFs)
  {
    // get the name of the DOF
    std::string DOF(_combobox->currentText().ascii());

    // DOF index is zero by default, but change it if present between []
    unsigned dof_idx = 0;
    if (DOF.find('[') != std::string::npos)
    {
      dof_idx = DOF[DOF.find('[')+1] - '0';
      DOF = DOF.substr(0,DOF.find('[')-1);
    }

    // look for one of the joints, if it is an articulated body
    if (_abody)
    {
      const std::vector<JointPtr>& joints = _abody->get_joints();
      for (unsigned i=0; i< joints.size(); i++)
        if (joints[i]->id == DOF)
        {
          // set the joint value
          VectorN& q = joints[i]->q;
          q[dof_idx] = value;

          // update the display
          _mutex->lock();
          _abody->update_link_transforms();
          _mutex->unlock();
        }
    }

    // look for free-body control 
    RigidBodyPtr rbody = (_rbody) ? _rbody : RigidBodyPtr(_abody->get_base_link());
    if (DOF == "x-translation" || DOF == "y-translation" || DOF == "z-translation")
    {
      Vector3 x = rbody->get_position();
      if (DOF[0] == 'x')
        x[X] = value;    
      else if (DOF[0] == 'y')
        x[Y] = value;  
      else
        x[Z] = value;

      // update the body
      _mutex->lock();
      rbody->set_position(x);
      if (_abody)
        _abody->update_link_transforms();
      _mutex->unlock();
    }
    else if (DOF == "x-rotation" || DOF == "y-rotation" || DOF == "z-rotation")
    {
      if (DOF[0] == 'x')
        _alpha = value;
      else if (DOF[0] == 'y')
        _beta = value;
      else
        _gamma = value;
      Matrix3 R = ea_to_matrix(_alpha, _beta, _gamma);
      Quat q(&R);

      // update the body
      _mutex->lock();
      rbody->set_orientation(q);
      if (_abody)
        _abody->update_link_transforms();
      _mutex->unlock();
    }
  }

  // call the collision detectors
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, _coldet)
  {
    // set the mode to report all contacts
    cd->mode = CollisionDetection::eAllContacts;

    // see whether there is a collision at the current configurations
    if (cd->is_collision())
    {
      // get the pairs in collision
      const std::set<sorted_pair<CollisionGeometryPtr> >& cpairs = cd->colliding_pairs;
      _table->setNumRows((int) cpairs.size());

      int j = 0;
      for (std::set<sorted_pair<CollisionGeometryPtr> >::const_iterator i = cpairs.begin(); i != cpairs.end(); i++)
      {
        // get the two single bodies
        SingleBodyPtr sb1(i->first->get_single_body());
        SingleBodyPtr sb2(i->second->get_single_body());

        // set the elements in the table
        _table->setText(j,0,QString(sb1->id.c_str()));
        _table->setText(j,1,QString(sb2->id.c_str()));
        j++;
      }
    }
    else
      _table->setNumRows(0);
  }

  // draw colliding triangles
  _oform->drawCollidingTris();
}

