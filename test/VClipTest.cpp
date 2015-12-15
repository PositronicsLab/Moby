#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>
#include <Moby/Polyhedron.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CompGeom.h>
#include <Moby/TessellatedPolyhedron.h>


using namespace Ravelin;
using namespace Moby;

  double get_random (double r_min, double r_max)
  {
    return (r_max-r_min) * ((double) rand() / (double) RAND_MAX) + r_min;
  }

  TEST(Vclip, Apart_BB_VClip)
  {
    const double TOL = 1e-6;
    const double TRANS_RND_MAX = 10.0;
    const double TRANS_RND_MIN = 2.5;
    const int MAX_ITERATION = 1000;

    std::vector<Origin3d> v(8);
      v[0] = Origin3d(-1.0, -1.0, -1.0);
      v[1] = Origin3d(-1.0, -1.0, +1.0);
      v[2] = Origin3d(-1.0, +1.0, -1.0);
      v[3] = Origin3d(-1.0, +1.0, +1.0);
      v[4] = Origin3d(+1.0, -1.0, -1.0);
      v[5] = Origin3d(+1.0, -1.0, +1.0);
      v[6] = Origin3d(+1.0, +1.0, -1.0);
      v[7] = Origin3d(+1.0, +1.0, +1.0);


      boost::shared_ptr<Moby::BoxPrimitive> p(new BoxPrimitive(2,2,2));
    boost::shared_ptr<Moby::BoxPrimitive> q(new BoxPrimitive(2,2,2));

    TessellatedPolyhedronPtr p_tess = CompGeom::calc_convex_hull(v.begin(), v.end());
    TessellatedPolyhedronPtr q_tess = CompGeom::calc_convex_hull(v.begin(), v.end());

    double trans_q_x, trans_q_y, trans_q_z, quat_q_x, quat_q_y, quat_q_z, quat_q_w;

    for(int i = 0; i< MAX_ITERATION; i++)
    {
      trans_q_x = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_y = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_z = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
        quat_q_x = get_random(-1.0, 1.0);
        quat_q_y = get_random(-1.0, 1.0);
         quat_q_z = get_random(-1.0, 1.0);
        quat_q_w = get_random(-1.0, 1.0);     

        Origin3d q_trans(trans_q_x, trans_q_y, trans_q_z);
        Quatd q_quat(quat_q_x, quat_q_y, quat_q_z, quat_q_w);
        q_quat.normalize();


        boost::shared_ptr<Ravelin::Pose3d> p_pose (new Pose3d(Origin3d(0, 0, 0)));
        boost::shared_ptr<Ravelin::Pose3d> q_pose (new Pose3d(q_quat,q_trans));

        boost::shared_ptr<const Polyhedron::Feature> closestP = *(p->get_polyhedron().get_vertices().begin());
        boost::shared_ptr<const Polyhedron::Feature> closestQ = *(q->get_polyhedron().get_vertices().begin());
       
        //calculating distance using vclip
        Point3d pointp(p_pose);
        Point3d pointq(q_pose);
        double dist_vclip = p->calc_signed_dist(q, pointp, pointq);

        TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

        Ravelin::Origin3d o(0,0,0);
        double dist_m = m_diff->calc_signed_distance(o);

        EXPECT_NEAR(dist_vclip, dist_m, TOL);
      }
  }


  TEST(Vclip, Penetrating_BB_Vclip)
  {
    const double TOL = 1e-6;
    const double TRANS_RND_MAX = 1.0;
    const double TRANS_RND_MIN = -1.0;
    const int MAX_ITERATION = 1000;

    std::vector<Origin3d> v(8);
      v[0] = Origin3d(-1.0, -1.0, -1.0);
      v[1] = Origin3d(-1.0, -1.0, +1.0);
      v[2] = Origin3d(-1.0, +1.0, -1.0);
      v[3] = Origin3d(-1.0, +1.0, +1.0);
      v[4] = Origin3d(+1.0, -1.0, -1.0);
      v[5] = Origin3d(+1.0, -1.0, +1.0);
      v[6] = Origin3d(+1.0, +1.0, -1.0);
      v[7] = Origin3d(+1.0, +1.0, +1.0);


      boost::shared_ptr<Moby::BoxPrimitive> p(new BoxPrimitive(2,2,2));
    boost::shared_ptr<Moby::BoxPrimitive> q(new BoxPrimitive(2,2,2));

    TessellatedPolyhedronPtr p_tess = CompGeom::calc_convex_hull(v.begin(), v.end());
    TessellatedPolyhedronPtr q_tess = CompGeom::calc_convex_hull(v.begin(), v.end());

    double trans_q_x, trans_q_y, trans_q_z, quat_q_x, quat_q_y, quat_q_z, quat_q_w;

    for(int i = 0; i< MAX_ITERATION; i++)
    {
      trans_q_x = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_y = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_z = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
        quat_q_x = get_random(-1.0, 1.0);
        quat_q_y = get_random(-1.0, 1.0);
         quat_q_z = get_random(-1.0, 1.0);
        quat_q_w = get_random(-1.0, 1.0);     

        Origin3d q_trans(trans_q_x, trans_q_y, trans_q_z);
        Quatd q_quat(quat_q_x, quat_q_y, quat_q_z, quat_q_w);
        q_quat.normalize();

        boost::shared_ptr<Ravelin::Pose3d> p_pose (new Pose3d(Origin3d(0, 0, 0)));
        boost::shared_ptr<Ravelin::Pose3d> q_pose (new Pose3d(q_quat,q_trans));

        boost::shared_ptr<const Polyhedron::Feature> closestP = *(p->get_polyhedron().get_vertices().begin());
        boost::shared_ptr<const Polyhedron::Feature> closestQ = *(q->get_polyhedron().get_vertices().begin());
       
       //calculating distance using vclip
        Point3d pointp(p_pose);
        Point3d pointq(q_pose);
        double dist_vclip = p->calc_signed_dist(q, pointp, pointq);

        TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

        Ravelin::Origin3d o(0,0,0);
        double dist_m = m_diff->calc_signed_distance(o);

        EXPECT_NEAR(dist_vclip, dist_m, TOL);
      }
  }
