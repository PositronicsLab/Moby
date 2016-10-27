#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>
#include <Moby/Polyhedron.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/CompGeom.h>
#include <Moby/TessellatedPolyhedron.h>
#include <Moby/Log.h>
#include <gperftools/profiler.h> 
#include <sys/times.h>
#ifdef USE_LIBCCD
#include <ccd/ccd.h>
#include <ccd/quat.h> // for work with quaternions
#endif
#include "support.h"


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
    const int MAX_ITERATION = 1;

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
    Polyhedron p_test;
    p_tess->to_polyhedron(p_test);
std::ofstream out("box.wrl");
    Polyhedron::to_vrml(out, p_test);
out.close();

    double trans_q_x, trans_q_y, trans_q_z, quat_q_x, quat_q_y, quat_q_z, quat_q_w;
    double total_vc, total_m;
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
        
        tms vcstart;  
        clock_t v_start_c = times(&vcstart);
        
        double dist_vclip = p->calc_signed_dist(q, pointp, pointq);
        
        tms vcstop;  
        clock_t v_end_c = times(&vcstart);

        total_vc += (v_end_c-v_start_c);

        tms mdstart;  
        clock_t m_start_c = times(&mdstart);

        TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

        tms mdstop;  
        clock_t m_end_c = times(&mdstop);
        total_m += m_end_c - m_start_c;


        Ravelin::Origin3d o(0,0,0);
        double dist_m = m_diff->calc_signed_distance(o);

        EXPECT_NEAR(dist_vclip, dist_m, TOL);
    }
    // std::cout << "average vclip time: " << total_vc/(double) MAX_ITERATION <<std::endl;
    // std::cout<< "average minkowski time " << total_m/(double) MAX_ITERATION << std::endl;
  }


/*  TEST(Vclip, Kissing_BB_Vclip)
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


      boost::shared_ptr<Ravelin::Pose3d> p_pose (new Pose3d(Origin3d(0, 0, 0)));
      boost::shared_ptr<Ravelin::Pose3d> q_pose (new Pose3d(q_quat,q_trans));

      boost::shared_ptr<const Polyhedron::Feature> closestP = *(p->get_polyhedron().get_vertices().begin());
      boost::shared_ptr<const Polyhedron::Feature> closestQ = *(q->get_polyhedron().get_vertices().begin());

      TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

      
       
        //calculating distance using vclip
      Point3d pointp(p_pose);
      Point3d pointq(q_pose);
//      double dist_vclip = p->calc_signed_dist(q, pointp, pointq);
//    double dist_vclip = Polyhedron::vclip(p, q, p_pose, q_pose, closestP, closestQ);


        Ravelin::Origin3d o(0,0,0);
        unsigned & closest_facet;
        double dist_m = m_diff->calc_signed_distance(o, closest_facet);


        EXPECT_NEAR(dist_vclip, dist_m, TOL);
      }
  }
*/

  TEST(Vclip, Penetrating_Poly_Vclip)
  {

    Moby::Log<Moby::OutputToFile>::reporting_level = (LOG_COLDET);
    Moby::OutputToFile::stream.open("logging.out");
    const double TOL = 1e-6;
    const double TRANS_RND_MAX = 1.0;
    const double TRANS_RND_MIN = -1.0;
    const int MAX_ITERATION = 10000;
    const int VERTEX_NUMBER = 8;

    // set up ccd_t struct
    #ifdef USE_LIBCCD
    ccd_t ccd;
    CCD_BOX(box1);
    CCD_BOX(box2);
    int res;
    ccd_vec3_t axis;
    ccd_quat_t rot;
    ccd_real_t depth;
    ccd_vec3_t dir, pos;
    box1.x = box1.y = box1.z = 2.;
    box2.x = 2;
    box2.y = 2.;
    box2.z = 2;
    fprintf(stderr, "\n\n\n---- boxboxPenetration ----\n\n\n");

    CCD_INIT(&ccd);
    ccd.support1 = ccdSupport;
    ccd.support2 = ccdSupport;
    #endif
    
      // now intersect holds 0 if obj1 and obj2 intersect, -1 otherwise
      // in depth, dir and pos is stored penetration depth, direction of
      // separation vector and position in global coordinate system
  
    std::vector<Origin3d> v(VERTEX_NUMBER);
      v[0] = Origin3d(-1.0, -1.0, -1.0);
      v[1] = Origin3d(-1.0, -1.0, +1.0);
      v[2] = Origin3d(-1.0, +1.0, -1.0);
      v[3] = Origin3d(-1.0, +1.0, +1.0);
      v[4] = Origin3d(+1.0, -1.0, -1.0);
      v[5] = Origin3d(+1.0, -1.0, +1.0);
      v[6] = Origin3d(+1.0, +1.0, -1.0);
      v[7] = Origin3d(+1.0, +1.0, +1.0);
      // for (int i=8; i<VERTEX_NUMBER; i++)
      // {
      //   Origin3d o = Origin3d(get_random(-1.0, 1.0),get_random(-1.0, 1.0),get_random(-1.0, 1.0));
      //   o.normalize();
      //   o = o * sqrt(3.0);
      //   v[i] = o;
      // }



  // boost::shared_ptr<Moby::BoxPrimitive> p(new BoxPrimitive(2,2,2));
  //   boost::shared_ptr<Moby::BoxPrimitive> q(new BoxPrimitive(2,2,2));

    TessellatedPolyhedronPtr p_tess = CompGeom::calc_convex_hull(v.begin(), v.end());
    TessellatedPolyhedronPtr q_tess = CompGeom::calc_convex_hull(v.begin(), v.end());

    Polyhedron p_poly,q_poly;
    p_tess->to_polyhedron(p_poly);
    q_tess->to_polyhedron(q_poly);
    boost::shared_ptr<Moby::PolyhedralPrimitive> p(new PolyhedralPrimitive());
    p->set_polyhedron(p_poly);
    boost::shared_ptr<Moby::PolyhedralPrimitive> q(new PolyhedralPrimitive());
    q->set_polyhedron(q_poly);
    boost::shared_ptr<const Moby::Primitive> qconst = q;

    double trans_q_x, trans_q_y, trans_q_z, quat_q_x, quat_q_y, quat_q_z, quat_q_w;
    double total_vc = 0,total_m = 0,total_ccd = 0, total_ccd2 = 0;

    for(int i = 0; i< MAX_ITERATION; i++)
    {
      trans_q_x = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_y = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      trans_q_z = get_random(TRANS_RND_MIN,TRANS_RND_MAX);
      quat_q_x = get_random(-1.0, 1.0);
      quat_q_y = get_random(-1.0, 1.0);
      quat_q_z = get_random(-1.0, 1.0);
      quat_q_w = get_random(-1.0, 1.0);     
        // quat_q_w=quat_q_z=quat_q_y=quat_q_x=0;
        Origin3d q_trans(trans_q_x, trans_q_y, trans_q_z);
        Quatd q_quat(quat_q_x, quat_q_y, quat_q_z, quat_q_w);
        q_quat.normalize();

        boost::shared_ptr<Ravelin::Pose3d> p_pose (new Pose3d(Origin3d(0, 0, 0)));
        boost::shared_ptr<Ravelin::Pose3d> q_pose (new Pose3d(q_quat,q_trans));

        boost::shared_ptr<const Polyhedron::Feature> closestP = *(p->get_polyhedron().get_faces().begin());
        boost::shared_ptr<const Polyhedron::Feature> closestQ = *(q->get_polyhedron().get_faces().begin());



       //calculating distance using vclip
        Point3d pointp(p_pose);
        Point3d pointq(q_pose);

        tms vcstart;  
        clock_t v_start_c = times(&vcstart);
        
        double dist_vclip = p->calc_signed_dist(qconst, pointp, pointq);
        
        tms vcstop;  
        clock_t v_end_c = times(&vcstart);

        total_vc += (v_end_c-v_start_c);

        //minkowskii
        tms mdstart;  
        clock_t m_start_c = times(&mdstart);

        TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

        Ravelin::Origin3d o(0,0,0);
        double dist_m = m_diff->calc_signed_distance(o);

        tms mdstop;  
        clock_t m_end_c = times(&mdstop);
        total_m += m_end_c - m_start_c;

        std::cout<< i << std::endl;
        EXPECT_NEAR(dist_vclip, dist_m, TOL) << *q_pose <<std::endl;

        // ccd
        tms ccdstart;  
        clock_t ccd_start_c = times(&ccdstart);

        #ifdef USE_LIBCCD        
        ccdVec3Set(&box2.pos, trans_q_x, trans_q_y, trans_q_z);
        ccdQuatSet(&box2.quat, quat_q_x, quat_q_y, quat_q_z, quat_q_w);
        ccdQuatNormalize(&box2.quat);
        res = ccdGJKPenetration(&box1, &box2, &ccd, &depth, &dir, &pos);
        EXPECT_NEAR(depth,dist_m,TOL);
        clock_t ccd_end_c = times(&ccdstart);
        total_ccd += (ccd_end_c-ccd_start_c);

        ccd_start_c = times(&ccdstart);
        int intersect = ccdGJKIntersect(&p, &q, &ccd);
        ccd_end_c = times(&ccdstart);
        total_ccd2 += (ccd_end_c-ccd_start_c);
        depth = 0;
        #endif
    }
   // ProfilerStop();

    for(int i = 0; i < VERTEX_NUMBER; i++){
      std::cout << v[i] <<std::endl;
    }
    std::cout << "average vclip time: " << total_vc/(double) MAX_ITERATION <<std::endl;
    std::cout<< "average minkowski time " << total_m/(double) MAX_ITERATION << std::endl;
    std::cout<< "average ccd time " << total_ccd/(double) MAX_ITERATION << std::endl;
    std::cout<< "average ccd intersection check time " << total_ccd2/(double) MAX_ITERATION << std::endl;

  }

  TEST(Vclip, Penetrating_BB_Vclip)
  {
    //ProfilerStart("prof.out");
    Moby::Log<Moby::OutputToFile>::reporting_level = (LOG_COLDET);
    Moby::OutputToFile::stream.open("logging.out");
    const double TOL = 1e-6;
    const double TRANS_RND_MAX = 1.0;
    const double TRANS_RND_MIN = -1.0;
    const int MAX_ITERATION = 100;
    const int VERTEX_NUMBER = 8;
    const double WIDTH = 2, HEIGHT = 2, DEPTH = 2;

    fprintf(stderr, "\n\n\n---- boxboxPenetration ----\n\n\n");

    // set up ccd_t struct
    #ifdef USE_LIBCCD
    ccd_t ccd;
    CCD_BOX(box1);
    CCD_BOX(box2);
    int res;
    ccd_vec3_t axis;
    ccd_quat_t rot;
    ccd_real_t depth;
    ccd_vec3_t dir, pos;

    // setup box dimensions 
    box1.x = box2.x = WIDTH;
    box1.y = box2.y = HEIGHT;
    box1.z = box2.z = DEPTH;

    // setup libccd
    CCD_INIT(&ccd);
    ccd.support1 = ccdSupport;
    ccd.support2 = ccdSupport;
    #endif

    // create box primitives    
    boost::shared_ptr<Moby::BoxPrimitive> p(new BoxPrimitive(WIDTH,HEIGHT,DEPTH));
    boost::shared_ptr<Moby::BoxPrimitive> q(new BoxPrimitive(WIDTH,HEIGHT,DEPTH));
    std::vector<Origin3d> v(VERTEX_NUMBER);
    v[0] = Origin3d(-WIDTH/2.0, -HEIGHT/2, -DEPTH/2);
    v[1] = Origin3d(-WIDTH/2.0, -HEIGHT/2, +DEPTH/2);
    v[2] = Origin3d(-WIDTH/2.0, +HEIGHT/2, -DEPTH/2);
    v[3] = Origin3d(-WIDTH/2.0, +HEIGHT/2, +DEPTH/2);
    v[4] = Origin3d(+WIDTH/2.0, -HEIGHT/2, -DEPTH/2);
    v[5] = Origin3d(+WIDTH/2.0, -HEIGHT/2, +DEPTH/2);
    v[6] = Origin3d(+WIDTH/2.0, +HEIGHT/2, -DEPTH/2);
    v[7] = Origin3d(+WIDTH/2.0, +HEIGHT/2, +DEPTH/2);
    TessellatedPolyhedronPtr p_tess = CompGeom::calc_convex_hull(v.begin(), v.end());
    TessellatedPolyhedronPtr q_tess = CompGeom::calc_convex_hull(v.begin(), v.end());

    const Polyhedron& p_pol = p->get_polyhedron();
    const Polyhedron& q_pol = q->get_polyhedron();
    boost::shared_ptr<const Moby::Primitive> qconst = q;

    double trans_q_x, trans_q_y, trans_q_z, quat_q_x, quat_q_y, quat_q_z, quat_q_w;
    double total_vc = 0,total_m = 0,total_ccd = 0, total_ccd2 = 0;

    for(int i = 0; i< MAX_ITERATION; i++)
    {
      // setup random positions and orientations for q
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
      boost::shared_ptr<Ravelin::Pose3d> p_pose(new Pose3d(Origin3d(0, 0, 0)));
      boost::shared_ptr<Ravelin::Pose3d> q_pose(new Pose3d(q_quat,q_trans));

      // setup arbitrary closest features
      boost::shared_ptr<const Polyhedron::Feature> closestP = p_pol.get_faces().front();
      boost::shared_ptr<const Polyhedron::Feature> closestQ = q_pol.get_faces().front();

      //calculate distance using vclip
      Point3d pointp(p_pose);
      Point3d pointq(q_pose);

      tms vcstart;  
      clock_t v_start_c = times(&vcstart);
      double dist_vclip = p->calc_signed_dist(qconst, pointp, pointq);
      tms vcstop;  
      clock_t v_end_c = times(&vcstart);
      total_vc += (v_end_c-v_start_c);

      // Minkowskii
      tms mdstart;  
      clock_t m_start_c = times(&mdstart);

      TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

      Ravelin::Origin3d o(0,0,0);
      double dist_m = m_diff->calc_signed_distance(o);

      tms mdstop;  
      clock_t m_end_c = times(&mdstop);
      total_m += m_end_c - m_start_c;

      std::cout<< i << std::endl;
      EXPECT_NEAR(dist_vclip, dist_m, TOL) << *q_pose <<std::endl;

      #ifdef USE_LIBCCD
      // try CCD
      tms ccdstart;  
      clock_t ccd_start_c = times(&ccdstart);            
      ccdVec3Set(&box2.pos, trans_q_x, trans_q_y, trans_q_z);
      ccdQuatSet(&box2.quat, q_quat.x, q_quat.y, q_quat.z, q_quat.w);
      res = ccdGJKPenetration(&box1, &box2, &ccd, &depth, &dir, &pos);
      clock_t ccd_end_c = times(&ccdstart);
      if (dist_m >= 0.0)
        EXPECT_NEAR(depth,0.0,TOL);
      else
        EXPECT_NEAR(-depth,dist_m,TOL);

      // try libccd using the (faster) intersection test
      ccd_start_c = times(&ccdstart);            
      int intersect = ccdGJKIntersect(&p, &q, &ccd);
      depth = 0;
      ccd_end_c = times(&ccdstart);
      total_ccd2 += (ccd_end_c-ccd_start_c);
      #endif
    }
   // ProfilerStop();

    for(int i = 0; i < VERTEX_NUMBER; i++){
      std::cout << v[i] <<std::endl;
    }
    std::cout << "average vclip time: " << total_vc/(double) MAX_ITERATION <<std::endl;
    std::cout<< "average minkowski time " << total_m/(double) MAX_ITERATION << std::endl;
    #ifdef USE_LIBCCD
    std::cout<< "average ccd time " << total_ccd/(double) MAX_ITERATION << std::endl;
    std::cout<< "average ccd intersection check time " << total_ccd2/(double) MAX_ITERATION << std::endl;
    #endif
  }

