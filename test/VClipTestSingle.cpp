#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <Moby/Polyhedron.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/CompGeom.h>
#include <Moby/TessellatedPolyhedron.h>
#include <Moby/Log.h>

#ifdef __cplusplus
extern "C"
{
#endif

  #include <ccd/ccd.h>
  #include <ccd/quat.h> // for work with quaternions
  #include <ccd/ccd.h>
  #include <ccd/quat.h> // for work with quaternions
#include "support.h"



#ifdef __cplusplus
}
#endif




// using namespace Ravelin;
// using namespace Moby;


  double get_random (double r_min, double r_max)
  {
    return (r_max-r_min) * ((double) rand() / (double) RAND_MAX) + r_min;
  }


int main(){

	// Moby::Log<Moby::OutputToFile>::reporting_level = 32;
	// Moby::OutputToFile::stream.open("logfile.txt", std::ofstream::out | std::ofstream::app);

	const double TOL = 1e-6;
    const double TRANS_RND_MAX = 1.0;
    const double TRANS_RND_MIN = -1.0;
    const int MAX_ITERATION = 10000;
    const int VERTEX_NUMBER = 20;

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

    CCD_INIT(&ccd);
    ccd.support1 = ccdSupport;
    ccd.support2 = ccdSupport;

    std::vector<Origin3d> v(8);
      v[0] = Origin3d(-1.0, -1.0, -1.0);
      v[1] = Origin3d(-1.0, -1.0, +1.0);
      v[2] = Origin3d(-1.0, +1.0, -1.0);
      v[3] = Origin3d(-1.0, +1.0, +1.0);
      v[4] = Origin3d(+1.0, -1.0, -1.0);
      v[5] = Origin3d(+1.0, -1.0, +1.0);
      v[6] = Origin3d(+1.0, +1.0, -1.0);
      v[7] = Origin3d(+1.0, +1.0, +1.0);


      for (int i=8; i<VERTEX_NUMBER; i++)
      {
        Origin3d o = Origin3d(get_random(-1.0, 1.0),get_random(-1.0, 1.0),get_random(-1.0, 1.0));
        o.normalize();
        o = o * sqrt(3.0);
        v[i] = o;
      }


  boost::shared_ptr<Moby::BoxPrimitive> p(new BoxPrimitive(2,2,2));
    boost::shared_ptr<Moby::BoxPrimitive> q(new BoxPrimitive(2,2,2));

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
	  

	  trans_q_x = 0;
      trans_q_y = 0.5;
      trans_q_z = 0.5;
      quat_q_x = 0;
      quat_q_y = 0;
      quat_q_z = 0;
      quat_q_w = 1;    
        quat_q_w=quat_q_z=quat_q_y=quat_q_x=0;
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
        
        double dist_vclip = p->calc_signed_dist(qconst, pointp, pointq);

		ccdVec3Set(&box2.pos, trans_q_x, trans_q_y, trans_q_z);
        ccdQuatSet(&box2.quat, quat_q_x, quat_q_y, quat_q_z, quat_q_w);
        ccdQuatNormalize(&box2.quat);
        res = ccdGJKPenetration(&box1, &box2, &ccd, &depth, &dir, &pos);    
        std::cout << "ccd_depth: " << depth <<std::endl <<"dir: " << ccdVec3X(&dir) <<", " << ccdVec3Y(&dir)  <<", "<< ccdVec3Z(&dir) << std::endl;
        std::cout << "quat: " << box2.quat.q[0]<< ", "<< box2.quat.q[1] << ", "<< box2.quat.q[2] << ", "<< box2.quat.q[3] << std::endl;
        std::cout << "pos: "  << ccdVec3X(&box2.pos) <<", " << ccdVec3Y(&box2.pos)  <<", "<< ccdVec3Z(&box2.pos) << std::endl;

        TessellatedPolyhedronPtr m_diff = TessellatedPolyhedron::minkowski(*p_tess, p_pose, *q_tess, q_pose);

        Ravelin::Origin3d o(0,0,0);
        double dist_m = m_diff->calc_signed_distance(o);


        std::cout<<*q_pose<<std::endl;
        std::cout<<"vclip: " << dist_vclip << std::endl;
        std::cout<<"minkowski: " << dist_m << std::endl;
        //EXPECT_NEAR(dist_vclip, dist_m, TOL);
    

}