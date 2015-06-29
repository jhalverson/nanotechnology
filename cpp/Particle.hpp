#ifndef _NANOBD_PARTICLE_HPP
#define _NANOBD_PARTICLE_HPP

#include "SiteCollection.hpp"

class Particle {

  public:
    // basic properties
    // do not confuse id with the index of the particle in the vector
    // id is not used so far but may be needed later for special cases or parallel code
    int id;
    int type;
    int bonds;
    double mass;
    double radius;
    double b;

    // center-of-mass position vector
    double x;
    double y;
    double z;

    // initial center-of-mass position vector
    double x_initial;
    double y_initial;
    double z_initial;

    // flags for unwrapping coordinates with periodic boundary conditions
    int xflag;
    int yflag;
    int zflag;

    // position when PairList is built
    double x_skin;
    double y_skin;
    double z_skin;

    // angular velocity vector
    double wx;
    double wy;
    double wz;

    // force vector
    double fx;
    double fy;
    double fz;

    // torque vector
    double tx;
    double ty;
    double tz;

    // rotation matrix
    double A[3][3];

    SiteCollection sites;

    /*
    double q0;
    double q1;
    double q2;
    double q3;
    double phi_euler;
    double theta_euler;
    double psi_euler;
    */

    // constructor
    Particle(int id_, int type_, double radius_, double x_, double y_, double z_, const double a0_);

    // destructor
    ~Particle() {};

    // update rotation matrix and displacement vectors and positions of sites
    void updateSites(double nwx_, double nwy_, double nwz_, double dtheta_);

    // update position at PairList build
    void updatePositionAtPairListBuild(double side__);

    // scale quaternions to have a magnitude of unity
    //void normalizeQuaternions();

};

#endif
