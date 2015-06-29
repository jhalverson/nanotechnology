#include <cmath>
#include <iostream>
#include "Particle.hpp"

/* This routine assumes particles coordinates are in the central box. */
Particle::Particle(int id_, int type_, double radius_, double x_,
                   double y_, double z_, const double a0_) {
  id = id_;
  type = type_;
  bonds = 0;
  radius = radius_;
  b = a0_ / radius_;

  x = x_;
  y = y_;
  z = z_;

  x_initial = x_;
  y_initial = y_;
  z_initial = z_;

  x_skin = x_;
  y_skin = y_;
  z_skin = z_;

  xflag = 0;
  yflag = 0;
  zflag = 0;

  fx = 0.0;
  fy = 0.0;
  fz = 0.0;

  tx = 0.0;
  ty = 0.0;
  tz = 0.0;

  A[0][0] = 1.0;
  A[1][0] = 0.0;
  A[2][0] = 0.0;
  A[0][1] = 0.0;
  A[1][1] = 1.0;
  A[2][1] = 0.0;
  A[0][2] = 0.0;
  A[1][2] = 0.0;
  A[2][2] = 1.0;

  /*
  phi_euler =   2.0 * atan(1.0) * drand48();
  theta_euler = 2.0 * atan(1.0) * drand48();
  psi_euler =   2.0 * atan(1.0) * drand48();
  q0 = cos(0.5 * theta_euler) * cos(0.5 * (phi_euler + psi_euler));
  q1 = sin(0.5 * theta_euler) * cos(0.5 * (phi_euler - psi_euler));
  q2 = sin(0.5 * theta_euler) * sin(0.5 * (phi_euler - psi_euler));
  q3 = cos(0.5 * theta_euler) * sin(0.5 * (phi_euler + psi_euler));
  normalizeQuaternions();
  updateRotationMatrix();
  */
}

void Particle::updatePositionAtPairListBuild(double side__) {
  x_skin = x + xflag * side__;
  y_skin = y + yflag * side__;
  z_skin = z + zflag * side__;
}

void Particle::updateSites(double nwx_, double nwy_, double nwz_, double dtheta_) {
  // update the rotation matrix (note that nw is the normalized angular velocity)
  A[0][0] = cos(dtheta_) + nwx_ * nwx_ * (1.0 - cos(dtheta_));
  A[1][0] = nwy_ * nwx_ * (1.0 - cos(dtheta_)) + nwz_ * sin(dtheta_);
  A[2][0] = nwz_ * nwx_ * (1.0 - cos(dtheta_)) - nwy_ * sin(dtheta_);
  A[0][1] = nwx_ * nwy_ * (1.0 - cos(dtheta_)) - nwz_ * sin(dtheta_);
  A[1][1] = cos(dtheta_) + nwy_ * nwy_ * (1.0 - cos(dtheta_));
  A[2][1] = nwz_ * nwy_ * (1.0 - cos(dtheta_)) + nwx_ * sin(dtheta_);
  A[0][2] = nwx_ * nwz_ * (1.0 - cos(dtheta_)) + nwy_ * sin(dtheta_);
  A[1][2] = nwy_ * nwz_ * (1.0 - cos(dtheta_)) - nwx_ * sin(dtheta_);
  A[2][2] = cos(dtheta_) + nwz_ * nwz_ * (1.0 - cos(dtheta_));

  // apply the rotation matrix to the displacement vector of each site
  for (SiteCollection::iterator it = sites.begin(); it < sites.end(); it++) {
    Site& s = *it;
    double xtmp = s.dx;
    double ytmp = s.dy;
    s.dx = A[0][0] * s.dx + A[0][1] * s.dy + A[0][2] * s.dz;
    s.dy = A[1][0] * xtmp + A[1][1] * s.dy + A[1][2] * s.dz;
    s.dz = A[2][0] * xtmp + A[2][1] * ytmp + A[2][2] * s.dz;

    // normalize displacement vector due to numerical roundoff
    double fac = radius / sqrt(s.dx * s.dx + s.dy * s.dy + s.dz * s.dz);
    s.dx *= fac;
    s.dy *= fac;
    s.dz *= fac;

    // update the position of the site
    s.x = x + s.dx;
    s.y = y + s.dy;
    s.z = z + s.dz;
  }
}

/*
void Particle::normalizeQuaternions() {
  double q2sqrt = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
  q0 = q0 / q2sqrt;
  q1 = q1 / q2sqrt;
  q2 = q2 / q2sqrt;
  q3 = q3 / q2sqrt;
}
*/
