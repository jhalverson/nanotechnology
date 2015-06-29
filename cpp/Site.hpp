#ifndef _NANOBD_SITE_HPP
#define _NANOBD_SITE_HPP

/*
A site is a vector with magnitude equal to the particle radius.
It points from the particle COM to the surface position. It
needs to be updated every step after the particle orientation
at time t+dt has been found. The position of the of the surface
site is found by adding the particle COM to the displacement
vector of the site.
*/

class Site {

  public:
    // basic properties
    int id;
    int type;

    // displacement vector (magnitude is the particle radius)
    double dx;
    double dy;
    double dz;

    // position
    double x;
    double y;
    double z;

    // flag to describe whether or not the site is bonded to another site
    bool bonded;

    // life time of bond
    long long life_time;

    // constructor
    Site(int id_, int type_, int bonded_, double x_, double y_, double z_,
                                          double dx_, double dy_, double dz_) {
      id = id_;
      type = type_;
      bonded = bonded_;
      x = x_ + dx_;
      y = y_ + dy_;
      z = z_ + dz_;
      dx = dx_;
      dy = dy_;
      dz = dz_;
    }

    // destructor
    ~Site() {};

};

#endif
