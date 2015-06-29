#ifndef _NANOBD_QUAD_HPP
#define _NANOBD_QUAD_HPP

class Quad {

  public:
    int first;
    int second;
    int third;
    int fourth;

    Quad(int w, int x, int y, int z) : first(w),
                                       second(x),
                                       third(y),
                                       fourth(z)
                                       {}

    bool operator==(Quad v) const {
      return ((first == v.first) &&
              (second == v.second) &&
              (third == v.third) &&
              (fourth == v.fourth)) ||
             ((first == v.third) &&
              (second == v.fourth) &&
              (third == v.first) &&
              (fourth == v.second));
    }

    void operator=(Quad v) {
      first = v.first;
      second = v.second;
      third = v.third;
      fourth = v.fourth;
    }


    Quad (const Quad& v) : first(v.first), second(v.second), third(v.third), fourth(v.fourth) {}

};

#endif
