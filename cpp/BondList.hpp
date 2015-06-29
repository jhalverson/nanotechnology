#ifndef _NANOBD_BONDLIST_HPP
#define _NANOBD_BONDLIST_HPP

#include <vector>
#include "Quad.hpp"

class BondList : public std::vector<Quad> {
  public:
    int bondsFormed;
    int bondsBroken;
    BondList() : std::vector<Quad>() {
      bondsFormed = 0;
      bondsBroken = 0;
  }
};

#endif
