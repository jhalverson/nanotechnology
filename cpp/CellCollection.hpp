#ifndef _NANOBD_CELLCOLLECTION_HPP
#define _NANOBD_CELLCOLLECTION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include "Cell.hpp"

class CellCollection : public std::vector<Cell> {
  public:
    CellCollection() : std::vector<Cell>() {};
    int cells_per_dimension;
    double cell_side;
    void init(double side_, double rc_, double rs_, int verlet_);
};

#endif
