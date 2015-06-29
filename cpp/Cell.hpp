#ifndef _NANOBD_CELL_HPP
#define _NANOBD_CELL_HPP

#include <vector>

class Cell : public std::vector<unsigned int> {
  public:
    Cell() : std::vector<unsigned int>() {}
    long neighbor_ids[13];
};

#endif
