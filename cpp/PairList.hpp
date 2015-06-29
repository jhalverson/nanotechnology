#ifndef _NANOBD_PAIRLIST_HPP
#define _NANOBD_PAIRLIST_HPP

#include <vector>

class PairList : public std::vector<std::pair<size_t, size_t> > {
  public:
    int rebuilds;
    PairList() : std::vector<std::pair<size_t, size_t> >() {
      rebuilds = 0;
  }
};

#endif
