#ifndef _NANOBD_SITECOLLECTION_HPP
#define _NANOBD_SITECOLLECTION_HPP

#include <vector>
#include "Site.hpp"

class SiteCollection : public std::vector<Site> {

  public:
    SiteCollection() : std::vector<Site>() {}
    int getNumberOfSites();
};

#endif
