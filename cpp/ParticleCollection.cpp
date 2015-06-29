#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "ParticleCollection.hpp"
#include "SiteCollection.hpp"
#include <gsl/gsl_linalg.h>

void mic(double& xij_, double& yij_, double& zij_, const double side_) {
  double side_half = 0.5 * side_;
  if(xij_ >  side_half) xij_ -= side_;
  if(xij_ < -side_half) xij_ += side_;
  if(yij_ >  side_half) yij_ -= side_;
  if(yij_ < -side_half) yij_ += side_;
  if(zij_ >  side_half) zij_ -= side_;
  if(zij_ < -side_half) zij_ += side_;
}

long neighbor2(double x_, double y_, double z_, double side_,
               double cell_side_, int cells_per_dimension_) {
  if (x_ > side_) x_ -= side_;
  if (x_ < 0.0)   x_ += side_;
  if (y_ > side_) y_ -= side_;
  if (y_ < 0.0)   y_ += side_;
  if (z_ > side_) z_ -= side_;
  if (z_ < 0.0)   z_ += side_;
  return int(x_ / cell_side_) +
         int(y_ / cell_side_) * cells_per_dimension_ +
         int(z_ / cell_side_) * cells_per_dimension_ * cells_per_dimension_;
}

/* In building Verlet lists rs is the surface-to-surface distance while for cells it
   is the cell size. In the case of cells one must be careful to make sure the size of
   the cell is at least the surface-to-distance plus rc. Note that cells only strictly
   work when all particles have a radius of 1. Verlet list is fine for all cases. */
void ParticleCollection::buildPairs(PairList& pairs_, CellCollection& cells_, const double side_,
                                    const int verlet_, double& rs_, const double rc_,
                                    bool& skin_violation_, const double rs_next_) const {
  pairs_.clear();
  unsigned int num = this->size();

  if (verlet_ == 2) {
    rs_ = rs_next_;
    if (int(side_ / (rc_ + rs_)) != cells_.cells_per_dimension)
      cells_.init(side_, rc_, rs_, verlet_);
    else
      for (long k = 0; k < cells_.size(); k++) cells_[k].clear();

    // assign particles to cells
    for (unsigned int i = 0; i < num; i++) {
      const Particle pi = this->at(i);
      long cell_id = neighbor2(pi.x , pi.y, pi.z, side_, cells_.cell_side, cells_.cells_per_dimension);
      cells_[cell_id].push_back(i);
    }

    // loop over cells looking for pairs
    for (long k = 0; k < cells_.size(); k++) {
      unsigned int numk = cells_[k].size();
      if (numk > 0) {
	double xi;
	double yi;
	double zi;
	double xij;
	double yij;
	double zij;
	double rijsq;
	double radiusi;
	double radiusj;
        if (numk > 1) {
	  // loop over pairs of particles in cell k
	  for (unsigned int i = 0; i < numk - 1; i++) {
	    const Particle pi = this->at(cells_[k].at(i));
	    xi = pi.x;
	    yi = pi.y;
	    zi = pi.z;
	    radiusi = pi.radius;
	    for (unsigned int j = i + 1; j < numk; j++) {
	      const Particle pj = this->at(cells_[k].at(j));
	      xij = xi - pj.x;
	      yij = yi - pj.y;
	      zij = zi - pj.z;
	      radiusj = pj.radius;
	      mic(xij, yij, zij, side_);
	      rijsq = xij * xij + yij * yij + zij * zij;
	      if(sqrt(rijsq) - radiusi - radiusj < rs_)
                pairs_.push_back(std::make_pair(cells_[k].at(i), cells_[k].at(j)));
	    }
	  }
        }
	// for each particle in cell k check against particles in neighboring cells
	for (unsigned int i = 0; i < numk; i++) {
	  const Particle pi = this->at(cells_[k].at(i));
	  xi = pi.x;
	  yi = pi.y;
	  zi = pi.z;
	  radiusi = pi.radius;
	  for (unsigned int n = 0; n < 13; n++) {
	    long ngh = cells_[k].neighbor_ids[n];
	    unsigned int numl = cells_[ngh].size();
	    for (unsigned int j = 0; j < numl; j++) {
	      const Particle pj = this->at(cells_[ngh].at(j));
	      xij = xi - pj.x;
	      yij = yi - pj.y;
	      zij = zi - pj.z;
	      radiusj = pj.radius;
	      mic(xij, yij, zij, side_);
	      rijsq = xij * xij + yij * yij + zij * zij;
	      if(sqrt(rijsq) - radiusi - radiusj < rs_)
                pairs_.push_back(std::make_pair(cells_[k].at(i), cells_[ngh].at(j)));
	    }
	  }
	}
      }
    }
  }

  if (verlet_ == 1) {
    rs_ = rs_next_; // automatic skin thickness updating
    double xi;
    double yi;
    double zi;
    double xij;
    double yij;
    double zij;
    double rijsq;
    double radiusi;
    double radiusj;
    for (size_t i = 0; i < num - 1; i++) {
      const Particle pi = this->at(i);
      xi = pi.x;
      yi = pi.y;
      zi = pi.z;
      radiusi = pi.radius;
      for (size_t j = i + 1; j < num; j++) {
	const Particle pj = this->at(j);
	xij = xi - pj.x;
	yij = yi - pj.y;
	zij = zi - pj.z;
	radiusj = pj.radius;
        mic(xij, yij, zij, side_);
	rijsq = xij * xij + yij * yij + zij * zij;
	if(sqrt(rijsq) - radiusi - radiusj < rs_) pairs_.push_back(std::make_pair(i, j));
      }
    }
  }

  if (verlet_ == 0) {
    for (size_t i = 0; i < num - 1; i++) {
      for (size_t j = i + 1; j < num; j++) {
        pairs_.push_back(std::make_pair(i, j));
      }
    }
  }

  skin_violation_ = false;
  pairs_.rebuilds++;
}

int poisson(double mu_) {
  double L = exp(-mu_);
  double p = 1.0;
  int k = 0;
  do {
    k++;
    p *= drand48();
  } while (p > L);
  return k - 1;
}

void ParticleCollection::updateBondList(const PairList pairs_, double side_, double dt_,
                                        BondList& bonds_, const BondList unbreakable_,
                                        double xcp_[], double Lcp_[], long long t_,
                                        VectorLL& life_times_, size_t max_site_type_) {
  size_t num = this->size();
  double xab;
  double yab;
  double zab;
  double rabsq;
  BondList possible;

  const double pi_const = 4.0 * atan(1.0);
  const double assoc_rate = 10.0;

  std::vector<std::pair<size_t, size_t> >::const_iterator it;
  for (it = pairs_.begin(); it < pairs_.end(); it++) {
    Particle& pi = this->at(it->first);
    SiteCollection& si = pi.sites;
    Particle& pj = this->at(it->second);
    SiteCollection& sj = pj.sites;
    for (size_t a = 0; a < si.size(); a++) {
      Site& sa = si.at(a);
      if (!sa.bonded) {
	for (size_t b = 0; b < sj.size(); b++) {
	  Site& sb = sj.at(b);
	  if (!sb.bonded && sa.type + sb.type == 0) {
	    xab = sa.x - sb.x;
	    yab = sa.y - sb.y;
	    zab = sa.z - sb.z;
            mic(xab, yab, zab, side_);
	    rabsq = xab * xab + yab * yab + zab * zab;
	    if (rabsq < pow(Lcp_[abs(sa.type)], 2)) {
	      Quad q(it->first, a, it->second, b);
	      possible.push_back(q);
	    }
	  }
	}
      }
    }
  }

  // try a bond forming event
  if (!possible.empty()) {
    double half_sphere = (2.0 / 3.0) * pi_const * pow(Lcp_[1], 3); // TODO: which value of Lcp to use
    double P_assoc = dt_ * assoc_rate / half_sphere;
    int number_to_form = poisson(P_assoc * possible.size());
    if (number_to_form > possible.size()) number_to_form = possible.size();
    int number_formed = 0;
    std::vector<int> already_picked;
    while (number_to_form - number_formed > 0) {
      int k = possible.size();
      while (k == possible.size() || std::find(already_picked.begin(), already_picked.end(), k) != already_picked.end())
        k = int((double(rand()) / RAND_MAX) * possible.size());
      already_picked.push_back(k);
      Quad d = possible.at(k);
      Particle& pi = this->at(d.first);
      Site& sa = pi.sites.at(d.second);
      Particle& pj = this->at(d.third);
      Site& sb = pj.sites.at(d.fourth);
      pi.bonds++;
      pj.bonds++;
      sa.bonded = true;
      sb.bonded = true;
      sa.life_time = t_;
      sb.life_time = t_;
      bonds_.push_back(d);
      bonds_.bondsFormed++;
      number_formed++;
    }
  }

  // could store bonds along with their types in a pair and then count the types above
  // and then when that type is selected then choose a random number between 0 and number of that type

  // bond_type_count should be computed once at beginning unless restarting
  int bond_type_count[512];
  for (int i = 0; i < 512; i++) bond_type_count[i] = 0;

  BondList::const_iterator its;
  for (its = bonds_.begin(); its < bonds_.end(); its++) {
    Particle& pi = this->at(its->first);
    Site sa = pi.sites.at(its->second);
    //Particle& pj = this->at(it->third);
    //Site sb = pj.sites.at(it->fourth);
    bond_type_count[abs(sa.type)] += 1;
  }

  // if all bonds have same x and L then why bother with dealing by types?

  // max_site_type_ assumes counting from 0 up to abs(max site type) inclusive
  // this means any MACRO define should be number of colors plus 1 or may go out of bounds
  // xcp[0] and Lcp[0] are overwritten in readScript
  double weight_sum = 0.0;
  double weight[512];
  for (int i = 0; i <= max_site_type_; i++) {
    //weight[i] = (2.0 / 3.0) * pi_const * pow(Lcp_[i], 3) * exp(xcp_[i]) * bond_type_count[i];
    weight[i] = exp(xcp_[i]) * bond_type_count[i];
    weight_sum += weight[i];
    //std::cout << i << " " << xcp_[i] << " " << exp(xcp_[i]) << " " << bond_type_count[i] << " " << weight[i] << " " << weight_sum << std::endl;
  }
  double w = 0.0;
  double bounds[512];
  for (int i = 0; i <= max_site_type_; i++) {
    w += weight[i];
    bounds[i] = w / weight_sum;
    //std::cout << i << " " << bounds[i] << std::endl;
  }
  // int number_to_break = poisson(2.0 * this->size() / pow(side_, 3) * weight_sum * dt_ * assoc_rate);
  // generate RN betweeon 0 and 1. If 0.23 then find the bond type between these bounds and randomly
  // destroy a bond of that type

  // bonds types must be sequential beginning at ; Add a check for this in readScript
  // what about type 0

  // does V_R really cancel when computing individual weights where each bond type has its own L

  // after breaking a bond should the weights be recomputed?

  // with the while loops seems like should be a problem if only bond left is low probability and it must be broken
  // maybe if try 100 times then manually break it

  int number_to_break = poisson(2.0 * weight_sum * dt_ * assoc_rate * this->size() / pow(side_, 3));
  // std::cout << "number_to_break = " << number_to_break << ", " << 2.0 * weight_sum * dt_ * assoc_rate * this->size() / pow(side_, 3) << std::endl;
  int number_broken = 0;
 
  // number_to_break = 1;

  // could be problem when have 20 bonds and need to break 3. If only 2 of one type and that type
  // comes up 3 times the algorithm doesn't know that those bonds are gone; it will break due to ct count

  while (bonds_.size() - unbreakable_.size() > 0 && number_to_break - number_broken > 0) {
    double r = drand48();
    int type_to_break = 0;
    while (r > bounds[type_to_break]) type_to_break++;
    if (bond_type_count[type_to_break] == 0) std::cout << "ERROR: bond type to break\n" << std::endl;
    //std::cout << "type: " << type_to_break << std::endl;
    int ct = 0;
    int stype;
    int k = bonds_.size();
    while (k == bonds_.size() || stype != type_to_break) {
      k = int((double(rand()) / RAND_MAX) * bonds_.size());
      Quad d = bonds_.at(k);
      Particle& pi = this->at(d.first);
      Site& sa = pi.sites.at(d.second);
      stype = abs(sa.type);
      ct++;
      //if (ct > 1000) { std::cout << "ct > 1000, type = " << type_to_break << std::endl; break; }
      if (ct > 1000) break;
    }
    if (std::find(unbreakable_.begin(), unbreakable_.end(), bonds_.at(k)) == unbreakable_.end() && ct < 1000) {
      Quad d = bonds_.at(k);
      Particle& pi = this->at(d.first);
      Site& sa = pi.sites.at(d.second);
      Particle& pj = this->at(d.third);
      Site& sb = pj.sites.at(d.fourth);
      pi.bonds--;
      pj.bonds--;
      sa.bonded = false;
      sb.bonded = false;
      life_times_.push_back(std::make_pair(t_ - sa.life_time, abs(sa.type)));
      bond_type_count[abs(sa.type)]--;
      bonds_.erase(bonds_.begin() + k);
      bonds_.bondsBroken++;
    }
    number_broken++;
  }
}

void ParticleCollection::computeSpringForcesAndTorques(double side_, BondList bonds_,
                                                       double Lcp_[], double k_spring_) {
  //double side_half = 0.5 * side_;
  double xab;
  double yab;
  double zab;
  double rabsq;
  double rab;
  double fmagn_over_rab;
  double fabx;
  double faby;
  double fabz;

  BondList::const_iterator it;
  for (it = bonds_.begin(); it < bonds_.end(); it++) {
    Particle& pi = this->at(it->first);
    Site sa = pi.sites.at(it->second);
    Particle& pj = this->at(it->third);
    Site sb = pj.sites.at(it->fourth);

    xab = sa.x - sb.x;
    yab = sa.y - sb.y;
    zab = sa.z - sb.z;

    mic(xab, yab, zab, side_);

    rabsq = xab * xab + yab * yab + zab * zab;
    rab = sqrt(rabsq);
    fmagn_over_rab = k_spring_ * (Lcp_[abs(sa.type)] - rab) / rab;
    fabx = fmagn_over_rab * xab;
    faby = fmagn_over_rab * yab;
    fabz = fmagn_over_rab * zab;

    // update forces on particles i and j
    pi.fx += fabx;
    pi.fy += faby;
    pi.fz += fabz;
    pj.fx -= fabx;
    pj.fy -= faby;
    pj.fz -= fabz;

    // update torques on particles i and j
    pi.tx += sa.dy * fabz - sa.dz * faby;
    pi.ty += sa.dz * fabx - sa.dx * fabz;
    pi.tz += sa.dx * faby - sa.dy * fabx;
    pj.tx -= sb.dy * fabz - sb.dz * faby;
    pj.ty -= sb.dz * fabx - sb.dx * fabz;
    pj.tz -= sb.dx * faby - sb.dy * fabx;
  }
}

/* all particles should have same size (assume a_i = 1.0) to use hydrodynamics, otherwise must modify Dij */
void ParticleCollection::updateDiffusionTensor(const PairList pairs_, double side_,
                                               const double a0_, gsl_matrix* Dij_) {
  double xij;
  double yij;
  double zij;
  double rij;
  double rijsq;
  double radiusi;
  double radiusj;
  double Dab;

  double I[3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      (a == b) ? I[a][b] = 1.0 : 0.0;

  // fill in upper triangle
  for (PairList::const_iterator it = pairs_.begin(); it < pairs_.end(); it++) {
    int i = it->first;
    int j = it->second;
    if (i > j) std::cout << "i j problem" << std::endl;
    /* swap is needed if use Verlet since must create upper half of the matrix here if (i > j) std::swap(i, j);*/

    Particle& pi = this->at(i);
    radiusi = pi.radius;
    Particle& pj = this->at(j);
    radiusj = pj.radius;

    xij = pi.x - pj.x;
    yij = pi.y - pj.y;
    zij = pi.z - pj.z;

    mic(xij, yij, zij, side_);

    rijsq = xij * xij + yij * yij + zij * zij;
    rij = sqrt(rijsq);

    double rab[3] = {xij / rij, yij / rij, zij / rij};
    // ignore symmetry at the moment
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        Dab = (6.0 * a0_ / (8.0 * rij)) * (I[a][b] + rab[a] * rab[b] + (2.0 * pow(1.0, 2) / (3.0 * rijsq)) * (I[a][b] - 3.0 * rab[a] * rab[b]));
        gsl_matrix_set(Dij_, 3 * i + a, 3 * j + b, Dab);
      }
    }
  }
  // fill in diagonals
  for (int i = 0; i < this->size(); i++) {
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        Dab = 0.0;
        if (a == b) Dab = 1.0; // general case is Dab = a_0 / a_i = b_i
        gsl_matrix_set(Dij_, 3 * i + a, 3 * i + b, Dab);
      }
    }
  }
  // fill in the lower half
  for (int i = 0; i < this->size() - 1; i++) {
    for (int j = i + 1; j < this->size(); j++) {
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          Dab = gsl_matrix_get(Dij_, 3 * i + a, 3 * j + b);
          gsl_matrix_set(Dij_, 3 * j + b, 3 * i + a, Dab);
        }
      }
    }
  }
}

void ParticleCollection::computeExcludedVolumeForces(const PairList pairs_, double side_,
                                                     bool remove_overlap_) {
  double xij;
  double yij;
  double zij;
  double rij;
  double rijsq;
  double surface_to_surface;
  double rnx, rny, rnz;
  double fmagn_over_rij;
  double A = 50.0;
  double B = 10.0;

  for (PairList::const_iterator it = pairs_.begin(); it < pairs_.end(); it++) {
    Particle& pi = this->at(it->first);
    Particle& pj = this->at(it->second);

    xij = pi.x - pj.x;
    yij = pi.y - pj.y;
    zij = pi.z - pj.z;

    mic(xij, yij, zij, side_);

    rijsq = xij * xij + yij * yij + zij * zij;
    rij = sqrt(rijsq);
    surface_to_surface = rij - pi.radius - pj.radius;
    if (remove_overlap_ && surface_to_surface < 0.0) {
      rnx = xij / rij;
      rny = yij / rij;
      rnz = zij / rij;
      //std::cout << surface_to_surface << std::endl;
      //std::cout << xij << " " << yij << " " << zij << std::endl;
      //std::cout << pi.x << " " << pi.y << " " << pi.z << std::endl;
      //std::cout << pj.x << " " << pj.y << " " << pj.z << std::endl;
      pi.x -= 0.5 * surface_to_surface * rnx;
      pi.y -= 0.5 * surface_to_surface * rny;
      pi.z -= 0.5 * surface_to_surface * rnz;
      pj.x += 0.5 * surface_to_surface * rnx;
      pj.y += 0.5 * surface_to_surface * rny;
      pj.z += 0.5 * surface_to_surface * rnz;
      //std::cout << pi.x << " " << pi.y << " " << pi.z << std::endl;
      //std::cout << pj.x << " " << pj.y << " " << pj.z << std::endl;
      pi.updateSites(0.0, 0.0, 0.0, 0.0);
      pj.updateSites(0.0, 0.0, 0.0, 0.0);
      xij = (pi.radius + pj.radius) * rnx;
      yij = (pi.radius + pj.radius) * rny;
      zij = (pi.radius + pj.radius) * rnz;
      //std::cout << xij << " " << yij << " " << zij << std::endl;
      //std::cout << "--" << std::endl;
      rij = pi.radius + pj.radius;
      surface_to_surface = 0.0;
    }
    fmagn_over_rij = A * exp(-B * surface_to_surface) / rij;
    pi.fx += fmagn_over_rij * xij;
    pi.fy += fmagn_over_rij * yij;
    pi.fz += fmagn_over_rij * zij;
    pj.fx -= fmagn_over_rij * xij;
    pj.fy -= fmagn_over_rij * yij;
    pj.fz -= fmagn_over_rij * zij;
  }
  //if(rijmin < 0.0) std::cout << rijmin << std::endl;

  // compute force with boundaries
  /*
  ParticleCollection::iterator its;
  for (its = this->begin(); its < this->end(); its++) {
    Particle& p = *its;
    p.fx += pow(1.0 / (p.x - p.radius), 12);
    p.fx -= pow(1.0 / (side_ - p.x - p.radius), 12);
    p.fy += pow(1.0 / (p.y - p.radius), 12);
    p.fy -= pow(1.0 / (side_ - p.y - p.radius), 12);
    p.fz += pow(1.0 / (p.z - p.radius), 12);
    p.fz -= pow(1.0 / (side_ - p.z - p.radius), 12);
  }
  */
}

double gauss() {
  // Box-Muller algorithm
  double x1 = log(drand48());
  double x2 = 2.0 * (4.0 * atan(1.0)) * drand48();
  double y1 = sqrt(-2.0 * x1) * cos(x2);
  double y2 = sqrt(-2.0 * x1) * sin(x2);
  return (drand48() > 0.5) ? y1 : y2;
}

/* rc is the minimum surface-to-surface separation where the particles first see each other.
   Skin is violated when 0.5 * (rs - rc) > rmsd of any particle. Particles are added to Verlet
   list if the surface-to-suface distance is less than rs. */
void ParticleCollection::updatePositionsAndOrientations(double side_, double dt_, int verlet_, double rc_,
                                                        double rs_, bool& skin_violation_, const double a0_,
                                                        gsl_matrix* Dij_, bool hydrodynamics_) {
  std::vector<double> R_uncorrelated;
  std::vector<double> R_correlated;
  std::vector<double> junk;
  if (hydrodynamics_) {
    for (int i = 0; i < this->size(); i++) {
      Particle& pi = this->at(i);
      double DijFjx = 0.0;
      double DijFjy = 0.0;
      double DijFjz = 0.0;
      for (int j = 0; j < this->size(); j++) {
        DijFjx += gsl_matrix_get(Dij_, 3 * i + 0, 3 * j + 0) * this->at(j).fx +
                  gsl_matrix_get(Dij_, 3 * i + 0, 3 * j + 1) * this->at(j).fy +
                  gsl_matrix_get(Dij_, 3 * i + 0, 3 * j + 2) * this->at(j).fz;
        DijFjy += gsl_matrix_get(Dij_, 3 * i + 1, 3 * j + 0) * this->at(j).fx +
                  gsl_matrix_get(Dij_, 3 * i + 1, 3 * j + 1) * this->at(j).fy +
                  gsl_matrix_get(Dij_, 3 * i + 1, 3 * j + 2) * this->at(j).fz;
        DijFjz += gsl_matrix_get(Dij_, 3 * i + 2, 3 * j + 0) * this->at(j).fx +
                  gsl_matrix_get(Dij_, 3 * i + 2, 3 * j + 1) * this->at(j).fy +
                  gsl_matrix_get(Dij_, 3 * i + 2, 3 * j + 2) * this->at(j).fz;
      }
      pi.x += DijFjx * dt_;
      pi.y += DijFjy * dt_;
      pi.z += DijFjz * dt_;
      junk.push_back(DijFjx * dt_);
      junk.push_back(DijFjy * dt_);
      junk.push_back(DijFjz * dt_);
    }

    /* These functions factorize the symmetric, positive-definite square matrix A into
       the Cholesky decomposition A = L L^T (or A = L L^H for the complex case). On input,
       the values from the diagonal and lower-triangular part of the matrix A are used
       (the upper triangular part is ignored). On output the diagonal and lower triangular
       part of the input matrix A contain the matrix L, while the upper triangular part of
       the input matrix is overwritten with L^T (the diagonal terms being identical for
       both L and L^T). If the matrix is not positive-definite then the decomposition will
       fail, returning the error code GSL_EDOM. */

    gsl_linalg_cholesky_decomp(Dij_);
    // create and store random numbers in vector R
    double b_i = 1.0;
    for (int i = 0; i < 3 * this->size(); i++) {
      R_uncorrelated.push_back(sqrt(24.0 * b_i * dt_) * (drand48() - 0.5));
      double Rsum = 0.0;
      for (int j = 0; j <= i; j++)
        Rsum += gsl_matrix_get(Dij_, i, j) * R_uncorrelated[j];
      R_correlated.push_back(Rsum);
    }
  }

  double wmagn;
  double nwx;
  double nwy;
  double nwz;
  double dtheta;
  double rsq;

  double ek = 0.0;
  int i = 0;
  for (ParticleCollection::iterator it = this->begin(); it < this->end(); it++) {
    Particle& p = *it;
    if (hydrodynamics_) {
      p.x += R_correlated[3 * i + 0];
      p.y += R_correlated[3 * i + 1];
      p.z += R_correlated[3 * i + 2];
      ek += (pow(junk[3 + i + 0] + R_correlated[3 * i + 0], 2) +
             pow(junk[3 + i + 1] + R_correlated[3 * i + 1], 2) +
             pow(junk[3 + i + 2] + R_correlated[3 * i + 2], 2)) / pow(dt_, 2);
    }
    else {
      double rndx = sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      double rndy = sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      double rndz = sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      //p.x += p.fx * p.b * dt_ + sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      //p.y += p.fy * p.b * dt_ + sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      //p.z += p.fz * p.b * dt_ + sqrt(24.0 * p.b * dt_) * (drand48() - 0.5);
      p.x += p.fx * p.b * dt_ + rndx;
      p.y += p.fy * p.b * dt_ + rndy;
      p.z += p.fz * p.b * dt_ + rndz;
      ek += (pow(p.fx * p.b * dt_ + rndx, 2) +
             pow(p.fy * p.b * dt_ + rndy, 2) +
             pow(p.fz * p.b * dt_ + rndz, 2)) / pow(dt_, 2);
    }

    if(p.x < 0.0) {p.x += side_; p.xflag -= 1;}
    if(p.y < 0.0) {p.y += side_; p.yflag -= 1;}
    if(p.z < 0.0) {p.z += side_; p.zflag -= 1;}
    if(p.x > side_) {p.x -= side_; p.xflag += 1;}
    if(p.y > side_) {p.y -= side_; p.yflag += 1;}
    if(p.z > side_) {p.z -= side_; p.zflag += 1;}

    // check for a skin violation
    if (verlet_ > 0 && !skin_violation_) {
      rsq = pow(p.x + p.xflag * side_ - p.x_skin, 2) +
            pow(p.y + p.yflag * side_ - p.y_skin, 2) +
            pow(p.z + p.zflag * side_ - p.z_skin, 2);
      if (rsq > pow(0.5 * (rs_ - rc_), 2)) skin_violation_ = true;
    }

    // compute angular velocity
    double fac1 = (6.0 / 8.0) * pow(p.b, 3);
    double fac2 = sqrt(2.0 * (6.0 / 8.0) * pow(p.b, 3) / dt_);
    p.wx = fac1 * p.tx + fac2 * (drand48() - 0.5);
    p.wy = fac1 * p.ty + fac2 * (drand48() - 0.5);
    p.wz = fac1 * p.tz + fac2 * (drand48() - 0.5);

    // rotate the particle by an angle dtheta about the axis <nwx, nwy, nwz>
    wmagn = sqrt(p.wx * p.wx + p.wy * p.wy + p.wz * p.wz);
    nwx = p.wx / wmagn;
    nwy = p.wy / wmagn;
    nwz = p.wz / wmagn;
    dtheta = wmagn * dt_;

    // update orientation of site vectors and the positions of the sites
    p.updateSites(nwx, nwy, nwz, dtheta);

    // zero forces and torques for next time step
    p.fx = 0.0;
    p.fy = 0.0;
    p.fz = 0.0;
    p.tx = 0.0;
    p.ty = 0.0;
    p.tz = 0.0;

    /*
    p.q0 += 0.5 * dt_ * (-p.q1 * p.wxBody - p.q2 * p.wyBody - p.q3 * p.wzBody);
    p.q1 += 0.5 * dt_ * ( p.q0 * p.wxBody - p.q3 * p.wyBody + p.q2 * p.wzBody);
    p.q2 += 0.5 * dt_ * ( p.q3 * p.wxBody + p.q0 * p.wyBody - p.q1 * p.wzBody);
    p.q3 += 0.5 * dt_ * (-p.q2 * p.wxBody + p.q1 * p.wyBody + p.q0 * p.wzBody);
    p.normalizeQuaternions();
    p.updateRotationMatrix();
    */
    i++;
  }
  if (verlet_ > 0 && skin_violation_) {
    for (ParticleCollection::iterator it = this->begin(); it < this->end(); it++) {
      Particle& p = *it;
      p.updatePositionAtPairListBuild(side_);
    }
  }

  /*
  std::ofstream ekin;
  ekin.open("kinetic.dat", std::ios::app);
  ekin << ek / this->size() << std::endl;
  ekin.close();
  */

}

void ParticleCollection::writeBondStats(long long t_, double dt_, const BondList bonds_) {
  std::ofstream ratio;
  ratio.open("bondStats.dat", std::ios::app);
  ratio << t_ * dt_ << " " << bonds_.size() << " " << bonds_.bondsFormed << " " <<
                                                      bonds_.bondsBroken << std::endl;
  ratio.close();
}

void ParticleCollection::writeClusterSizes(long long t_, double dt_, const BondList bonds_,
                                           bool composition_) const {

  unsigned int tetrahedrons = 0;
  unsigned int cubes = 0;
  unsigned int dodecahedrons = 0;
  unsigned int truncatedIcosahedrons = 0;
  unsigned int pyramids24 = 0;
  unsigned int pyramids32 = 0;
  unsigned int pyramids60 = 0;
  unsigned int pyramids62 = 0;
  unsigned int pyramids64 = 0;
  unsigned int pyramids66 = 0;
  unsigned int box = 0;

  std::vector<int> assigned;
  unsigned int count[1024]; // allocate dynamically
  unsigned int ave_bonds[1024]; // allocate dynamically

  std::ofstream comp;
  if (composition_) {
    comp.open("clusterCompositions.dat", std::ios::app);
    comp << "t_ * dt_ = " << t_ * dt_ << std::endl;
    comp.close();
  }

  for (unsigned int i = 0; i < 1024; i++) { count[i] = 0; ave_bonds[i] = 0; }
 
  for (int i = 0; i < this->size(); i++) {
    if (std::find(assigned.begin(), assigned.end(), i) == assigned.end()) {
      std::vector<int> home;
      std::vector<int> cluster;
      assigned.push_back(i);
      home.push_back(i);
      cluster.push_back(i);
      bool more_neighbors = true;
      while (more_neighbors) {
        std::vector<int> neighbor;
        more_neighbors = false;
	for (std::vector<int>::iterator itr = home.begin(); itr < home.end(); itr++) {
	  int home_index = *itr;
	  BondList::const_iterator it;
	  for (it = bonds_.begin(); it < bonds_.end(); it++) {
	    if (it->first == home_index) neighbor.push_back(it->third);
	    if (it->third == home_index) neighbor.push_back(it->first);
          }
	}
        home.clear();
        for (int j = 0; j < neighbor.size(); j++) {
          if (std::find(cluster.begin(), cluster.end(), neighbor[j]) == cluster.end()) {
            cluster.push_back(neighbor[j]);
            assigned.push_back(neighbor[j]);
            home.push_back(neighbor[j]);
            more_neighbors = true;
          }
        }
      }
      std::sort(cluster.begin(), cluster.end());
      std::vector<int>::iterator its = std::unique(cluster.begin(), cluster.end());
      cluster.resize(its - cluster.begin());
      if (cluster.size() < 1024) count[cluster.size()]++;
      else std::cout << "ERROR: cluster size exceeded memory" << std::endl;

      // count specific structures and average number of bonds per cluster
      int bonds_per_cluster = 0;
      bool allThree = true;
      for (std::vector<int>::iterator itr = cluster.begin(); itr < cluster.end(); itr++) {
        if (this->at(*itr).bonds != 3) allThree = false;
        bonds_per_cluster += this->at(*itr).bonds;
      }
      bonds_per_cluster /= 2;

      if (cluster.size() == 4  && allThree) tetrahedrons++;
      if (cluster.size() == 8  && allThree) cubes++;
      if (cluster.size() == 20 && allThree) dodecahedrons++;
      if (cluster.size() == 60 && allThree) truncatedIcosahedrons++;
      if (cluster.size() == 14 && bonds_per_cluster == 24) pyramids24++;
      if (cluster.size() == 14 && bonds_per_cluster == 32) pyramids32++;
      if (cluster.size() == 30 && bonds_per_cluster == 60) pyramids60++;
      if (cluster.size() == 30 && bonds_per_cluster == 62) pyramids62++;
      if (cluster.size() == 30 && bonds_per_cluster == 64) pyramids64++;
      if (cluster.size() == 30 && bonds_per_cluster == 66) pyramids66++;
      if (cluster.size() == 63 && bonds_per_cluster == 114) box++;
      ave_bonds[cluster.size()] += bonds_per_cluster;

      // get composition by writing out the bonds of the particles in the cluster
      if (composition_ && cluster.size() > 1) {
        std::vector<std::pair<int, int> > bonded_pairs;
	for (std::vector<int>::iterator itr = cluster.begin(); itr < cluster.end(); itr++) {
	  BondList::const_iterator it;
	  for (it = bonds_.begin(); it < bonds_.end(); it++)
	    if (it->first == *itr) bonded_pairs.push_back(std::make_pair(it->first, it->third));
	}
	comp.open("clusterCompositions.dat", std::ios::app);
	std::vector<std::pair<int, int> >::const_iterator it;
	for (it = bonded_pairs.begin(); it < bonded_pairs.end(); it++)
	  comp << "(" << it->first << " " << it->second << ") ";
	comp << std::endl;
	comp.close();
      }
    }
  }

  // consistency check
  int psum = 0;
  for (int i = 0; i < 1024; i++) psum += count[i] * i;
  if (psum != this->size()) std::cout << "WARNING: cluster mismatch: " << this->size() << " " << psum << std::endl;

  // count plus size clusters
  int sumPlusSizes = 0;
  int max_cluster_size = 63;
  for (int i = max_cluster_size + 1; i < 1024; i++)
    sumPlusSizes += count[i];

  std::ofstream agg;
  agg.open("clusterSizes.dat", std::ios::app);
  agg << t_ * dt_ << " ";
  for (int i = 1; i <= max_cluster_size; i++)
    agg << count[i] << " ";
  agg << sumPlusSizes << std::endl;
  agg.close();

  agg.open("clusterStats.dat", std::ios::app);
  agg << t_ * dt_ << " " << tetrahedrons << " "
                         << cubes << " "
                         << dodecahedrons << " "
                         << truncatedIcosahedrons << " "
                         << pyramids24 + pyramids32 << " "
                         << pyramids60 << " "
                         << pyramids62 << " "
                         << pyramids64 << " "
                         << pyramids66 << " "
                         << box << std::endl;
  agg.close();
}

void ParticleCollection::writeClusterStats(long long t_, double dt_, const BondList bonds_) const {

  int unbound = 0;
  int dimers = 0;
  int linear_trimers = 0;
  int triangles = 0;
  int squares = 0;
  int corners = 0;
  int tetrahedrons = 0;
  int cubes = 0;

  // unbound
  ParticleCollection::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    const Particle p = *it;
    if (p.bonds == 0) unbound++;
  }

  // dimers
  BondList::const_iterator itr;
  for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
    const Particle pi = this->at(itr->first);
    const Particle pj = this->at(itr->third);
    if (pi.bonds == 1 && pj.bonds == 1) dimers++;
  }

  // linear trimers, triangles and squares
  for (int i = 0; i < this->size(); i++) {
    const Particle p = this->at(i);
    int partner[2];
    if (p.bonds == 2) {
      int ct = 0;
      for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
        if (itr->first == i) { partner[ct] = itr->third; ct++; }
        if (itr->third == i) { partner[ct] = itr->first; ct++; }
      }
      const Particle pj = this->at(partner[0]);
      const Particle pk = this->at(partner[1]);
      int other_j;
      int other_k;
      if (pj.bonds == 2 && pk.bonds == 2) {
	for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
	  if (itr->first == partner[0] && itr->third != i) other_j = itr->third;
	  if (itr->third == partner[0] && itr->first != i) other_j = itr->first;
	  if (itr->first == partner[1] && itr->third != i) other_k = itr->third;
	  if (itr->third == partner[1] && itr->first != i) other_k = itr->first;
	}
	if (other_j == partner[1] && other_k == partner[0]) {
	  triangles++;
	}
	else {
	  if (other_j == other_k && this->at(other_j).bonds == 2 && this->at(other_k).bonds == 2) squares++;
	}
      }
      else {
	if (pj.bonds == 1 && pk.bonds == 1) linear_trimers++;
      }
    }
  }

  // tetrahedron
  for (int i = 0; i < this->size(); i++) {
    const Particle p = this->at(i);
    int partner[3];
    if (p.bonds == 3) {
      int ct = 0;
      for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
        if (itr->first == i) { partner[ct] = itr->third; ct++; }
        if (itr->third == i) { partner[ct] = itr->first; ct++; }
      }
      const Particle pj = this->at(partner[0]);
      const Particle pk = this->at(partner[1]);
      const Particle pl = this->at(partner[2]);
      if (pj.bonds == 3 && pk.bonds == 3 && pl.bonds == 3) {
	std::vector<int> ids;
	for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
	  if (itr->first == partner[0] || itr->third == partner[0]) { ids.push_back(itr->first); ids.push_back(itr->third); }
	  if (itr->first == partner[1] || itr->third == partner[1]) { ids.push_back(itr->first); ids.push_back(itr->third); }
	  if (itr->first == partner[2] || itr->third == partner[2]) { ids.push_back(itr->first); ids.push_back(itr->third); }
	}
        std::sort(ids.begin(), ids.end());
	std::vector<int>::iterator its = std::unique(ids.begin(), ids.end());
	ids.resize(its - ids.begin());
	if (ids.size() == 4) tetrahedrons++;
      }
      else {
        if (pj.bonds == 1 && pk.bonds == 1 && pl.bonds == 1) corners++;
      }
    }
  }

  // cubes
  for (int i = 0; i < this->size(); i++) {
    const Particle p = this->at(i);
    int neigh[3];
    if (p.bonds == 3) {
      int ct = 0;
      for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
        if (itr->first == i) { neigh[ct] = itr->third; ct++; }
        if (itr->third == i) { neigh[ct] = itr->first; ct++; }
      }
      const Particle pj = this->at(neigh[0]);
      const Particle pk = this->at(neigh[1]);
      const Particle pl = this->at(neigh[2]);
      int neigh_of_neigh[6];
      if (pj.bonds == 3 && pk.bonds == 3 && pl.bonds == 3) {
	std::vector<int> ids;
	int ct = 0;
	for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
          for (int m = 0; m < 3; m++) {
	    if (itr->first == neigh[m] && itr->third != i) { neigh_of_neigh[ct] = itr->third; ct++; }
	    if (itr->third == neigh[m] && itr->first != i) { neigh_of_neigh[ct] = itr->first; ct++; }
	    if (itr->first == neigh[m] || itr->third == neigh[m]) { ids.push_back(itr->first); ids.push_back(itr->third); }
          }
	}
	const Particle p0 = this->at(neigh_of_neigh[0]);
	const Particle p1 = this->at(neigh_of_neigh[1]);
	const Particle p2 = this->at(neigh_of_neigh[2]);
	const Particle p3 = this->at(neigh_of_neigh[3]);
	const Particle p4 = this->at(neigh_of_neigh[4]);
	const Particle p5 = this->at(neigh_of_neigh[5]);
	if(p0.bonds == 3 && p1.bonds == 3 && p2.bonds == 3 && p3.bonds == 3 && p4.bonds == 3 && p5.bonds == 3) {
	  int ct = 0;
	  for (itr = bonds_.begin(); itr < bonds_.end(); itr++) {
	    for (int m = 0; m < 6; m++)
	      if (itr->first == neigh_of_neigh[m] || itr->third == neigh_of_neigh[m]) { ids.push_back(itr->first); ids.push_back(itr->third); }
	  }
          std::sort(ids.begin(), ids.end());
	  std::vector<int>::iterator its = std::unique(ids.begin(), ids.end());
	  ids.resize(its - ids.begin());
	  if (ids.size() == 8) {
	    bool threeBonds = true;
	    for (int n = 0; n < 8; n++) {
	      const Particle pn = this->at(ids.at(n));
	      if (pn.bonds != 3) threeBonds = false;
	    }
	    if (threeBonds) cubes++;
          }
	}
      }
    }
  }

  std::ofstream agg;
  agg.open("clusterStats.dat", std::ios::app);
  agg << t_ * dt_ << " " << unbound << " "
                         << dimers << " "
                         << linear_trimers << " "
                         << triangles / 3 << " "
                         << squares / 4 << " "
                         << corners << " "
                         << tetrahedrons / 4 << " "
                         << cubes / 8 << std::endl;
  agg.close();
}

void ParticleCollection::writeConfiguration(long long t_, long long steps_, long long freqWrite_,
                                            double dt_, double side_, double k_spring_,
                                            BondList bonds_, BondList unbreakable_, const double a0_) const {
  // find number of particle and site types
  size_t num_particle_types;
  size_t num_site_types;
  std::vector<size_t> ptypes;
  std::vector<size_t> stypes;
  ParticleCollection::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    Particle p = *it;
    if (std::count(ptypes.begin(), ptypes.end(), p.type) == 0) ptypes.push_back(p.type);
    SiteCollection::const_iterator its;
    SiteCollection sc = p.sites;
    for (its = sc.begin(); its < sc.end(); its++) {
      Site s = *its;
      if (std::count(stypes.begin(), stypes.end(), s.type) == 0) stypes.push_back(s.type);
    }
  }
  num_particle_types = ptypes.size();
  num_site_types = stypes.size();

  // convert int to string and make file name
  std::stringstream ss;
  ss << "N" << this->size() << "P" << num_particle_types << "S" << num_site_types << "." << t_;
  std::string s(ss.str());
  const char* fname = s.c_str();
 
  // write out the file
  bool explosion = false;
  std::ofstream config;
  config.open(fname);
  config << "<step>" << t_ << "</step>" << std::endl;
  config << "<steps>" << steps_ << "</steps>" << std::endl;
  config << "<timestep>" << dt_ << "</timestep>" << std::endl;
  config << "<numberOfParticles>" << this->size() << "</numberOfParticles>" << std::endl;
  config << "<numberOfParticleTypes>" << num_particle_types << "</numberOfParticleTypes>" << std::endl;
  config << "<numberOfSiteTypes>" << num_site_types << "</numberOfSiteTypes>" << std::endl;
  config << "<a0>" << a0_ << "</a0>" << std::endl;
  config << "<boxSide>" << side_ << "</boxSide>" << std::endl;
  config << std::endl;
  for (it = this->begin(); it < this->end(); it++) {
    Particle p = *it;
    config << "<particle>" << std::endl;
    config << "  <id>" << p.id << "</id>" << std::endl;
    config << "  <type>" << p.type << "</type>" << std::endl;
    config << "  <radius>" << p.radius << "</radius>" << std::endl;
    config << "  <x>" << p.x << "</x>" << std::endl;
    config << "  <y>" << p.y << "</y>" << std::endl;
    config << "  <z>" << p.z << "</z>" << std::endl;

    // test for numerical explosion
    if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) explosion = true;

    SiteCollection sc = p.sites;
    config << "  <numberOfSites>" << sc.size() << "</numberOfSites>" << std::endl;
    SiteCollection::const_iterator its;
    for (its = sc.begin(); its < sc.end(); its++) {
      Site s = *its;
      config << "  <site>" << std::endl;
      config << "    <id>" << s.id << "</id>" << std::endl;
      config << "    <type>" << s.type << "</type>" << std::endl;
      config << "    <bonded>" << s.bonded << "</bonded>" << std::endl;
      config << "    <x>" << p.x + s.dx << "</x>" << std::endl;
      config << "    <y>" << p.y + s.dy << "</y>" << std::endl;
      config << "    <z>" << p.z + s.dz << "</z>" << std::endl;
      config << "  </site>" << std::endl;
    }
    config << "</particle>" << std::endl;
  }
  config << std::endl;
  if (bonds_.size() != 0) config << "<numberOfBonds>" << bonds_.size() << "</numberOfBonds>" << std::endl;
  if (unbreakable_.size() != 0) config << "<numberOfExcludedBonds>" << unbreakable_.size() << "</numberOfExcludedBonds>" << std::endl;
  if (bonds_.size() != 0 || unbreakable_.size() != 0) config << "<bonds>" << std::endl;
  if (bonds_.size() != 0) {
    for (BondList::const_iterator it = bonds_.begin(); it < bonds_.end(); it++) {
      Quad q = *it;
      config << "  <bond>" << std::endl;
      config << "    <i>" << q.first << "</i>" << std::endl;
      config << "    <a>" << q.second << "</a>" << std::endl;
      config << "    <j>" << q.third << "</j>" << std::endl;
      config << "    <b>" << q.fourth << "</b>" << std::endl;
      config << "  </bond>" << std::endl;
    }
  }
  if (unbreakable_.size() != 0) {
    for (BondList::const_iterator it = unbreakable_.begin(); it < unbreakable_.end(); it++) {
      Quad q = *it;
      config << "  <exclude>" << std::endl;
      config << "    <i>" << q.first << "</i>" << std::endl;
      config << "    <a>" << q.second << "</a>" << std::endl;
      config << "    <j>" << q.third << "</j>" << std::endl;
      config << "    <b>" << q.fourth << "</b>" << std::endl;
      config << "  </exclude>" << std::endl;
    }
  }
  if (bonds_.size() != 0 || unbreakable_.size() != 0) config << "</bonds>" << std::endl;
  config.close();
  if (explosion) { std::cout << "NaN detected. Exiting ..." << std::endl; exit(1); }
}

void ParticleCollection::writeBondLifeTimes(VectorLL& life_times_, double dt_) const {
  // bond life times in units of tau_t
  if (!life_times_.empty()) {
    std::ofstream blt;
    blt.open("lifeTimes.dat", std::ios::app);
    for (VectorLL::iterator it = life_times_.begin(); it < life_times_.end(); it++)
      blt << (*it).first * dt_ << " " << (*it).second << std::endl;
    blt.close();
    life_times_.clear();
  }
}

template <typename T>
T StringToNumber (const std::string& Text) {
  std::stringstream ss(Text);
  T result;
  return ss >> result ? result : 0;
}

std::string stripValue(const std::string& line_) {
  size_t start = line_.find_first_of(">");
  size_t end = line_.find_last_of("<");
  return line_.substr(start + 1, end - start);
}

std::string strip2ndValue(std::string line_) {
  std::string::iterator end_pos = std::remove(line_.begin(), line_.end(), ' ');
  line_.erase(end_pos, line_.end());
  return line_;
}

void ParticleCollection::readScript(const char* f_, long long& step_, long long& steps_,
                                    long long& freqWrite_,
                                    double& dt_, double& side_, double& k_spring_,
                                    int& verlet_, double& rc_, double& rs_,
                                    BondList& bonds_, unsigned int& seed_, bool& restart_,
                                    BondList& unbreakable_, double xcp_[], double Lcp_[],
                                    size_t& max_site_type_, const double a0_) {
  double L_;
  double expArg_;
  std::string line;
  std::string value;

  std::ifstream fname(f_);
  if (fname.is_open()) {
    std::cout << std::endl;
    std::cout << std::endl;
    while (fname.good()) {
      getline(fname, line);
      if (line.size() > 0) {
        if (strip2ndValue(line).compare(0, 1, "#") == 0) {
          std::cout << "Ignoring comment: " << line << std::endl; }
        else if (line.find("step") != std::string::npos && line.find("time") == std::string::npos) {
          value = strip2ndValue(line.substr(line.find("p") + 1, line.length() - line.find("p") - 1));
          step_ = StringToNumber<long long>(value); }
        else if (line.find("timestep") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("p") + 1, line.length() - line.find("p") - 1));
          dt_ = StringToNumber<double>(value); }
        else if (line.find("run") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("n") + 1, line.length() - line.find("n") - 1));
          steps_ = StringToNumber<long long>(value); }
        else if (line.find("seed") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("d") + 1, line.length() - line.find("d") - 1));
          seed_ = StringToNumber<int>(value); }
        else if (line.find("Lij") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("j") + 1, line.length() - line.find("j") - 1));
          L_ = StringToNumber<double>(value);
          for (int i = 0; i < 512; i++) Lcp_[i] = L_; }
        else if (line.find("Lcp") != std::string::npos) {
          std::istringstream iss(line);
          for (int i = 0; i < 3; i++) {
            std::string sub;
            iss >> sub;
            if (i == 1) value = sub;
            if (i == 2) Lcp_[StringToNumber<int>(value)] = StringToNumber<double>(sub); } }
        else if (line.find("config") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("g") + 1, line.length() - line.find("g") - 1));
          readConfiguration(value.c_str(), step_, side_, bonds_, unbreakable_, a0_); }
        else if (line.find("kspring") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("g") + 1, line.length() - line.find("g") - 1));
          k_spring_ = StringToNumber<double>(value); }
        else if (line.find("rc") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("c") + 1, line.length() - line.find("c") - 1));
          rc_ = StringToNumber<double>(value); }
        else if (line.find("rs") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("s") + 1, line.length() - line.find("s") - 1));
          rs_ = StringToNumber<double>(value); }
        else if (line.find("xij") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("j") + 1, line.length() - line.find("j") - 1));
          expArg_ = StringToNumber<double>(value);
          for (int i = 0; i < 512; i++) xcp_[i] = expArg_; }
        else if (line.find("bonding") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("g") + 1, line.length() - line.find("g") - 1));
          expArg_ = StringToNumber<double>(value);
          for (int i = 0; i < 512; i++) xcp_[i] = expArg_; }
        else if (line.find("write") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("e") + 1, line.length() - line.find("e") - 1));
          freqWrite_ = StringToNumber<long long>(value); }
        else if (line.find("restart") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("rt") + 2, line.length() - line.find("rt") - 2));
          if (value.compare("no") == 0) restart_ = false;
          if (value.compare("yes") == 0) restart_ = true; }
        else if (line.find("verlet") != std::string::npos) {
          value = strip2ndValue(line.substr(line.find("t") + 1, line.length() - line.find("t") - 1));
          verlet_ = StringToNumber<int>(value); }
        else if (line.find("xcp") != std::string::npos) {
          std::istringstream iss(line);
          for (int i = 0; i < 3; i++) {
            std::string sub;
	    iss >> sub;
            if (i == 1) value = sub;
            if (i == 2) xcp_[StringToNumber<int>(value)] = StringToNumber<double>(sub); } }
        else
          std::cout << "WARNING: " << line << " was ignored." << std::endl;
      }
    }
    fname.close();

    // write parameters to standard output
    std::cout << std::endl;
    std::cout << "step = " << step_ << std::endl;
    std::cout << "steps = " << steps_ << std::endl;
    if (step_ != 0) std::cout << "end step = " << step_ + steps_ << std::endl;
    std::cout << "freqWrite = " << freqWrite_ << std::endl;
    std::cout << "dt = " << dt_ << std::endl;
    std::cout << "k_spring = " << k_spring_ << std::endl;
    std::cout << "rc = " << rc_ << std::endl;
    std::cout << "rs = " << rs_ << std::endl;
    std::cout << "xij = " << expArg_ << std::endl;
    std::cout << "Lij = " << L_ << std::endl;
    std::cout << "N = " << this->size() << std::endl;
    std::cout << "verlet = " << verlet_ << std::endl;
    std::cout << "restart = " << restart_ << std::endl;
    std::cout << "seed = " << seed_ << std::endl;
    std::cout << "max Brownian step = " << 0.5 * sqrt(24.0 * dt_) << std::endl;
    std::cout << "max spring step = " << k_spring_ * sqrt(24.0 * dt_) * dt_ << std::endl;
    std::cout << "number of bonds = " << bonds_.size() << std::endl;
    std::cout << "number of exclusions = " << unbreakable_.size() << std::endl;

    bool nonstandard = false;
    for (int i = 0; i < 512; i++) if (xcp_[i] != expArg_) nonstandard = true;
    if (nonstandard) {
      std::cout << "xcp:" << std::endl;
      for (int i = 0; i < 512; i++)
        if (xcp_[i] != expArg_) std::cout << " " << i << " " << xcp_[i] << std::endl;
    }
    nonstandard = false;
    for (int i = 0; i < 512; i++) if (Lcp_[i] != L_) nonstandard = true;
    if (nonstandard) {
      std::cout << "Lcp:" << std::endl;
      for (int i = 0; i < 512; i++)
      if (Lcp_[i] != L_) std::cout << " " << i << " " << Lcp_[i] << std::endl;
    }

    // find maximum site type and check for consecutive types
    int zero_type = 0;
    int site_types_sum = 0;
    std::vector<int> stypes;
    ParticleCollection::const_iterator it;
    for (it = this->begin(); it < this->end(); it++) {
      Particle p = *it;
      SiteCollection::const_iterator its;
      SiteCollection sc = p.sites;
      for (its = sc.begin(); its < sc.end(); its++) {
        Site s = *its;
        site_types_sum += s.type;
        if (std::count(stypes.begin(), stypes.end(), s.type) == 0) { stypes.push_back(s.type); }
        if (s.type == 0) { zero_type = 1; std::cout << "WARNING: site type of 0 encountered." << std::endl; }
      }
    }
    if (site_types_sum != 0) std::cout << "WARNING: site types do not sum to 0 but " << site_types_sum << std::endl;
    max_site_type_ = stypes.size() / 2;
    std::cout << "max_site_type = " << max_site_type_ << std::endl;

    std::sort(stypes.begin(), stypes.end());
    std::cout << "site types = ";
    for (int i = 0; i < stypes.size(); i++)
      std::cout << stypes.at(i) << " ";
    std::cout << std::endl;

    try {
      if (zero_type == 1) {
	int j = 0;
	for (int i = stypes.front(); i <= stypes.back(); i++) {
	  if (stypes.at(j) != i) std::cout << "ERROR: site types not consecutive: " << i << " " << stypes.at(j) << std::endl;
	  j++;
	}
      }
      else {
	int j = 0;
	for (int i = stypes.front(); i <= stypes.back(); i++) {
	  if (i != 0) {
	    if (stypes.at(j) != i) std::cout << "ERROR: site types not consecutive: " << i << " " << stypes.at(j) << std::endl;
	    j++;
	  }
	}
      }
    }
    catch (...) {
      std::cout << "ERRORHANDLER: site types are not continuous." << std::endl;
      exit(1);
    }
  }
  else {
    std::cout << "ERROR: unable to open script file." << std::endl;
    exit(1);
  }
  if (!restart_) {
    // delete files
    std::remove("bondStats.dat");
    std::remove("clusterStats.dat");
    std::remove("clusterSizes.dat");
    std::remove("clusterCompositions.dat");
    std::remove("averageBondsPerCluster.dat");
    std::remove("lifeTimes.dat");
  }

  // write out particle indices and types
  std::ofstream comp;
  comp.open("clusterCompositions.dat", std::ios::app);
  comp << "# particle_index  particle_type" << std::endl;
  int i = 0;
  ParticleCollection::const_iterator itr;
  for (itr = this->begin(); itr < this->end(); itr++) {
    comp << i << " " << itr->type << std::endl;
    i++;
  }
  comp.close();
}

void ParticleCollection::readConfiguration(const char* f_, long long& step_, double& side_, BondList& bonds_,
                                           BondList& unbreakable_, const double a0_) {
  std::string line;
  std::string tag;
  std::string value;
  size_t num_particles = 0;
  size_t num_particle_types = 0;
  size_t num_site_types = 0;
  size_t num_bonds = 0;
  int id;
  int type;
  double radius;
  double x;
  double y;
  double z;
  int pi;
  int sa;
  int pj;
  int sb;
  int bonded;
  SiteCollection sites;

  std::ifstream fname(f_);
  if (fname.is_open()) {
    while (fname.good()) {
      getline(fname, line);
      if (line.size() > 0) {
        tag = line.substr(line.find("<"), line.find(">") - line.find("<") + 1);
        if (tag == "<step>") {
          value = stripValue(line);
          step_ = StringToNumber<long long>(value); }
        else if (tag == "<steps>") {
          value = stripValue(line);}
          //steps_ = StringToNumber<int>(value); }
        else if (tag == "<freqWrite>") {
          value = stripValue(line);}
          //freqWrite_ = StringToNumber<int>(value); }
        else if (tag == "<timestep>") {
          value = stripValue(line);}
          //dt_ = StringToNumber<double>(value); }
        else if (tag == "<boxSide>") {
          value = stripValue(line);
          side_ = StringToNumber<double>(value); }
        else if (tag == "<averageBondLength>") {
          value = stripValue(line);}
          //L_ = StringToNumber<double>(value); }
        else if (tag == "<springConstant>") {
          value = stripValue(line);}
          //k_spring_ = StringToNumber<double>(value); }
        else if (tag == "<expArg>") {
          value = stripValue(line);}
          //expArg_ = StringToNumber<double>(value); }
        else if (tag == "<numberOfParticles>") {
          value = stripValue(line);
          num_particles = StringToNumber<size_t>(value); }
        else if (tag == "<numberOfParticleTypes>") {
          value = stripValue(line);
          num_particle_types = StringToNumber<size_t>(value); }
        else if (tag == "<numberOfSiteTypes>") {
          value = stripValue(line);
          num_site_types= StringToNumber<size_t>(value); }
        else if (tag == "<id>") {
          value = stripValue(line);
          id = StringToNumber<int>(value); }
        else if (tag == "<type>") {
          value = stripValue(line);
          type = StringToNumber<int>(value); }
        else if (tag == "<radius>") {
          value = stripValue(line);
          radius = StringToNumber<double>(value); }
        else if (tag == "<x>") {
          value = stripValue(line);
          x = StringToNumber<double>(value); }
        else if (tag == "<y>") {
          value = stripValue(line);
          y = StringToNumber<double>(value); }
        else if (tag == "<z>") {
          value = stripValue(line);
          z = StringToNumber<double>(value); }
        else if (tag == "<numberOfSites>") {
          Particle p(id, type, radius, x, y, z, a0_);
          this->push_back(p); }
        else if (tag == "<bonded>") {
          value = stripValue(line);
          bonded = StringToNumber<int>(value); }
        else if (tag == "</site>") {
          Site s(id, type, bonded, this->back().x,
                                   this->back().y,
                                   this->back().z,
                                   x - this->back().x,
                                   y - this->back().y,
                                   z - this->back().z);
          sites.push_back(s); }
        else if (tag == "</particle>") {
          this->back().sites = sites;
          sites.clear(); }
        else if (tag == "<numberOfBonds>") {
          value = stripValue(line);
          num_bonds = StringToNumber<size_t>(value); }
        else if (tag == "<i>") {
          value = stripValue(line);
          pi = StringToNumber<int>(value); }
        else if (tag == "<a>") {
          value = stripValue(line);
          sa = StringToNumber<int>(value); }
        else if (tag == "<j>") {
          value = stripValue(line);
          pj = StringToNumber<int>(value); }
        else if (tag == "<b>") {
          value = stripValue(line);
          sb = StringToNumber<int>(value); }
        else if (tag == "</bond>") {
          Quad q(pi, sa, pj, sb);
          this->at(pi).bonds++;
          this->at(pj).bonds++;
          this->at(pi).sites.at(sa).bonded = true;
          this->at(pj).sites.at(sb).bonded = true;
          bonds_.push_back(q); }
        else if (tag == "</exclude>") {
          Quad q(pi, sa, pj, sb);
          unbreakable_.push_back(q); }
        else if (tag == "<particle>" || tag == "<site>" || tag == "<bond>" ||
                 tag == "<bonds>" || tag == "</bonds>" || tag == "<exclude>" ||
                 tag == "<numberOfExcludedBonds>") {
          ;}
        else
          std::cout << "WARNING: " << tag << ": tag in configuration file was ignored." << std::endl;

      }
    }
    fname.close();
    if(num_particles != this->size())
      std::cout << "WARNING: number of particles does not match:" << num_particles << " vs. " << this->size() << std::endl;
    if(num_bonds != bonds_.size())
      std::cout << "WARNING: number of bonds does not match." << std::endl;
  }
  else {
    std::cout << "ERROR: unable to open configuration file." << std::endl;
    exit(1);
  }
}

// create a simple system with 27 particles and two sites per particle
void ParticleCollection::initialize(double& side_, const double a0_) {
  side_ = 10.0;
  int ct = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        Particle p(ct, 0, 1.0, 3.0 * i + 2.0, 3.0 * j + 2.0, 3.0 * k + 2.0, a0_);
        this->push_back(p);
        ct++;
      }
    }
  }
  double nx;
  double ny;
  double nz;
  ParticleCollection::iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    Particle& p = *it;
    SiteCollection sites;
    Site s0(0, 0, 0, p.x, p.y, p.z,  p.radius, 0.0, 0.0);
    Site s1(1, 0, 0, p.x, p.y, p.z, -p.radius, 0.0, 0.0);
    sites.push_back(s0);
    sites.push_back(s1);
    p.sites = sites;
    nx = 0.5 - drand48();
    ny = 0.5 - drand48();
    nz = 0.5 - drand48();
    nx = nx / sqrt(nx * nx + ny * ny + nz * nz);
    ny = ny / sqrt(nx * nx + ny * ny + nz * nz);
    nz = nz / sqrt(nx * nx + ny * ny + nz * nz);
    p.updateSites(nx, ny, nz, 0.5 * 3.1415 * drand48());
  }
}

void ParticleCollection::printPositions(long long t_) const {
  std::cout << " step = " << t_ << std::endl;
  std::cout << " id  x  y  z" << std::endl;
  std::vector<Particle>::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    std::cout << it->id << " " << it->x << " " << it->y << " " << it->z << std::endl;
  }
}

void ParticleCollection::printSiteCollection(long long t_) const {
  std::cout << " step = " << t_ << std::endl;
  std::cout << " id  x  y  z" << std::endl;
  std::vector<Particle>::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    std::cout << "particle x y z" << std::endl;
    std::cout << it->id << " " << it->x << " " << it->y << " " << it->z << std::endl;
    std::cout << "sites dx dy dz" << std::endl;
    std::vector<Site>::const_iterator its;
    for (its = it->sites.begin(); its < it->sites.end(); its++) {
      std::cout << its->id << " " << its->dx << " " << its->dy << " " << its->dz << std::endl;
    }
  }
}

void ParticleCollection::printMSD(long long t_, double side_, double dt_) const {
  double x_final;
  double y_final;
  double z_final;
  double msd;
  double msd_ave = 0.0;
  std::vector<Particle>::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    x_final = it->x + it->xflag * side_;
    y_final = it->y + it->yflag * side_;
    z_final = it->z + it->zflag * side_;
    msd = pow(it->x_initial - x_final, 2) + pow(it->y_initial - y_final, 2) + pow(it->z_initial - z_final, 2);
    msd_ave += msd;
  }
  std::cout << t_ * dt_ << " " << msd_ave / this->size() << std::endl;
}

void ParticleCollection::printCosTheta(long long t_, double dt_) const {
  double thetaSq = 0.0;
  double costheta = 0.0;
  double costhetaSq = 0.0;
  std::vector<Particle>::const_iterator it;
  for (it = this->begin(); it < this->end(); it++) {
    Particle p = *it;
    thetaSq += pow(acos(p.sites[0].dz / p.radius), 2);
    costheta += p.sites[0].dz / p.radius;
    costhetaSq += 0.5 * (3.0 * pow(p.sites[0].dz / p.radius, 2) - 1.0);
  }
  std::cout << t_ * dt_ << " " << costheta / this->size() << " " <<
               thetaSq / this->size() << " " << costhetaSq / this->size() << std::endl;
}

void ParticleCollection::printBonds(long long t_, double dt_, const BondList& bonds_) const {
  std::cout << "step = " << t_ << std::endl;
  BondList::const_iterator it;
  for (BondList::const_iterator it = bonds_.begin(); it < bonds_.end(); it++) {
    Quad q = *it;
    std::cout << q.first << " " << q.second << " " << q.third << " " << q.fourth << std::endl;
  }
  std::cout << std::endl;
}
