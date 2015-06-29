#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "ParticleCollection.hpp"
#include "BondList.hpp"
#include "PairList.hpp"
#include "CellCollection.hpp"
#include "config.hpp"
#include <gsl/gsl_matrix.h>

int main(int argc, char* argv[]) {

  // check for input file
  if (argc == 1) {
    std::cout << "Pass the input script as a command line argument. Exiting ...\n" << std::endl;
    exit(1);
  }

  // input file name as command-line argument
  const char* f = argv[1];

  // default system parameters
  long long step = 0;
  long long steps = 50000000;
  long long freqWrite = 1000000;
  unsigned int seed = 123456789;
  double a0 = 1.0; // radius with dimensions
  double dt = 0.001;
  double side = 10.0;
  double k_spring = 50.0;
  double rc = 1.0;
  double rs = 5.0;
  int verlet = 1; // 0 all pairs, 1 verlet, 2 verlet and cells
  bool restart = false;
  bool skin_violation = true;
  bool remove_overlap = false;
  bool hydrodynamics = false;
  bool composition = false;

  ParticleCollection particles;
  PairList pairs;
  BondList bonds;
  BondList unbreakable;
  VectorLL life_times;
  CellCollection cells;
  double xcp[512];
  double Lcp[512];
  size_t max_site_type;

  // automatic skin thickness updating
  clock_t start = clock();
  clock_t end;
  double elapsed;
  double elapsed_prev = 0;
  double rs_sign;
  double elapsed_sign;
  double ave_elapsed_time = 0.0;
  long intervals = 0;

  gsl_matrix* Dij;
  Dij = gsl_matrix_alloc(9, 9); // 3N x 3N

  particles.readScript(f, step, steps, freqWrite, dt, side, k_spring, verlet, rc, rs,
                       bonds, seed, restart, unbreakable, xcp, Lcp, max_site_type, a0);
  (seed == 123456789) ? srand(time(0)) : srand(seed);
  long long end_step = step + steps;
  particles.writeConfiguration(step, end_step, freqWrite, dt, side, k_spring, bonds, unbreakable, a0);
  cells.init(side, rc, rs, verlet);
  double rs_next = rs;
  double rs_prev = rs - 0.05;
  particles.buildPairs(pairs, cells, side, verlet, rs, rc, skin_violation, rs_next);
  particles.updateBondList(pairs, side, dt, bonds, unbreakable, xcp, Lcp, step, life_times, max_site_type);

  long long t;
  for (t = step + 1; t <= end_step; t++) {
    particles.computeExcludedVolumeForces(pairs, side, remove_overlap);
    particles.computeSpringForcesAndTorques(side, bonds, Lcp, k_spring);
    if (hydrodynamics) particles.updateDiffusionTensor(pairs, side, a0, Dij);
    /*
    for (int i=0;i<9;i++) {
      for (int j=0;j<9;j++)
        printf("%f ",gsl_matrix_get(Dij,i,j));
      printf("\n"); }
    printf("\n");
    */
    particles.updatePositionsAndOrientations(side, dt, verlet, rc, rs, skin_violation, a0, Dij, hydrodynamics);
    if (verlet > 0 && skin_violation) particles.buildPairs(pairs, cells, side, verlet, rs, rc, skin_violation, rs_next);
    particles.updateBondList(pairs, side, dt, bonds, unbreakable, xcp, Lcp, t, life_times, max_site_type);

    if (t % freqWrite == 0) particles.writeConfiguration(t, end_step, freqWrite, dt, side, k_spring, bonds, unbreakable, a0);
    if (t % 100000 == 0) particles.writeBondStats(t, dt, bonds);
    if (t <= 100000000 && t % 10000  == 0) particles.writeClusterSizes(t, dt, bonds, composition);
    if (t >  100000000 && t % 100000 == 0) particles.writeClusterSizes(t, dt, bonds, composition);
    if (t % 10000 == 0) particles.writeBondLifeTimes(life_times, dt);
    if (t % 10000 == 0) {
      end = clock();
      elapsed = double(end - start);
      if (elapsed == elapsed_prev) elapsed_sign = 1.0;
      else elapsed_sign = (elapsed - elapsed_prev) / fabs(elapsed - elapsed_prev);
      if (rs == rs_prev) rs_sign = 1.0;
      else rs_sign = (rs - rs_prev) / fabs(rs - rs_prev);
      rs_prev = rs;
      rs_next = rs - 0.05 * rs_sign * elapsed_sign;
      if (rs_next < rc + 1.0) rs_next = rc + 1.0;
      if (rs_next > 0.5 * side) rs_next = 0.5 * side;
      elapsed_prev = elapsed;
      start = clock();
      ave_elapsed_time += elapsed;
      intervals++;
    }
  }
  particles.writeConfiguration(t - 1, end_step, freqWrite, dt, side, k_spring, bonds, unbreakable, a0);
  particles.writeBondLifeTimes(life_times, dt);
 
  gsl_matrix_free(Dij);

  std::cout << std::endl;
  std::cout << "number of bonds = " << bonds.size() << std::endl;
  std::cout << "bonds formed = " << bonds.bondsFormed << ", broken = " << bonds.bondsBroken << std::endl;
  std::cout << "number of builds = " << pairs.rebuilds << std::endl;
  std::cout << "steps per build = " << steps / double(pairs.rebuilds) << std::endl;
  std::cout << "neighbors per particle (final list) = " << pairs.size() / double(particles.size()) << std::endl;
  std::cout << "average elapsed time for 10000 steps = " << ave_elapsed_time / (intervals * CLOCKS_PER_SEC) << " s" << std::endl;
  return 0;
}
