#include "CellCollection.hpp"

long neighbor(double x_, double y_, double z_, double side_,
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

void CellCollection::init(double side_, double rc_, double rs_, int verlet_) {
  this->clear();

  if (verlet_ < 2) cells_per_dimension = 1;
  else cells_per_dimension = int(side_ / (rc_ + rs_));
  cell_side = side_ / cells_per_dimension;
  /*std::cout << "Cells: " << cells_per_dimension << " x " <<
                            cells_per_dimension << " x " <<
                            cells_per_dimension << std::endl; */
  // initialize CellCollection with the right number of elements
  for (long k = 0; k < pow(cells_per_dimension, 3); k++) {
    Cell c;
    this->push_back(c);
  }
  // apply stencil to each cell to get neighbors
  for (int i = 0; i < cells_per_dimension; i++) {
    for (int j = 0; j < cells_per_dimension; j++) {
      for (int k = 0; k < cells_per_dimension; k++) {
	double x = 0.5 * cell_side + cell_side * i;
	double y = 0.5 * cell_side + cell_side * j;
	double z = 0.5 * cell_side + cell_side * k;
	Cell c;
	c.neighbor_ids[0]  = neighbor(x + cell_side, y            , z            , side_, cell_side, cells_per_dimension);
	c.neighbor_ids[1]  = neighbor(x            , y + cell_side, z            , side_, cell_side, cells_per_dimension);
	c.neighbor_ids[2]  = neighbor(x + cell_side, y + cell_side, z            , side_, cell_side, cells_per_dimension);
	c.neighbor_ids[3]  = neighbor(x + cell_side, y - cell_side, z            , side_, cell_side, cells_per_dimension);
	c.neighbor_ids[4]  = neighbor(x + cell_side, y            , z - cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[5]  = neighbor(x            , y + cell_side, z - cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[6]  = neighbor(x + cell_side, y + cell_side, z - cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[7]  = neighbor(x + cell_side, y - cell_side, z - cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[8]  = neighbor(x + cell_side, y            , z + cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[9]  = neighbor(x            , y + cell_side, z + cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[10] = neighbor(x + cell_side, y + cell_side, z + cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[11] = neighbor(x + cell_side, y - cell_side, z + cell_side, side_, cell_side, cells_per_dimension);
	c.neighbor_ids[12] = neighbor(x            , y            , z + cell_side, side_, cell_side, cells_per_dimension);
        long cell_id = neighbor(x, y, z, side_, cell_side, cells_per_dimension);
        this->at(cell_id) = c;
      }
    }
  }
}
