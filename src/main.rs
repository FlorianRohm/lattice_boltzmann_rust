extern crate lattice_boltzmann;

use lattice_boltzmann::Cell;
use lattice_boltzmann::Config;
use lattice_boltzmann::getGrid;
use lattice_boltzmann::iterateTimes;

pub fn main() {
  let config = Config {
    size_x: 2,
    size_y: 3,
  };

  let initialGrid = getGrid(&config);
  let grid = getGrid(&config);

  let iterated = iterateTimes(initialGrid, grid, 4);

  println!("{:?}", iterated.0);

}
