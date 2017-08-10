#[macro_use]
extern crate ndarray;

use ndarray::prelude::*;

#[derive(Debug)]
pub struct Config {
  pub size_x: usize,
  pub size_y: usize,
}

#[derive(Debug)]
#[derive(Clone)]
pub struct Cell {
  u: f64,
  v: f64,
  rho: f64,
  f_00: f64,
  f_p0: f64,
  f_n0: f64,
  f_0p: f64,
  f_pp: f64,
  f_np: f64,
  f_0n: f64,
  f_pn: f64,
  f_nn: f64,
}
impl Default for Cell {
  fn default() -> Cell {
    Cell {
      u: 0.,
      v: 0.,
      rho: 0.,
      f_00: 0.,
      f_p0: 0.,
      f_n0: 0.,
      f_0p: 0.,
      f_pp: 0.,
      f_np: 0.,
      f_0n: 0.,
      f_pn: 0.,
      f_nn: 0.,
    }
  }
}

pub type Grid = Array2<Cell>;

pub fn getGrid(config: &Config) -> Grid {
  let cell = Cell {
    ..Default::default()
  };
  Grid::from_elem((config.size_x, config.size_y), cell)
}

pub fn update(current: Grid, mut next: Grid) -> (Grid, Grid) {
  (next, current)
}

pub fn iterateTimes(mut current: Grid, mut next: Grid, times: usize) -> (Grid, Grid) {
  let mut tupel = update(current, next);

  next = tupel.0;
  current = tupel.1;

  println!("The value of x is: {:?}", next);
  println!("The value of y is: {:?}", current);;

  for _ in 1..times {
    tupel = update(current, next);

    next = tupel.0;
    current = tupel.1;

    println!("The value of x is: {:?}", next);
    println!("The value of y is: {:?}", current);
  }
  (next, current)
}
