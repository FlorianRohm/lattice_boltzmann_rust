#![allow(dead_code)]
pub mod boundaries;
pub mod constants;

use crate::lbm::boundaries::BoundaryType;
use crate::lbm::constants::{C, C_F64, C_SRC, T};
use ndarray::prelude::{s, Axis};
use ndarray::{Array2, Zip};

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum CollisionType {
    BGK { omega: f64 },
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum NodeType {
    Boundary(BoundaryType),
    Inward(CollisionType),
    Ghost,
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Node {
    pub velocities: [f64; 9],
    pub node_type: NodeType,
}

impl Node {
    pub fn collide(mut self) -> Self {
        self.velocities = match self.node_type {
            NodeType::Boundary(boundary) => boundary.collide(self.velocities),
            NodeType::Inward(collision_type) => collision_type.collide(self.velocities),
            NodeType::Ghost => self.velocities,
        };

        self
    }
}

impl CollisionType {
    pub fn collide(self, velocities: [f64; 9]) -> [f64; 9] {
        match self {
            Self::BGK { omega } => bgk(velocities, omega),
        }
    }
}

// bgk collision term
pub fn bgk(mut f_pop: [f64; 9], omega: f64) -> [f64; 9] {
    let (rho, ux, uy) = compute_macros(&f_pop);
    let u_sqr = ux * ux + uy * uy;

    f_pop[0] = f_pop[0] * (1. - omega) + omega * compute_equilibrium(0, rho, ux, uy, u_sqr);
    f_pop[1] = f_pop[1] * (1. - omega) + omega * compute_equilibrium(1, rho, ux, uy, u_sqr);
    f_pop[2] = f_pop[2] * (1. - omega) + omega * compute_equilibrium(2, rho, ux, uy, u_sqr);
    f_pop[3] = f_pop[3] * (1. - omega) + omega * compute_equilibrium(3, rho, ux, uy, u_sqr);
    f_pop[4] = f_pop[4] * (1. - omega) + omega * compute_equilibrium(4, rho, ux, uy, u_sqr);
    f_pop[5] = f_pop[5] * (1. - omega) + omega * compute_equilibrium(5, rho, ux, uy, u_sqr);
    f_pop[6] = f_pop[6] * (1. - omega) + omega * compute_equilibrium(6, rho, ux, uy, u_sqr);
    f_pop[7] = f_pop[7] * (1. - omega) + omega * compute_equilibrium(7, rho, ux, uy, u_sqr);
    f_pop[8] = f_pop[8] * (1. - omega) + omega * compute_equilibrium(8, rho, ux, uy, u_sqr);

    f_pop
}

pub fn compute_equilibrium(i_pop: usize, rho: f64, ux: f64, uy: f64, u_sqr: f64) -> f64 {
    let c_u = C_F64[i_pop][0] * ux + C_F64[i_pop][1] * uy;
    rho * T[i_pop] * (1. + 3. * c_u + 4.5 * c_u * c_u - 1.5 * u_sqr)
}

pub fn ini_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u_sqr = ux * ux + uy * uy;
    [
        compute_equilibrium(0, rho, ux, uy, u_sqr),
        compute_equilibrium(1, rho, ux, uy, u_sqr),
        compute_equilibrium(2, rho, ux, uy, u_sqr),
        compute_equilibrium(3, rho, ux, uy, u_sqr),
        compute_equilibrium(4, rho, ux, uy, u_sqr),
        compute_equilibrium(5, rho, ux, uy, u_sqr),
        compute_equilibrium(6, rho, ux, uy, u_sqr),
        compute_equilibrium(7, rho, ux, uy, u_sqr),
        compute_equilibrium(8, rho, ux, uy, u_sqr),
    ]
}

pub fn compute_macros(f: &[f64; 9]) -> (f64, f64, f64) {
    let upper_line = f[2] + f[5] + f[6];
    let medium_line = f[0] + f[1] + f[3];
    let lower_line = f[4] + f[7] + f[8];
    let rho = upper_line + medium_line + lower_line;
    let ux = (f[1] + f[5] + f[8] - (f[3] + f[6] + f[7])) / rho;
    let uy = (upper_line - lower_line) / rho;

    (rho, ux, uy)
}

#[derive(PartialEq, Debug)]
pub struct Simulation {
    pub l_x: i32,
    pub l_y: i32,
    pub lattice: Array2<Node>,
    pub tmp_lattice: Array2<Node>,
}

macro_rules! inner_slice {
    ( $x:expr ) => {
        s!(1..=$x.l_x, 1..=$x.l_y)
    };
}

impl Simulation {
    pub fn collide(mut self) -> Self {
        self.lattice
            // .slice_mut(inner_slice!(self))
            .map_inplace(|node| *node = node.collide());
        self
    }

    pub fn propagate_push(self) -> Self {
        self.propagate_push_inner().make_periodic_push()
    }

    pub fn propagate_pull(self) -> Self {
        self.make_periodic_pull().propagate_pull_inner()
    }

    fn propagate_pull_inner(mut self) -> Self {
        let source = self.lattice.windows((3, 3));
        let mut target = self.tmp_lattice.slice_mut(inner_slice!(self));
        Zip::from(&mut target).and(source).for_each(|tar, src| {
            for i_pop in 0..9 {
                let c = C_SRC[i_pop];
                let vel = src[[c[0], c[1]]].velocities[i_pop];
                tar.velocities[i_pop] = vel;
            }
        });

        std::mem::swap(&mut self.tmp_lattice, &mut self.lattice);

        self
    }

    pub fn propagate_push_inner(mut self) -> Self {
        for i_x in 1..=self.l_x as usize {
            for i_y in 1..=self.l_y as usize {
                for i_pop in 0..9 {
                    let next_x = (i_x as i32 + C[i_pop][0]) as usize;
                    let next_y = (i_y as i32 + C[i_pop][1]) as usize;
                    self.tmp_lattice[[next_x, next_y]].velocities[i_pop] =
                        self.lattice[[i_x, i_y]].velocities[i_pop];
                }
            }
        }

        std::mem::swap(&mut self.tmp_lattice, &mut self.lattice);

        self
    }

    // implement periodic boundary conditions (to be called after
    //   the propagation step)
    fn make_periodic_push(mut self) -> Self {
        let l_x: usize = self.l_x as usize;
        let l_y: usize = self.l_y as usize;
        let lat = &mut self.lattice;

        for i_x in 1..=l_x {
            lat[[i_x, l_y]].velocities[4] = lat[[i_x, 0]].velocities[4];
            lat[[i_x, l_y]].velocities[7] = lat[[i_x, 0]].velocities[7];
            lat[[i_x, l_y]].velocities[8] = lat[[i_x, 0]].velocities[8];
            lat[[i_x, 1]].velocities[2] = lat[[i_x, l_y + 1]].velocities[2];
            lat[[i_x, 1]].velocities[5] = lat[[i_x, l_y + 1]].velocities[5];
            lat[[i_x, 1]].velocities[6] = lat[[i_x, l_y + 1]].velocities[6];
        }

        for i_y in 1..=l_y {
            lat[[1, i_y]].velocities[1] = lat[[l_x + 1, i_y]].velocities[1];
            lat[[1, i_y]].velocities[5] = lat[[l_x + 1, i_y]].velocities[5];
            lat[[1, i_y]].velocities[8] = lat[[l_x + 1, i_y]].velocities[8];
            lat[[l_x, i_y]].velocities[3] = lat[[0, i_y]].velocities[3];
            lat[[l_x, i_y]].velocities[6] = lat[[0, i_y]].velocities[6];
            lat[[l_x, i_y]].velocities[7] = lat[[0, i_y]].velocities[7];
        }

        lat[[1, 1]].velocities[5] = lat[[l_x + 1, l_y + 1]].velocities[5];
        lat[[l_x, 1]].velocities[6] = lat[[0, l_y + 1]].velocities[6];
        lat[[l_x, l_y]].velocities[7] = lat[[0, 0]].velocities[7];
        lat[[1, l_y]].velocities[8] = lat[[l_x + 1, 0]].velocities[8];

        self
    }

    fn make_periodic_pull(mut self) -> Self {
        let l_x: usize = self.l_x as usize;
        let l_y: usize = self.l_y as usize;
        let lat = &mut self.lattice;

        for i_x in 1..=l_x {
            lat[[i_x, 0]].velocities[2] = lat[[i_x, l_y]].velocities[2];
            lat[[i_x, 0]].velocities[5] = lat[[i_x, l_y]].velocities[5];
            lat[[i_x, 0]].velocities[6] = lat[[i_x, l_y]].velocities[6];
            lat[[i_x, l_y + 1]].velocities[4] = lat[[i_x, 1]].velocities[4];
            lat[[i_x, l_y + 1]].velocities[7] = lat[[i_x, 1]].velocities[7];
            lat[[i_x, l_y + 1]].velocities[8] = lat[[i_x, 1]].velocities[8];
        }

        for i_y in 1..=l_y {
            lat[[l_x + 1, i_y]].velocities[3] = lat[[1, i_y]].velocities[3];
            lat[[l_x + 1, i_y]].velocities[6] = lat[[1, i_y]].velocities[6];
            lat[[l_x + 1, i_y]].velocities[7] = lat[[1, i_y]].velocities[7];
            lat[[0, i_y]].velocities[1] = lat[[l_x, i_y]].velocities[1];
            lat[[0, i_y]].velocities[5] = lat[[l_x, i_y]].velocities[5];
            lat[[0, i_y]].velocities[8] = lat[[l_x, i_y]].velocities[8];
        }

        lat[[l_x + 1, l_y + 1]].velocities[7] = lat[[1, 1]].velocities[7];
        lat[[0, l_y + 1]].velocities[8] = lat[[l_x, 1]].velocities[8];
        lat[[0, 0]].velocities[5] = lat[[l_x, l_y]].velocities[5];
        lat[[l_x + 1, 0]].velocities[6] = lat[[1, l_y]].velocities[6];
        self
    }
}

#[allow(dead_code)]
pub fn save_vel(sim: &Simulation, f_name: &str) {
    use std::io::Write;
    let mut buffer = std::fs::File::create(f_name).expect("should not fail");

    sim.lattice
        .slice(inner_slice!(sim))
        .lanes(Axis(0))
        .into_iter()
        .for_each(|lane| {
            lane.iter().for_each(|f| {
                let (_, ux, uy) = compute_macros(&f.velocities);
                let u_norm = (ux * ux + uy * uy).sqrt();
                write!(buffer, "{:.6} ", u_norm).expect("writing should not fail");
            });
            writeln!(buffer).expect("writing should not fail");
        });
}

#[allow(dead_code)]
pub fn save_field(sim: &Array2<Node>, f_name: &str) {
    use std::io::Write;
    let mut buffer = std::fs::File::create(f_name).expect("should not fail");

    sim.lanes(Axis(0)).into_iter().for_each(|lane| {
        lane.iter().for_each(|f| {
            let vel = &f.velocities;
            write!(
                buffer,
                "[{:.9}\n {:.9}\n {:.9}\n {:.9}\n {:.9}\n {:.9}\n {:.9}\n {:.9}\n {:.9}]\n\n",
                vel[0], vel[1], vel[2], vel[3], vel[4], vel[5], vel[6], vel[7], vel[8]
            )
            .expect("writing should not fail");
        });
        write!(buffer, "\n\n").expect("writing should not fail");
    });
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ndarray::arr2;

    #[test]
    pub fn test_propagate_pushpull_1x1() {
        let node_type = NodeType::Ghost;
        let node = Node {
            velocities: [0.4; 9],
            node_type,
        };
        let node_zero = Node {
            velocities: [0.; 9],
            node_type,
        };
        let lattice = arr2(&[[node, node, node], [node, node, node], [node, node, node]]);
        let lattice2 = arr2(&[[node, node, node], [node, node, node], [node, node, node]]);

        let tmp_lattice = arr2(&[
            [node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero],
        ]);
        let tmp_lattice2 = arr2(&[
            [node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero],
        ]);

        let mut sim_pull = Simulation {
            l_x: 1,
            l_y: 1,
            lattice,
            tmp_lattice,
        };
        let mut sim_push = Simulation {
            l_x: 1,
            l_y: 1,
            lattice: lattice2,
            tmp_lattice: tmp_lattice2,
        };

        sim_pull = sim_pull.propagate_pull();
        sim_push = sim_push.propagate_push();

        assert_eq!(
            sim_pull.lattice.slice(inner_slice!(sim_pull)),
            sim_push.lattice.slice(inner_slice!(sim_push))
        );
    }

    #[test]
    pub fn test_propagate_pushpull_3x3_complete() {
        let node_type = NodeType::Ghost;

        let nodes: Vec<Node> = (0..9)
            .into_iter()
            .map(|a: i32| Node {
                velocities: [a as f64 / 10.; 9],
                node_type,
            })
            .collect();

        let node_zero = Node {
            velocities: [0.; 9],
            node_type,
        };
        let lattice = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, nodes[0], nodes[1], nodes[2], node_zero],
            [node_zero, nodes[3], nodes[4], nodes[5], node_zero],
            [node_zero, nodes[6], nodes[7], nodes[8], node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ]);
        let lattice2 = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, nodes[0], nodes[1], nodes[2], node_zero],
            [node_zero, nodes[3], nodes[4], nodes[5], node_zero],
            [node_zero, nodes[6], nodes[7], nodes[8], node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ]);

        let tmp_lattice = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ]);
        let tmp_lattice2 = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ]);

        let mut sim_pull = Simulation {
            l_x: 3,
            l_y: 3,
            lattice,
            tmp_lattice,
        };
        let mut sim_push = Simulation {
            l_x: 3,
            l_y: 3,
            lattice: lattice2,
            tmp_lattice: tmp_lattice2,
        };

        sim_pull = sim_pull.propagate_pull().propagate_pull().propagate_pull();
        sim_push = sim_push.propagate_push().propagate_push().propagate_push();

        print_slice("pull", sim_pull.lattice.slice(inner_slice!(sim_pull)));
        print_slice(
            "push, correct",
            sim_push.lattice.slice(inner_slice!(sim_push)),
        );

        assert_eq!(
            sim_pull.lattice.slice(inner_slice!(sim_pull)),
            sim_push.lattice.slice(inner_slice!(sim_push))
        );
    }

    #[test]
    pub fn test_propagate_pushpull_3x3_single_step() {
        let node_type = NodeType::Ghost;

        let nodes: Vec<Node> = (0..9)
            .into_iter()
            .map(|a: i32| Node {
                velocities: [a as f64 / 10.; 9],
                node_type,
            })
            .collect();

        let node_zero = Node {
            velocities: [0.; 9],
            node_type,
        };
        let lattice = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, nodes[0], nodes[1], nodes[2], node_zero],
            [node_zero, nodes[3], nodes[4], nodes[5], node_zero],
            [node_zero, nodes[6], nodes[7], nodes[8], node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ])
        .reversed_axes();
        let lattice2 = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, nodes[0], nodes[1], nodes[2], node_zero],
            [node_zero, nodes[3], nodes[4], nodes[5], node_zero],
            [node_zero, nodes[6], nodes[7], nodes[8], node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ])
        .reversed_axes();

        let tmp_lattice = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ])
        .reversed_axes();
        let tmp_lattice2 = arr2(&[
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
            [node_zero, node_zero, node_zero, node_zero, node_zero],
        ])
        .reversed_axes();

        let mut sim_pull = Simulation {
            l_x: 3,
            l_y: 3,
            lattice,
            tmp_lattice,
        };
        let mut sim_push = Simulation {
            l_x: 3,
            l_y: 3,
            lattice: lattice2,
            tmp_lattice: tmp_lattice2,
        };

        sim_pull = sim_pull.propagate_pull();
        sim_push = sim_push.propagate_push();

        print_slice("pull", sim_pull.lattice.slice(inner_slice!(sim_pull)));
        print_slice(
            "push, correct",
            sim_push.lattice.slice(inner_slice!(sim_push)),
        );

        assert_eq!(
            sim_pull.lattice.slice(inner_slice!(sim_pull)),
            sim_push.lattice.slice(inner_slice!(sim_push))
        );
    }

    #[allow(dead_code)]
    fn print_lattice(title: &str, sim: &Simulation) {
        println!("{}", title);
        sim.lattice.lanes(Axis(0)).into_iter().for_each(|lane| {
            lane.iter().for_each(|f| {
                let vel = &f.velocities;
                print!(
                    "[{:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1}] ",
                    vel[0], vel[1], vel[2], vel[3], vel[4], vel[5], vel[6], vel[7], vel[8]
                )
            });
            print!("\n");
        });
    }

    fn print_slice(title: &str, slice: ndarray::ArrayBase<ndarray::ViewRepr<&Node>, ndarray::Ix2>) {
        println!("{}", title);
        slice.lanes(Axis(0)).into_iter().for_each(|lane| {
            lane.iter().for_each(|f| {
                let vel = &f.velocities;
                print!(
                    "[{:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1} {:.1}] ",
                    vel[0], vel[1], vel[2], vel[3], vel[4], vel[5], vel[6], vel[7], vel[8]
                )
            });
            print!("\n");
        });
    }
}
