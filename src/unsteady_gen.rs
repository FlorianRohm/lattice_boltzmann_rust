#![warn(clippy::all, clippy::nursery)]

mod generic;
mod lbm;

use crate::generic::RectangleSimulation;
use crate::generic::{
    compute_macros, BounceBack, Collidable, D2Q9Node, Ghost, LeftRegularized, LowerRegularized,
    RightPressureRegularized, Simulation, UpperRegularized, BGK, D2Q9, save_vel
};
use crate::lbm::ini_equilibrium;

pub struct Settings {
    l_x: usize,
    l_y: usize,
    obst_x: usize,
    obst_y: usize,
    obst_r: usize,
    u_max: f64,
    omega: f64,
    max_t: i32,
}

// compute parabolic Poiseuille profile
pub fn compute_poiseuille(i: usize, u_max: f64, l_y: usize) -> f64 {
    let y: f64 = (i - 1) as f64;
    let l: f64 = (l_y - 1) as f64;
    4. * u_max / l.powi(2) * (l * y - y.powi(2))
}

pub enum NodeTypes {
    BGK(BGK),
    LeftRegularized(LeftRegularized<BGK>),
    RightPressureRegularized(RightPressureRegularized<BGK>),
    UpperRegularized(UpperRegularized<BGK>),
    LowerRegularized(LowerRegularized<BGK>),
    Ghost(Ghost),
    BounceBack(BounceBack),
}

impl Collidable for NodeTypes {
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        match self {
            NodeTypes::BGK(x) => x.collide(vel),
            NodeTypes::LeftRegularized(x) => x.collide(vel),
            NodeTypes::RightPressureRegularized(x) => x.collide(vel),
            NodeTypes::UpperRegularized(x) => x.collide(vel),
            NodeTypes::LowerRegularized(x) => x.collide(vel),
            NodeTypes::Ghost(x) => x.collide(vel),
            NodeTypes::BounceBack(x) => x.collide(vel),
        }
    }
}

pub type Node = D2Q9Node<NodeTypes>;

pub fn construct_sim_gen(settings: &Settings) -> RectangleSimulation<Node> {
    let generation_function = |(i, j): (usize, usize)| {
        let (l_x, l_y) = (settings.l_x, settings.l_y);

        let outer_edge =
            |(i, j): (usize, usize)| (i == 0 || i == l_x + 1 || j == 0 || j == l_y + 1);
        let obstacle = |(i, j): (usize, usize)| {
            (i - settings.obst_x).pow(2) + (j - settings.obst_y).pow(2) <= settings.obst_r.pow(2)
        };

        let u_poiseuille: f64 = compute_poiseuille(j, settings.u_max, l_y);
        let equibrilium = ini_equilibrium(1., u_poiseuille, 0.);
        let collision_type = BGK {
            omega: settings.omega,
        };      
        let collision_node = NodeTypes::BGK(BGK {
            omega: settings.omega,
        });

        match (i, j) {
            coords if outer_edge(coords) => Node {
                velocities: [0.; 9],
                node_type: NodeTypes::Ghost(Ghost()),
            },
            (_, y) if y == 1 => Node {
                velocities: equibrilium,
                node_type: NodeTypes::LowerRegularized(LowerRegularized {
                    ux: 0.,
                    uy: 0.,
                    collision_type,
                }),
            },
            (_, y) if y == l_y => Node {
                velocities: equibrilium,
                node_type: NodeTypes::UpperRegularized(UpperRegularized {
                    ux: 0.,
                    uy: 0.,
                    collision_type,
                }),
            },
            (x, _) if x == 1 => Node {
                velocities: equibrilium,
                node_type: NodeTypes::LeftRegularized(LeftRegularized {
                    ux: u_poiseuille,
                    uy: 0.,
                    collision_type,
                }),
            },
            (x, _) if x == l_x => Node {
                velocities: equibrilium,
                node_type: NodeTypes::RightPressureRegularized(RightPressureRegularized {
                    rho: 1.,
                    u_par: 0.,
                    collision_type,
                }),
            },
            coords if obstacle(coords) => Node {
                velocities: equibrilium,
                node_type: NodeTypes::BounceBack(BounceBack()),
            },
            _ => Node {
                velocities: equibrilium,
                node_type: collision_node,
            },
        }
    };

    RectangleSimulation::from_function(settings.l_x, settings.l_y, generation_function)
}

fn update_zero_gradient_boundary(settings: &Settings, sim: &mut RectangleSimulation<Node>) {
    let l_x = settings.l_x as usize;
    let l_y = settings.l_y as usize;
    for i_y in 2..l_y as usize {
        let third_to_last = &sim.lattice[[(l_x - 1), i_y]];
        let second_to_last = &sim.lattice[[(l_x - 2), i_y]];

        let (rho1, _ux1, _uy1) = compute_macros(&third_to_last.velocities);
        let (rho2, _ux2, _uy2) = compute_macros(&second_to_last.velocities);

        let last: &mut Node = &mut sim.lattice[[l_x, i_y]];

        if let NodeTypes::RightPressureRegularized(RightPressureRegularized {
            ref mut rho,
            ref mut u_par,
            ..
        }) = last.node_type
        {
            *rho = 4. / 3. * rho1 - 1. / 3. * rho2;
            *u_par = 0.;
        } else {
            panic!("this function only works if the last column are RightPressureRegularized")
        }
    }
}
pub fn main() {
    // channel length
    let l_x = 250;

    // channel height
    let l_y = 50;

    // maximum velocity of the Poiseuille inflow
    let u_max = 0.02;

    // position of the cylinder, the cylinder is
    let obst_x = l_x / 5;

    // offset from the center to break symmetry
    let obst_y = l_y / 2;

    // radius of the cylinder
    let obst_r = l_y / 10 + 1;

    // Reynolds number
    let re = 100.;

    // kinematic fluid viscosity
    let nu = u_max * 2.0f64 * obst_r as f64 / re;

    // total number of iterations
    let max_t = 10_000;

    // relaxation parameter
    let omega = 1. / (3. * nu + 1. / 2.);

    let settings = Settings {
        l_x,
        l_y,
        obst_x,
        obst_y,
        obst_r,
        u_max,
        omega,
        max_t,
    };

    let mut sim = construct_sim_gen(&settings);
    println!(
        "\nlx={}, ly={}, omega={:.6}\n\n",
        settings.l_x, settings.l_y, settings.omega
    );

    for i_t in 0..settings.max_t {
        if i_t % 1000 == 0 {
            println!("t={}", i_t);
        }
        update_zero_gradient_boundary(&settings, &mut sim);

        // on the right boundary, outlet condition grad_x u = 0

        sim = sim.collide();

        // By default: periodic boundary conditions. In this case,
        //   this is important, because the upper and lower
        //   boundaries are horizontally periodic, so that no
        //   special corner nodes are needed.
        sim = sim.stream();
    }
    save_vel(&sim, "rust_result_gen.dat");

}
