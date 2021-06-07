#![warn(
    clippy::all,
    clippy::nursery,
)]

mod lbm;
mod generic;

use crate::lbm::boundaries::BoundaryType;
use crate::lbm::{ save_vel, Node, NodeType};
use crate::lbm::{
    compute_macros, ini_equilibrium, CollisionType, Simulation
};
use ndarray::Array;

pub struct Settings {
    l_x: i32,
    l_y: i32,
    obst_x: i32,
    obst_y: i32,
    obst_r: i32,
    u_max: f64,
    omega: f64,
    max_t: i32,
}

// compute parabolic Poiseuille profile
pub fn compute_poiseuille(i: i32, u_max: f64, l_y: i32) -> f64 {
    let y: f64 = (i - 1) as f64;
    let l: f64 = (l_y - 1) as f64;
    4. * u_max / l.powi(2) * (l * y - y.powi(2))
}

pub fn construct_sim(settings: &Settings) -> Simulation {
    let shape = (250 + 2, 50 + 2);

    let lattice = Array::from_shape_fn(shape, |(i, j)| {
        let x = i as i32;
        let y = j as i32;
        let (l_x, l_y) = (settings.l_x, settings.l_y);

        let outer_edge = |(i, j): (i32, i32)| (i == 0 || i == l_x + 1 || j == 0 || j == l_y + 1);
        let obstacle = |(i, j): (i32, i32)| (i - settings.obst_x).pow(2) + (j - settings.obst_y).pow(2) <= settings.obst_r.pow(2);

        let u_poiseuille: f64 = compute_poiseuille(y, settings.u_max, l_y);
        let equibrilium = ini_equilibrium(1., u_poiseuille, 0.);
        let collision_type = CollisionType::BGK {
            omega: settings.omega,
        };

        let node_type = match (x, y) {
            coords if outer_edge(coords) => NodeType::Ghost,
            (_, y) if y == 1 => NodeType::Boundary(BoundaryType::LowerRegularized { ux: 0., uy: 0., collision_type, }),
            (_, y) if y == l_y => NodeType::Boundary(BoundaryType::UpperRegularized { ux: 0., uy: 0., collision_type, }),
            (x, _) if x == 1 => NodeType::Boundary(BoundaryType::LeftRegularized { ux: u_poiseuille, uy: 0., collision_type, }),
            (x, _) if x == l_x => NodeType::Boundary(BoundaryType::RightPressureRegularized { rho: 1., u_par: 0., collision_type, }),
            coords if obstacle(coords) => NodeType::Boundary(BoundaryType::BounceBack), 
            _ => NodeType::Inward(collision_type), 
        };
        
        match node_type {
            node_type @ NodeType:: Ghost => Node{velocities: [0.;9], node_type},
            node_type => Node { velocities: equibrilium, node_type,   }
        }
    });

    let tmp_lattice = lattice.clone();
    
    Simulation {
        tmp_lattice,
        lattice,
        l_x: settings.l_x,
        l_y: settings.l_y,
    }
}

// Compute a second order extrapolation on the right boundary to
// ensure a zero-gradient boundary condition on the pressure.
// This must be recomputed at every time step. The velocity is
// constrained to be perpendicular to the outflow surface.
pub fn update_zero_gradient_boundary(settings: &Settings, sim: &mut Simulation) {
    let l_x = settings.l_x as usize;
    let l_y = settings.l_y as usize;
    for i_y in 2..l_y as usize {
        let third_to_last = &sim.lattice[[(l_x - 1), i_y]];
        let second_to_last = &sim.lattice[[(l_x - 2), i_y]];

        let (rho1, _ux1, _uy1) = compute_macros(&third_to_last.velocities);
        let (rho2, _ux2, _uy2) = compute_macros(&second_to_last.velocities);

        let last = &mut sim.lattice[[l_x, i_y]];

        if let NodeType::Boundary(BoundaryType::RightPressureRegularized {
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
    let omega =  1. / (3. * nu + 1. / 2.); 


    let settings = Settings {l_x, l_y, obst_x, obst_y, obst_r, u_max, omega, max_t};

    let mut sim = construct_sim(&settings);
    println!("\nlx={}, ly={}, omega={:.6}\n\n", settings.l_x, settings.l_y, settings.omega);

    for i_t in 0..settings.max_t {
        if i_t % 1000 == 0 {
            println!("t={}", i_t);
        }

        // on the right boundary, outlet condition grad_x u = 0
        update_zero_gradient_boundary(&settings, &mut sim);

        sim = sim.collide();

        // By default: periodic boundary conditions. In this case,
        //   this is important, because the upper and lower
        //   boundaries are horizontally periodic, so that no
        //   special corner nodes are needed.
        sim = sim.propagate_pull();
    }
    save_vel(&sim, "rust_result.dat");
}
