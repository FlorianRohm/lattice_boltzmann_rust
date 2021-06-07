#![allow(dead_code)]

use ndarray::{s, Array2, Zip, Axis};

pub trait Collidable: Sized {
    fn collide(&self, vel: D2Q9) -> D2Q9;
}

pub trait BulkCollide: Collidable + Sized {}
pub const C_F64: [[f64; 2]; 9] = [
    [0., 0.],
    [1., 0.],
    [0., 1.],
    [-1., 0.],
    [0., -1.],
    [1., 1.],
    [-1., 1.],
    [-1., -1.],
    [1., -1.],
];
pub const C_INDEX: [[usize; 2]; 9] = [
    [1, 1],
    [0, 1],
    [1, 0],
    [2, 1],
    [1, 2],
    [0, 0],
    [2, 0],
    [2, 2],
    [0, 2],
];
pub const OPPOSITE_OF: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
pub const C: [[i32; 2]; 9] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];
pub const T: [f64; 9] = [
    4. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
];

#[rustfmt::skip]
pub struct BGK { pub omega: f64}
impl Collidable for BGK {
    fn collide(&self, mut vel: D2Q9) -> D2Q9 {
        let (rho, ux, uy) = compute_macros(&vel);
        let omega = self.omega;
        let u_sqr = ux * ux + uy * uy;

        vel[0] = vel[0] * (1. - omega) + omega * compute_equilibrium(0, rho, ux, uy, u_sqr);
        vel[1] = vel[1] * (1. - omega) + omega * compute_equilibrium(1, rho, ux, uy, u_sqr);
        vel[2] = vel[2] * (1. - omega) + omega * compute_equilibrium(2, rho, ux, uy, u_sqr);
        vel[3] = vel[3] * (1. - omega) + omega * compute_equilibrium(3, rho, ux, uy, u_sqr);
        vel[4] = vel[4] * (1. - omega) + omega * compute_equilibrium(4, rho, ux, uy, u_sqr);
        vel[5] = vel[5] * (1. - omega) + omega * compute_equilibrium(5, rho, ux, uy, u_sqr);
        vel[6] = vel[6] * (1. - omega) + omega * compute_equilibrium(6, rho, ux, uy, u_sqr);
        vel[7] = vel[7] * (1. - omega) + omega * compute_equilibrium(7, rho, ux, uy, u_sqr);
        vel[8] = vel[8] * (1. - omega) + omega * compute_equilibrium(8, rho, ux, uy, u_sqr);

        vel
    }
}
impl BulkCollide for BGK {}

pub struct BounceBack();
impl Collidable for BounceBack {
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        bounce_back(vel)
    }
}
pub struct Ghost();
impl Collidable for Ghost {
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        vel
    }
}

pub struct LeftRegularized<T>
where
    T: BulkCollide + Sized,
{
    pub ux: f64,
    pub uy: f64,
    pub collision_type: T,
}

impl<T> Collidable for LeftRegularized<T>
where
    T: BulkCollide + Sized,
{
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        left_regularized(vel, self.ux, self.uy, &self.collision_type)
    }
}
pub struct LowerRegularized<T>
where
    T: BulkCollide + Sized,
{
    pub ux: f64,
    pub uy: f64,
    pub collision_type: T,
}

impl<T> Collidable for LowerRegularized<T>
where
    T: BulkCollide + Sized,
{
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        lower_regularized(vel, self.ux, self.uy, &self.collision_type)
    }
}

pub struct UpperRegularized<T>
where
    T: BulkCollide + Sized,
{
    pub ux: f64,
    pub uy: f64,
    pub collision_type: T,
}

impl<T> Collidable for UpperRegularized<T>
where
    T: BulkCollide + Sized,
{
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        upper_regularized(vel, self.ux, self.uy, &self.collision_type)
    }
}
pub struct RightPressureRegularized<T>
where
    T: BulkCollide + Sized,
{
    pub rho: f64,
    pub u_par: f64,
    pub collision_type: T,
}

impl<T> Collidable for RightPressureRegularized<T>
where
    T: BulkCollide + Sized,
{
    fn collide(&self, vel: D2Q9) -> D2Q9 {
        right_pressure_regularized(vel, self.rho, self.u_par, &self.collision_type)
    }
}
pub fn left_regularized<T>(f_pop: D2Q9, ux: f64, uy: f64, collision_type: &T) -> D2Q9
where
    T: BulkCollide,
{
    let rho: f64 = left_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = left_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

fn left_rho(f_pop: &D2Q9, ux: f64) -> f64 {
    return 1. / (1. - ux)
        * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[3] + f_pop[6] + f_pop[7]));
}

fn left_neq_p(f_neq: &D2Q9) -> P {
    let neq_p_xx = 2. * (f_neq[3] + f_neq[6] + f_neq[7]);
    let neq_p_yy = f_neq[2] + f_neq[4] + 2. * (f_neq[6] + f_neq[7]);
    let neq_p_xy = 2. * (f_neq[7] - f_neq[6]);
    P {
        xx: neq_p_xx,
        yy: neq_p_yy,
        xy: neq_p_xy,
    }
}
pub fn lower_regularized<T>(f_pop: D2Q9, ux: f64, uy: f64, collision_type: &T) -> D2Q9
where
    T: BulkCollide,
{
    let rho: f64 = lower_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = lower_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

fn lower_rho(f_pop: &D2Q9, uy: f64) -> f64 {
    return 1. / (1. - uy)
        * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[8] + f_pop[4] + f_pop[7]));
}

fn lower_neq_p(f_neq: &D2Q9) -> P {
    let neq_p_xx = f_neq[1] + f_neq[3] + 2. * (f_neq[7] + f_neq[8]);
    let neq_p_yy = 2. * (f_neq[4] + f_neq[7] + f_neq[8]);
    let neq_p_xy = 2. * (f_neq[7] - f_neq[8]);
    P {
        xx: neq_p_xx,
        yy: neq_p_yy,
        xy: neq_p_xy,
    }
}

pub fn upper_regularized<T>(f_pop: D2Q9, ux: f64, uy: f64, collision_type: &T) -> D2Q9
where
    T: BulkCollide,
{
    let rho: f64 = upper_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = upper_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

fn upper_rho(f_pop: &D2Q9, uy: f64) -> f64 {
    return 1. / (1. + uy)
        * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[6] + f_pop[2] + f_pop[5]));
}

fn upper_neq_p(f_neq: &D2Q9) -> P {
    let neq_p_xx = f_neq[1] + f_neq[3] + 2. * (f_neq[5] + f_neq[6]);
    let neq_p_yy = 2. * (f_neq[6] + f_neq[2] + f_neq[5]);
    let neq_p_xy = 2. * (f_neq[5] - f_neq[6]);
    P {
        xx: neq_p_xx,
        yy: neq_p_yy,
        xy: neq_p_xy,
    }
}

pub fn right_pressure_regularized<T>(f_pop: D2Q9, rho: f64, u_par: f64, collision_type: &T) -> D2Q9
where
    T: BulkCollide,
{
    let ux = right_u(&f_pop, rho);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, u_par);

    let p_neq = right_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

fn right_u(f_pop: &D2Q9, rho: f64) -> f64 {
    return -1.
        + (1. / rho) * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[1] + f_pop[5] + f_pop[8]));
}

fn regularize_f(mut f: D2Q9, f_eq: D2Q9, p_neq: P) -> D2Q9 {
    for i_pop in 0..9 {
        f[i_pop] = f_eq[i_pop]
            + 9. / 2.
                * T[i_pop]
                * (((C[i_pop][0] * C[i_pop][0]) as f64 - 1. / 3.) * p_neq.xx
                    + ((C[i_pop][1] * C[i_pop][1]) as f64 - 1. / 3.) * p_neq.yy
                    + 2. * C[i_pop][0] as f64 * C[i_pop][1] as f64 * p_neq.xy);
    }
    f
}

fn split_eq_neq(f: &D2Q9, rho: f64, ux: f64, uy: f64) -> (D2Q9, D2Q9) {
    let mut f_eq = [0.; 9];
    let mut f_neq = [0.; 9];

    let u_sqr = ux * ux + uy * uy;

    for i_pop in 0..9 {
        f_eq[i_pop] = compute_equilibrium(i_pop, rho, ux, uy, u_sqr);
        f_neq[i_pop] = f[i_pop] - f_eq[i_pop]
    }
    (f_eq, f_neq)
}

fn right_neq_p(f_neq: &D2Q9) -> P {
    let neq_p_xx = 2. * (f_neq[1] + f_neq[5] + f_neq[8]);
    let neq_p_yy = f_neq[2] + f_neq[4] + 2. * (f_neq[5] + f_neq[8]);
    let neq_p_xy = 2. * (f_neq[5] - f_neq[8]);
    P {
        xx: neq_p_xx,
        yy: neq_p_yy,
        xy: neq_p_xy,
    }
}
struct P {
    xx: f64,
    yy: f64,
    xy: f64,
}

pub fn bounce_back(mut f_pop: D2Q9) -> D2Q9 {
    let mut f_tmp: D2Q9 = [0.; 9];
    for i_pop in 0..9 {
        f_tmp[i_pop] = f_pop[OPPOSITE_OF[i_pop]];
    }
    for i_pop in 0..9 {
        f_pop[i_pop] = f_tmp[i_pop];
    }
    f_pop
}

pub fn compute_macros(f: &D2Q9) -> (f64, f64, f64) {
    let upper_line = f[2] + f[5] + f[6];
    let medium_line = f[0] + f[1] + f[3];
    let lower_line = f[4] + f[7] + f[8];
    let rho = upper_line + medium_line + lower_line;
    let ux = (f[1] + f[5] + f[8] - (f[3] + f[6] + f[7])) / rho;
    let uy = (upper_line - lower_line) / rho;

    (rho, ux, uy)
}

pub fn compute_equilibrium(i_pop: usize, rho: f64, ux: f64, uy: f64, u_sqr: f64) -> f64 {
    let c_u = C_F64[i_pop][0] * ux + C_F64[i_pop][1] * uy;
    rho * T[i_pop] * (1. + 3. * c_u + 4.5 * c_u * c_u - 1.5 * u_sqr)
}
pub trait Node {
    fn collide(&mut self) -> ();
}

pub trait Simulation {
    fn collide(self) -> Self;
    fn stream(self) -> Self;
}

pub type D2Q9 = [f64; 9];

pub struct D2Q9Node<T>
where
    T: Collidable,
{
    pub velocities: D2Q9,
    pub node_type: T,
}

impl<T> Node for D2Q9Node<T>
where
    T: Collidable,
{
    fn collide(&mut self) -> () {
        self.velocities = self.node_type.collide(self.velocities);
    }
}

pub struct RectangleSimulation<N>
where
    N: Node + Sized,
{
    pub l_x: usize,
    pub l_y: usize,
    pub lattice: Array2<N>,
    tmp_lattice: Array2<N>,
}

impl<N: Node> RectangleSimulation<N> {
    pub fn from_function<F>(l_x: usize, l_y: usize, func: F) -> Self
    where
        F: FnMut((usize, usize)) -> N + Copy,
    {
        let shape = (l_x + 2, l_y + 2);
        let lattice = Array2::from_shape_fn(shape, func.clone());
        let tmp_lattice = Array2::from_shape_fn(shape, func);

        RectangleSimulation {
            lattice,
            tmp_lattice,
            l_x,
            l_y,
        }
    }
}

impl<NodeType: Collidable> Simulation for RectangleSimulation<D2Q9Node<NodeType>> {
    fn collide(mut self) -> Self {
        self.lattice.map_inplace(|node| node.collide());
        self
    }
    fn stream(self) -> Self {
        self.propagate_pull()
    }
}

macro_rules! inner_slice {
    ( $x:expr ) => {
        s!(1..=$x.l_x, 1..=$x.l_y)
    };
}

impl<T> RectangleSimulation<D2Q9Node<T>>
where
    T: Collidable,
{
    fn propagate_pull(self) -> Self {
        self.make_periodic_pull().propagate_pull_inner()
    }

    fn propagate_pull_inner(mut self) -> Self {
        let source = self.lattice.windows((3, 3));
        let mut target = self.tmp_lattice.slice_mut(inner_slice!(self));
        Zip::from(&mut target).and(source).for_each(|tar, src| {
            for i_pop in 0..9 {
                let c = C_INDEX[i_pop];
                let vel = src[[c[0], c[1]]].velocities[i_pop];
                tar.velocities[i_pop] = vel;
            }
        });

        std::mem::swap(&mut self.tmp_lattice, &mut self.lattice);

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
pub fn save_vel<T:Collidable>(sim: &RectangleSimulation<D2Q9Node<T>>, f_name: &str) {
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