use crate::lbm::constants::{C, OPPOSITE_OF, T};
use crate::lbm::{compute_equilibrium, CollisionType};

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum BoundaryType {
    BounceBack,
    LowerRegularized {
        ux: f64,
        uy: f64,
        collision_type: CollisionType,
    },
    UpperRegularized {
        ux: f64,
        uy: f64,
        collision_type: CollisionType,
    },
    LeftRegularized {
        ux: f64,
        uy: f64,
        collision_type: CollisionType,
    },
    RightPressureRegularized {
        rho: f64,
        u_par: f64,
        collision_type: CollisionType,
    },
    ZouHeVelocityLeft {
        vx: f64,
        vy: f64,
        collision_type: CollisionType,
    },
    ZouHeVelocityRight {
        vx: f64,
        vy: f64,
        collision_type: CollisionType,
    },
    ZouHeVelocityTop {
        vx: f64,
        vy: f64,
        collision_type: CollisionType,
    },
    ZouHeVelocityBot {
        vx: f64,
        vy: f64,
        collision_type: CollisionType,
    },
    ZouHePressureLeft {
        rho: f64,
        u_par: f64,
        collision_type: CollisionType,
    },
    ZouHePressureRight {
        rho: f64,
        u_par: f64,
        collision_type: CollisionType,
    },
    ZouHePressureTop {
        rho: f64,
        u_par: f64,
        collision_type: CollisionType,
    },
    ZouHePressureBot {
        rho: f64,
        u_par: f64,
        collision_type: CollisionType,
    },
}

impl BoundaryType {
    pub fn collide(&self, velocities: [f64; 9]) -> [f64; 9] {
        match self {
            BoundaryType::BounceBack => bounce_back(velocities),
            BoundaryType::LeftRegularized {
                ux,
                uy,
                collision_type,
            } => left_regularized(velocities, *ux, *uy, collision_type),
            BoundaryType::LowerRegularized {
                ux,
                uy,
                collision_type,
            } => lower_regularized(velocities, *ux, *uy, collision_type),
            BoundaryType::UpperRegularized {
                ux,
                uy,
                collision_type,
            } => upper_regularized(velocities, *ux, *uy, collision_type),
            BoundaryType::RightPressureRegularized {
                rho,
                u_par,
                collision_type,
            } => right_pressure_regularized(velocities, *rho, *u_par, collision_type),
            _ => unimplemented!(),
        }
    }
}

struct P {
    xx: f64,
    yy: f64,
    xy: f64
}
pub fn bounce_back(mut f_pop: [f64; 9]) -> [f64; 9] {
    let mut f_tmp: [f64; 9] = [0.; 9];
    for i_pop in 0..9 {
        f_tmp[i_pop] = f_pop[OPPOSITE_OF[i_pop]];
    }
    for i_pop in 0..9 {
        f_pop[i_pop] = f_tmp[i_pop];
    }
    f_pop
}

pub fn left_regularized(f_pop: [f64; 9], ux: f64, uy: f64, collision_type: &CollisionType) ->  [f64; 9] {
    let rho: f64 = left_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = left_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

pub fn lower_regularized( f_pop: [f64; 9], ux: f64, uy: f64, collision_type: &CollisionType) ->  [f64; 9] {
    let rho: f64 = lower_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = lower_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

pub fn upper_regularized(f_pop: [f64; 9], ux: f64, uy: f64, collision_type: &CollisionType) ->  [f64; 9] {
    let rho: f64 = upper_rho(&f_pop, ux);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, uy);
    let p_neq = upper_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

#[no_mangle]
pub fn right_pressure_regularized(f_pop: [f64; 9], rho: f64, u_par: f64, collision_type: &CollisionType, ) -> [f64; 9] {
    let ux = right_u(&f_pop, rho);
    let (f_eq, f_neq) = split_eq_neq(&f_pop, rho, ux, u_par);

    let p_neq = right_neq_p(&f_neq);
    let f_reg = regularize_f(f_pop, f_eq, p_neq);
    collision_type.collide(f_reg)
}

// compute non-equilibrium stress tensors on left boundary
#[inline]
fn left_neq_p(f_neq: &[f64; 9]) -> P {
    let neq_p_xx = 2. * (f_neq[3] + f_neq[6] + f_neq[7]);
    let neq_p_yy = f_neq[2] + f_neq[4] + 2. * (f_neq[6] + f_neq[7]);
    let neq_p_xy = 2. * (f_neq[7] - f_neq[6]);
    P{ xx: neq_p_xx, yy: neq_p_yy, xy: neq_p_xy}
}
// compute non-equilibrium stress tensor on upper boundary
#[inline]
fn upper_neq_p(f_neq: &[f64; 9]) -> P {
    let neq_p_xx = f_neq[1] + f_neq[3] + 2. * (f_neq[5] + f_neq[6]);
    let neq_p_yy = 2. * (f_neq[6] + f_neq[2] + f_neq[5]);
    let neq_p_xy = 2. * (f_neq[5] - f_neq[6]);
    P{ xx: neq_p_xx, yy: neq_p_yy, xy: neq_p_xy}
}

// compute non-equilibrium stress tensor on lower boundary
#[inline]
fn lower_neq_p(f_neq: &[f64; 9]) -> P {
    let neq_p_xx = f_neq[1] + f_neq[3] + 2. * (f_neq[7] + f_neq[8]);
    let neq_p_yy = 2. * (f_neq[4] + f_neq[7] + f_neq[8]);
    let neq_p_xy = 2. * (f_neq[7] - f_neq[8]);
    P{ xx: neq_p_xx, yy: neq_p_yy, xy: neq_p_xy}
}

// compute non-equilibrium stress tensor on right boundary
#[inline]
fn right_neq_p(f_neq: &[f64; 9]) -> P {
    let neq_p_xx = 2. * (f_neq[1] + f_neq[5] + f_neq[8]);
    let neq_p_yy = f_neq[2] + f_neq[4] + 2. * (f_neq[5] + f_neq[8]);
    let neq_p_xy = 2. * (f_neq[5] - f_neq[8]);
    P{ xx: neq_p_xx, yy: neq_p_yy, xy: neq_p_xy}
}

#[inline]
fn split_eq_neq(f: &[f64; 9], rho: f64, ux: f64, uy: f64) -> ([f64; 9], [f64; 9]) {
    let mut f_eq = [0.; 9];
    let mut f_neq = [0.; 9];

    let u_sqr = ux * ux + uy * uy;

    for i_pop in 0..9 {
        f_eq[i_pop] = compute_equilibrium(i_pop, rho, ux, uy, u_sqr);
        f_neq[i_pop] = f[i_pop] - f_eq[i_pop]
    }
    (f_eq, f_neq)
}

#[inline]
fn regularize_f(mut f: [f64; 9], f_eq: [f64; 9], p_neq: P) -> [f64; 9] {
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

/* Compute density on wall from bulk information on upper boundary. */
#[inline]
fn upper_rho(f_pop: &[f64; 9], uy: f64) -> f64 {
    return 1. / (1. + uy)
        * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[6] + f_pop[2] + f_pop[5]));
}
/* Compute uy on wall from bulk information on upper boundary. */
#[inline]
#[allow(dead_code)]
fn upper_u(f_pop: &[f64; 9], rho: f64) -> f64 {
    return -1.
        + 1.0 / rho * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[6] + f_pop[2] + f_pop[5]));
}
/* Compute density on wall from bulk information on lower boundary. */
#[inline]
fn lower_rho(f_pop: &[f64; 9], uy: f64) -> f64 {
    return 1. / (1. - uy)
        * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[8] + f_pop[4] + f_pop[7]));
}
/* Compute uy on wall from bulk information on lower boundary. */
#[inline]
#[allow(dead_code)]
fn lower_u(f_pop: &[f64; 9], rho: f64) -> f64 {
    return 1.
        - 1. / rho * (f_pop[0] + f_pop[3] + f_pop[1] + 2. * (f_pop[8] + f_pop[4] + f_pop[7]));
}
/* Compute density on wall from bulk information on right boundary. */
#[inline]
#[allow(dead_code)]
fn right_rho(f_pop: &[f64; 9], ux: f64) -> f64 {
    return 1. / (1. + ux)
        * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[1] + f_pop[5] + f_pop[8]));
}
/* Compute ux on wall from bulk information on right boundary. */
#[inline]
fn right_u(f_pop: &[f64; 9], rho: f64) -> f64 {
    return -1. + (1. / rho) * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[1] + f_pop[5] + f_pop[8]));
}
/* Compute density on wall from bulk information on left boundary. */
#[inline]
fn left_rho(f_pop: &[f64; 9], ux: f64) -> f64 {
    return 1. / (1. - ux)
        * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[3] + f_pop[6] + f_pop[7]));
}
/* Compute ux on wall from bulk information on left boundary. */
#[inline]
#[allow(dead_code)]
fn left_u(f_pop: &[f64; 9], rho: f64) -> f64 {
    return 1.
        - 1. / rho * (f_pop[0] + f_pop[2] + f_pop[4] + 2. * (f_pop[3] + f_pop[6] + f_pop[7]));
}

/* Zou/He helper functions and boundary implmenetations          */
/* ****************************************************************/
#[inline]
#[allow(dead_code)]
fn complete_upper(mut f_pop: [f64; 9], ux: f64, uy: f64, rho: f64) -> [f64; 9] {
    f_pop[7] = f_pop[5] + 0.5 * (f_pop[1] - f_pop[3]) - rho * uy / 6. - rho * ux / 2.;
    f_pop[8] = f_pop[6] + 0.5 * (f_pop[3] - f_pop[1]) - rho * uy / 6. + rho * ux / 2.;
    f_pop[4] = f_pop[2] - 2. / 3. * rho * uy;
    f_pop
}
#[inline]
#[allow(dead_code)]
fn complete_lower(mut f_pop: [f64; 9], ux: f64, uy: f64, rho: f64) -> [f64; 9] {
    f_pop[6] = f_pop[8] + 0.5 * (f_pop[1] - f_pop[3]) + rho * uy / 6. - rho * ux / 2.;
    f_pop[5] = f_pop[7] + 0.5 * (f_pop[3] - f_pop[1]) + rho * uy / 6. + rho * ux / 2.;
    f_pop[2] = f_pop[4] + 2. / 3. * rho * uy;
    f_pop
}
#[inline]
#[allow(dead_code)]
fn complete_right(mut f_pop: [f64; 9], ux: f64, uy: f64, rho: f64) -> [f64; 9] {
    f_pop[6] = f_pop[8] + 0.5 * (f_pop[4] - f_pop[2]) - rho * ux / 6. + rho * uy / 2.;
    f_pop[7] = f_pop[5] + 0.5 * (f_pop[2] - f_pop[4]) - rho * ux / 6. - rho * uy / 2.;
    f_pop[3] = f_pop[1] - 2. / 3. * rho * ux;
    f_pop
}
#[inline]
#[allow(dead_code)]
fn complete_left(mut f_pop: [f64; 9], ux: f64, uy: f64, rho: f64) -> [f64; 9] {
    f_pop[5] = f_pop[7] + 0.5 * (f_pop[4] - f_pop[2]) + rho * ux / 6. + rho * uy / 2.;
    f_pop[8] = f_pop[6] + 0.5 * (f_pop[2] - f_pop[4]) + rho * ux / 6. - rho * uy / 2.;
    f_pop[1] = f_pop[3] + 2. / 3. * rho * ux;
    f_pop
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bounce_back() {
        let mut sample = [0., 1., 2., 3., 4., 5., 6., 7., 8.];

        sample = bounce_back( sample);

        assert_eq!(sample, [0., 3., 4., 1., 2., 7., 8., 5., 6.]);
    }
}