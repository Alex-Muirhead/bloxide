/*!
    Core routines for bloxide: A compressible boundary layer analysis code written in Rust.

    @author: Nick Gibbons
*/

#![allow(non_snake_case)]
#![allow(unused_variables)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::{Add, Mul};
use std::path::Path;

pub mod config;
pub mod parameters;
#[cfg(feature = "python")]
pub mod python;
pub mod state;
pub mod viscosity;
use num_dual::{DualSVec64, first_derivative};

use crate::parameters::Parameters;
use crate::state::{Abs, Number, State};
use crate::viscosity::sutherland_mu;

pub fn rkf45_step<S>(f: impl Fn(f64, S) -> S, t: f64, y: S, h: f64) -> (f64, S, S)
where
    S: Copy + Add<Output = S> + Mul<f64, Output = S> + Abs,
{
    // Build up the sample point information as per the text book descriptions.
    let k1 = f(t, y);
    let k2 = f(t + h / 4.0, y + k1 * 0.25 * h);
    let k3 = f(
        t + h * (3.0 / 8.0),
        y + k1 * h * (3.0 / 32.0) + k2 * h * (9.0 / 32.0),
    );
    let k4 = f(
        t + h * (12.0 / 13.0),
        y + k1 * h * (1932.0 / 2197.0) + k2 * h * (-7200.0 / 2197.0) + k3 * h * (7296.0 / 2197.0),
    );
    let k5 = f(
        t + h,
        y + k1 * h * (439.0 / 216.0)
            + k2 * h * -8.0
            + k3 * h * (3680.0 / 513.0)
            + k4 * h * (-845.0 / 4104.0),
    );
    let k6 = f(
        t + h / 2.0,
        y + k1 * h * (-8.0 / 27.0)
            + k2 * h * 2.0
            + k3 * h * (-3544.0 / 2565.0)
            + k4 * h * (1859.0 / 4104.0)
            + k5 * h * (-11.0 / 40.0),
    );
    // Now, do the integration as a weighting of the sampled data.
    let y1 = y
        + k1 * h * (16.0 / 135.0)
        + k3 * h * (6656.0 / 12825.0)
        + k4 * h * (28561.0 / 56430.0)
        + k5 * h * (-9.0 / 50.0)
        + k6 * h * (2.0 / 55.0);
    let err = k1 * h * (1.0 / 360.0)
        + k3 * h * (-128.0 / 4275.0)
        + k4 * h * (-2197.0 / 75240.0)
        + k5 * h * (1.0 / 50.0)
        + k6 * h * (2.0 / 55.0);
    (t + h, y1, err.abs())
}

pub fn soft_max<T: Number>(a: T, b: f64) -> T {
    let da = a - b;
    let scale = (a + b) * 0.5;
    let eps = scale * 1e-3 + 1e-3;

    scale + (da * da + eps * eps).sqrt() / 2.0
}

pub fn density_viscosity_product<T: Number>(g: T, pm: &Parameters) -> T {
    /*!
        Ratio of density x viscosity product at a given point in the boundary layer
    */
    let Temp = g * pm.h_e / pm.C_p;
    let softTemp = soft_max(Temp, 60.0);
    let rho = (softTemp * pm.R).inv() * pm.p_e;
    let mu = sutherland_mu(softTemp);
    rho * mu / (pm.rho_e * pm.mu_e)
}

pub fn self_similar_ode<T: Number>(_t: f64, z: State<T>, pm: &Parameters) -> State<T> {
    let State { f, fd, fdd, g, gd, y } = z;

    let (C, dCdg) = first_derivative(|g| density_viscosity_product(g, pm), g);
    let Cd = dCdg * gd; // oops it's dCdeta = dCdg*dgdeta

    let fddd = C.inv() * (-f * fdd - Cd * fdd);
    let gdd =
        C.inv() * pm.Pr * (-gd * (Cd / pm.Pr + f) - C * pm.u_e.powi(2) / pm.h_e * fdd.powi(2));
    let yd = g * f64::sqrt(2.0 * pm.xi) / pm.u_e * pm.h_e / pm.p_e * (pm.gamma - 1.0) / pm.gamma;

    //println!("        Called ODE: g {:#} dCdg {:#} dCdg2 {:#}", g.re, dCdg.re, dCdg2.re);

    State { f: fd, fd: fdd, fdd: fddd, g: gd, gd: gdd, y: yd }
}

const NSTEPS: usize = 500;
pub fn integrate_through_bl<T: Number>(state0: State<T>, pm: &Parameters) -> Vec<State<T>> {
    /*!
        Fixed-step integrator through eta = [0, 5]
    */
    let eta0 = 0.0;
    let eta_final = 5.0;
    let d_eta = (eta_final - eta0) / (NSTEPS as f64);

    let mut zs = Vec::with_capacity(NSTEPS + 1);
    zs.push(state0);

    let mut eta = eta0;
    let mut z = state0;
    for _ in 0..NSTEPS {
        (eta, z, _) = rkf45_step(|t, y| self_similar_ode(t, y, pm), eta, z, d_eta);
        zs.push(z);
    }

    zs
}

pub fn skin_friction<T: Number>(z: State<T>, pm: &Parameters) -> T {
    /*!
        Return tau, using equations 6.71 and 6.59 from Anderson
    */
    let rhomuw_on_rhomue = density_viscosity_product(z.g, pm);
    let fddw = z.fdd;
    let Rex = pm.rho_e * pm.u_e / pm.mu_e * pm.x;
    let cf = rhomuw_on_rhomue * fddw * f64::sqrt(2.0) / f64::sqrt(Rex);

    cf * 0.5 * pm.rho_e * pm.u_e * pm.u_e
}

pub fn heat_transfer<T: Number>(z: State<T>, pm: &Parameters) -> T {
    /*!
        Return q, using equations 6.79 and ??? from Anderson.
    */
    let rhomuw_on_rhomue = density_viscosity_product(z.g, pm);
    let gdw = z.gd;
    let Rex = pm.rho_e * pm.u_e / pm.mu_e * pm.x;

    rhomuw_on_rhomue * gdw * pm.u_e * pm.rho_e / f64::sqrt(2.0) / pm.Pr * pm.h_e / f64::sqrt(Rex)
}

pub fn recovery_enthalpy<T: Number>(z: State<T>, pm: &Parameters) -> f64 {
    /*!
        The enthlpy the gas reaches after being stagnated in the boundary
        layer. Anderson calls this h_aw, and this expression is equation 6.88.
    */
    pm.h_e + 0.5 * pm.u_e * pm.u_e * f64::sqrt(pm.Pr)
}

pub fn heat_transfer_coefficient<T: Number>(z: State<T>, pm: &Parameters) -> T {
    /*!
        Return pieces needed for equation 6.88 from Anderson. This may
        be useful for material response codes that need to model a
        changing wall temperature.
    */
    let qw = heat_transfer(z, pm);
    let hr = recovery_enthalpy(z, pm);

    qw / ((hr - pm.h_wall) * pm.rho_e * pm.u_e)
}

pub fn boundary_layer_size<T: Number>(states: &[State<T>]) -> Option<T> {
    /*!
        Use 99.9% of the freestream velocity to get the BL size.
    */
    states.iter().find(|z| z.fd > 0.999).map(|z| z.y)
}

pub fn reynolds_number(rho: f64, vel: f64, mu: f64, x: f64) -> f64 {
    rho * vel * x / mu
}

pub fn solve_boundary_layer(pm: &Parameters) -> Vec<State<f64>> {
    let mut error = 1e99;
    let tol = 1e-10;
    let mut iterations = 0;
    let mut fdd = 0.5;
    let mut gd = 1.0;

    while error > tol {
        let dual_wall_state = {
            // Tag the derivatives we care about
            let mut state: State<DualSVec64<2>> =
                State::wall_state(fdd, gd, pm.h_wall, pm.h_e).cast();
            state.fdd = state.fdd.derivative(0);
            state.gd = state.gd.derivative(1);
            state
        };

        let (state, [partial_fdd, partial_gd]) = integrate_through_bl(dual_wall_state, pm)
            .last()
            .expect("Integration didn't return any values")
            .split();

        let fd_err = state.fd - 1.0; // f1 == fd_err
        let g_err = state.g - 1.0; // f2 == g_err
        error = f64::sqrt(fd_err * fd_err + g_err * g_err);

        // Two equation Newton's Method has a 2x2 jacobian that can be
        // inverted analytically. Do this here to get fdd and gd corrections.
        let determinant = partial_fdd.g * partial_gd.fd - partial_gd.g * partial_fdd.fd;
        let diff_fdd = (fd_err * partial_gd.g - g_err * partial_gd.fd) / determinant;
        let diff_gd = (g_err * partial_fdd.fd - fd_err * partial_fdd.g) / determinant;
        fdd += diff_fdd;
        gd += diff_gd;

        iterations += 1;
        if iterations > 100 {
            panic!("Too many iterations of newton solve");
        }
    }

    println!("Solved boundary layer in {:?} iters", iterations);
    let wall_state_final = State::wall_state(fdd, gd, pm.h_wall, pm.h_e);

    integrate_through_bl(wall_state_final, pm)
}

pub fn solve_adiabatic_boundary_layer(pm: &Parameters) -> Vec<State<f64>> {
    let mut error = 1e99;
    let tol = 1e-10;
    let mut iterations = 0;
    let mut fdd = 0.5;
    let mut g = 1.0;
    let eps = 1e-16;

    while error > tol {
        let dual_wall_state = {
            // Tag the derivatives we care about
            let mut state: State<DualSVec64<2>> = State::adiabatic_wall_state(fdd, g).cast();
            state.fdd = state.fdd.derivative(0);
            state.g = state.g.derivative(1);
            state
        };

        let (state, [partial_fdd, partial_g]) = integrate_through_bl(dual_wall_state, pm)
            .last()
            .expect("Integration didn't return any values")
            .split();

        let fd_err = state.fd - 1.0; // f1 == fd_err
        let g_err = state.g - 1.0; // f2 == g_err
        error = f64::sqrt(fd_err * fd_err + g_err * g_err);

        // Two equation Newton's Method has a 2x2 jacobian that can be
        // inverted analytically. Do this here to get fdd and gd corrections.
        let determinant = partial_fdd.g * partial_g.fd - partial_g.g * partial_fdd.fd;
        let diff_fdd = (fd_err * partial_g.g - g_err * partial_g.fd) / determinant;
        let diff_g = (g_err * partial_fdd.fd - fd_err * partial_fdd.g) / determinant;
        fdd += diff_fdd;
        g += diff_g;

        iterations += 1;
        if iterations > 100 {
            panic!("Too many iterations of newton solve");
        }
    }

    println!("Solved adiabatic boundary layer in {:?} iters", iterations);
    let wall_state_final = State::adiabatic_wall_state(fdd, g);

    integrate_through_bl(wall_state_final, pm)
}

pub fn write_dat_file(states: Vec<State<f64>>, filename: impl AsRef<Path>, pm: &Parameters) {
    let file = File::create(filename).expect("Unable to open for writing");
    let mut buf = BufWriter::new(file);
    buf.write_all(b"# y vel T rho p\n")
        .expect("Unable to write to file");

    for zi in states {
        let h = zi.g * pm.h_e;
        let T = h / pm.C_p;
        let rho = pm.p_e / (pm.R * T);
        let u = zi.fd * pm.u_e;
        let y = zi.y;
        let p = pm.p_e;

        writeln!(
            buf,
            "{:16.16e} {:16.16e} {:16.16e} {:16.16e} {:16.16e}",
            y, u, T, rho, p
        )
        .expect("Unable to write line to file");
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn get_heat_transfer(x: f64) -> f64 {
    1.65 * x
}
