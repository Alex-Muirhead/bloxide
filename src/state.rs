/*!
    Let's use structs with named fields instead of array. I can't
    get the complex numbers to work properly with those.

    @author: Nick Gibbons
*/

use derive_more::{Add, Div, Mul, Neg, Sub};
use num_dual::{Dual64, DualNum, DualSVec64};

pub trait Number: DualNum<f64> + Copy {}
impl<T: DualNum<f64> + Copy> Number for T {}

// TODO: Maybe we don't want copy here???
#[derive(Clone, Copy, Debug)]
#[derive(Add, Sub, Mul, Div, Neg)]
pub struct State<T> {
    pub f: T,
    pub fd: T,
    pub fdd: T,
    pub g: T,
    pub gd: T,
    pub y: T,
}

pub trait Abs {
    fn abs(self) -> Self;
}

impl<T> State<T> {
    fn map<U>(self, mut f: impl FnMut(T) -> U) -> State<U> {
        State {
            f: f(self.f),
            fd: f(self.fd),
            fdd: f(self.fdd),
            g: f(self.g),
            gd: f(self.gd),
            y: f(self.y),
        }
    }

    fn map_ref<U>(self, mut f: impl FnMut(&T) -> U) -> State<U> {
        State {
            f: f(&self.f),
            fd: f(&self.fd),
            fdd: f(&self.fdd),
            g: f(&self.g),
            gd: f(&self.gd),
            y: f(&self.y),
        }
    }
}

impl<T: Number> std::ops::Mul<State<T>> for f64 {
    type Output = State<T>;

    fn mul(self, rhs: State<T>) -> State<T> {
        rhs * self
    }
}

impl<T: Number> Abs for State<T> {
    fn abs(self) -> Self {
        self.map_ref(T::abs)
    }
}

impl State<f64> {
    pub fn cast<T: Number>(self) -> State<T> {
        self.map(T::from)
    }

    pub fn wall_state(fdd: f64, gd: f64, h_wall: f64, h_e: f64) -> Self {
        State {
            f: 0.0,
            fd: 0.0,
            fdd,
            g: (h_wall / h_e),
            gd,
            y: 0.0,
        }
    }

    pub fn adiabatic_wall_state(fdd: f64, g: f64) -> Self {
        State { f: 0.0, fd: 0.0, fdd, g, gd: 0.0, y: 0.0 }
    }
}

impl State<Dual64> {
    pub fn split(self) -> (State<f64>, State<f64>) {
        let value = self.map(|d| d.re);
        let deriv = self.map(|d| d.eps);
        (value, deriv)
    }
}

impl<const N: usize> State<DualSVec64<N>> {
    /// Split into the value and one State per derivative
    pub fn split(self) -> (State<f64>, [State<f64>; N]) {
        let value = self.map(|d| d.re);
        let derivs = std::array::from_fn(|i| {
            // A bit of icky "plumbing" code to get the specific element
            self.map(|d| d.eps.unwrap_generic(nalgebra::Const::<N>, nalgebra::Const)[i])
        });
        (value, derivs)
    }
}
