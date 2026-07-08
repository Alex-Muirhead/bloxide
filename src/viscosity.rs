/*!
    Code for computing viscosity with the Sutherland expression

    @author: Nick Gibbons
*/

use crate::state::Number;

const MU_REF: f64 = 1.716e-05;
const T_REF: f64 = 273.0;
const S: f64 = 111.0;

// T here is a generic type, not the temperature!
pub fn sutherland_mu<N: Number>(temp: N) -> N {
    (temp / T_REF).sqrt() * (temp / T_REF) * (T_REF + S) / (temp + S) * MU_REF
}

pub fn sutherland_mu_derivative<N: Number>(temp: N) -> N {
    (temp / T_REF).sqrt() * (temp + 3.0 * S) * MU_REF * (T_REF + S)
        / ((temp + S).powi(2) * 2.0 * T_REF)
}
