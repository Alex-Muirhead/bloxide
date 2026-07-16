/*!
    Python bindings for bloxide.

    All Python-facing adaptation lives here: core functions keep their
    Rust signatures, and thin #[pyfunction] wrappers handle type
    conversion and error mapping for Python callers.
*/

use pyo3::prelude::*;

use crate::config::Config;

#[pymethods]
impl Config {
    #[new]
    fn py_new(
        R: f64,
        gamma: f64,
        Pr: f64,
        p_e: f64,
        u_e: f64,
        T_e: f64,
        T_wall: f64,
        x: f64,
    ) -> Self {
        Config::new(R, gamma, Pr, p_e, u_e, T_e, T_wall, x)
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self)
    }
}

#[pyo3::pymodule(name = "_core")]
mod _core {
    use std::path::PathBuf;

    use pyo3::exceptions::{PyFileNotFoundError, PyValueError};
    use pyo3::prelude::*;

    #[pymodule_export]
    use crate::config::Config;

    #[pyfunction]
    fn get_heat_transfer(x: f64) -> f64 {
        1.65 * x
    }

    #[pyfunction]
    fn read_config_file(path: PathBuf) -> PyResult<Config> {
        if !path.exists() {
            return Err(PyFileNotFoundError::new_err(format!(
                "No such config file: {}",
                path.display()
            )));
        }
        let filename = path
            .to_str()
            .ok_or_else(|| PyValueError::new_err("config path is not valid UTF-8"))?;
        Ok(crate::config::read_config_file(filename))
    }
}
