/*!
    Config file reading and associated data structures.

    @author: Nick Gibbons
*/

use std::fs;
use std::path::Path;

use serde::Deserialize;

#[derive(Debug, thiserror::Error)]
pub enum ConfigError {
    #[error("unable to read config file: {0}")]
    Io(#[from] std::io::Error),
    #[error("failed to parse config file: {0}")]
    Parse(#[from] toml::de::Error),
}

#[derive(Debug, PartialEq, Deserialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, eq, module = "bloxidepy")
)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub R: f64,
    pub gamma: f64,
    pub Pr: f64,
    pub p_e: f64,
    pub u_e: f64,
    pub T_e: f64,
    pub T_wall: f64,
    pub x: f64,
}

impl Config {
    /// Parse a Config from a TOML string. Pure parsing, no I/O.
    pub fn from_toml(s: &str) -> Result<Self, toml::de::Error> {
        toml::from_str(s)
    }

    /// Read and parse a Config from a TOML file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, ConfigError> {
        let buffer = fs::read_to_string(path)?;
        Ok(Self::from_toml(&buffer)?)
    }
}
