/*!
    bloxide: A compressible boundary layer analysis code.

    References:
     - "Hypersonic and High Temperature Gas Dyanmics", John D. Anderson

    @author: Nick Gibbons and Peter Jacobs
*/

#![allow(non_snake_case)]
#![allow(unused_variables)]

use std::path::PathBuf;

use clap::Parser;

use bloxide::config::Config;
use bloxide::parameters::Parameters;
use bloxide::*;

#[derive(Parser)]
#[command(name = "bloxide")]
#[command(about = "A compressible boundary layer analysis code.", long_about = None)]
struct Cli {
    /// Path to a config file
    config_path: PathBuf,
}

fn main() {
    let cli = Cli::parse();
    let config_file_name = cli.config_path;

    let config = Config::from_file(&config_file_name).unwrap_or_else(|err| {
        panic!(
            "Failed to load config file {:?}: {}",
            &config_file_name, err
        )
    });
    let pm = Parameters::new(&config);
    println!("{:#?}", config);

    let states = solve_boundary_layer(&pm);
    let state_initial = states[0];
    let state_final = states.last().unwrap();

    let tauw = skin_friction(state_initial, &pm);
    let qw = heat_transfer(state_initial, &pm);
    let ybl = boundary_layer_size(&states).expect("Cannot find bl size");
    let Rex = reynolds_number(pm.rho_e, pm.u_e, pm.mu_e, pm.x);
    let Ret = reynolds_number(pm.rho_e, pm.u_e, pm.mu_e, ybl);

    println!("Skin Friction:   {:5.5} N/m2", tauw);
    println!("Heat Transfer:   {:5.5} W/cm2", qw / 1e4);
    println!("99.9% BL size:   {:5.5} mm", ybl * 1000.0);
    println!("Reynolds Number: {:5.5} million (x)", Rex / 1e6);
    println!("Reynolds Number: {:5.2} (BL thickness)\n", Ret);

    let hr = recovery_enthalpy(state_initial, &pm);
    let CH = heat_transfer_coefficient(state_initial, &pm);
    println!("Freestream cp:     {:7.2} J/kg", pm.C_p);
    println!("Wall rhoe*ue*CH:   {:7.7e} kg/m2/s", pm.rho_e * pm.u_e * CH);
    println!("Heat Transfer CH:  {:7.7e}", CH);
    println!("Recovery Enthalpy: {:7.7e} J/kg", hr);
    println!();

    let adiabatic_states = solve_adiabatic_boundary_layer(&pm);
    let adiabatic_state_initial = adiabatic_states[0];
    let hwall = adiabatic_state_initial.g * pm.h_e;
    let Twall = hwall / pm.C_p;
    println!("Adiabatic Wall Temp: {:5.5} K", Twall);

    let filename = config_file_name.with_extension("dat");
    write_dat_file(states, filename, &pm);

    println!("Done.");
}
