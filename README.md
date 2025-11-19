# bloxide -- A compressible boundary layer analysis code

bloxide is a tool for analysing self-similar compressible boundary layers, written in Rust.

## Build Instructions
Install dependancies:

    curl https://sh.rustup.rs -sSf | sh

Clone the repository:

     git clone https://github.com/uqngibbo/bloxide.git

Build the code:

    cd bloxide
    cargo install --path .

## Example Use

    cd examples
    bloxide cold.yaml

This should output some details about the boundary layer as well as a file called cold.dat which has the entire profile saved in it.

    Config {
        R: 287.1,
        gamma: 1.4,
        Pr: 0.71,
        p_e: 2303.0,
        u_e: 604.5,
        T_e: 108.1,
        T_wall: 269.5,
        x: 0.5,
    }
    Solved boundary layer in 5 iters
    Skin Friction:   5.05226 N/m2
    Heat Transfer:   -0.00927 W/cm2
    99.9% BL size:   2.77701 mm
    Reynolds Number: 2.99296 million (x)
    Reynolds Number: 16622.98 (BL thickness)
    
    Freestream cp:     1004.85 J/kg
    Wall rhoe*ue*CH:   1.1271221e-2 kg/m2/s
    Heat Transfer CH:  2.5126940e-4
    Recovery Enthalpy: 2.6257857e-1 MJ/kg
    
    Solved adiabatic boundary layer in 5 iters
    Adiabatic Wall Temp: 260.68950 K
    Done.

## Author:
Nick Gibbons (n.gibbons@uq.edu.au) and Peter Jacobs

## License:
bloxide use is governed by the GNU General Public License 3. This is a copyleft license, which means that any code based on bloxide must also be made freely available under a similar license.
See the file gpl.txt for further details.
