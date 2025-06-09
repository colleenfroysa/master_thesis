# CO₂ Absorption Column Simulation in Julia

This repository contains a dynamic simulation model of a packed absorption column for capturing carbon dioxide (CO₂) using an aqueous monoethanolamine (MEA) solution. The model is implemented in Julia using `ModelingToolkit.jl` and solves coupled mass and energy balance equations with chemical reaction and mass transfer kinetics. The simulation is validated against experimental data from Tobiesen (2007).

---

## Features

- Rate-based, non-equilibrium (NEQ) simulation of CO₂ absorption
- Two-film theory for gas–liquid mass transfer
- Reaction kinetics for MEA–CO₂ chemistry
- Heat transfer with species- and temperature-dependent properties
- Hydrodynamics of structured packing (Mellapak-250Y)
- Validation against pilot-scale experimental data

---

## Dependencies

This model uses the Julia language and the SciML ecosystem. Install required packages using Julia's package manager:

"""julia
using Pkg
Pkg.add([
    "DifferentialEquations",
    "ModelingToolkit",
    "MethodOfLines",
    "OrdinaryDiffEq",
    "DomainSets",
    "Plots"
])
"""

---

## How to Run

1. Ensure the following project structure:

"""
.
├── src/
│   ├── simulation.jl         # Main simulation script
│   └── exp_data.jl           # Experimental data and input parameters
├── Project.toml              # Julia project environment
├── Manifest.toml             # Package version lock file
├── README.md                 # This file
└── LICENSE                   # License information
"""

2. Make sure `exp_data.jl` defines experimental input for Tobiesen Run 1:
   - `Q_v_m3h_run1`, `T_v_in_C_run1`, `P_in_kPa_run1`, `yCO2_molar_dry_run1`
   - `Q_l_Lpm_run1`, `T_l_in_C_run1`, `loading_CO2_run1`, `rho_l_val_run1`
   - `CO2_absorbed_exp_run1`, `rich_loading_exp_run1`, `pCO2_exp_run1`

2. In `src/exp_data.jl`, you can:
   - Select experimental input from Tobiesen Runs 1 to 20 by changing inputs:
        - `Q_v_m3h_run1`, `T_v_in_C_run1`, `P_in_kPa_run1`, `yCO2_molar_dry_run1`
        - `Q_l_Lpm_run1`, `T_l_in_C_run1`, `loading_CO2_run1`, `rho_l_val_run1`
        - `CO2_absorbed_exp_run1`, `rich_loading_exp_run1`, `pCO2_exp_run1`
   - **Or** define your own data manually

   Ensure that all custom input values are in the correct units:
   - Volumetric flow rates in m³/h for vapour phase and L/min for liquid phase 
   - Pressures in kPa
   - Temperatures in °C
   - Loadings in mol CO₂ / mol MEA

3. Run the simulation from the Julia REPL or a script:

"""julia
include("exp_data.jl")
include("simulation.jl")
"""

4. The script will print:
   - Simulated vs. experimental CO₂ absorbed [kg/h]
   - Outlet partial pressure of CO₂ [kPa]
   - Rich CO₂ loading in the solvent [mol/mol]
   - Outlet liquid and vapor temperatures
   - Relative errors vs. experiment

---

## Output Metrics

- CO₂ capture efficiency [%]
- Rich loading [mol CO₂ / mol MEA]
- CO₂ partial pressure at gas outlet [kPa]
- Temperature profiles of vapor and liquid streams
- Mass transfer fluxes [mol/m²·s]

---

## Reference

> Finn Andrew Tobiesen, Hallvard F. Svendsen, and Olav Juliussen. “Experimental validation of a rigorous absorber model for CO2 postcombustion capture”. In: AIChE Journal 53.4 (2007), pp. 846–865. doi: https://doi. org/10.1002/aic.11133. url: https://aiche.onlinelibrary.wiley.com/doi/abs/10.1002/aic.11133

---

## License

This project is provided for academic and educational purposes. If used in publications or derivative work, please credit the original source appropriately.