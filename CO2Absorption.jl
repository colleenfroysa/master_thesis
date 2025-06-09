include("exp_data.jl")

using DifferentialEquations, ModelingToolkit, MethodOfLines, Plots, DomainSets, OrdinaryDiffEq

# --- PDE Definition ---
@parameters z t
@variables cCO2l(..) cH2Ol(..) cMEAl(..) Tl(..)
@variables cCO2v(..) cH2Ov(..) cMEAv(..) Tv(..)
@variables cCO2b(..)
@variables nCO2t(..) nH2Ot(..) nMEAt(..)
@variables nCO2g(..) nH2Og(..) nMEAg(..)
@variables p_CO2(..) p_H2O(..)
@variables cCO2_v_star(..) cH2O_v_star(..)
@variables E(..) K_vov(..) cN2v(..)

# --- Parameters ---
@parameters Q_l Q_v Q_v_mol u_l u_v A_c P_inlet T_v_inlet yCO2_inlet
@parameters h_L A_w h_ov rho_l rho_v 
@parameters ΔH_abs ΔH_vap K_v K_l C R C_L C_V D_l D_v a a_ph ε μ_v ν_v
@parameters A_H2O B_H2O C_H2O A_MEA EA_MEA C1 C2 C3 C4 M_CO2 M_H2O
@parameters θ_deg θ FSE S_p σ cosγ μ_L C_l C_e

# --- Constants ---
const Z = 4.36
const y_N2, y_CO2, y_H2O = 0.70, 0.12, 0.18

# --- Utility Functions ---
function H_CO2_MEA(T, cCO2l_val, cCO2b_val, cMEAl_val)
    α = (cCO2l_val + cCO2b_val) / (cMEAl_val + 2 * cCO2b_val)
    term1 = 496.563 + 34169 * (α / T)
    term2 = 1.69131 * α^2 - 1472.25 / T - (128338 * α) / T^2
    H = term1 * exp(term2)  # mol/m³/Pa
    return H  
end

k2_expr(T) = 4.61 * 1e6 * exp(-4931/T) 

Cp_MEA_expr(T) = 2.5749 + 6.612e-3 * (T - 273.15) - 1.9e-5 * (T - 273.15)^2
Cp_H2O_expr(T) = 4.1908 - 6.62e-4 * (T - 273.15) + 9.14e-6 * (T - 273.15)^2
Cp_l_mass_expr(T) = ((1 - C) * Cp_H2O_expr(T) + C * Cp_MEA_expr(T) + C * (1 - C) * (-0.9198 + 0.013969 * (T - 273.15) + 69.643 * C / (T - 273.15)^1.5859)) *1e3

function Cp_N2(T)
    a1, a2, a3, a4, a5 = 28.98641, 1.853978e-3, -9.647459e-7, 1.471623e-9, -7.682396e-13
    return a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^(-2)
end

function Cp_CO2(T)
    a1, a2, a3, a4, a5 = 24.99735, 55.18696e-3, -33.69137e-6, 7.948387e-9, -0.136638e-11
    return a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^(-2)
end

function Cp_H2O(T)
    a1, a2, a3, a4, a5 = 30.09200, 6.832514e-3, 6.793435e-6, -2.534480e-9, 0.082139e-11
    return a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^(-2)
end

function Cp_v_mol_expr(T)
    return y_N2 * Cp_N2(T) + y_CO2 * Cp_CO2(T) + y_H2O * Cp_H2O(T)
end

# --- Gas and Liquid Feed Setup ---
"""
    gas_feed_conditions(Q_v_m3h, T_v_in_C, P_in_kPa, yCO2_vol)

Compute and return a NamedTuple with
  - Q_v       :: m³/s
  - Q_v_mol   :: mol/s (dry‐basis)
  - T_v_inlet :: K
  - P_inlet   :: Pa
  - cH2Ov0    :: mol/m³ (inlet vapor H₂O)
  - cCO2v0    :: mol/m³ (inlet vapor CO₂)
  - cN2v0     :: mol/m³ (inlet vapor N₂)
  - A_c       :: m² (column cross‐section for u_v_target)
  - u_v       :: m/s (vapor superficial velocity)
"""
function gas_feed_conditions_castor(
    Q_v_m3h::Float64,
    T_v_in_C::Float64,
    P_in_kPa::Float64,
    yCO2_vol::Float64;
    RH::Float64 = 0.80,
    D_col::Float64 = 0.15,
    A_H2O::Float64 = 16.3872, B_H2O::Float64 = 3885.70, C_H2O::Float64 = 230.170,
    R::Float64 = 8.314, μ_v::Float64 = 1.93e-5
)
     # Convert to SI
     Q_v       = Q_v_m3h / 3600.0        # m³/s
     T_v_inlet = T_v_in_C + 273.15       # K
     P_inlet   = P_in_kPa * 1e3          # Pa

     # Column area
     A_c   = π*(D_col/2)^2
     u_v = Q_v / A_c

     # Antoine gives sat‐pressure in kPa
     pH2O_sat_kPa = exp(A_H2O - B_H2O/(T_v_in_C + C_H2O))

     # Convert kPa to Pa for ideal‐gas
     pH2O_sat = pH2O_sat_kPa * 1e3      # Pa

     # pH2O_in  = pH2O_sat
     pH2O_in  = RH*pH2O_sat
     cH2Ov_init   = pH2O_in/(R*T_v_inlet)

     # Dry‐gas molar flow
     c_tot_in = (P_inlet - pH2O_in)/(R*T_v_inlet)
     Q_v_mol  = c_tot_in*Q_v

     # Species inlet concentrations
     cCO2v_init = yCO2_vol*c_tot_in
     cN2v_init  = (1.0-yCO2_vol)*c_tot_in

     M_flue_gas_combustion = 27.7e-3 # kg/mol

    # Density of vapour 
    rho_v = (P_inlet * M_flue_gas_combustion) / (R * T_v_inlet)

    return (
        Q_v       = Q_v,
        Q_v_mol   = Q_v_mol,
        T_v_inlet = T_v_inlet,
        P_inlet   = P_inlet,
        cH2Ov_init= cH2Ov_init,
        cCO2v_init= cCO2v_init,
        cN2v_init = cN2v_init,
        A_c       = A_c,
        u_v       = u_v, 
        rho_v     = rho_v,
        μ_v       = μ_v)
end


"""
    liquid_feed_conditions(loading_CO2, cMEA_kmol_m3, Q_l_Lpm, T_l_in_C; rho_l=1060.0)

Compute and return a NamedTuple with
  - Q_l        :: m³/s
  - T_l_inlet  :: K
  - cMEAl0     :: mol/m³
  - cCO2l0     :: mol/m³ (free CO₂—inlet)
  - cH2Ol0     :: mol/m³ (inlet H₂O)
"""
function liquid_feed_conditions_castor(
    loading_CO2::Float64,
    Q_l_Lpm::Float64,
    T_l_in_C::Float64,
    rho_l::Float64;
    wt_MEA::Float64 = 0.30,
    M_MEA::Float64 = 61.08e-3,     # kg/mol
    M_CO2::Float64 = 44.01e-3,     # kg/mol
    M_H2O::Float64 = 18.015e-3     # kg/mol
)
    # Convert to SI
    Q_l       = Q_l_Lpm / 1000 / 60           # m³/s
    T_l_inlet = T_l_in_C + 273.15             # K

    # Calculate MEA concentration from wt%
    m_MEA = wt_MEA * rho_l                    # kg MEA / m³
    cMEAl = m_MEA / M_MEA                     # mol MEA / m³
    cMEAl_mol_L = cMEAl / 1000                # mol/L

    # CO₂ concentration (mol/m³)
    cCO2l = loading_CO2 * cMEAl

    # Mass of CO₂ and MEA per m³
    m_CO2 = cCO2l * M_CO2                     # kg/m³
    m_H2O = rho_l - m_MEA - m_CO2             # kg/m³
    cH2Ol = m_H2O / M_H2O                     # mol/m³

    return (
        Q_l        = Q_l,
        T_l_inlet  = T_l_inlet,
        cMEAl_init = cMEAl,
        cCO2l_init = cCO2l,
        cH2Ol_init = cH2Ol,
        cMEAl_mol_L = cMEAl_mol_L
    )
end

Q_v_m3h = Q_v_m3h_run1
T_v_in_C = T_v_in_C_run1
P_in_kPa = P_in_kPa_run1
yCO2_molar_dry = yCO2_molar_dry_run1

gas = gas_feed_conditions_castor(Q_v_m3h, T_v_in_C, P_in_kPa, yCO2_molar_dry)

Q_l_Lpm = Q_l_Lpm_run1
T_l_in_C = T_l_in_C_run1
loading_CO2 = loading_CO2_run1
rho_l_val = rho_l_val_run1

liq = liquid_feed_conditions_castor(loading_CO2, Q_l_Lpm, T_l_in_C, rho_l_val)

# --- PDE Equations ---
Dt = Differential(t)
Dz = Differential(z)

# Liquid phase (downward flow)
eq1 = Dt(cCO2l(z,t)) ~ (u_l/h_L)*Dz(cCO2l(z,t)) + (1/h_L)*nCO2t(z,t) + (1/h_L)*nCO2g(z,t)
eq2 = Dt(cH2Ol(z,t)) ~ (u_l/h_L)*Dz(cH2Ol(z,t)) + (1/h_L)*nH2Ot(z,t) + (1/h_L)*nH2Og(z,t)
eq3 = Dt(cMEAl(z,t)) ~ (u_l/h_L)*Dz(cMEAl(z,t)) + (1/h_L)*nMEAt(z,t) + (1/h_L)*nMEAg(z,t)
eq4 = Dt(cCO2b(z,t)) ~ (u_l/h_L)*Dz(cCO2b(z,t)) - (1/h_L)nCO2g(z,t)

# Vapor phase (upward flow)
eq5 = Dt(cCO2v(z,t)) ~ - (u_v/(1-h_L))*Dz(cCO2v(z,t)) - (1/(1-h_L))*nCO2t(z,t)
eq6 = Dt(cH2Ov(z,t)) ~ - (u_v/(1-h_L))*Dz(cH2Ov(z,t)) - (1/(1-h_L))*nH2Ot(z,t)
eq7 = Dt(cMEAv(z,t)) ~ - (u_v/(1-h_L))*Dz(cMEAv(z,t)) - (1/(1-h_L))*nMEAt(z,t)
eq8 = Dt(cN2v(z,t)) ~ - (u_v/(1-h_L))*Dz(cN2v(z,t))

# Energy balances 
eq9 = Dt(Tl(z,t)) ~
    (u_l/h_L)*Dz(Tl(z,t)) +
    (h_ov*A_w/(h_L*rho_l*Cp_l_mass_expr(Tl(z,t))))*(Tv(z,t)-Tl(z,t)) +
    (nCO2t(z,t) * M_CO2 /(h_L*rho_l*Cp_l_mass_expr(Tl(z,t))))*(-ΔH_abs) +
    (nH2Ot(z,t) * M_H2O /(h_L*rho_l*Cp_l_mass_expr(Tl(z,t))))*(-ΔH_vap) 

eq10 = Dt(Tv(z,t)) ~ - (u_v/(1-h_L))*Dz(Tv(z,t)) - (h_ov*A_w)/ ((1-h_L) *
        (cCO2v(z,t) + cH2Ov(z,t) + cN2v(z,t))* Cp_v_mol_expr(Tv(z,t))) * (Tv(z,t) - Tl(z,t))

# Algebraic equations
eq11 = E(z,t) ~ max(sqrt(k2_expr(Tl(z,t)) * cMEAl(z,t) * D_l) / K_l, 1e-12)
eq12 = K_vov(z,t) ~  1 / (1/K_v + 1 / (H_CO2_MEA(Tl(z,t), cCO2l(z,t), cCO2b(z,t), cMEAl(z,t))*K_l*E(z,t)))
eq13 = nCO2t(z,t) ~ K_vov(z,t)*a*(cCO2v(z,t) - cCO2_v_star(z,t)) 
eq14 = nH2Ot(z,t) ~ K_v*a*(cH2Ov(z,t) - cH2O_v_star(z,t))
eq15 = nMEAt(z,t) ~ 0
eq16 = nCO2g(z,t) ~ - k2_expr(Tl(z,t)) * cCO2l(z,t) * cMEAl(z,t)
eq17 = nMEAg(z,t) ~ - 2*k2_expr(Tl(z,t)) * cCO2l(z,t) * cMEAl(z,t)
eq18 = nH2Og(z,t) ~ 0
eq19 = p_CO2(z,t) ~ cCO2l(z,t) / H_CO2_MEA(Tl(z,t), cCO2l(z,t), cCO2b(z,t), cMEAl(z,t))
eq20 = p_H2O(z,t) ~ exp(A_H2O - B_H2O/((Tv(z,t) - 273.15) + C_H2O)) * 1e3
eq21 = cCO2_v_star(z,t) ~ p_CO2(z,t) / (R * Tl(z,t))
eq22 = cH2O_v_star(z,t) ~ p_H2O(z,t) / (R*Tv(z,t))

eqs = [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq15,eq16,eq17,eq18,eq19,eq20,eq21,eq22]

# --- Parameter Assignments ---
Q_v_val       = gas.Q_v
Q_v_mol_val   = gas.Q_v_mol
P_inlet_val   = gas.P_inlet
T_v_inlet_val = gas.T_v_inlet
A_c_val       = gas.A_c
u_v_val       = gas.u_v
rho_v_val     = gas.rho_v
μ_v_val       = gas.μ_v
cCO2v0_val    = gas.cCO2v_init
cH2Ov0_val    = gas.cH2Ov_init
cN2v0_val     = gas.cN2v_init
yCO2_inlet_val= yCO2_molar_dry

Q_l_val       = liq.Q_l
T_l_inlet_val = liq.T_l_inlet
cMEAl0_val    = liq.cMEAl_init
cCO2l0_val    = liq.cCO2l_init
cH2Ol0_val    = liq.cH2Ol_init

R_val = 8.314
g = 9.18

#── packing geometry for Mellapak-250Y ────────────────────────────
a_ph_val = 250.0           # projected area [m²/m³]
FSE_val = 0.35             # surface enhancement factor
ε_val = 0.95               # void fraction
θ_deg_val = 45.0           # corrugation angle [degrees]
θ_val = θ_deg_val * π/180  # [radians]
a_val = a_ph_val * FSE_val
S_p_val = 0.018
a_pack_val = a_ph_val * FSE_val 

num_params = Dict(
    Q_l       => Q_l_val,
    Q_v       => Q_v_val,
    Q_v_mol   => Q_v_mol_val,
    P_inlet   => P_inlet_val,
    yCO2_inlet=> yCO2_inlet_val,
    T_v_inlet => T_v_inlet_val,
    A_c       => A_c_val,
    u_v       => u_v_val, 
    u_l       => Q_l_val / A_c_val,
    R         => R_val,
    A_w       => 15,
    h_ov      => 430.0,
    rho_l     => rho_l_val,      
    rho_v     => rho_v_val,  
    μ_v       => μ_v_val,       
    ν_v       => μ_v_val / rho_v_val,
    ΔH_vap    => 48000,
    ΔH_abs    => -82000,
    C         => 0.3,
    M_CO2     => 44.01e-3,
    M_H2O     => 18.015e-3,
    a_ph      => a_ph_val,
    FSE       => FSE_val,
    ε         => ε_val,
    a         => a_ph_val * FSE_val,
    S_p       => S_p_val,
    σ         => 55e-3,
    cosγ      => 0.9,
    μ_L       => 1.8426e-3,
    C_l       => 0.614 + 71.35 * S_p_val,
    C_e       => 0.9,
    θ         => θ_val,
    A_H2O     => 16.3872, B_H2O => 3885.70, C_H2O => 230.170,
    C_L       => 1.334, C_V => 0.385,
    D_l       => 2.0398e-9, D_v => 1.933e-5,
)

# Initial guess for h_L
h_L_guess = 0.038 * Z

# effective film velocities
u_ve = num_params[u_v] /( num_params[ε]*(1-h_L_guess)*sin(num_params[θ]) )
u_Le = num_params[u_l] /( num_params[ε]* h_L_guess *sin(num_params[θ]) )

Ft = (29.12 * num_params[u_l]^0.4 * num_params[u_l]^0.8 * num_params[S_p]^0.359)/((1 - 0.93*num_params[cosγ]) * (sin(num_params[θ]))^0.3 * (num_params[ε])^0.2) * (rho_l_val/(num_params[σ]*g))^0.15

a_eff = num_params[a_ph] * num_params[FSE] * Ft

ΔP_Z_dry = (0.177 * rho_l_val) / (num_params[S_p] * num_params[ε]^2 * (sin(num_params[θ]))^2) * u_v_val^2 + (88.774 * num_params[μ_v]) / (num_params[S_p]^2 * num_params[ε] * sin(num_params[θ])) * u_v_val
ΔP_Z_flood = 1025 # Pa/m

ΔP_Z = min(ΔP_Z_dry * (1 / (1 - num_params[C_l] * h_L_guess))^5, 0.9 * ΔP_Z_flood)

g_eff = g * ( ( (rho_l_val - rho_v_val) / rho_l_val ) * (1 - ΔP_Z/ΔP_Z_flood))

num_params[h_L] = (4 * Ft/num_params[S_p])^(2/3) * ((3*num_params[μ_L]*num_params[u_l])/(rho_l_val*num_params[ε]*g_eff*sin(num_params[θ])))^(1/3)
num_params[K_v] = 0.054 * ((num_params[D_v]/num_params[S_p])*(u_ve + u_Le)*rho_v_val*num_params[S_p]/num_params[μ_v])^0.8 * (num_params[μ_v]/(rho_v_val*num_params[D_v]))^0.33
num_params[K_l] = 2 * sqrt(num_params[D_l] * num_params[C_e] * u_Le / (π*num_params[S_p]))

#--------------------------------------------------------------------------
# Define domains
z_min = 0.0; z_max = Z
t_min = 0.0; t_max = 3000.0
domains = [z ∈ Interval(z_min, z_max), t ∈ Interval(t_min, t_max)]

#--------------------------------------------------------------------------
# Specify initial conditions (at t = 0)
# Assume uniform initial conditions in z
# Liquid-phase initial conditions:
const cMEA_feed = cMEAl0_val   # total MEA in lean solvent, mol/m³

# Liquid‐phase ICs:
cH2Ol0(z) = cH2Ol0_val
cCO2l0(z) = 0.0
cCO2b0(z) = loading_CO2 * cMEA_feed
cMEAl0(z) = cMEA_feed - 2.0 * cCO2b0(z)
Tl0(z)    = T_l_inlet_val

# Vapor-phase initial conditions:
cCO2v0(z) = cCO2v0_val
cH2Ov0(z) = cH2Ov0_val
cMEAv0(z) = 0.0

cN2v0(z) = cN2v0_val
Tv0(z)   = T_v_inlet_val

p_CO2_0(z) =  cCO2l0(z)  / H_CO2_MEA(Tl0(z), cCO2l0(z), cCO2b0(z), cMEAl0(z))

p_H2O_0(z) = exp(num_params[A_H2O] - num_params[B_H2O]/((Tv(z,t) - 273.15) + num_params[C_H2O])) *1e3
cH2O_v_star_0(z) = p_H2O_0(z) / (num_params[R] * Tv0(z))
cCO2_v_star_0(z) = p_CO2_0(z) / (R * Tl0(z))

# Algebraic variable initial conditions:
E0(z) =  sqrt(k2_expr(Tl0(z)) * cMEAl0(z) * D_l) / K_l
K_vov0(z) =  1 / (1/K_v + 1 / (H_CO2_MEA(Tl0(z), cCO2l0(z), cCO2b0(z), cMEAl0(z))*K_l*E0(z)))
nCO2t0(z) = K_vov0(z)*a*(cCO2v0(z) - cCO2_v_star_0(z))
nH2Ot0(z) = K_v*a*(cH2Ov0(z) - cH2O_v_star_0(z))
nMEAt0(z) = 0.0
nCO2g0(z) = - k2_expr(Tl0(z))*cCO2l0(z)*cMEAl0(z)
nMEAg0(z) = - 2*k2_expr(Tl0(z))*cCO2l0(z)*cMEAl0(z)
nH2Og0(z) = 0.0

#--------------------------------------------------------------------------
# Specify boundary conditions
bcs = [
    # Initial conditions at t = 0:
    cCO2l(z,0) ~ cCO2l0(z),
    cH2Ol(z,0) ~ cH2Ol0(z),
    cMEAl(z,0) ~ cMEAl0(z),
    cCO2b(z,0) ~ cCO2b0(z),   
    cN2v(z,0)  ~ cN2v0(z),
    Tl(z,0)    ~ Tl0(z),
    cCO2v(z,0) ~ cCO2v0(z),
    cH2Ov(z,0) ~ cH2Ov0(z),
    cMEAv(z,0) ~ cMEAv0(z),
    Tv(z,0)    ~ Tv0(z),
    nCO2t(z,0) ~ nCO2t0(z),
    nH2Ot(z,0) ~ nH2Ot0(z),
    nMEAt(z,0) ~ nMEAt0(z),
    nCO2g(z,0) ~ nCO2g0(z),
    nMEAg(z,0) ~ nMEAg0(z),
    nH2Og(z,0) ~ nH2Og0(z),
    p_CO2(z,0) ~ p_CO2_0(z),
    p_H2O(z,0) ~ p_H2O_0(z),
    cCO2_v_star(z,0) ~ cCO2_v_star_0(z),
    cH2O_v_star(z,0) ~ cH2O_v_star_0(z),
    E(z,0) ~ E0(z), 
    K_vov(z,0) ~ K_vov0(z), 
    
    # --- Liquid phase BCs ---
    # Liquid enters at the top (z = z_max)
    cCO2l(z_max,t) ~ cCO2l0(z_max),
    cH2Ol(z_max,t) ~ cH2Ol0(z_max),
    cMEAl(z_max,t) ~ cMEAl0(z_max),
    cCO2b(z_max,t) ~ cCO2b0(z_max),  
    Tl(z_max,t)    ~ Tl0(z_max),
    
    # --- Vapor phase BCs ---
    # Vapor enters at the bottom (z = z_min)
    cCO2v(z_min,t) ~ cCO2v0(z_min),
    cH2Ov(z_min,t) ~ cH2Ov0(z_min),
    cMEAv(z_min,t) ~ cMEAv0(z_min),
    cN2v(z_min,t) ~ cN2v0(z_min),
    Tv(z_min,t)    ~ Tv0(z_min)
]

# --- Substitution and PDESystem construction ---
num_eqs = substitute.(eqs,  Ref(num_params))
num_bcs = substitute.(bcs,  Ref(num_params))

# Construct the PDESystem
@named pdesys_num = PDESystem(num_eqs, num_bcs, domains, [z,t],
    [ cCO2l(z,t), cH2Ol(z,t), cMEAl(z,t), Tl(z,t), 
      cCO2b(z,t), cN2v(z,t),
      cCO2v(z,t), cH2Ov(z,t), cMEAv(z,t), Tv(z,t),
      p_CO2(z,t), p_H2O(z,t), cCO2_v_star(z,t), cH2O_v_star(z,t),
      E(z,t), K_vov(z,t), 
      nCO2t(z,t), nH2Ot(z,t), nMEAt(z,t),
      nCO2g(z,t), nMEAg(z,t), nH2Og(z,t)
      ]
)

# --- Discretization ---
N = 64 # number of spatial grid points
discretization = MOLFiniteDifference([z => N], t, approx_order=4, advection_scheme = UpwindScheme())

prob = discretize(pdesys_num, discretization)

# --- Solve ODE ---
sol = solve(prob, FBDF(); reltol=1e-6, abstol=1e-8)

# --- Post-processing ---
# Inlet CO₂ molar flow (dry basis)
CO2_in = Q_v_mol_val * yCO2_molar_dry

# Read exit concentrations at z = 0 (first grid point), t = final
Nz = size(sol[cCO2v(z,t)], 1)         # number of z‐points
cCO2v_out = sol[cCO2v(z,t)][Nz, end]  # gas exit at z = Z
cH2Ov_out = sol[cH2Ov(z,t)][Nz, end]
cCO2b_sol = sol[cCO2b(z,t)]

# Actual outlet CO₂ mole fraction in the gas
cN2v_out = sol[cN2v(z,t)][Nz,end]
yCO2_out = cCO2v_out / (cCO2v_out + cH2Ov_out + cN2v_out)

# Outlet CO₂ molar flow & capture fraction
CO2_out      = Q_v_mol_val * yCO2_out
CO2_captured = CO2_in - CO2_out
capture_eff  = CO2_captured / CO2_in * 100  # in %

# Rich‐solvent loading
MEA_flow     = Q_l_val * cMEAl0_val
rich_loading = loading_CO2 + CO2_captured / (2 * MEA_flow)

# --- Extract outlet temperatures ---
# Number of spatial grid points
Nz = size(sol[Tl(z, t)], 1)

# Liquid outlet at z = 0 → first row in the solution matrix
Tl_outlet_K = sol[Tl(z, t)][1, end]
# Vapor outlet at z = L → last row
Tv_outlet_K = sol[Tv(z, t)][Nz, end]

# Convert to °C
Tl_outlet_C = Tl_outlet_K - 273.15
Tv_outlet_C = Tv_outlet_K - 273.15

# Spatial grid
z_vals    = sol[z]

# Extract final‐time profiles
nCO2t_sol = sol[nCO2t(z, t)]
Tl_sol    = sol[Tl(z, t)]
Tv_sol    = sol[Tv(z, t)]
cCO2l_sol = sol[cCO2l(z, t)]
cCO2v_sol = sol[cCO2v(z, t)]
cCO2b_sol = sol[cCO2b(z, t)]
cMEAl_sol = sol[cMEAl(z, t)]
cMEAv_sol = sol[cMEAv(z, t)]

# Compute loading
total_MEA   = cMEAl_sol .+ 2 .* cCO2b_sol
loading_sol = (cCO2l_sol .+ cCO2b_sol) ./ total_MEA

# Simulated outlet pCO₂
pCO2_out_Pa = yCO2_out * P_inlet_val
pCO2_out_kPa = pCO2_out_Pa / 1000

# --- Experimental results from Tobiesen Run 1 ---
rich_loading_exp = rich_loading_exp_run1
CO2_absorbed_exp = CO2_absorbed_exp_run1
pCO2_exp = pCO2_exp_run1

CO2_absorbed_kgph = CO2_captured * num_params[M_CO2] * 3600  # [mol/s * kg/mol * s/h]

# --- Error Reporting ---
error_loading = (rich_loading - rich_loading_exp) / rich_loading_exp * 100
error_pCO2 = (pCO2_out_kPa - pCO2_exp) /  pCO2_exp * 100
error_absorbed = (CO2_absorbed_kgph - CO2_absorbed_exp) / CO2_absorbed_exp * 100

println("\n--- Comparison to Experiment (Tobiesen 2007, Run 1) ---")
println("Simulated CO₂ absorbed:  $(round(CO2_absorbed_kgph, sigdigits=4)) kg/h")
println("Experimental CO₂ absorbed:  $(CO2_absorbed_exp) kg/h")
println("Relative error:  $(round(error_absorbed, sigdigits=4)) %")

println("Simulated outlet pCO₂:   $(round(pCO2_out_kPa, sigdigits=4)) kPa")
println("Experimental pCO₂: $(round(pCO2_exp, sigdigits=4)) kPa")
println("Relative Error in pCO₂:  $(round(error_pCO2, sigdigits=4)) %")

println("Simulated Rich Loading: $(round(rich_loading, sigdigits=4)) mol/mol")
println("Experimental Rich Loading: $(rich_loading_exp) mol/mol")
println("Relative Error in Loading: $(round(error_loading, sigdigits=4)) %")

println("\n--- Outlet Temperatures ---")
println("Liquid outlet temperature:  $(round(Tl_outlet_K, sigdigits=5)) K  ($(round(Tl_outlet_C, sigdigits=4)) °C)")
println("Vapor  outlet temperature:  $(round(Tv_outlet_K, sigdigits=5)) K  ($(round(Tv_outlet_C, sigdigits=4)) °C)")