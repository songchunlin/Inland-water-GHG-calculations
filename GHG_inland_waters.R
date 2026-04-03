# =============================================================================
# GHG Concentrations and Air-Water Fluxes for Inland Waters
# =============================================================================


########################################################################
# GHG_inland_waters.R
#
# Calculates dissolved GHG concentrations (pCO2, pCH4, pN2O) and
# air-water fluxes (FCO2, FCH4, FN2O) from raw GC headspace data for
# inland freshwater systems (streams, rivers, lakes, reservoirs).
#
# Author: Chunlin Song, Sichuan University
# Date:   2026-04-02
# Contact: songchunlin@scu.edu.cn
#
# ── Inputs ────────────────────────────────────────────────────────────
#   ghg_input.csv   Raw field and GC measurements (see column guide below)
#
# ── Outputs ───────────────────────────────────────────────────────────
#   ghg_results.csv  Concentrations, Schmidt numbers, k600, k_gas, Kh,
#                    and air-water fluxes for CO2, CH4, and N2O
#
# ── Code provenance ───────────────────────────────────────────────────
# The Rheadspace_GHG() function in this script is derived from two sources:
#
# (1) Rheadspace.R — complete headspace CO2 method with carbonate equilibrium
#     Authors: R. Marcé (ICRA), J. Kim (UQAM), Y.T. Prairie (UQAM)
#     Date: December 2020; bug fixes July 2023 and April 2025 (R. Marcé)
#     Reference: Koschorreck M., Prairie Y.T., Kim J., Marcé R. (2021).
#       Technical note: CO2 is not like CH4 – limits of the headspace method
#       to analyse pCO2 in water. Biogeosciences, 18, 1619–1627.
#       https://doi.org/10.5194/bg-18-1619-2021
#     License: GNU General Public License v3
#
# (2) headspace_calcs — simple headspace mass-balance for CH4 and N2O
#     Author: Kelly S. Aho
#     Repository: https://github.com/kellyaho/headspace_calcs/
#
# Modifications relative to the above sources:
#   - Freshwater carbonate equilibrium constants (Millero 1979) used throughout
#   - Henry's Law constants for CO2, CH4, and N2O updated to Sander (2015)
#   - CH4 and N2O processing added alongside the CO2 complete headspace routine
#   - Air-water flux calculations added using gas transfer velocity k600
#     (Raymond et al. 2012 energy-dissipation model, two-regime) and Henry's
#     solubility from marelac::gas_solubility()
#
# ── Key references ────────────────────────────────────────────────────
#   Millero F. (1979). The thermodynamics of the carbonate system in seawater.
#     Geochim. Cosmochim. Acta, 43(10), 1651–1661.
#   Sander R. (2015). Compilation of Henry's law constants (version 4.0)
#     for water as solvent. Atmos. Chem. Phys., 15, 4399–4981.
#     https://doi.org/10.5194/acp-15-4399-2015
#   Raymond P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic
#     geometry in streams and small rivers. Limnol. Oceanogr. Fluids Environ.,
#     2, 41–53. https://doi.org/10.1215/21573689-1597669
#   Wanninkhof R. (1992). Relationship between wind speed and gas exchange over
#     the ocean. J. Geophys. Res., 97(C5), 7373–7382.
#
# ── Column guide for ghg_input.csv ───────────────────────────────────
#   site                  Sample site ID
#   datetime.EST          Sample date/time [YYYY-MM-DD HH:MM]
#   HS.mCH4/CO2/N2O.before  Headspace mole fraction before equilibration [ppmv]
#   HS.mCH4/CO2/N2O.after   Headspace mole fraction after equilibration  [ppmv]
#   Temp.insitu           In-situ water temperature [°C]
#   Temp.equil            Water temperature at equilibration [°C]
#   Alkalinity.measured   Total alkalinity [µeq/L]
#   Volume.gas            Headspace gas volume [mL]
#   Volume.water          Sample water volume [mL]
#   Bar.pressure          Barometric pressure [kPa]  (101.325 kPa = 1 atm)
#   Tw                    In-situ water temperature for Schmidt number [°C]
#   v_ms                  Mean flow velocity [m/s]  (streams/rivers only)
#   slope                 Channel slope [m/m]       (streams/rivers only)
#   SAL                   Salinity [PSU]
#   CO2_air_ppm           Atmospheric CO2 [ppm]     (global mean 2024: 423)
#   CH4_air_ppm           Atmospheric CH4 [ppm]     (global mean 2024: 1.93)
#   N2O_air_ppm           Atmospheric N2O [ppm]     (global mean 2024: 0.338)
########################################################################

# ── 0. Libraries ──────────────────────────────────────────────────────────────
if (!requireNamespace("marelac", quietly = TRUE)) install.packages("marelac")
library(marelac)   # gas_solubility()
library(readr)
library(dplyr)

# ── 1. Headspace GHG concentration function ───────────────────────────────────
# Returns 16-column data frame:
#   cols 1-2  : Site, Timestamp
#   cols 3-6  : CO2 complete headspace (mCO2, pCO2, [CO2], pH)
#   cols 7-10 : CO2 simple headspace   (mCO2, pCO2, [CO2], % error vs complete)
#   cols 11-13: CH4 simple headspace   (mCH4, pCH4, [CH4])
#   cols 14-16: N2O simple headspace   (mN2O, pN2O, [N2O])
#
# Units: mGAS = ppmv, pGAS = µatm, [GAS] = µmol/L
#
# Bar.pressure input must be in kPa (e.g., 57.75 kPa = 0.57 atm at ~4500 m)

Rheadspace_GHG <- function(input.table) {
  
  if (!is.data.frame(input.table) || ncol(input.table) < 14)
    stop("Input must be a data frame with at least 14 columns.", call. = FALSE)
  
  Site           <- as.character(input.table$site)
  Timestamp      <- input.table$datetime.EST
  mCH4_headspace <- as.numeric(input.table$HS.mCH4.before)
  mCO2_headspace <- as.numeric(input.table$HS.mCO2.before)
  mN2O_headspace <- as.numeric(input.table$HS.mN2O.before)
  mCH4_eq        <- as.numeric(input.table$HS.mCH4.after)
  mCO2_eq        <- as.numeric(input.table$HS.mCO2.after)
  mN2O_eq        <- as.numeric(input.table$HS.mN2O.after)
  temp_insitu    <- as.numeric(input.table$Temp.insitu)
  temp_eq        <- as.numeric(input.table$Temp.equil)
  alk            <- as.numeric(input.table$Alkalinity.measured)
  vol_gas        <- as.numeric(input.table$Volume.gas)
  vol_water      <- as.numeric(input.table$Volume.water)
  Bar.pressure   <- as.numeric(input.table$Bar.pressure)   # kPa
  
  R_const <- 0.082057338  # L·atm·K⁻¹·mol⁻¹
  
  out <- data.frame(matrix(NA_real_, nrow = length(mCO2_eq), ncol = 16))
  names(out) <- c(
    "Site", "Timestamp",
    "mCO2_complete_ppmv", "pCO2_complete_uatm", "CO2_complete_umolL", "pH",
    "mCO2_simple_ppmv",   "pCO2_simple_uatm",   "CO2_simple_umolL",  "pct_error",
    "mCH4_simple_ppmv",   "pCH4_simple_uatm",   "CH4_simple_umolL",
    "mN2O_simple_ppmv",   "pN2O_simple_uatm",   "N2O_simple_umolL"
  )
  
  for (i in seq_along(mCO2_eq)) {
    
    AT   <- alk[i] * 1e-6           # mol/L
    Teq  <- temp_eq[i] + 273.15     # K
    Tis  <- temp_insitu[i] + 273.15 # K
    BP   <- Bar.pressure[i]         # kPa
    
    # Equilibrium constants (Lueker et al. 2000 / Millero 1995)
    K1  <- 10^-(-126.34048 + 6320.813/Teq + 19.568224*log(Teq))
    K2  <- 10^-(-90.18333  + 5143.692/Teq + 14.613358*log(Teq))
    Kw  <- exp(148.9652 - 13847.26/Teq - 23.6521*log(Teq))
    
    # Henry's constants (Sander 2015) – at equilibration temp and in-situ temp
    Kh_CO2_eq  <- 0.00033  * exp(2400*(1/Teq - 1/298.15)) * 101325/1000
    Kh_CO2_is  <- 0.00033  * exp(2400*(1/Tis - 1/298.15)) * 101325/1000
    Kh_CH4_eq  <- 0.000014 * exp(1900*(1/Teq - 1/298.15)) * 101325/1000
    Kh_CH4_is  <- 0.000014 * exp(1900*(1/Tis - 1/298.15)) * 101325/1000
    Kh_N2O_eq  <- 0.00024  * exp(2700*(1/Teq - 1/298.15)) * 101325/1000
    Kh_N2O_is  <- 0.00024  * exp(2700*(1/Tis - 1/298.15)) * 101325/1000
    
    HS <- vol_gas[i] / vol_water[i]  # headspace ratio
    
    # ── Complete headspace CO2 (carbonate equilibrium) ──────────────────────
    co2_eq_val <- Kh_CO2_eq * mCO2_eq[i] / 1e6
    h_all      <- polyroot(c(-(2*K1*K2*co2_eq_val), -(co2_eq_val*K1 + Kw), AT, 1))
    h          <- Re(h_all)[Re(h_all) > 0]
    
    DIC_eq  <- co2_eq_val * (1 + K1/h + K1*K2/(h^2))
    DIC_ori <- DIC_eq + (mCO2_eq[i] - mCO2_headspace[i]) / 1e6 /
      (R_const * Teq) * HS
    
    h_all_ori <- polyroot(c(
      -(K1*K2*Kw),
      K1*K2*AT - K1*Kw - 2*DIC_ori*K1*K2,
      AT*K1 - Kw + K1*K2 - DIC_ori*K1,
      AT + K1, 1))
    h_ori     <- Re(h_all_ori)[Re(h_all_ori) > 0]
    co2_final <- h_ori * (DIC_ori * h_ori * K1 / (h_ori^2 + K1*h_ori + K1*K2)) / K1
    
    out[i, 1]  <- Site[i]
    out[i, 2]  <- Timestamp[i]
    out[i, 3]  <- co2_final / Kh_CO2_is * 1e6            # mCO2 ppmv
    out[i, 4]  <- co2_final / Kh_CO2_is * 1e6 * BP/101.325  # pCO2 µatm
    out[i, 5]  <- co2_final * 1e6                         # [CO2] µmol/L
    out[i, 6]  <- -log10(h_ori)                           # pH
    
    # ── Simple headspace CO2 (mass balance) ─────────────────────────────────
    CO2_sol  <- mCO2_eq[i] / 1e6 * Kh_CO2_eq
    CO2_smass <- CO2_sol * vol_water[i] / 1000
    CO2_gmass <- mCO2_eq[i]       / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    CO2_gmass0 <- mCO2_headspace[i] / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    CO2_conc <- (CO2_smass + CO2_gmass - CO2_gmass0) / (vol_water[i]/1000)
    
    out[i, 7]  <- CO2_conc / Kh_CO2_is * 1e6
    out[i, 8]  <- CO2_conc / Kh_CO2_is * 1e6 * BP/101.325
    out[i, 9]  <- CO2_conc * 1e6
    out[i, 10] <- (out[i, 7] - out[i, 3]) / out[i, 3] * 100  # % error simple vs complete
    
    # ── Simple headspace CH4 ─────────────────────────────────────────────────
    CH4_sol   <- mCH4_eq[i] / 1e6 * Kh_CH4_eq
    CH4_smass <- CH4_sol * vol_water[i] / 1000
    CH4_gmass <- mCH4_eq[i]       / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    CH4_gmass0 <- mCH4_headspace[i] / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    CH4_conc  <- (CH4_smass + CH4_gmass - CH4_gmass0) / (vol_water[i]/1000)
    
    out[i, 11] <- CH4_conc / Kh_CH4_is * 1e6
    out[i, 12] <- CH4_conc / Kh_CH4_is * 1e6 * BP/101.325
    out[i, 13] <- CH4_conc * 1e6
    
    # ── Simple headspace N2O ─────────────────────────────────────────────────
    N2O_sol   <- mN2O_eq[i] / 1e6 * Kh_N2O_eq
    N2O_smass <- N2O_sol * vol_water[i] / 1000
    N2O_gmass <- mN2O_eq[i]       / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    N2O_gmass0 <- mN2O_headspace[i] / 1e6 * (vol_gas[i]/1000) / (R_const * Teq)
    N2O_conc  <- (N2O_smass + N2O_gmass - N2O_gmass0) / (vol_water[i]/1000)
    
    out[i, 14] <- N2O_conc / Kh_N2O_is * 1e6
    out[i, 15] <- N2O_conc / Kh_N2O_is * 1e6 * BP/101.325
    out[i, 16] <- N2O_conc * 1e6
  }
  
  out$Timestamp <- as.character(out$Timestamp)
  return(out)
}

# ── 2. Load input data ────────────────────────────────────────────────────────
# skip = 2 jumps over the units row (row 2) and notes row (row 3) in the template;
# remove it if your CSV has only a single header row with no descriptor rows below it.
raw <- read_csv("ghg_input.csv", skip = 2, col_names = TRUE, show_col_types = FALSE)

# Re-attach column names from row 1 (read_csv with skip=2 loses them; read separately)
col_names <- names(read_csv("ghg_input.csv", n_max = 0, show_col_types = FALSE))
names(raw) <- col_names

# Remove rows with missing Site
raw <- raw[!is.na(raw$site), ]

# Coerce all measurement columns to numeric (guards against stray text in the CSV)
numeric_cols <- setdiff(names(raw), c("site", "datetime.EST"))
raw[numeric_cols] <- lapply(raw[numeric_cols], as.numeric)

# ── 3. GHG concentrations (Rheadspace_GHG) ───────────────────────────────────
hs_cols <- c("site", "datetime.EST",
             "HS.mCH4.before", "HS.mCO2.before", "HS.mN2O.before",
             "HS.mCH4.after",  "HS.mCO2.after",  "HS.mN2O.after",
             "Temp.insitu", "Temp.equil",
             "Alkalinity.measured", "Volume.gas", "Volume.water", "Bar.pressure")

pGHG <- Rheadspace_GHG(raw[, hs_cols])

# ── 4. Henry's solubility constants (marelac) ─────────────────────────────────
# gas_solubility() returns mmol/m³/bar; Tw = in-situ water temp (°C), SAL = salinity (PSU)
Kh_CO2 <- gas_solubility(S = raw$SAL, t = raw$Tw, species = "CO2")
Kh_CH4 <- gas_solubility(S = raw$SAL, t = raw$Tw, species = "CH4")
Kh_N2O <- gas_solubility(S = raw$SAL, t = raw$Tw, species = "N2O")

# ── 5. Gas transfer velocity ────────
# ── 5a. Schmidt numbers ───────────────────────────────────────────────────────
# Schmidt number (Sc) = kinematic viscosity / molecular diffusivity of the gas.
# Used to scale k600 to gas-specific piston velocity: k_gas = k600 × (Sc/600)^n
# where n = -0.5 for rough surfaces (streams, wind-mixed lakes) or -2/3 for
# smooth surfaces (calm lakes). n = -0.5 is used here.
#
# ACTIVE — Wanninkhof (1992), freshwater, 3rd-order polynomials:
Sc_CO2 <- 1911.1 - 118.11*raw$Tw + 3.4527*raw$Tw^2 - 0.04132*raw$Tw^3
Sc_CH4 <- 1897.8 - 114.28*raw$Tw + 3.2902*raw$Tw^2 - 0.03906*raw$Tw^3
Sc_N2O <- 2055.6 - 137.11*raw$Tw + 4.3173*raw$Tw^2 - 0.05435*raw$Tw^3

# ALTERNATIVE — Wanninkhof (2014), freshwater, 4th-order polynomials (more accurate,
# especially above 20°C; recommended for studies spanning a wide temperature range):
#   Sc_CO2 <- 1923.6 - 125.06*raw$Tw + 4.3773*raw$Tw^2 - 0.085681*raw$Tw^3 + 0.00070284*raw$Tw^4
#   Sc_CH4 <- 1909.4 - 120.78*raw$Tw + 4.1555*raw$Tw^2 - 0.080578*raw$Tw^3 + 0.00065777*raw$Tw^4
#   Sc_N2O <- 2141.2 - 152.56*raw$Tw + 5.8963*raw$Tw^2 - 0.12411*raw$Tw^3  + 0.0010655*raw$Tw^4
#   Reference: Wanninkhof R. (2014). Relationship between wind speed and gas exchange
#   over the ocean revisited. Limnol. Oceanogr. Methods, 12, 351-362.
#   https://doi.org/10.4319/lom.2014.12.351

# ── 5b. k600 — streams and rivers ────────────────────────────────────────────
# ACTIVE — Ulseth et al. (2019) two-regime energy-dissipation model.
# eD [m²/s³] = g × v [m/s] × S [m/m]  (gravitational acceleration × velocity × slope)
# The two-regime breakpoint (eD = 0.02 m²/s³) separates smooth low-gradient rivers
# from turbulent high-gradient streams.  Output is k600 [m/d].
# Reference: Ulseth A.J. et al. (2019). Distinct air-water gas exchange regimes in
#   low- and high-energy streams. Nature Geoscience, 12, 259–263.
#   https://doi.org/10.1038/s41561-019-0324-8
elevition <- 500  # m; ADJUST this to your own elevation
g_accel <- 9.81 * (1 - 2*elevition/637100)          # gravitational acceleration (m/s²)
eD      <- g_accel * raw$v_ms * raw$slope       # energy dissipation rate (m²/s³)

k600 <- ifelse(eD > 0.02,
               exp(6.43 + 1.18 * log(eD)),      # high-energy / turbulent streams
               exp(3.10 + 0.35 * log(eD)))       # low-energy  / smooth rivers

# ALTERNATIVE stream/river k600 models (replace the active block above):
# All equations below assume v = velocity [m/s], d = mean depth [m], S = slope [m/m].
# Conversion: cm/h × 24/100 = m/d
#
# Raymond et al. (2012) — single power-law, does not distinguish energy regimes:
#   k600 <- v_ms*slope*2841 + 2.02  # V*S*2841±107 + 2.02±0.209
#   Reference: Raymond P.A. et al. (2012). Scaling the gas transfer velocity and
#   hydraulic geometry in streams and small rivers. Limnol. Oceanogr. Fluids Environ.,
#   2, 41-53. https://doi.org/10.1215/21573689-1597669


# ── 5c. k600 alternatives for LAKES and RESERVOIRS ───────────────────────────
# Replace the entire 5b block above with one of the following for lentic systems.
# U10 = wind speed at 10 m height [m/s]; add a U10 column to ghg_input.csv.
# Lake_area = surface area [km²]; add a Lake_area column if using Read et al. (2012).
#
# Cole & Caraco (1998) — widely used for lakes, low-to-moderate wind:
#   k600 <- 2.07 + 0.215 * raw$U10^1.7                     # cm/h
#   k600 <- k600 * 24 / 100                                 # m/d
#   Reference: Cole J.J. & Caraco N.F. (1998). Atmospheric exchange of carbon dioxide
#   in a low-wind oligotrophic lake measured by the addition of SF6. Limnol. Oceanogr.,
#   43(4), 647-656. https://doi.org/10.4319/lo.1998.43.4.0647
#
# Wanninkhof (1992) — originally oceanic, commonly applied to large lakes:
#   k600 <- 0.31 * raw$U10^2                                # cm/h (Sc normalised to 660)
#   k600 <- k600 * (600/660)^(-0.5) * 24 / 100             # scale to Sc=600, convert m/d
#   Reference: Wanninkhof R. (1992). Relationship between wind speed and gas exchange
#   over the ocean. J. Geophys. Res., 97(C5), 7373-7382.
#
# Crusius & Wanninkhof (2003) — two-regime wind model for lakes:
#   k600 <- ifelse(raw$U10 < 3.7,
#                  0.72 * raw$U10,
#                  4.33 * raw$U10 - 13.3)                   # cm/h
#   k600 <- k600 * 24 / 100                                 # m/d
#   Reference: Crusius J. & Wanninkhof R. (2003). Gas transfer velocities measured at
#   low wind speed over a lake. Limnol. Oceanogr., 48(3), 1010-1017.
#   https://doi.org/10.4319/lo.2003.48.3.1010
#
# Read et al. (2012) — wind + lake area; accounts for fetch in small lakes:
#   k600 <- 2.51 + 1.48 * raw$U10 + 0.39 * raw$U10 * log10(raw$Lake_area) # cm/h
#   k600 <- k600 * 24 / 100                                 # m/d
#   Reference: Read J.S. et al. (2012). Lake-size dependency of wind shear and convection
#   as controls on gas exchange. Geophys. Res. Lett., 39, L09405.
#   https://doi.org/10.1029/2012GL051886
#
# MacIntyre et al. (2010) — wind + buoyancy flux; recommended for stratified reservoirs:
#   k600 <- 0.0277 * raw$U10^2 + 0.216                      # m/d (direct)
#   Reference: MacIntyre S. et al. (2010). Buoyancy flux, turbulence, and the gas
#   transfer coefficient in a stratified lake. Geophys. Res. Lett., 37, L24604.
#   https://doi.org/10.1029/2010GL044164
# ─────────────────────────────────────────────────────────────────────────────

kCO2 <- k600 * (Sc_CO2 / 600)^(-0.5)           # gas-specific piston velocity (m/d)
kCH4 <- k600 * (Sc_CH4 / 600)^(-0.5)
kN2O <- k600 * (Sc_N2O / 600)^(-0.5)

# ── 6. Air-water fluxes ───────────────────────────────────────────────────────
# F [mmol/m²/d] = k [m/d] × (pGAS_water [µatm] - GAS_air [ppm] × P [atm])
#                          × Kh [mmol/m³/bar] / 1e6 [µatm→atm] / 0.986923 [atm→bar]
# Positive = outgassing; Negative = uptake
#
# pCO2/pCH4/pN2O from complete headspace (cols 4, 12, 15 of pGHG)
# Atmospheric concentrations default to global means if not provided;
#   override by supplying CO2_air_ppm, CH4_air_ppm, N2O_air_ppm in ghg_input.csv

P_atm <- raw$Bar.pressure / 101.325            # kPa → atm

FCO2 <- kCO2 * (as.numeric(pGHG$pCO2_complete_uatm) - raw$CO2_air_ppm * P_atm) *
  Kh_CO2 / 1e6 / 0.986923
FCH4 <- kCH4 * (as.numeric(pGHG$pCH4_simple_uatm)   - raw$CH4_air_ppm * P_atm) *
  Kh_CH4 / 1e6 / 0.986923
FN2O <- kN2O * (as.numeric(pGHG$pN2O_simple_uatm)   - raw$N2O_air_ppm * P_atm) *
  Kh_N2O / 1e6 / 0.986923

# ── 7. Assemble and export results ────────────────────────────────────────────
results <- bind_cols(
  pGHG,
  tibble(
    Tw_C           = raw$Tw,
    SAL            = raw$SAL,
    k600_md        = k600,
    Sc_CO2         = Sc_CO2,
    Sc_CH4         = Sc_CH4,
    Sc_N2O         = Sc_N2O,
    kCO2_md        = kCO2,
    kCH4_md        = kCH4,
    kN2O_md        = kN2O,
    CO2_air_ppm    = raw$CO2_air_ppm,
    CH4_air_ppm    = raw$CH4_air_ppm,
    N2O_air_ppm    = raw$N2O_air_ppm,
    Kh_CO2_mmolm3bar = Kh_CO2,
    Kh_CH4_mmolm3bar = Kh_CH4,
    Kh_N2O_mmolm3bar = Kh_N2O,
    FCO2_mmolm2d   = FCO2,
    FCH4_mmolm2d   = FCH4,
    FN2O_mmolm2d   = FN2O
  )
)

write_csv(results, "ghg_results.csv")
message("Done. Results written to ghg_results.csv")
