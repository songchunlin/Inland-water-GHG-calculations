# Inland water GHG calculations

R script to calculate dissolved greenhouse gas concentrations and air-water fluxes (CO<sub>2</sub>, CH<sub>4</sub>, N<sub>2</sub>O) from raw GC headspace data collected in inland freshwater systems.

## Overview

`GHG_inland_waters.R` takes raw field measurements and gas chromatograph headspace data as input and returns a single CSV file containing dissolved GHG concentrations, gas transfer velocities, Henry's solubility constants, and air-water fluxes. It is designed for streams and rivers but includes commented-out wind-based alternatives for lakes and reservoirs.

The concentration calculations are built around the complete headspace method for CO<sub>2</sub>, which accounts for carbonate equilibrium in the equilibration vessel and is more accurate than the simple headspace calculation approach, particularly at high alkalinity or high CO<sub>2</sub> concentrations. CH<sub>4</sub> and N<sub>2</sub>O are calculated using a simple headspace method, which is appropriate for these gases because they do not participate in acid-base equilibria in water.

## Repository contents

| File | Description |
|---|---|
| `GHG_inland_waters.R` | Main script — concentrations, k<sub>600</sub>, K<sub>h</sub>, and fluxes |
| `ghg_input.csv` | Input template with column names, units, and an example row |
| `README.md` | This file |

## Dependencies

The script uses three R packages. `marelac` must be installed separately; `readr` and `dplyr` are part of the tidyverse.

```r
install.packages(c("marelac", "readr", "dplyr"))
```

Developed and tested on R ≥ 4.1.

## Usage

1. Clone or download this repository.
2. Fill in `ghg_input.csv` with your data (see [Input format](#input-format) below). Delete the units row (row 2) and notes row (row 3) if you prefer a plain single-header CSV, and remove the `skip = 2` argument from `read_csv()` in the script accordingly.
3. Place `GHG_inland_waters.R` and your completed `ghg_input.csv` in the same working directory.
4. Run the script:

```r
source("GHG_inland_waters.R")
```

Results are written to `ghg_results.csv` in the same directory.

## Input format

`ghg_input.csv` must contain the following 21 columns. Column names are case-sensitive.

| Column | Units | Description |
|---|---|---|
| `site` | — | Site identifier |
| `datetime.EST` | YYYY-MM-DD HH:MM | Sample date and time |
| `HS.mCH4.before` | ppmv | CH<sub>4</sub> mole fraction in headspace **before** equilibration |
| `HS.mCO2.before` | ppmv | CO<sub>2</sub> mole fraction in headspace **before** equilibration |
| `HS.mN2O.before` | ppmv | N<sub>2</sub>O mole fraction in headspace **before** equilibration |
| `HS.mCH4.after` | ppmv | CH<sub>4</sub> mole fraction in headspace **after** equilibration |
| `HS.mCO2.after` | ppmv | CO<sub>2</sub> mole fraction in headspace **after** equilibration |
| `HS.mN2O.after` | ppmv | N<sub>2</sub>O mole fraction in headspace **after** equilibration |
| `Temp.insitu` | °C | In-situ water temperature at sampling |
| `Temp.equil` | °C | Water temperature during headspace equilibration |
| `Alkalinity.measured` | µeq/L | Total alkalinity |
| `Volume.gas` | mL | Volume of gas in the headspace vessel |
| `Volume.water` | mL | Volume of water in the headspace vessel |
| `Bar.pressure` | kPa | Barometric pressure at field conditions (101.325 kPa = 1 atm) |
| `Tw` | °C | In-situ water temperature for Schmidt number calculation (typically = `Temp.insitu`) |
| `v_ms` | m/s | Mean flow velocity (streams/rivers; see note for lakes) |
| `slope` | m/m | Channel slope (streams/rivers; see note for lakes) |
| `SAL` | PSU | Salinity (set to 0 for freshwater) |
| `CO2_air_ppm` | ppm | Atmospheric CO<sub>2</sub> concentration (global mean in 2024: 423) |
| `CH4_air_ppm` | ppm | Atmospheric CH<sub>4</sub> concentration (global mean in 2024: 1.93) |
| `N2O_air_ppm` | ppm | Atmospheric N<sub>2</sub>O concentration (global mean in 2024: 0.338) |

## Output format

`ghg_results.csv` contains one row per sample and the following column groups:

**Dissolved concentrations** (from `Rheadspace_GHG()`)
- CO<sub>2</sub> complete headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L), pH
- CO<sub>2</sub> simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L), % error relative to complete method
- CH<sub>4</sub> simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L)
- N<sub>2</sub>O simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L)

**Gas transfer velocity**
- `k600_md`: piston velocity normalised to Sc = 600 (m/d)
- `Sc_CO2`, `Sc_CH4`, `Sc_N2O`: Schmidt numbers at in-situ temperature
- `kCO2_md`, `kCH4_md`, `kN2O_md`: gas-specific piston velocities (m/d)

**Henry's solubility**
- `Kh_CO2_mmolm3bar`, `Kh_CH4_mmolm3bar`, `Kh_N2O_mmolm3bar`: solubility constants (mmol/m<sup>3</sup>/bar) from `marelac::gas_solubility()`

**Air-water fluxes**
- `FCO2_mmolm2d`, `FCH4_mmolm2d`, `FN2O_mmolm2d`: fluxes in mmol/m<sup>2</sup>/d. Positive values indicate outgassing to the atmosphere; negative values indicate uptake.

## Methods summary

### Dissolved CO<sub>2</sub> — complete headspace method
CO<sub>2</sub> concentrations are calculated following the complete headspace method of Koschorreck et al. (2021), which reconstructs the original dissolved inorganic carbon (DIC) and pH of the water sample by solving the full carbonate equilibrium system. This accounts for CO<sub>2</sub> that shifted between the gas and liquid phases during equilibration. Freshwater carbonate equilibrium constants (Millero 1979) are used throughout. Henry's Law constants are from Sander (2015).

The script also outputs CO<sub>2</sub> from the simple mass-balance method alongside a % error column to allow comparison between approaches.

### Dissolved CH<sub>4</sub> and N<sub>2</sub>O — simple headspace mass balance
CH<sub>4</sub> and N<sub>2</sub>O do not participate in carbonate equilibria, so a straightforward mass balance between the gas and liquid phases is used. Henry's Law constants are from Sander (2015).

### Henry's solubility constants (K<sub>h</sub>)
`marelac::gas_solubility()` returns the solubility coefficient K<sub>0</sub> [mmol/m<sup>3</sup>/bar], which is equivalent to the Henry's volatility constant K<sub>h</sub> in the form [gas]_aq / p_gas — the correct form for the flux equation F = k × Δp × K<sub>h</sub>. Using K<sub>0</sub> as K<sub>h</sub> is standard practice in inland water GHG studies. The underlying equations are Weiss (1974) for CO<sub>2</sub>, Wiesenburg & Guinasso (1979) for CH<sub>4</sub>, and Weiss & Price (1980) for N<sub>2</sub>O, all of which are valid for freshwater (SAL = 0) and brackish systems across the 0–30°C range.

Note that K<sub>0</sub> from `marelac` differs from the Sander (2015) K<sub>h</sub> used inside `Rheadspace_GHG()` for the headspace equilibrium step. Both formulations are appropriate for their respective purposes; the partial pressures produced by the headspace calculation are independent of which K<sub>h</sub> is used for the flux step.

### Schmidt numbers
Schmidt numbers (Sc = kinematic viscosity / molecular diffusivity) are used to scale k<sub>600</sub> to a gas-specific piston velocity: k_gas = k<sub>600</sub> × (Sc/600)^−0.5, where the exponent −0.5 applies to wind- or turbulence-mixed surfaces. The script uses 3rd-order freshwater polynomials from Wanninkhof (1992) by default. A commented-out alternative using 4th-order polynomials from Wanninkhof (2014) is also provided and is more accurate across a wider temperature range, particularly above 20°C.

### Gas transfer velocity (k<sub>600</sub>)

**Active — streams and rivers:** The two-regime energy-dissipation model of Ulseth et al. (2019) is used as the default and recommended method:

- **eD > 0.02 m<sup>2</sup>/s<sup>3</sup>** (turbulent, high-gradient streams): `k600 = exp(6.43 + 1.18 × ln(eD))`
- **eD ≤ 0.02 m<sup>2</sup>/s<sup>3</sup>** (smooth, low-gradient rivers): `k600 = exp(3.10 + 0.35 × ln(eD))`

where eD = g × v × S (gravitational acceleration × flow velocity × channel slope). This two-regime structure captures the disproportionately high gas exchange in mountain and steep headwater streams that a single power law underestimates.

**Alternative stream/river models** (commented out in section 5b of the script):

| Model | Equation | Notes |
|---|---|---|
| Raymond et al. (2012) | k<sub>600</sub> = 2841 × V × S + 2.02 (m/d) | Linear model using velocity and slope, k600 = VS × 2841 ± 107+2.02 ± 0.209; also widely used |

**Alternative lake/reservoir models** (commented out in section 5c of the script):

| Model | Equation | Notes |
|---|---|---|
| Cole & Caraco (1998) | k<sub>600</sub> = 2.07 + 0.215 × U<sub>10</sub><sup>1.7</sup> (cm/h) | Default for small/medium lakes |
| Wanninkhof (1992) | k<sub>600</sub> = 0.31 × U<sub>10</sub><sup>2</sup> (cm/h) | Large lakes, cross-study comparisons |
| Crusius & Wanninkhof (2003) | Two-regime wind model (cm/h) | Smooth-to-rough wind transition |
| Read et al. (2012) | k<sub>600</sub> = 2.51 + 1.48×U<sub>10</sub> + 0.39×U<sub>10</sub>×log<sub>10</sub>(A) (cm/h) | Accounts for lake fetch/area |
| MacIntyre et al. (2010) | k<sub>600</sub> = 0.0277 × U<sub>10</sub><sup>2</sup> + 0.216 (m/d) | Stratified reservoirs with buoyancy flux |

### Air-water fluxes
Fluxes are calculated as:

F = k_gas × (pGAS_water − pGAS_air) × K<sub>h</sub>

where pGAS_water is the dissolved partial pressure (µatm) from the headspace calculation, pGAS_air is the atmospheric partial pressure (atmospheric concentration × barometric pressure), and K<sub>h</sub> is the Henry's solubility constant (mmol/m<sup>3</sup>/bar) from `marelac::gas_solubility()`. Pressure units are converted from µatm to bar internally.

## Code provenance

`Rheadspace_GHG()` is derived from two sources:

- **Rheadspace.R** (Marcé, Kim & Prairie 2020) — complete headspace CO<sub>2</sub> method. Licensed under GNU GPL v3. Reference: Koschorreck et al. (2021), *Biogeosciences*, 18, 1619–1627. https://doi.org/10.5194/bg-18-1619-2021
- **headspace_calcs** (K.S. Aho) — simple headspace mass balance for CH<sub>4</sub> and N<sub>2</sub>O. https://github.com/kellyaho/headspace_calcs/

Modifications: freshwater-only carbonate constants; Henry's Law constants updated to Sander (2015); CH<sub>4</sub> and N<sub>2</sub>O processing added; air-water flux calculations added.

## References

Cole, J.J. & Caraco, N.F. (1998). Atmospheric exchange of carbon dioxide in a low-wind oligotrophic lake measured by the addition of SF<sub>6</sub>. *Limnology and Oceanography*, 43(4), 647–656.

Crusius, J. & Wanninkhof, R. (2003). Gas transfer velocities measured at low wind speed over a lake. *Limnology and Oceanography*, 48(3), 1010–1017.

Koschorreck, M., Prairie, Y.T., Kim, J. & Marcé, R. (2021). Technical note: CO<sub>2</sub> is not like CH<sub>4</sub> – limits of the headspace method to analyse pCO<sub>2</sub> in water. *Biogeosciences*, 18, 1619–1627. https://doi.org/10.5194/bg-18-1619-2021

MacIntyre, S., Jonsson, A., Jansson, M., Aberg, J., Turney, D.E. & Miller, S.D. (2010). Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake. *Geophysical Research Letters*, 37, L24604.

Millero, F. (1979). The thermodynamics of the carbonate system in seawater. *Geochimica et Cosmochimica Acta*, 43(10), 1651–1661.

O'Connor, D.J. & Dobbins, W.E. (1958). Mechanism of reaeration in natural streams. *Transactions of the American Society of Civil Engineers*, 123, 641–666.

Owens, M., Edwards, R.W. & Gibbs, J.W. (1964). Some reaeration studies in streams. *International Journal of Air and Water Pollution*, 8, 469–486.

Raymond, P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers. *Limnology and Oceanography: Fluids and Environments*, 2, 41–53. https://doi.org/10.1215/21573689-1597669

Read, J.S. et al. (2012). Lake-size dependency of wind shear and convection as controls on gas exchange. *Geophysical Research Letters*, 39, L09405. https://doi.org/10.1029/2012GL051886

Sander, R. (2015). Compilation of Henry's law constants (version 4.0) for water as solvent. *Atmospheric Chemistry and Physics*, 15, 4399–4981. https://doi.org/10.5194/acp-15-4399-2015

Ulseth, A.J. et al. (2019). Distinct air-water gas exchange regimes in low- and high-energy streams. *Nature Geoscience*, 12, 259–263. https://doi.org/10.1038/s41561-019-0324-8

Wanninkhof, R. (1992). Relationship between wind speed and gas exchange over the ocean. *Journal of Geophysical Research*, 97(C5), 7373–7382.

Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. *Limnology and Oceanography: Methods*, 12, 351–362. https://doi.org/10.4319/lom.2014.12.351

Weiss, R.F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. *Marine Chemistry*, 2, 203–215.

Weiss, R.F. & Price, B.A. (1980). Nitrous oxide solubility in water and seawater. *Marine Chemistry*, 8(5), 347–359.

Wiesenburg, D.A. & Guinasso, N.L. (1979). Equilibrium solubilities of methane, carbon monoxide, and hydrogen in water and sea water. *Journal of Chemical and Engineering Data*, 24(4), 356–360.

## License

This script is shared under the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html), consistent with the license of the Rheadspace.R code from which it is partly derived.
