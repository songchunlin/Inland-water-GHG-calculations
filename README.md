# GHG_inland_waters

R script to calculate dissolved greenhouse gas concentrations and air-water fluxes (CO₂, CH₄, N₂O) from raw GC headspace data collected in inland freshwater systems.

## Overview

`GHG_inland_waters.R` takes raw field measurements and gas chromatograph headspace data as input and returns a single CSV file containing dissolved GHG concentrations, gas transfer velocities, Henry's solubility constants, and air-water fluxes. It is designed for streams and rivers but includes commented-out wind-based alternatives for lakes and reservoirs.

The concentration calculations are built around the complete headspace method for CO₂, which accounts for carbonate equilibrium in the equilibration vessel and is more accurate than the simple mass-balance approach, particularly at high alkalinity or high CO₂ concentrations. CH₄ and N₂O are calculated using a simple headspace mass balance, which is appropriate for these gases because they do not participate in acid-base equilibria in water.

## Repository contents

| File | Description |
|---|---|
| `GHG_inland_waters.R` | Main script — concentrations, k₆₀₀, Kₕ, and fluxes |
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
| `HS.mCH4.before` | ppmv | CH₄ mole fraction in headspace **before** equilibration |
| `HS.mCO2.before` | ppmv | CO₂ mole fraction in headspace **before** equilibration |
| `HS.mN2O.before` | ppmv | N₂O mole fraction in headspace **before** equilibration |
| `HS.mCH4.after` | ppmv | CH₄ mole fraction in headspace **after** equilibration |
| `HS.mCO2.after` | ppmv | CO₂ mole fraction in headspace **after** equilibration |
| `HS.mN2O.after` | ppmv | N₂O mole fraction in headspace **after** equilibration |
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
| `CO2_air_ppm` | ppm | Atmospheric CO₂ concentration (global mean: 423) |
| `CH4_air_ppm` | ppm | Atmospheric CH₄ concentration (global mean: 1.93) |
| `N2O_air_ppm` | ppm | Atmospheric N₂O concentration (global mean: 0.338) |

## Output format

`ghg_results.csv` contains one row per sample and the following column groups:

**Dissolved concentrations** (from `Rheadspace_GHG()`)
- CO₂ complete headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L), pH
- CO₂ simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L), % error relative to complete method
- CH₄ simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L)
- N₂O simple headspace: mole fraction (ppmv), partial pressure (µatm), concentration (µmol/L)

**Gas transfer velocity**
- `k600_md`: piston velocity normalised to Sc = 600 (m/d)
- `Sc_CO2`, `Sc_CH4`, `Sc_N2O`: Schmidt numbers at in-situ temperature
- `kCO2_md`, `kCH4_md`, `kN2O_md`: gas-specific piston velocities (m/d)

**Henry's solubility**
- `Kh_CO2_mmolm3bar`, `Kh_CH4_mmolm3bar`, `Kh_N2O_mmolm3bar`: solubility constants (mmol/m³/bar) from `marelac::gas_solubility()`

**Air-water fluxes**
- `FCO2_mmolm2d`, `FCH4_mmolm2d`, `FN2O_mmolm2d`: fluxes in mmol/m²/d. Positive values indicate outgassing to the atmosphere; negative values indicate uptake.

## Methods summary

### Dissolved CO₂ — complete headspace method
CO₂ concentrations are calculated following the complete headspace method of Koschorreck et al. (2021), which reconstructs the original dissolved inorganic carbon (DIC) and pH of the water sample by solving the full carbonate equilibrium system. This accounts for CO₂ that shifted between the gas and liquid phases during equilibration. Freshwater carbonate equilibrium constants (Millero 1979) are used throughout. Henry's Law constants are from Sander (2015).

The script also outputs CO₂ from the simple mass-balance method alongside a % error column to allow comparison between approaches.

### Dissolved CH₄ and N₂O — simple headspace mass balance
CH₄ and N₂O do not participate in carbonate equilibria, so a straightforward mass balance between the gas and liquid phases is used. Henry's Law constants are from Sander (2015).

### Gas transfer velocity (k₆₀₀)
k₆₀₀ is estimated from the energy dissipation model of Raymond et al. (2012) using two regimes:

- **eD > 0.02 m²/s³** (turbulent): `k600 = exp(6.43 + 1.18 × ln(eD))`
- **eD ≤ 0.02 m²/s³** (smooth): `k600 = exp(3.10 + 0.35 × ln(eD))`

where eD = g × v × S (gravitational acceleration × flow velocity × channel slope).

Gas-specific piston velocities are calculated by scaling k₆₀₀ by the ratio of the gas Schmidt number to 600: k_gas = k₆₀₀ × (Sc/600)^−0.5. Schmidt number polynomials follow Wanninkhof (1992).

**For lakes and reservoirs**, wind-based k₆₀₀ parameterisations (Cole & Caraco 1998; Wanninkhof 1992; MacIntyre et al. 2010) are available as commented-out alternatives in section 5 of the script.

### Air-water fluxes
Fluxes are calculated as:

F = k_gas × (pGAS_water − pGAS_air) × Kₕ

where pGAS_water is the dissolved partial pressure (µatm) from the headspace calculation, pGAS_air is the atmospheric partial pressure (atmospheric concentration × barometric pressure), and Kₕ is the Henry's solubility constant (mmol/m³/bar) from `marelac::gas_solubility()`. Pressure units are converted from µatm to bar internally.

## Code provenance

`Rheadspace_GHG()` is derived from two sources:

- **Rheadspace.R** (Marcé, Kim & Prairie 2020) — complete headspace CO₂ method. Licensed under GNU GPL v3. Reference: Koschorreck et al. (2021), *Biogeosciences*, 18, 1619–1627. https://doi.org/10.5194/bg-18-1619-2021
- **headspace_calcs** (K.S. Aho) — simple headspace mass balance for CH₄ and N₂O. https://github.com/kellyaho/headspace_calcs/

Modifications: freshwater-only carbonate constants; Henry's Law constants updated to Sander (2015); CH₄ and N₂O processing added; air-water flux calculations added.

## References

Cole, J.J. & Caraco, N.F. (1998). Atmospheric exchange of carbon dioxide in a low-wind oligotrophic lake measured by the addition of SF₆. *Limnology and Oceanography*, 43(4), 647–656.

Koschorreck M., Prairie Y.T., Kim J. & Marcé R. (2021). Technical note: CO₂ is not like CH₄ – limits of the headspace method to analyse pCO₂ in water. *Biogeosciences*, 18, 1619–1627. https://doi.org/10.5194/bg-18-1619-2021

MacIntyre S., Jonsson A., Jansson M., Aberg J., Turney D.E. & Miller S.D. (2010). Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake. *Geophysical Research Letters*, 37, L24604.

Millero F. (1979). The thermodynamics of the carbonate system in seawater. *Geochimica et Cosmochimica Acta*, 43(10), 1651–1661.

Raymond P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers. *Limnology and Oceanography: Fluids and Environments*, 2, 41–53. https://doi.org/10.1215/21573689-1597669

Sander R. (2015). Compilation of Henry's law constants (version 4.0) for water as solvent. *Atmospheric Chemistry and Physics*, 15, 4399–4981. https://doi.org/10.5194/acp-15-4399-2015

Wanninkhof R. (1992). Relationship between wind speed and gas exchange over the ocean. *Journal of Geophysical Research*, 97(C5), 7373–7382.

## License

This script is shared under the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html), consistent with the license of the Rheadspace.R code from which it is partly derived.
