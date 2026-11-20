# Where dimensionality matters and how to quantify it

R code and data for the paper, written in R 4.3.3 with base R and the package `parallel`.


## 1. Layout

```
code/
  01_exact_simul_figures.R      main scripts 01 to 11, run in this order
    ...                         (see 3. Scripts)
  11_supp_note11_third_moment.R

  04_realdata_core.R            functions used by 04 (not run on its own)
  99_verify_all.R               implementation checks (run separately)
  paths.R                       folders and the FULL switch

  run_local.bat                 runs 01 to 11 on Windows
  run                           (Code Ocean launcher)

data/
  quakes.csv, magic.csv, wine_white.csv, letter.csv, supercond.csv
  sim_ch6.rds (simulation cache)

results/
  figures/  tables/  source_data/  checks/  intermediate/  logs/
  organize_results.bat    sorts a fresh run into the folders above (run after all scripts)
```


## 2. How to run

On Windows, run `code/run_local.bat`. It finds the folders by itself, writes a log for each script in `results/`, and records one line per script (exit code, whether the expected output file was written, elapsed time) in `results/run_summary.txt`. A failing script does not stop the run; check `run_summary.txt` before using the results.

Otherwise, set `LOCAL_DATA`, `LOCAL_OUT` and `CODE_DIR` (or edit the defaults in `code/paths.R`) and run the scripts with Rscript in numerical order.

By default the simulation is read from the cache and the real-data part uses reduced settings (8 subsamples × 100 trees). This takes about an hour, most of it in 07 and 10. With `FULL=1` the simulation is rerun, the cache is rewritten to `results/`, and the real-data part uses the settings of the paper (20 × 200); this takes several hours.


## 3. Scripts

| No. | Description |
|---|---|
| 01 | exact recursion, simulation, Figs. 2–3, Extended Data Figs. 1–3, Source Data |
| 02 | closed forms against the recursion |
| 03 | Supplementary Table 1 |
| 04 | real data: five data sets, original and rank-transformed |
| 05 | Table 2, Extended Data Table 2 and Fig. 4 (drawn from Extended Data Table 2) |
| 06 | Supplementary Table 2 |
| 07 | Supplementary Table 3 (copulas) |
| 08 | Supplementary Tables 4 and 5 |
| 09 | Supplementary Table 6 (zeta limits) |
| 10 | Supplementary Tables 7–9 (maximal state) |
| 11 | third moment of R, closed form against numerical value |
| 99 | implementation checks (closed forms against the recursion, Lemmas 5–6, kernel row sums); not part of `run_local.bat`, run separately with Rscript; writes `verify_closed_forms.csv` |

Scripts 04, 06, 07 and 10 are Monte Carlo with fixed seeds. 07 distributes jobs dynamically, so its last digits can differ between runs within the standard errors.


## 4. Output

The scripts write everything to `results/`. After all scripts have finished without error, running `results/organize_results.bat` sorts the files into the folders below. Do not run it between scripts: 02, 03 and 05 read files from the top of `results/`.

| | |
|---|---|
| `figures/` | `Fig2`–`Fig4`, `ExtData_Fig1`–`ExtData_Fig3` (pdf, png) |
| `tables/` | `ExtDataTable1_exact_moments.csv` (Table 1 is the finite-d rows), `Table2_*`, `ExtDataTable2_*`, `SuppTable1_*` to `SuppTable9_*` |
| `source_data/` | one csv per figure |
| `checks/` | `closed_form_check.csv` (02), `third_moment_check.csv` (11) |
| `intermediate/` | `realdata_*.csv` (04, read by 05), `copula_*` (07) |
| `logs/` | one log per script and `run_summary.txt` |


## 5. Data

`quakes.csv` (1,000 × 5) is from the R datasets package. `magic.csv` (19,020 × 10), `wine_white.csv` (4,898 × 11), `letter.csv` (20,000 × 16) and `supercond.csv` (21,263 × 81) are from the UCI Machine Learning Repository. Only the explanatory variables are kept.

`sim_ch6.rds` holds the simulation of the model for n = 256, 512 and d = 1–10, 20, 50, 100, with 20 replicates × 500 trees and seed 20260928.


## 6. License

Code: MIT (see LICENSE).

Data:

| File | Source | License |
|---|---|---|
| quakes.csv | R datasets package _ base | GPL-3 |
| magic.csv | UCI, MAGIC Gamma Telescope, doi:10.24432/C52C8B | CC BY 4.0 |
| wine_white.csv | UCI, Wine Quality, doi:10.24432/C56S3T | CC BY 4.0 |
| letter.csv | UCI, Letter Recognition, doi:10.24432/C5ZP40 | CC BY 4.0 |
| supercond.csv | UCI, Superconductivity Data, doi:10.24432/C53P47 | CC BY 4.0 |
