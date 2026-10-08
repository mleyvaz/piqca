# PI-QCA: partial-identification QCA

Replication code for:

> Leyva-Vazquez, M. Y. (2026). *What QCA Cannot See: Identification Bounds for Configurational Solutions with Missing
> Set Memberships*. Preprint.

PI-QCA represents each missing set-membership score by an interval, considers every completion of the data, and
reports exact bounds for the consistency and coverage of any term, the **sure solution** (terms present in the
solution of every completion), the **possible region** of solutions (an inner bound from sampled completions and an
outer bound from enumeration) and the **identification gap**. It is implemented in R on top of the `QCA` package.

## Structure
| Path | Content |
|---|---|
| `R/cotas_consistencia.R` | Exact bounds of consistency and coverage of a term over all completions (Dinkelbach) |
| `R/region_solucion.R` | Row statuses, inner bound, outer bound and sure solution |
| `R/test_fuerza_bruta.R` | Check of the exact bounds against exhaustive grid search |
| `R/simulacion_v1.R`, `R/combinar_sim_v1.R`, `R/resumen_simulacion_v1.R` | Simulation study (900 replicates) and its summary table |
| `R/simulacion_v1_ruido004.R` | Robustness variant with lower noise |
| `R/c01_*.R` | Reanalysis of Ettensperger and Schleutker (2025) |
| `R/c02_*.R` | Reanalysis of Kurz and Ettensperger (2024) |
| `R/aplicacion_wdr2026.R` | Illustration with 131 economies (data on request) |
| `outputs/` | Outputs reported in the article |

Comments inside the scripts are in Spanish.

## Reproducing the results
Requirements: R (tested with 4.5.1) and the packages `QCA` (tested with 3.23) and `parallel`. Run every script from
the `R/` folder, for example:

```
Rscript simulacion_v1.R 100 6      # 100 replicates per scenario, 6 cores; resumes if interrupted
Rscript combinar_sim_v1.R
Rscript resumen_simulacion_v1.R
Rscript c02_preparar.R && Rscript c02_piqca.R
```

For the reanalysis of Ettensperger and Schleutker (2025), first download their data (see `data/C01_ettensperger/`).

## Licence
Code: MIT. Data in `data/C02_kurz`: CC BY 4.0 (see its README).

## Use of generative AI
Generative AI tools (Claude, Anthropic; Codex, OpenAI) were used to write and check parts of this code; the author
verified the results against brute-force checks and the outputs reported here.
