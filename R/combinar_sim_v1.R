## Combina las replicas de salidas/sim_v1 (las primeras 145 no tienen exterior_timeout ni segundos:
## se guardaron antes de anadir el limite de tiempo el 7-oct-2026 a las 23:20) y escribe
## salidas/simulacion_v1_replicas.csv, que lee resumen_simulacion_v1.R.
## Uso: Rscript combinar_sim_v1.R [carpeta_replicas] [csv_salida]   (desde R/)
args <- commandArgs(trailingOnly = TRUE)
DIR <- if (length(args) >= 1) args[1] else "../outputs/sim_v1"
OUTF <- if (length(args) >= 2) args[2] else "../outputs/simulacion_v1_replicas.csv"
res <- lapply(list.files(DIR, pattern = "^t[0-9]+[.]csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE)
ok <- vapply(res, function(d) !"error" %in% names(d), logical(1))
cat("replicas:", length(res), " con error:", sum(!ok), "\n")
cols <- unique(unlist(lapply(res[ok], names)))
R <- do.call(rbind, lapply(res[ok], function(d) {
  for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
  d[, cols]
}))
write.csv(R, OUTF, row.names = FALSE)
cat("escrito:", OUTF, nrow(R), "filas\n")
