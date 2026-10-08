## Diagnostico (8-oct-2026): ¿cuantas minimizaciones fallidas ("(sin solucion)") hay en la
## simulacion v1? nucleo() trata una minimizacion fallida como solucion sin terminos y
## devuelve una solucion cierta VACIA; region_exterior() la cuenta como una solucion mas.
## Reproduce las replicas indicadas de simulacion_v1.R (mismas semillas, N_REP = 100) y
## cuenta los fallos en la cota interior y en la exterior.
## Uso: Rscript diag_sin_solucion.R [reps_por_escenario] [nucleos] [ruido]   (desde R/)
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 5
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 4
RUIDO <- if (length(args) >= 3) as.numeric(args[3]) else 0.08
N_REP <- 100
library(parallel); WD <- getwd()
ESC <- expand.grid(p = c(0.05, 0.10, 0.20), mecanismo = c("MCAR", "MNAR_alto", "MNAR_bajo"), stringsAsFactors = FALSE)
tareas <- list(); s <- 0
for (e in seq_len(nrow(ESC))) for (r in seq_len(N_REP)) {
  s <- s + 1
  if (r <= REPS) tareas[[length(tareas) + 1]] <- list(p = ESC$p[e], mecanismo = ESC$mecanismo[e], rep = r, semilla = 20261007 + s)
}
trabajo <- function(t) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  source("cotas_consistencia.R"); source("region_solucion.R")
  set.seed(t$semilla); CONDS <- c("A", "B", "C", "D"); INCL <- 0.85; n <- 60
  X <- matrix(rbeta(n * 4, 0.6, 0.6), n, 4, dimnames = list(NULL, CONDS))
  Y <- pmin(1, pmax(0, pmax(pmin(X[, "A"], X[, "B"]), X[, "C"]) + rnorm(n, 0, RUIDO)))
  prob <- switch(t$mecanismo, MCAR = matrix(t$p, n, 4), MNAR_alto = matrix(pmin(0.95, 2 * t$p * Y), n, 4),
                 MNAR_bajo = matrix(pmin(0.95, 2 * t$p * (1 - Y)), n, 4))
  M <- matrix(runif(n * 4) < prob, n, 4)
  Xl <- X; Xu <- X; Xl[M] <- 0; Xu[M] <- 1; colnames(Xl) <- colnames(Xu) <- CONDS
  int <- region_interior(Xl, Xu, Y, Y, INCL, n_sim = 40)
  E <- estado_filas(Xl, Xu, Y, Y, INCL)
  ext <- tryCatch(region_exterior(E, CONDS, max_ind = 8), error = function(e) NULL)
  ss <- function(tab) if (is.null(tab)) NA else if ("(sin solucion)" %in% names(tab)) tab[["(sin solucion)"]] else 0
  ok <- int[names(int) != "(sin solucion)"]
  data.frame(p = t$p, mecanismo = t$mecanismo, rep = t$rep, int_sin_sol = ss(int), int_vacia = if ("(vacia)" %in% names(int)) int[["(vacia)"]] else 0,
             ext_sin_sol = ss(ext), nucleo_vacio = length(nucleo(int)) == 0,
             nucleo_vacio_sin_fallos = length(nucleo(ok)) == 0, stringsAsFactors = FALSE)
}
cl <- makeCluster(N_CORES); clusterExport(cl, c("trabajo", "WD", "RUIDO")); invisible(clusterEvalQ(cl, setwd(WD)))
R <- do.call(rbind, parLapplyLB(cl, tareas, function(t) tryCatch(trabajo(t), error = function(e)
  data.frame(p = t$p, mecanismo = t$mecanismo, rep = t$rep, int_sin_sol = NA, int_vacia = NA, ext_sin_sol = NA,
             nucleo_vacio = NA, nucleo_vacio_sin_fallos = NA))))
stopCluster(cl)
f <- sprintf("../salidas/diag_sin_solucion_ruido%s.csv", sub("\\.", "", format(RUIDO)))
write.csv(R, f, row.names = FALSE)
cat("Replicas:", nrow(R), "| con fallos en la interior:", sum(R$int_sin_sol > 0, na.rm = TRUE),
    "| fallos totales interior:", sum(R$int_sin_sol, na.rm = TRUE), "| replicas con fallos en la exterior:", sum(R$ext_sin_sol > 0, na.rm = TRUE),
    "| nucleo vacio:", sum(R$nucleo_vacio, na.rm = TRUE), "-> sin fallos:", sum(R$nucleo_vacio_sin_fallos, na.rm = TRUE), "\n")
