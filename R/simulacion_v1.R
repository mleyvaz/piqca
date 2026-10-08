## =============================================================================
## simulacion_v1.R  --  Eliminacion de casos, imputacion y cotas (v1)
## =============================================================================
## Correcciones respecto de la v0:
##  (1) linea base con datos completos en cada replica (mismo QCA, sin faltantes);
##  (2) prueba EXACTA de imposibilidad: la solucion por eliminacion se busca en
##      la region EXTERIOR (que contiene a la verdadera); si no esta, ninguna
##      completacion la produce. Solo cuando la region exterior es enumerable;
##      si no, NA;
##  (3) 200 replicas por escenario, en paralelo;
##  (4) "termino falso": implicante que no implica el modelo verdadero A*B + C,
##      es decir, cuyo conjunto de literales no contiene {A, B} ni {C}.
## Uso: Rscript simulacion_v1.R [n_rep] [n_nucleos]   (desde la carpeta R/)
## =============================================================================
args <- commandArgs(trailingOnly = TRUE)
N_REP <- if (length(args) >= 1) as.integer(args[1]) else 200
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 3
OUT <- "../outputs"; dir.create(OUT, showWarnings = FALSE)
LIM_EXT <- 120; LIM_INT <- 900   # segundos por etapa
library(parallel)
WD <- getwd()

ESCENARIOS <- expand.grid(p = c(0.05, 0.10, 0.20),
                          mecanismo = c("MCAR", "MNAR_alto", "MNAR_bajo"),
                          stringsAsFactors = FALSE)

trabajo <- function(tarea) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  source("cotas_consistencia.R"); source("region_solucion.R")
  set.seed(tarea$semilla)
  CONDS <- c("A", "B", "C", "D"); INCL <- 0.85
  n <- 60
  X <- matrix(rbeta(n * 4, 0.6, 0.6), n, 4, dimnames = list(NULL, CONDS))
  Y <- pmin(1, pmax(0, pmax(pmin(X[, "A"], X[, "B"]), X[, "C"]) + rnorm(n, 0, 0.08)))
  prob <- switch(tarea$mecanismo,
    MCAR      = matrix(tarea$p, n, 4),
    MNAR_alto = matrix(pmin(0.95, 2 * tarea$p * Y), n, 4),
    MNAR_bajo = matrix(pmin(0.95, 2 * tarea$p * (1 - Y)), n, 4))
  M <- matrix(runif(n * 4) < prob, n, 4)
  Xl <- X; Xu <- X; Xl[M] <- 0; Xu[M] <- 1; colnames(Xl) <- colnames(Xu) <- CONDS

  falso <- function(term) {
    lit <- strsplit(term, "\\*")[[1]]
    !(("C" %in% lit) || all(c("A", "B") %in% lit))
  }
  n_falsos <- function(sol) sum(vapply(.terminos(sol), falso, logical(1)))
  exacta <- function(sol) identical(sort(.terminos(sol)), c("A*B", "C"))

  sol_full <- solucion_completacion(X, Y, INCL)
  completo <- rowSums(M) == 0
  sol_del <- if (sum(completo) >= 5)
    solucion_completacion(X[completo, , drop = FALSE], Y[completo], INCL) else "(sin solucion)"
  Xm <- X; for (j in 1:4) if (any(M[, j])) Xm[M[, j], j] <- mean(X[!M[, j], j])
  sol_imp <- solucion_completacion(Xm, Y, INCL)

  ## Limite de tiempo por etapa (anadido 7-oct-2026, 23:20): dos replicas dejaron la corrida
  ## parada ~3 h. Si la region exterior no se enumera en LIM_EXT segundos se registra como no
  ## enumerada (NA), como ya preveia el diseno; si la interior no termina en LIM_INT, la replica
  ## se guarda como error (se cuenta y se informa).
  con_limite <- function(expr, s) {
    setTimeLimit(elapsed = s, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
    expr
  }
  t0 <- proc.time()[["elapsed"]]
  int <- con_limite(region_interior(Xl, Xu, Y, Y, INCL, n_sim = 40), LIM_INT)
  nuc <- nucleo(int)
  E <- estado_filas(Xl, Xu, Y, Y, INCL)
  ext_timeout <- FALSE
  ext <- tryCatch(con_limite(region_exterior(E, CONDS, max_ind = 8), LIM_EXT), error = function(e) {
    if (grepl("time limit|limite de tiempo|elapsed", conditionMessage(e), ignore.case = TRUE)) ext_timeout <<- TRUE
    NULL })

  data.frame(
    p = tarea$p, mecanismo = tarea$mecanismo, rep = tarea$rep,
    n_completos = sum(completo), celdas_faltantes = sum(M),
    full_exacta = exacta(sol_full), full_falsos = n_falsos(sol_full),
    del_exacta = exacta(sol_del), del_falsos = n_falsos(sol_del),
    del_igual_full = sol_del == sol_full,
    imp_exacta = exacta(sol_imp), imp_falsos = n_falsos(sol_imp),
    imp_igual_full = sol_imp == sol_full,
    nucleo = paste(nuc, collapse = " + "), nucleo_vacio = length(nuc) == 0,
    nucleo_falsos = sum(vapply(nuc, falso, logical(1))),
    nucleo_en_full = all(nuc %in% .terminos(sol_full)),
    full_en_interior = sol_full %in% names(int),
    del_en_interior = sol_del %in% names(int),
    exterior_enumerada = !is.null(ext),
    del_imposible = if (is.null(ext)) NA else !(sol_del %in% names(ext)),
    n_sol_interior = length(int),
    exterior_timeout = ext_timeout,
    segundos = round(proc.time()[["elapsed"]] - t0, 1),
    stringsAsFactors = FALSE)
}

tareas <- list(); s <- 0
for (e in seq_len(nrow(ESCENARIOS))) for (r in seq_len(N_REP)) {
  s <- s + 1
  tareas[[s]] <- list(p = ESCENARIOS$p[e], mecanismo = ESCENARIOS$mecanismo[e],
                      rep = r, semilla = 20261007 + s)
}
cat(format(Sys.time(), "%H:%M:%S"), "tareas:", length(tareas), "nucleos:", N_CORES, "\n")
cl <- makeCluster(N_CORES); clusterSetRNGStream(cl, 20261007)
clusterExport(cl, c("trabajo", "WD", "LIM_EXT", "LIM_INT"))
invisible(clusterEvalQ(cl, setwd(WD)))
## Guardado por replica: cada tarea escribe su fila y se salta si ya existe, de
## modo que una interrupcion no pierde el trabajo hecho (relanzar = reanudar).
DIRREP <- file.path(OUT, "sim_v1"); dir.create(DIRREP, showWarnings = FALSE)
clusterExport(cl, "DIRREP")
pendientes <- Filter(function(t) !file.exists(file.path(DIRREP, sprintf("t%04d.csv", t$id))),
                     lapply(seq_along(tareas), function(i) c(tareas[[i]], id = i)))
cat("pendientes:", length(pendientes), "\n")
invisible(parLapplyLB(cl, pendientes, function(t) {
  d <- tryCatch(trabajo(t), error = function(e)
    data.frame(p = t$p, mecanismo = t$mecanismo, rep = t$rep, error = conditionMessage(e)))
  write.csv(d, file.path(DIRREP, sprintf("t%04d.csv", t$id)), row.names = FALSE)
  NULL }))
stopCluster(cl)
res <- lapply(list.files(DIRREP, pattern = "^t[0-9]+\\.csv$", full.names = TRUE),
              read.csv, stringsAsFactors = FALSE)
ok <- vapply(res, function(d) !"error" %in% names(d), logical(1))
cat("replicas con error:", sum(!ok), "\n")
if (any(!ok)) print(head(do.call(rbind, lapply(res[!ok], function(d) d[, c("p", "mecanismo", "error")]))))
R <- do.call(rbind, res[ok])
write.csv(R, file.path(OUT, "simulacion_v1_replicas.csv"), row.names = FALSE)
cat(format(Sys.time(), "%H:%M:%S"), "fin\n")
