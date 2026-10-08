## Aplicacion al WDR2026: eliminacion de casos frente a region de identificacion.
## Uso: Rscript aplicacion_wdr2026.R  (desde la carpeta R/)
`%||%` <- function(a, b) if (is.null(a)) b else a
source("cotas_consistencia.R"); source("region_solucion.R")
OUT <- "../outputs"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "aplicacion_wdr2026.txt"), split = TRUE)

D <- read.csv("../data/wdr2026/datos_calibrados_completos.csv")  # not distributed; see data/wdr2026/README.md
rownames(D) <- D$iso3
cat("Casos:", nrow(D), "\n")

## Intervalos de las condiciones:
##  - NA (sin dato)                 -> [0, 1]
##  - AGENCY/STRAT "indeterminada"  -> [0.05, 0.95] (hoy codificada como 0.5)
intervalo <- function(x, indet = NULL) {
  l <- x; u <- x
  l[is.na(x)] <- 0; u[is.na(x)] <- 1
  if (!is.null(indet)) { l[indet] <- 0.05; u[indet] <- 0.95 }
  cbind(l, u)
}
I_cond <- list(
  NET     = intervalo(D$fNET),
  LANG    = intervalo(D$fLANG_endo_live),
  GOVTECH = intervalo(D$fGOVTECH),
  AGENCY  = intervalo(D$fAGENCY, D$AGENCY_indeterminada),
  STRAT   = intervalo(D$fSTRAT,  D$STRAT_indeterminada))

correr <- function(etiqueta, conds, Yl, Yu, Ypunto, abstiene, incl.cut) {
  cat("\n=====================================================================\n")
  cat(etiqueta, " | condiciones:", paste(conds, collapse = ", "),
      " | corte:", incl.cut, "\n")
  Xl <- sapply(conds, function(c) I_cond[[c]][, 1])
  Xu <- sapply(conds, function(c) I_cond[[c]][, 2])
  rownames(Xl) <- rownames(Xu) <- rownames(D)

  ## (1) Practica convencional: eliminar todo caso con algun faltante,
  ##     indeterminacion o resultado suspendido.
  completo <- rowSums(Xu - Xl > 0) == 0 & !abstiene
  sol_conv <- solucion_completacion(Xl[completo, , drop = FALSE], Ypunto[completo],
                                    incl.cut)
  cat("\n(1) Eliminacion de casos: N =", sum(completo), "de", nrow(D),
      "\n    solucion parsimoniosa:", sol_conv, "\n")

  ## (2) Estado exacto de las filas
  E <- estado_filas(Xl, Xu, Yl, Yu, incl.cut)
  cat("\n(2) Estado de las filas (todas las completaciones):\n")
  print(cbind(E[, conds], E[, c("n_cierto", "n_posible")],
              incl_min = round(E$incl_min, 3), incl_max = round(E$incl_max, 3),
              existencia = E$existencia, OUT = E$OUT), row.names = FALSE)

  ## (3) Region exterior e interior
  ext <- region_exterior(E, conds)
  int <- region_interior(Xl, Xu, Yl, Yu, incl.cut, n_sim = 400)
  cat("\n(3a) Region EXTERIOR (", length(ext), "soluciones distintas ):\n"); print(ext)
  cat("\n(3b) Region INTERIOR (", length(int), "soluciones en 400 completaciones ):\n"); print(int)
  nuc <- nucleo(ext)
  cat("\n(4) Nucleo identificado:", if (length(nuc)) paste(nuc, collapse = " + ") else "VACIO", "\n")
  terminos_conv <- .terminos(sol_conv)
  cat("    Terminos de la solucion convencional dentro del nucleo:",
      if (length(terminos_conv)) paste0(sum(terminos_conv %in% nuc), "/", length(terminos_conv)) else "-",
      "\n    Solucion convencional alcanzable (region interior):", sol_conv %in% names(int), "\n")
  invisible(list(E = E, ext = ext, int = int, nucleo = nuc, conv = sol_conv))
}

## Resultado: intervalo [Y_inf, Y_sup] del WDR2026; ~Y = [1 - Y_sup, 1 - Y_inf]
Yl <- D$Y_inf; Yu <- D$Y_sup

## Modelo A: el hallazgo principal del borrador (~AGENCY + ~STRAT => ~Y, corte 0.9)
r1 <- correr("A. ~Y institucional", c("GOVTECH", "AGENCY", "STRAT"),
             1 - Yu, 1 - Yl, 1 - D$Y_fuzzy, D$Y_abstiene, 0.90)

## Modelo B: solo los faltantes de CONDICIONES (resultado puntual), para aislar
## cuanto aporta cada fuente de ignorancia
r2 <- correr("B. ~Y institucional, resultado puntual", c("GOVTECH", "AGENCY", "STRAT"),
             1 - D$Y_fuzzy, 1 - D$Y_fuzzy, 1 - D$Y_fuzzy, rep(FALSE, nrow(D)), 0.90)

## Modelo C: integrado para Y (corte 0.75)
r3 <- correr("C. Y integrado", c("NET", "LANG", "GOVTECH"),
             Yl, Yu, D$Y_fuzzy, D$Y_abstiene, 0.75)

saveRDS(list(A = r1, B = r2, C = r3), file.path(OUT, "aplicacion_wdr2026.rds"))
sink()
