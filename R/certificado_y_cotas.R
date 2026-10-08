## Version 2 (SMR): S5 (certificado sin enumerar) y S7 (cotas de termino en datos reales).
## Uso: Rscript certificado_y_cotas.R (desde R/)
##
## S5. Condicion necesaria para que una solucion S0 pertenezca a la region posible R.
## En toda completacion x, una fila OBSERVADA en toda completacion y con OUT = 1 estable
## es una fila positiva de la tabla de x, y una con OUT = 0 estable es negativa. La
## minimizacion cubre todas las filas positivas y ninguna negativa; con include = ""
## (conservadora) tampoco cubre remanentes, y una fila remanente en toda completacion es
## remanente en x. Si S0 viola alguna de estas condiciones, S0 no esta en R (ni en R+).
## No hace falta enumerar.
##
## S7. Cotas exactas [incl-, incl+] y [cov-, cov+] de cada termino publicado sobre
## todas las completaciones (Proposicion 2), junto al valor sobre los casos analizados.
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages(library(QCA))
source("cotas_consistencia.R"); source("region_solucion.R")
OUT <- "../salidas"; sink(file.path(OUT, "certificado_y_cotas.txt"), split = TRUE)

literales <- function(term) {
  l <- strsplit(term, "*", fixed = TRUE)[[1]]
  data.frame(cond = sub("^~", "", l), pres = !startsWith(l, "~"), stringsAsFactors = FALSE)
}
cubre <- function(term, fila) { L <- literales(term); all(as.numeric(fila[L$cond]) == as.numeric(L$pres)) }

certificar <- function(S0, E, nombres, include) {
  tt <- .terminos(S0)
  cub <- vapply(seq_len(nrow(E)), function(r) any(vapply(tt, cubre, logical(1), fila = E[r, nombres])), logical(1))
  fila_txt <- function(r) paste(paste0(ifelse(E[r, nombres] == 1, "", "~"), nombres), collapse = "*")
  v1 <- which(E$existencia == "observada" & E$OUT == "1" & !cub)
  v0 <- which(E$existencia == "observada" & E$OUT == "0" & cub)
  vr <- if (identical(include, "")) which(E$existencia == "remanente" & cub) else integer(0)
  cat("    Filas observadas y OUT=1 estable NO cubiertas:", length(v1), "\n"); for (r in v1) cat("      ", fila_txt(r), "\n")
  cat("    Filas observadas y OUT=0 estable CUBIERTAS:", length(v0), "\n")
  for (r in v0) cat("      ", fila_txt(r), " (n_cierto =", E$n_cierto[r], ", incl max =", round(E$incl_max[r], 3), ")\n")
  if (identical(include, "")) cat("    Remanentes en toda completacion CUBIERTAS (conservadora):", length(vr), "\n")
  cert <- length(v1) + length(v0) + length(vr) > 0
  cat("    ==> S0 fuera de R certificado sin enumerar:", cert, "\n")
  ## terminos que ninguna solucion de R puede contener: cubren una fila observada con OUT=0 estable
  imposibles <- tt[vapply(tt, function(t) any(vapply(which(E$existencia == "observada" & E$OUT == "0"),
                  function(r) cubre(t, E[r, nombres]), logical(1))), logical(1))]
  cat("    Terminos publicados que ninguna completacion admite:", if (length(imposibles)) paste(imposibles, collapse = ", ") else "ninguno", "\n")
  invisible(list(cert = cert, v1 = length(v1), v0 = length(v0), vr = length(vr), imposibles = imposibles))
}

cotas_terminos <- function(S0, Xl, Xu, Yl, Yu, pub) {
  tt <- .terminos(S0)
  do.call(rbind, lapply(tt, function(t) {
    L <- literales(t); j <- match(L$cond, colnames(Xl))
    cond <- lapply(j, function(jj) cbind(Xl[, jj], Xu[, jj]))
    cb <- cotas_ajuste(cond, L$pres, cbind(Yl, Yu))
    condp <- lapply(j, function(jj) cbind(Xl[pub, jj], Xl[pub, jj]))
    cp <- cotas_ajuste(condp, L$pres, cbind(Yl[pub], Yl[pub]))
    data.frame(termino = t, incl_pub = round(cp$suf_min, 3), incl_min = round(cb$suf_min, 3), incl_max = round(cb$suf_max, 3),
               cov_pub = round(cp$cov_min, 3), cov_min = round(cb$cov_min, 3), cov_max = round(cb$cov_max, 3))
  }))
}

analizar <- function(etq, Xl, Xu, Yl, Yu, pub, incl.cut, include) {
  cat("\n=====================================================================\n", etq, "| corte", incl.cut, "| include =", dQuote(include, FALSE), "\n")
  S0 <- solucion_completacion(Xl[pub, ], Yl[pub], incl.cut, include = include, row.dom = TRUE)
  cat("  S0 (eliminacion, N =", length(pub), "):", S0, "\n")
  E <- estado_filas(Xl, Xu, Yl, Yu, incl.cut)
  ce <- certificar(S0, E, colnames(Xl), include)
  ct <- cotas_terminos(S0, Xl, Xu, Yl, Yu, pub)
  cat("  Cotas exactas de cada termino publicado (todas las completaciones de N =", nrow(Xl), "):\n"); print(ct, row.names = FALSE)
  invisible(list(S0 = S0, cert = ce, cotas = ct))
}

## ---- C01 Ettensperger & Schleutker 2025 ----
D <- readRDS(file.path(OUT, "c01_intervalos.rds"))
conds <- c("qca.relfrag", "qca.gdp", "qca.srel1900", "qca.fog", "qca.supp.relgrp", "qca.supp.urbanwk")
cn <- c("RELFRAG", "GDP", "SR1900", "FOG", "SUPPREL", "SUPPURB")
Xl <- as.matrix(D$L[, conds]); Xu <- as.matrix(D$U[, conds]); colnames(Xl) <- colnames(Xu) <- cn
y <- D$L$qca.outcome2.nxx; pub <- seq_len(D$n_pub)
r <- list()
r$c01_Y  <- analizar("C01 Y: altas restricciones", Xl, Xu, y, y, pub, 0.75, "?")
r$c01_nY <- analizar("C01 ~Y: bajas restricciones", Xl, Xu, 1 - y, 1 - y, pub, 0.90, "?")

## ---- C02 Kurz & Ettensperger 2024 ----
D <- readRDS(file.path(OUT, "c02_intervalos.rds"))
cn <- c("YOUTH", "GDP", "ELEC", "DEC", "PROG", "YOUNG")
Xl <- as.matrix(D$L[, cn]); Xu <- as.matrix(D$U[, cn]); y <- D$Y; pub <- seq_len(D$n_pub)
r$c02_Y_p  <- analizar("C02 Y: muchos jovenes (parsimoniosa)", Xl, Xu, y, y, pub, 0.8, "?")
r$c02_nY_p <- analizar("C02 ~Y: pocos jovenes (parsimoniosa)", Xl, Xu, 1 - y, 1 - y, pub, 0.8, "?")
r$c02_Y_c  <- analizar("C02 Y: muchos jovenes (conservadora)", Xl, Xu, y, y, pub, 0.8, "")
r$c02_nY_c <- analizar("C02 ~Y: pocos jovenes (conservadora)", Xl, Xu, 1 - y, 1 - y, pub, 0.8, "")
saveRDS(r, file.path(OUT, "certificado_y_cotas.rds"))
sink()
