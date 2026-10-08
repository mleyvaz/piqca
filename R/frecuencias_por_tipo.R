## Version 2 (SMR), S9: frecuencias de terminos en la cota interior separadas por tipo de
## completacion (extremas dirigidas, vertices aleatorios, puntos uniformes del interior)
## y curva de saturacion del numero de soluciones distintas.
## Misma semilla y mismo orden de extracciones que region_interior() (region_solucion.R),
## de modo que la union reproduce las tablas publicadas en el preprint v1.
## Uso: Rscript frecuencias_por_tipo.R (desde R/)
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages(library(QCA))
source("cotas_consistencia.R"); source("region_solucion.R")
OUT <- "../salidas"; sink(file.path(OUT, "frecuencias_por_tipo.txt"), split = TRUE)

interior_etiquetado <- function(Xl, Xu, Yl, Yu, incl.cut, include = "?", n_sim, semilla = 20261006) {
  set.seed(semilla)
  k <- ncol(Xl); G <- .filas(k, colnames(Xl))
  cond <- lapply(seq_len(k), function(j) cbind(Xl[, j], Xu[, j]))
  ext <- character(0)
  for (r in seq_len(nrow(G))) for (sen in c("min", "max")) {
    ce <- completacion_extrema(cond, G[r, ] == 1, cbind(Yl, Yu), "suf", sen)
    X <- ce$X; dimnames(X) <- dimnames(Xl)
    ext <- c(ext, solucion_completacion(X, ce$Y, incl.cut, 1, include, TRUE))
  }
  sols <- character(n_sim); tipo <- character(n_sim)
  for (s in seq_len(n_sim)) {
    modo <- if (s %% 2 == 0) "vertice" else "interior"
    X <- matrix(completar(Xl, Xu, modo), nrow(Xl), dimnames = dimnames(Xl))
    Y <- completar(Yl, Yu, modo)
    sols[s] <- solucion_completacion(X, Y, incl.cut, 1, include, TRUE); tipo[s] <- modo
  }
  data.frame(sol = c(sols, ext), tipo = c(tipo, rep("extrema", length(ext))), orden = c(seq_len(n_sim), rep(NA, length(ext))),
             stringsAsFactors = FALSE)
}

resumir <- function(etq, d, terminos) {
  cat("\n=====================================================================\n", etq, "\n")
  cat("  Completaciones por tipo:\n"); print(table(d$tipo))
  fr <- sapply(c("extrema", "vertice", "interior"), function(tp) {
    s <- d$sol[d$tipo == tp]
    sapply(terminos, function(t) mean(vapply(s, function(x) t %in% .terminos(x), logical(1))))
  })
  fr <- cbind(fr, todas = sapply(terminos, function(t) mean(vapply(d$sol, function(x) t %in% .terminos(x), logical(1)))))
  cat("  Frecuencia de cada termino publicado por tipo de completacion:\n"); print(round(fr, 3))
  cat("  Soluciones distintas por tipo:", sapply(c("extrema", "vertice", "interior"), function(tp) length(unique(d$sol[d$tipo == tp]))),
      "| total:", length(unique(d$sol)), "\n")
  ## saturacion: soluciones distintas acumuladas tras las extremas + las primeras m aleatorias
  base <- unique(d$sol[d$tipo == "extrema"]); al <- d$sol[d$tipo != "extrema"]
  hitos <- unique(pmin(c(10, 25, 50, 100, 200, 300, 400), length(al)))
  sat <- sapply(hitos, function(m) length(unique(c(base, al[seq_len(m)]))))
  cat("  Saturacion (extremas + m aleatorias -> soluciones distintas):\n"); print(setNames(sat, hitos))
  sat_al <- sapply(hitos, function(m) length(unique(al[seq_len(m)])))
  cat("  Solo aleatorias (m -> soluciones distintas):\n"); print(setNames(sat_al, hitos))
  invisible(list(frec = fr, saturacion = setNames(sat, hitos), saturacion_aleatorias = setNames(sat_al, hitos)))
}

res <- list()
## C01 altas restricciones (528 = 128 extremas + 400 aleatorias)
D <- readRDS(file.path(OUT, "c01_intervalos.rds"))
conds <- c("qca.relfrag", "qca.gdp", "qca.srel1900", "qca.fog", "qca.supp.relgrp", "qca.supp.urbanwk")
cn <- c("RELFRAG", "GDP", "SR1900", "FOG", "SUPPREL", "SUPPURB")
Xl <- as.matrix(D$L[, conds]); Xu <- as.matrix(D$U[, conds]); colnames(Xl) <- colnames(Xu) <- cn
rownames(Xl) <- rownames(Xu) <- D$L$name
y <- D$L$qca.outcome2.nxx
pubY <- c("~RELFRAG*SR1900*~FOG", "GDP*SR1900*~FOG", "GDP*SR1900*~SUPPURB", "SR1900*~FOG*~SUPPREL*SUPPURB", "SR1900*SUPPREL*~SUPPURB")
d <- interior_etiquetado(Xl, Xu, y, y, 0.75, "?", 400); res$c01_Y <- resumir("C01 Y altas restricciones", d, pubY); res$c01_Y$datos <- d

## C02 muchos / pocos jovenes, parsimoniosa (208 = 128 extremas + 80 aleatorias)
D <- readRDS(file.path(OUT, "c02_intervalos.rds"))
cn <- c("YOUTH", "GDP", "ELEC", "DEC", "PROG", "YOUNG")
Xl <- as.matrix(D$L[, cn]); Xu <- as.matrix(D$U[, cn]); rownames(Xl) <- rownames(Xu) <- rownames(D$L); y <- D$Y
d <- interior_etiquetado(Xl, Xu, y, y, 0.8, "?", 80)
res$c02_Y <- resumir("C02 Y muchos jovenes (parsimoniosa)", d, c("~YOUTH*ELEC*PROG", "GDP*ELEC*~DEC*YOUNG", "GDP*ELEC*PROG*~YOUNG")); res$c02_Y$datos <- d
d <- interior_etiquetado(Xl, Xu, 1 - y, 1 - y, 0.8, "?", 80)
res$c02_nY <- resumir("C02 ~Y pocos jovenes (parsimoniosa)", d, c("~ELEC", "~GDP*~PROG*~YOUNG", "~YOUTH*DEC*~PROG", "GDP*DEC*YOUNG", "YOUTH*~DEC*~YOUNG"))
res$c02_nY$datos <- d
saveRDS(res, file.path(OUT, "frecuencias_por_tipo.rds"))
sink()
