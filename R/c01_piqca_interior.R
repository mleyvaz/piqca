## PI-QCA sobre C01 (Ettensperger & Schleutker 2025). Uso: Rscript c01_piqca.R (desde R/)
`%||%` <- function(a, b) if (is.null(a)) b else a
source("cotas_consistencia.R"); source("region_solucion.R")
OUT <- "../outputs"; sink(file.path(OUT, "c01_piqca_interior.txt"), split = TRUE)
D <- readRDS(file.path(OUT, "c01_intervalos.rds"))
conds <- c("qca.relfrag", "qca.gdp", "qca.srel1900", "qca.fog", "qca.supp.relgrp", "qca.supp.urbanwk")
cn <- c("RELFRAG", "GDP", "SR1900", "FOG", "SUPPREL", "SUPPURB")
Xl <- as.matrix(D$L[, conds]); Xu <- as.matrix(D$U[, conds]); colnames(Xl) <- colnames(Xu) <- cn
rownames(Xl) <- rownames(Xu) <- D$L$name
y <- D$L$qca.outcome2.nxx
pub <- seq_len(D$n_pub)
cat("Casos: publicados", D$n_pub, "+ eliminados", length(D$eliminados), "=", nrow(Xl), "\n")

correr <- function(etq, Yl, Yu, incl.cut) {
  cat("\n=====================================================================\n", etq, "| corte", incl.cut, "\n")
  conv <- solucion_completacion(Xl[pub, ], Yl[pub], incl.cut, include = "?", row.dom = TRUE)
  cat("(1) Publicada / eliminacion de casos (N =", length(pub), "):", conv, "\n")
  E <- estado_filas(Xl, Xu, Yl, Yu, incl.cut)
  cat("(2) Filas: OUT=1 estable", sum(E$OUT == "1"), "| OUT=0 estable", sum(E$OUT == "0"),
      "| indeterminadas", sum(E$OUT == "ind"), "| solo posibles (abstencion)", sum(E$existencia == "posible"),
      "| remanentes", sum(E$existencia == "remanente"), "\n")
  print(E[E$OUT == "ind" | E$existencia == "posible", ], row.names = FALSE)
  ext <- try(stop("omitida en la corrida rapida"), silent = TRUE)
  int <- region_interior(Xl, Xu, Yl, Yu, incl.cut, n_sim = 400, row.dom = TRUE)
  cat("\n(3b) Region posible, cota interior (", length(int), "soluciones alcanzadas ):\n"); print(head(int, 15))
  if (!inherits(ext, "try-error")) {
    cat("\n(3a) Region posible, cota exterior:", length(ext), "soluciones\n")
    nuc <- nucleo(ext)
  } else { cat("\n(3a) exterior no calculable:", ext, "\n"); nuc <- nucleo(int); cat("   (nucleo sobre la interior: cota SUPERIOR de la solucion cierta)\n") }
  cat("\n(4) Solucion cierta:", if (length(nuc)) paste(nuc, collapse = " + ") else "VACIA", "\n")
  tc <- .terminos(conv)
  cat("    Terminos publicados en la solucion cierta:", sum(tc %in% nuc), "/", length(tc),
      "\n    Solucion publicada alcanzada por alguna completacion:", conv %in% names(int),
      "\n    Frecuencia de la publicada en la interior:", if (conv %in% names(int)) int[[conv]] else 0, "/", sum(int), "\n")
  ## estabilidad termino a termino
  todos <- unique(unlist(lapply(names(int), .terminos)))
  frec <- sapply(todos, function(t) sum(int[sapply(names(int), function(s) t %in% .terminos(s))]) / sum(int))
  cat("    Frecuencia de cada termino en la region interior:\n"); print(round(sort(frec, decreasing = TRUE), 3))
  invisible(list(conv = conv, E = E, int = int, ext = if (!inherits(ext, "try-error")) ext, nucleo = nuc))
}
r_pos <- correr("Y: altas restricciones (parsimoniosa)", y, y, 0.75)
r_neg <- correr("~Y: bajas restricciones (parsimoniosa)", 1 - y, 1 - y, 0.90)
saveRDS(list(pos = r_pos, neg = r_neg), file.path(OUT, "c01_piqca_interior.rds"))
sink()
