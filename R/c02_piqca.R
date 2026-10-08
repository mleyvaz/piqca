## PI-QCA sobre C02 (Kurz & Ettensperger 2024, EPS 23: 349-397). Uso: Rscript c02_piqca.R (desde R/)
## Requiere ../outputs/c02_intervalos.rds (c02_preparar.R, que verifica la replica exacta).
## Faltante principal (peor caso): DEC en [0, 1] para los 22 partidos sin dato PPDB.
## Opciones de los autores: incl.cut = 0.8, n.cut = 1, row.dom = TRUE;
## parsimoniosa include = "1, ?"; conservadora include = "1". La negacion se calcula
## con 1 - Y y el mismo incl.cut (asi la recalcula minimize() con outcome = "~...").
`%||%` <- function(a, b) if (is.null(a)) b else a
source("cotas_consistencia.R"); source("region_solucion.R")
OUT <- "../outputs"; sink(file.path(OUT, "c02_piqca.txt"), split = TRUE)
t0 <- Sys.time()
D <- readRDS(file.path(OUT, "c02_intervalos.rds"))
cn <- c("YOUTH", "GDP", "ELEC", "DEC", "PROG", "YOUNG")
Xl <- as.matrix(D$L[, cn]); Xu <- as.matrix(D$U[, cn])
rownames(Xl) <- rownames(Xu) <- rownames(D$L)
y <- D$Y
pub <- seq_len(D$n_pub)
cat("Casos: publicados", D$n_pub, "+ eliminados", length(D$eliminados), "=", nrow(Xl),
    "| faltante: DEC en [0,1] para los eliminados\n")

## Variante local de region_exterior() con row.dom = TRUE (opcion de los autores).
## Copia la logica de region_solucion.R sin modificar ese archivo.
.minimizar_out_rd <- function(G, out_vec, include = "?") {
  if (!any(out_vec == "1")) return("(vacia)")
  dat <- G; dat$Y <- 0
  tt <- suppressWarnings(QCA::truthTable(dat, outcome = "Y", conditions = names(G), incl.cut = 0.5, complete = TRUE))
  o <- out_vec[match(apply(tt$tt[, names(G), drop = FALSE], 1, paste, collapse = ""), apply(G, 1, paste, collapse = ""))]
  tt$tt$OUT <- ifelse(o == "?", "?", o)
  res <- try(QCA::minimize(tt, include = include, row.dom = TRUE), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res$solution)) return("(sin solucion)")
  paste(sort(vapply(res$solution, function(s) paste(sort(s), collapse = " + "), character(1))), collapse = " || ")
}
region_exterior_rd <- function(E, nombres, include = "?", max_ind = 12) {
  G <- E[, nombres]; base <- E$OUT
  idx <- which(base == "ind" | E$existencia == "posible")
  if (length(idx) > max_ind) stop("Demasiadas filas indeterminadas: ", length(idx))
  opciones <- lapply(idx, function(i) {
    o <- if (base[i] == "ind") c("0", "1") else base[i]
    if (E$existencia[i] == "posible") unique(c(o, "?")) else o })
  combos <- expand.grid(opciones, stringsAsFactors = FALSE)
  sols <- vapply(seq_len(nrow(combos)), function(c) { v <- base; v[idx] <- unlist(combos[c, ])
    .minimizar_out_rd(G, v, include) }, character(1))
  sort(table(sols), decreasing = TRUE)
}

correr <- function(etq, Yl, Yu, incl.cut, include, publicada) {
  cat("\n=====================================================================\n", etq, "| corte", incl.cut,
      "| include =", dQuote(include, FALSE), "| row.dom = TRUE\n")
  conv <- solucion_completacion(Xl[pub, ], Yl[pub], incl.cut, include = include, row.dom = TRUE)
  cat("(1) Publicada / eliminacion de casos (N =", length(pub), "):", conv,
      "\n    identica a la formula del articulo:", identical(conv, publicada), "\n")
  E <- estado_filas(Xl, Xu, Yl, Yu, incl.cut)
  cat("(2) Filas: OUT=1 estable", sum(E$OUT == "1"), "| OUT=0 estable", sum(E$OUT == "0"),
      "| indeterminadas", sum(E$OUT == "ind"), "| solo posibles (abstencion)", sum(E$existencia == "posible"),
      "| remanentes", sum(E$existencia == "remanente"), "\n")
  print(E[E$OUT == "ind" | E$existencia == "posible", ], row.names = FALSE)
  t1 <- Sys.time()
  int <- region_interior(Xl, Xu, Yl, Yu, incl.cut, include = include, n_sim = 80, row.dom = TRUE)
  cat("\n(3b) Region posible, cota interior (", length(int), "soluciones alcanzadas en", sum(int), "completaciones ):\n")
  print(head(int, 15))
  ext <- try(region_exterior_rd(E, cn, include = include, max_ind = 12), silent = TRUE)
  ext0 <- try(region_exterior(E, cn, include = include, max_ind = 12), silent = TRUE)
  if (!inherits(ext, "try-error")) {
    cat("\n(3a) Region posible, cota exterior (row.dom = TRUE):", length(ext), "soluciones\n")
    if (!inherits(ext0, "try-error")) cat("     [control: region_exterior() de region_solucion.R, sin row.dom:", length(ext0), "soluciones]\n")
    cat("     Interior contenida en la exterior:", all(names(int) %in% names(ext)), "\n")
    nuc <- nucleo(ext)
    fuera <- !(conv %in% names(ext))
  } else {
    cat("\n(3a) exterior no calculable:", ext, "\n"); nuc <- nucleo(int)
    cat("   (nucleo sobre la interior: cota SUPERIOR de la solucion cierta)\n"); fuera <- NA
  }
  cat("\n(4) Solucion cierta:", if (length(nuc)) paste(nuc, collapse = " + ") else "VACIA", "\n")
  if (!inherits(ext, "try-error")) cat("    Nucleo de la interior (cota superior):",
      if (length(nucleo(int))) paste(nucleo(int), collapse = " + ") else "VACIO",
      "| coincide con el de la exterior:", setequal(nucleo(int), nuc), "\n")
  tc <- .terminos(conv)
  cat("    Terminos publicados en la solucion cierta:", sum(tc %in% nuc), "/", length(tc),
      "\n    Solucion publicada alcanzada por alguna completacion:", conv %in% names(int),
      "\n    Frecuencia de la publicada en la interior:", if (conv %in% names(int)) int[[conv]] else 0, "/", sum(int),
      "\n    Publicada FUERA de la exterior (incompatibilidad certificada):", fuera, "\n")
  todos <- unique(unlist(lapply(names(int), .terminos)))
  frec <- sapply(todos, function(t) sum(int[sapply(names(int), function(s) t %in% .terminos(s))]) / sum(int))
  cat("    Frecuencia de cada termino PUBLICADO en la region interior:\n")
  fp <- setNames(sapply(tc, function(t) if (t %in% names(frec)) frec[[t]] else 0), tc); print(round(fp, 3))
  cat("    Frecuencia de cada termino en la region interior:\n"); print(round(sort(frec, decreasing = TRUE), 3))
  cat("    [tiempo:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min]\n")
  invisible(list(conv = conv, E = E, int = int, ext = if (!inherits(ext, "try-error")) ext,
                 ext_sin_rowdom = if (!inherits(ext0, "try-error")) ext0, nucleo = nuc, frec = frec, fuera = fuera))
}
P <- D$publicadas
r_pos  <- correr("Y: alta representacion joven ARI35 (parsimoniosa)", y, y, 0.8, "?", P$pars_Y)
r_neg  <- correr("~Y: baja representacion joven ARI35 (parsimoniosa)", 1 - y, 1 - y, 0.8, "?", P$pars_nY)
c_pos  <- correr("Y: alta representacion joven ARI35 (conservadora, secundaria)", y, y, 0.8, "", P$cons_Y)
c_neg  <- correr("~Y: baja representacion joven ARI35 (conservadora, secundaria)", 1 - y, 1 - y, 0.8, "", P$cons_nY)
saveRDS(list(pos = r_pos, neg = r_neg, cons_pos = c_pos, cons_neg = c_neg), file.path(OUT, "c02_piqca.rds"))
cat("\nTiempo total:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
sink()
