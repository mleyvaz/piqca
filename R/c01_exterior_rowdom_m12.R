## C01: cota exterior con row.dom = TRUE (mismo protocolo que la interior y que los autores).
## Corrige c01_piqca.R, cuya exterior se calculo sin row.dom. Uso: Rscript c01_exterior_rowdom.R
`%||%` <- function(a, b) if (is.null(a)) b else a
source("cotas_consistencia.R"); source("region_solucion.R")
x <- readRDS("../outputs/c01_piqca.rds")
cn <- c("RELFRAG", "GDP", "SR1900", "FOG", "SUPPREL", "SUPPURB")
sink("../outputs/c01_exterior_rowdom_m12.txt", split = TRUE)
for (r in c("neg", "pos")) {
  E <- x[[r]]$E
  m <- sum(E$OUT == "ind" | E$existencia == "posible")
  cat("\n==== ", r, " | filas indeterminadas o posibles:", m, "\n")
  if (m > 12) { cat("exterior con row.dom no enumerada (", m, "> 12 ); se reporta solo la interior\n"); next }
  ext <- region_exterior(E, cn, row.dom = TRUE, max_ind = 12)
  cat("Soluciones en la exterior (row.dom = TRUE):", length(ext), "\n")
  nE <- nucleo(ext); nI <- nucleo(x[[r]]$int)
  cat("Nucleo exterior K(R+):", paste(nE, collapse = " + "), "\n")
  cat("Nucleo interior K(R-):", paste(nI, collapse = " + "), "\n")
  cat("Cotas iguales (solucion cierta exacta):", setequal(nE, nI), "\n")
  cat("Publicada en la exterior:", x[[r]]$conv %in% names(ext), "\n")
  cat("Interior contenida en exterior:", all(names(x[[r]]$int) %in% names(ext)), "\n")
  tc <- .terminos(x[[r]]$conv); cat("Terminos publicados en K(R+):", sum(tc %in% nE), "/", length(tc), "\n")
}
sink()
