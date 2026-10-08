## C01 Ettensperger & Schleutker 2025: reconstruye la calibracion y recupera los
## casos eliminados por faltantes como intervalos. Uso: Rscript c01_preparar.R (desde R/)
suppressPackageStartupMessages(library(QCA))
dir <- "../data/C01_ettensperger"
F <- read.csv(file.path(dir, "Data_QCA_full.csv"))
P <- read.csv(file.path(dir, "Data_QCA2_nxx.csv"))
raw <- c(qca.relfrag = "al_religion2000", qca.gdp = "mad_gdppc", qca.srel1900 = "X1900",
         qca.fog = "fh_fog", qca.supp.relgrp = "vdem_20ya_regsup_g7",
         qca.supp.urbanwk = "vdem_20ya_regsup_g9", qca.outcome2.nxx = "NXX2014X")
cat("Full:", nrow(F), " Publicado:", nrow(P), "\n")
## infiere los anclajes (e, c, i) de la calibracion directa logistica con los 155
anclas <- function(x, y) {
  f <- function(p) { a <- sort(p); z <- try(calibrate(x, type = "fuzzy", method = "direct",
         thresholds = a, logistic = TRUE), silent = TRUE)
         if (inherits(z, "try-error") || any(!is.finite(z))) return(1e9); sum((round(z, 3) - y)^2) }
  q <- quantile(x, c(.1, .5, .9), na.rm = TRUE)
  o <- optim(q, f, control = list(maxit = 4000)); list(a = sort(o$par), sse = o$value)
}
A <- list()
for (v in names(raw)) {
  if (v == "qca.srel1900") next
  r <- anclas(P[[raw[v]]], P[[v]]); A[[v]] <- r$a
  cat(sprintf("%-18s anclas %s  max|err| = %.4f\n", v, paste(signif(r$a, 4), collapse = ", "),
      max(abs(round(calibrate(P[[raw[v]]], type = "fuzzy", method = "direct", thresholds = r$a, logistic = TRUE), 3) - P[[v]]))))
}
## casos eliminados
pub <- P$name_vdem
E <- F[!F$name_vdem %in% pub, ]
cat("\nEliminados:", nrow(E), "\n"); print(E[, c("name_vdem", raw)], row.names = FALSE)
## intervalos para TODOS los casos (publicados + eliminados)
U <- rbind(P[, c("name_vdem", raw)], E[, c("name_vdem", raw)])
L <- Uu <- data.frame(name = U$name_vdem)
for (v in names(raw)) {
  x <- U[[raw[v]]]
  m <- if (v == "qca.srel1900") x else round(calibrate(x, type = "fuzzy", method = "direct", thresholds = A[[v]], logistic = TRUE), 3)
  ## en los publicados se usa el valor publicado tal cual
  m[seq_len(nrow(P))] <- P[[v]]
  L[[v]] <- ifelse(is.na(m), 0, m); Uu[[v]] <- ifelse(is.na(m), 1, m)
}
saveRDS(list(L = L, U = Uu, n_pub = nrow(P), eliminados = E$name_vdem), "../outputs/c01_intervalos.rds")
cat("\nCeldas faltantes por variable entre los eliminados:\n"); print(colSums(Uu[-(1:nrow(P)), -1] - L[-(1:nrow(P)), -1] > 0))
