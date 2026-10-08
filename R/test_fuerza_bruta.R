## Validacion: cotas exactas frente a fuerza bruta sobre una rejilla fina.
## Uso: Rscript test_fuerza_bruta.R  (desde la carpeta R/)
source("cotas_consistencia.R")
set.seed(20261006)

fuerza_bruta <- function(cond, signo, Yint, paso = 0.05) {
  n <- nrow(Yint)
  ## variables libres: cada celda con l < u se discretiza
  libres <- list()
  for (j in seq_along(cond)) for (i in seq_len(n))
    if (cond[[j]][i, 2] > cond[[j]][i, 1]) libres[[length(libres) + 1]] <- c(j, i)
  for (i in seq_len(n)) if (Yint[i, 2] > Yint[i, 1]) libres[[length(libres) + 1]] <- c(0, i)
  rej <- lapply(libres, function(v) {
    iv <- if (v[1] == 0) Yint[v[2], ] else cond[[v[1]]][v[2], ]
    unique(c(seq(iv[1], iv[2], by = paso), iv[2]))
  })
  grid <- if (length(rej)) as.matrix(expand.grid(rej)) else matrix(0, 1, 0)
  res <- matrix(NA, nrow(grid), 2)
  for (g in seq_len(nrow(grid))) {
    X <- lapply(cond, function(m) m[, 1]); Y <- Yint[, 1]
    for (k in seq_along(libres)) {
      v <- libres[[k]]
      if (v[1] == 0) Y[v[2]] <- grid[g, k] else X[[v[1]]][v[2]] <- grid[g, k]
    }
    T <- rep(1, n)
    for (j in seq_along(X)) T <- pmin(T, if (signo[j]) X[[j]] else 1 - X[[j]])
    res[g, ] <- c(sum(pmin(T, Y)) / sum(T), sum(pmin(T, Y)) / sum(Y))
  }
  c(suf_min = min(res[, 1]), suf_max = max(res[, 1]),
    cov_min = min(res[, 2]), cov_max = max(res[, 2]))
}

errores <- c()
for (rep in 1:40) {
  n <- 8; k <- 2
  cond <- lapply(1:k, function(j) { x <- round(runif(n), 2); cbind(x, x) })
  y <- round(runif(n), 2); Yint <- cbind(y, y)
  ## 2-3 faltantes: totales [0,1] o parciales
  for (m in 1:sample(2:3, 1)) {
    j <- sample(0:k, 1); i <- sample(n, 1)
    iv <- if (runif(1) < 0.5) c(0, 1) else sort(round(runif(2), 2))
    if (j == 0) Yint[i, ] <- iv else cond[[j]][i, ] <- iv
  }
  signo <- sample(c(TRUE, FALSE), k, replace = TRUE)
  ex <- unlist(cotas_ajuste(cond, signo, Yint)[1:4])
  fb <- fuerza_bruta(cond, signo, Yint, paso = 0.02)
  ## exacto debe contener a la rejilla y quedar muy cerca de ella
  errores <- rbind(errores, c(ex - fb))
}
cat("Max |exacto - rejilla| por parametro:\n")
print(round(apply(abs(errores), 2, max), 5))
cat("Violaciones (exacto mas estrecho que la rejilla):",
    sum(errores[, c(1, 3)] > 1e-9) + sum(errores[, c(2, 4)] < -1e-9), "\n")

## La completacion extrema debe reproducir la cota calculada
dif <- c()
for (rep in 1:30) {
  n <- 10; k <- 3
  cond <- lapply(1:k, function(j) { x <- runif(n); m <- sample(n, 2); l <- x; u <- x
    l[m] <- 0; u[m] <- 1; cbind(l, u) })
  y <- runif(n); Yint <- cbind(y, y); Yint[sample(n, 1), ] <- c(0.2, 0.9)
  signo <- sample(c(TRUE, FALSE), k, replace = TRUE)
  for (sen in c("min", "max")) {
    ce <- completacion_extrema(cond, signo, Yint, "suf", sen)
    T <- apply(sapply(1:k, function(j) if (signo[j]) ce$X[, j] else 1 - ce$X[, j]), 1, min)
    dif <- c(dif, abs(sum(pmin(T, ce$Y)) / sum(T) - ce$valor))
    ## y la completacion respeta los intervalos
    stopifnot(all(sapply(1:k, function(j) all(ce$X[, j] >= cond[[j]][, 1] - 1e-12 &
                                                  ce$X[, j] <= cond[[j]][, 2] + 1e-12))))
  }
}
cat("Completacion extrema: max |consistencia realizada - cota| =", signif(max(dif), 3), "\n")
