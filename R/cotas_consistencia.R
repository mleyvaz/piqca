## =============================================================================
## cotas_consistencia.R  --  Nivel 1: cotas exactas de los parametros de ajuste
## =============================================================================
## Cada pertenencia es un intervalo [l, u]; un dato completo tiene l = u y un
## faltante total es [0, 1]. Para un termino T (conjuncion de condiciones, con
## negaciones) y un resultado Y se calculan el minimo y el maximo exactos de
##   consistencia de suficiencia  S = sum(min(T,Y)) / sum(T)
##   consistencia de necesidad    N = sum(min(T,Y)) / sum(Y)   (= cobertura de T)
## sobre todas las completaciones de la caja de intervalos.
##
## Por que es exacto:
##  - min() es monotono, asi que el termino T_i tambien es un intervalo
##    [min_j l_ij, min_j u_ij] (con las negaciones invertidas: 1-u, 1-l).
##  - Para el cociente se usa Dinkelbach: min N/D <=> raiz de
##    g(lambda) = min_x sum_i [num_i(x_i) - lambda * den_i(x_i)], separable
##    por caso. En cada caso la funcion es lineal por tramos con un solo
##    quiebre (T = Y), asi que el optimo esta en un conjunto finito de
##    candidatos: extremos del intervalo y el punto de cruce.
## =============================================================================

## Intervalo de un termino. `cond` es una lista de matrices n x 2 (l, u) y
## `signo` es un vector logico: TRUE = condicion presente, FALSE = negada.
intervalo_termino <- function(cond, signo) {
  L <- U <- rep(1, nrow(cond[[1]]))
  for (j in seq_along(cond)) {
    lj <- if (signo[j]) cond[[j]][, 1] else 1 - cond[[j]][, 2]
    uj <- if (signo[j]) cond[[j]][, 2] else 1 - cond[[j]][, 1]
    L <- pmin(L, lj); U <- pmin(U, uj)
  }
  cbind(L, U)
}

## Candidatos por caso para el par (t, y) en la caja [tl,tu] x [yl,yu].
.candidatos <- function(tl, tu, yl, yu) {
  ts <- unique(c(tl, tu, min(max(yl, tl), tu), min(max(yu, tl), tu)))
  ys <- unique(c(yl, yu, min(max(tl, yl), yu), min(max(tu, yl), yu)))
  expand.grid(t = ts, y = ys)
}

## Cota (sentido = "min" o "max") del cociente sum(num)/sum(den) con Dinkelbach.
## tipo = "suf" (den = T) o "nec" (den = Y).
.cota_cociente <- function(Tint, Yint, tipo, sentido, tol = 1e-12, maxit = 200,
                          devolver = FALSE) {
  n <- nrow(Tint)
  cand <- lapply(seq_len(n), function(i)
    .candidatos(Tint[i, 1], Tint[i, 2], Yint[i, 1], Yint[i, 2]))
  num_f <- function(cc) pmin(cc$t, cc$y)
  den_f <- function(cc) if (tipo == "suf") cc$t else cc$y
  ## punto de partida: completacion por el punto medio
  sel <- vapply(cand, function(cc) 1L, integer(1))
  lam <- {
    nn <- sum(mapply(function(cc, k) num_f(cc)[k], cand, sel))
    dd <- sum(mapply(function(cc, k) den_f(cc)[k], cand, sel))
    if (dd <= 0) 0 else nn / dd
  }
  for (it in seq_len(maxit)) {
    g <- 0; nn <- 0; dd <- 0; tsel <- ysel <- numeric(n)
    for (i in seq_len(n)) {
      cc <- cand[[i]]; v <- num_f(cc) - lam * den_f(cc)
      k <- if (sentido == "min") which.min(v) else which.max(v)
      ## en empate se prefiere mayor denominador (evita 0/0)
      emp <- which(abs(v - v[k]) < 1e-15); k <- emp[which.max(den_f(cc)[emp])]
      g <- g + v[k]; nn <- nn + num_f(cc)[k]; dd <- dd + den_f(cc)[k]
      tsel[i] <- cc$t[k]; ysel[i] <- cc$y[k]
    }
    if (dd <= 0) return(NA_real_)
    nuevo <- nn / dd
    if (abs(g) < tol || abs(nuevo - lam) < tol)
      return(if (devolver) list(valor = nuevo, t = tsel, y = ysel) else nuevo)
    lam <- nuevo
  }
  warning("Dinkelbach no convergio"); lam
}

## Interfaz principal: devuelve las cotas de consistencia de suficiencia y de
## necesidad (= cobertura) del termino.
cotas_ajuste <- function(cond, signo, Yint) {
  Tint <- intervalo_termino(cond, signo)
  data.frame(
    suf_min = .cota_cociente(Tint, Yint, "suf", "min"),
    suf_max = .cota_cociente(Tint, Yint, "suf", "max"),
    cov_min = .cota_cociente(Tint, Yint, "nec", "min"),
    cov_max = .cota_cociente(Tint, Yint, "nec", "max"),
    n_casos_con_faltante = sum(Tint[, 2] > Tint[, 1] | Yint[, 2] > Yint[, 1]))
}

## Completacion que alcanza una cota: dado el valor t_i buscado para el termino,
## cada literal j se fija en max(a_j, min(t_i, b_j)); su minimo es exactamente
## t_i porque t_i >= min_j a_j. Devuelve la matriz de condiciones y el vector Y.
completacion_extrema <- function(cond, signo, Yint, tipo = "suf", sentido = "min") {
  Tint <- intervalo_termino(cond, signo)
  r <- .cota_cociente(Tint, Yint, tipo, sentido, devolver = TRUE)
  X <- sapply(seq_along(cond), function(j) {
    a <- if (signo[j]) cond[[j]][, 1] else 1 - cond[[j]][, 2]
    b <- if (signo[j]) cond[[j]][, 2] else 1 - cond[[j]][, 1]
    lit <- pmax(a, pmin(r$t, b))
    if (signo[j]) lit else 1 - lit
  })
  list(X = X, Y = r$y, valor = r$valor)
}
