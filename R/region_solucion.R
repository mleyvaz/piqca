## =============================================================================
## region_solucion.R  --  Nivel 2: region de identificacion de la solucion QCA
## =============================================================================
## Entrada: pertenencias como intervalos. Xl, Xu son matrices n x k (condiciones)
## y Yl, Yu vectores (resultado). Dato completo: l = u. Faltante total: [0, 1].
##
## Paso A (exacto, fila a fila): para cada fila de la tabla de verdad se calculan
##   - las cotas exactas de su consistencia (Nivel 1);
##   - la poblacion cierta (casos con pertenencia > 0.5 en TODA completacion) y
##     la posible (en ALGUNA completacion).
##   Estado de la fila: OUT = 1 estable, OUT = 0 estable o indeterminado; y
##   existencia: observada, posible (abstencion) o remanente.
##
## Paso B (region EXTERIOR): se enumeran todas las combinaciones de las filas
##   indeterminadas y se minimiza cada una. Como el paso A trata las filas por
##   separado e ignora que comparten casos, esta region CONTIENE a la verdadera.
##
## Paso C (region INTERIOR): se minimizan completaciones concretas (vertices
##   aleatorios de la caja y puntos interiores). Toda solucion encontrada asi es
##   alcanzable, luego esta region esta CONTENIDA en la verdadera.
##
##   interior  ⊆  region verdadera  ⊆  exterior
##
## Nucleo identificado: implicantes presentes en TODAS las soluciones de la
## region exterior (si estan en todas las de la exterior, estan en todas las de
## la verdadera).
## =============================================================================

suppressPackageStartupMessages(library(QCA))
## requiere cotas_consistencia.R cargado antes

`%||%` <- function(a, b) if (is.null(a)) b else a

.filas <- function(k, nombres) {
  g <- expand.grid(rep(list(c(0, 1)), k)); names(g) <- nombres; g
}

## Paso A ---------------------------------------------------------------------
estado_filas <- function(Xl, Xu, Yl, Yu, incl.cut, n.cut = 1) {
  k <- ncol(Xl); nm <- colnames(Xl); G <- .filas(k, nm)
  cond <- lapply(seq_len(k), function(j) cbind(Xl[, j], Xu[, j]))
  Yint <- cbind(Yl, Yu)
  out <- lapply(seq_len(nrow(G)), function(r) {
    signo <- G[r, ] == 1
    Tint <- intervalo_termino(cond, signo)
    cb <- cotas_ajuste(cond, signo, Yint)
    data.frame(G[r, , drop = FALSE],
               n_cierto = sum(Tint[, 1] > 0.5), n_posible = sum(Tint[, 2] > 0.5),
               incl_min = cb$suf_min, incl_max = cb$suf_max)
  })
  E <- do.call(rbind, out)
  E$existencia <- ifelse(E$n_cierto >= n.cut, "observada",
                  ifelse(E$n_posible >= n.cut, "posible", "remanente"))
  E$OUT <- ifelse(E$existencia == "remanente", "?",
           ifelse(E$incl_min >= incl.cut, "1",
           ifelse(E$incl_max < incl.cut, "0", "ind")))
  E
}

## Minimizacion a partir de un vector de OUT por fila --------------------------
.minimizar_out <- function(G, out_vec, include = "?", row.dom = FALSE) {
  ## Se construye un truthTable de una sola fila por configuracion con datos
  ## ficticios y luego se sustituye OUT; asi se reutiliza QCA::minimize.
  if (!any(out_vec == "1")) return("(vacia)")
  dat <- G; dat$Y <- 0
  tt <- suppressWarnings(QCA::truthTable(dat, outcome = "Y",
          conditions = names(G), incl.cut = 0.5, complete = TRUE))
  clave_tt <- apply(tt$tt[, names(G), drop = FALSE], 1, paste, collapse = "")
  clave_G  <- apply(G, 1, paste, collapse = "")
  o <- out_vec[match(clave_tt, clave_G)]
  tt$tt$OUT <- ifelse(o == "?", "?", o)
  res <- try(QCA::minimize(tt, include = include, row.dom = row.dom), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res$solution)) return("(sin solucion)")
  paste(sort(vapply(res$solution, function(s) paste(sort(s), collapse = " + "),
                    character(1))), collapse = " || ")
}

.terminos <- function(sol) {
  if (sol %in% c("(vacia)", "(sin solucion)")) return(character(0))
  modelos <- strsplit(sol, " \\|\\| ")[[1]]
  ## un termino esta en la solucion si aparece en TODOS los modelos equivalentes
  Reduce(intersect, lapply(modelos, function(m) strsplit(m, " \\+ ")[[1]]))
}

## Paso B ---------------------------------------------------------------------
region_exterior <- function(E, nombres, include = "?", max_ind = 14, row.dom = FALSE) {
  ## row.dom DEBE coincidir con el de la region interior y el protocolo publicado;
  ## si no, la cota exterior no contiene a la region verdadera (Prop. 3).
  G <- E[, nombres]
  base <- E$OUT
  ## filas posibles (abstencion): pueden ser remanente o tener OUT 0/1
  idx_ind <- which(base == "ind" | E$existencia == "posible")
  if (length(idx_ind) > max_ind) stop("Demasiadas filas indeterminadas: ", length(idx_ind))
  opciones <- lapply(idx_ind, function(i) {
    o <- if (base[i] == "ind") c("0", "1") else base[i]
    if (E$existencia[i] == "posible") unique(c(o, "?")) else o
  })
  combos <- expand.grid(opciones, stringsAsFactors = FALSE)
  sols <- character(nrow(combos))
  for (c in seq_len(nrow(combos))) {
    v <- base; v[idx_ind] <- unlist(combos[c, ])
    v[v == "ind"] <- "0"   # no deberia quedar ninguno
    sols[c] <- .minimizar_out(G, v, include, row.dom)
  }
  sort(table(sols), decreasing = TRUE)
}

## Paso C ---------------------------------------------------------------------
completar <- function(l, u, modo) {
  if (modo == "vertice") ifelse(runif(length(l)) < 0.5, l, u)
  else runif(length(l), l, u)
}

solucion_completacion <- function(X, Y, incl.cut, n.cut = 1, include = "?", row.dom = FALSE) {
  d <- as.data.frame(X); d$Y <- Y
  tt <- try(suppressWarnings(QCA::truthTable(d, outcome = "Y", conditions = colnames(X),
            incl.cut = incl.cut, n.cut = n.cut, complete = TRUE)), silent = TRUE)
  if (inherits(tt, "try-error")) return("(sin solucion)")
  if (!any(tt$tt$OUT == "1")) return("(vacia)")
  res <- try(QCA::minimize(tt, include = include, row.dom = row.dom), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res$solution)) return("(sin solucion)")
  paste(sort(vapply(res$solution, function(s) paste(sort(s), collapse = " + "),
                    character(1))), collapse = " || ")
}

region_interior <- function(Xl, Xu, Yl, Yu, incl.cut, n.cut = 1, include = "?",
                            n_sim = 500, semilla = 20261006, dirigidas = TRUE,
                            row.dom = FALSE) {
  set.seed(semilla)
  sols <- character(n_sim)
  ## Completaciones dirigidas: para cada fila, las que llevan su consistencia
  ## al minimo y al maximo exactos (Nivel 1). Son las que mas mueven la tabla.
  extra <- character(0)
  if (dirigidas) {
    k <- ncol(Xl); G <- .filas(k, colnames(Xl))
    cond <- lapply(seq_len(k), function(j) cbind(Xl[, j], Xu[, j]))
    for (r in seq_len(nrow(G))) for (sen in c("min", "max")) {
      ce <- completacion_extrema(cond, G[r, ] == 1, cbind(Yl, Yu), "suf", sen)
      X <- ce$X; dimnames(X) <- dimnames(Xl)
      extra <- c(extra, solucion_completacion(X, ce$Y, incl.cut, n.cut, include, row.dom))
    }
  }
  for (s in seq_len(n_sim)) {
    modo <- if (s %% 2 == 0) "vertice" else "interior"
    X <- matrix(completar(Xl, Xu, modo), nrow(Xl), dimnames = dimnames(Xl))
    Y <- completar(Yl, Yu, modo)
    sols[s] <- solucion_completacion(X, Y, incl.cut, n.cut, include, row.dom)
  }
  sort(table(c(sols, extra)), decreasing = TRUE)
}

## Resumen --------------------------------------------------------------------
nucleo <- function(tabla_sols) {
  terms <- lapply(names(tabla_sols), .terminos)
  if (any(lengths(terms) == 0)) return(character(0))
  Reduce(intersect, terms)
}
