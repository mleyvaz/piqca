## Resumen de la variante con ruido 0.04 (v2 SMR) -> salidas/simulacion_v1_ruido004_resumen.txt
## Uso: Rscript resumen_simulacion_v1.R  (desde la carpeta R/)
R <- read.csv("../salidas/simulacion_v1_ruido004_replicas.csv", stringsAsFactors = FALSE)
sink("../salidas/simulacion_v1_ruido004_resumen.txt", split = TRUE)
cat("Replicas:", nrow(R), "\n\n")
pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
agg <- do.call(rbind, lapply(split(R, list(R$mecanismo, R$p), drop = TRUE), function(d) data.frame(
  mecanismo = d$mecanismo[1], p = d$p[1], n = nrow(d),
  N_completos = round(mean(d$n_completos), 1),
  full_con_falsos = pct(d$full_falsos > 0),
  elim_con_falsos = pct(d$del_falsos > 0),
  imput_con_falsos = pct(d$imp_falsos > 0),
  elim_distinta_full = pct(!d$del_igual_full),
  imput_distinta_full = pct(!d$imp_igual_full),
  cierta_vacia = pct(d$nucleo_vacio),
  cierta_con_falsos = pct(d$nucleo_falsos > 0),
  cierta_en_full = pct(d$nucleo_en_full),
  full_alcanzada = pct(d$full_en_interior),
  exterior_enumerada = pct(d$exterior_enumerada),
  elim_incompatible_certificada = pct(d$del_imposible),        # sobre las enumeradas
  n_incompatible = sum(d$del_imposible, na.rm = TRUE),
  n_enumeradas = sum(d$exterior_enumerada),
  soluciones_interior = round(mean(d$n_sol_interior), 1))))
agg <- agg[order(agg$mecanismo, agg$p), ]
rownames(agg) <- NULL
print(agg)
cat("\nDefiniciones: '_con_falsos' = % de replicas cuya solucion contiene al menos un termino que no",
    "implica el modelo verdadero A*B + C; 'distinta_full' = % en que la solucion difiere de la",
    "obtenida con datos completos; 'cierta_*' = solucion cierta estimada por la cota interior",
    "(K(R-), cota superior de la verdadera); 'elim_incompatible_certificada' = % de las replicas con",
    "cota exterior enumerada en que la solucion por eliminacion queda fuera de R+.\n")
## Totales para el texto
tot <- function(m) {
  d <- R[R$p >= 0.10, ]
  if (m != "todos") d <- d[d$mecanismo == m, ]
  c(n = nrow(d), cierta_con_falsos = pct(d$nucleo_falsos > 0),
    elim_con_falsos = pct(d$del_falsos > 0), imput_con_falsos = pct(d$imp_falsos > 0),
    full_con_falsos = pct(d$full_falsos > 0))
}
cat("\nAgregado p >= 0.10:\n"); print(sapply(c("todos", "MCAR", "MNAR_alto", "MNAR_bajo"), tot))
## Falsos de la solucion cierta en relacion con los de datos completos
cat("\nReplicas con terminos falsos en la cierta pero no en la solucion con datos completos:",
    sum(R$nucleo_falsos > 0 & R$full_falsos == 0), "de", nrow(R), "\n")
cat("Replicas en que la cierta contiene solo terminos de la solucion con datos completos:",
    sum(R$nucleo_en_full), "de", nrow(R), "\n")
sink()
