## Figura 2 del artículo (v2): porcentaje de réplicas cuya solución contiene un término falso,
## por tratamiento y mecanismo, agrupando p = 0.10 y 0.20 (simulación v2, 600 réplicas).
## Legible en blanco y negro: tonos de gris + tramas. Uso: Rscript figura_simulacion.R (desde R/)
R <- read.csv("../salidas/simulacion_v2_replicas.csv", stringsAsFactors = FALSE)
s <- R[R$p >= 0.10, ]
mec <- c("MCAR", "MNAR_alto", "MNAR_bajo")
etq_mec <- c("MCAR", "MNAR-high", "MNAR-low")
pct <- function(x) 100 * mean(x)
M <- sapply(mec, function(m) {
  d <- s[s$mecanismo == m, ]
  c(Complete = pct(d$full_falsos > 0), Deletion = pct(d$del_falsos > 0),
    Imputation = pct(d$imp_falsos > 0), Sure = pct(d$nucleo_falsos > 0))
})
colnames(M) <- etq_mec
print(round(M, 1))
dir.create("../manuscrito/figures", showWarnings = FALSE)
dibujar <- function() {
par(mar = c(3.2, 4.2, 0.8, 0.5))
grises <- c("grey15", "grey45", "grey70", "white")
bp <- barplot(M, beside = TRUE, ylim = c(0, 112), col = grises, border = "black", las = 1, axes = FALSE,
              ylab = "Replicates with a false term (%)", names.arg = etq_mec)
barplot(M, beside = TRUE, ylim = c(0, 112), add = TRUE, col = "black",
        density = c(0, 0, 18, 0), angle = c(0, 0, 45, 0), border = NA, axes = FALSE, names.arg = rep("", 3))
text(bp, M + 3, sprintf("%.1f", M), cex = 0.72)
lg <- c("Complete data", "Case deletion", "Mean imputation", "Sure solution (inner bound)")
legend("topleft", ncol = 2, bty = "n", cex = 0.85, legend = lg, fill = grises, border = "black")
legend("topleft", ncol = 2, bty = "n", cex = 0.85, legend = lg, fill = "black", density = c(0, 0, 18, 0),
       angle = 45, border = NA, text.col = NA)
axis(2, at = seq(0, 100, 20), las = 1)
}
pdf("../manuscrito/figures/fig_simulation.pdf", width = 6.5, height = 3.6, pointsize = 10)
dibujar()
dev.off()
if (nzchar(Sys.getenv("PREVIEW_PNG"))) { png(Sys.getenv("PREVIEW_PNG"), width = 6.5, height = 3.6, units = "in", res = 130, pointsize = 10); dibujar(); dev.off() }
cat("escrito ../manuscrito/figures/fig_simulation.pdf\n")
