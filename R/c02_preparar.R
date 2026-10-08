## C02 Kurz & Ettensperger 2024 (European Political Science 23: 349-397;
## DOI 10.1057/s41304-023-00460-7). Datos: Harvard Dataverse doi:10.7910/DVN/CDZMPM.
## Reconstruye la calibracion de los autores (Script 2, agosto 2023) para TODOS los
## partidos del universo (127, nr.MPs >= 10), verifica que con los 105 analizados
## (drop_na(qca.decentral)) se reproduce exactamente la solucion publicada
## (Script 3, .html) y guarda las pertenencias como intervalos: los 22 partidos sin
## Decentralization (PPDB) reciben [0, 1] en DEC.
## Uso: Rscript c02_preparar.R (desde R/)
suppressPackageStartupMessages(library(QCA))
source("region_solucion.R")   # solucion_completacion
dir <- "../data/C02_kurz"
OUT <- "../outputs"
sink(file.path(OUT, "c02_preparar.txt"), split = TRUE)
d <- read.csv(file.path(dir, "c02_Youth_QCA_Data.csv"), stringsAsFactors = FALSE)
d <- d[d$nr.MPs >= 10, ]
rownames(d) <- paste(d$ISO3, d$party.abb, sep = "-")
cat("Universo (nr.MPs >= 10):", nrow(d), "partidos\n")

## ---- calibracion: copia literal de los anclajes y redondeos del Script 2 ----
cal <- function(x, th, dig) round(calibrate(x, type = "fuzzy", method = "direct", thresholds = th, logistic = TRUE), dig)
X <- data.frame(row.names = rownames(d))
X$YOUTH <- d$Youth.index                                            # qca.youthindex (sin calibrar)
X$GDP   <- cal(d$Gdp.per.cap, c(25000, 36000, 55000), 2)            # qca.gdpcap
X$ELEC  <- ifelse(d$Elec.system %in% c(9, 11), 1,                   # qca.elecsystem
           ifelse(d$Elec.system %in% c(1, 2, 3, 10, 12), 0, NA))
X$DEC   <- cal(d$Decentralization, c(0.1, 0.48, 0.9), 2)           # qca.decentral
X$PROG  <- cal(d$Manifesto.progcons, c(10, -8, -20), 3)            # qca.progressive
X$YOUNG <- cal(d$party.age, c(50, 19.9, 3), 3)                      # qca.youngparty
Y       <- cal(d$ARI35, c(0.05, 0.6, 0.95), 3)                      # qca.outcome.ARI35
cat("\nFaltantes por condicion / resultado en el universo:\n"); print(c(colSums(is.na(X)), Y = sum(is.na(Y))))

an <- !is.na(X$DEC)                                                 # drop_na(qca.decentral)
cat("\nAnalizados:", sum(an), " Eliminados:", sum(!an), "\n")
cat("Eliminados (solo les falta DEC; Inclusiveness tampoco esta, pero no es condicion):\n")
print(data.frame(partido = rownames(d)[!an], ARI35 = round(d$ARI35[!an], 3), Y = Y[!an],
                 X[!an, c("YOUTH", "GDP", "ELEC", "PROG", "YOUNG")]), row.names = FALSE)

## ---- verificacion de la replica ----
Xa <- as.matrix(X[an, ]); ya <- Y[an]
d_tt <- as.data.frame(Xa); d_tt$Y <- ya
tt <- truthTable(d_tt, outcome = "Y", conditions = colnames(Xa), incl.cut = 0.8, n.cut = 1, complete = TRUE)
obs <- tt$tt[tt$tt$n > 0, ]
cat("\nFilas observadas:", nrow(obs), "| OUT=1:", sum(obs$OUT == "1"), "| casos:", sum(obs$n), "\n")
## filas con OUT=1 publicadas en el .html del Script 3 (numeracion de QCA; orden YOUTH GDP ELEC DEC PROG YOUNG)
pub_out1 <- c(28, 16, 26, 15, 27, 11, 63, 31, 59)
pub_incl <- c(`28` = .917, `16` = .907, `26` = .892, `15` = .880, `27` = .860, `11` = .857, `63` = .829,
              `31` = .825, `59` = .822, `25` = .797, `64` = .790, `30` = .790, `61` = .662, `55` = .336)
cat("Filas OUT=1 identicas a las publicadas:", setequal(rownames(obs)[obs$OUT == "1"], pub_out1), "\n")
cat("Max |incl reproducida - publicada| en 14 filas de control:",
    max(abs(round(as.numeric(as.character(obs[names(pub_incl), "incl"])), 3) - pub_incl)), "\n")

## formulas publicadas (Script 3 .html), traducidas a las etiquetas cortas
publicadas <- list(
  pars_Y  = "~YOUTH*ELEC*PROG + GDP*ELEC*~DEC*YOUNG + GDP*ELEC*PROG*~YOUNG",
  pars_nY = "~ELEC + ~YOUTH*DEC*~PROG + YOUTH*~DEC*~YOUNG + ~GDP*~PROG*~YOUNG + GDP*DEC*YOUNG",
  cons_Y  = "~YOUTH*ELEC*PROG*~YOUNG + GDP*ELEC*PROG*~YOUNG + ~YOUTH*~GDP*ELEC*DEC*PROG + ~YOUTH*GDP*ELEC*~DEC*YOUNG",
  cons_nY = paste("GDP*~ELEC*~YOUNG + ~YOUTH*GDP*DEC*~PROG + ~YOUTH*~ELEC*~PROG*YOUNG + YOUTH*ELEC*~DEC*~YOUNG",
                  "+ ~ELEC*~DEC*PROG*~YOUNG + ~YOUTH*~GDP*~DEC*~PROG*~YOUNG + YOUTH*GDP*ELEC*DEC*PROG*YOUNG"))
canon <- function(s) paste(sort(strsplit(s, " \\+ ")[[1]]), collapse = " + ")
publicadas <- lapply(publicadas, canon)
rep <- list(
  pars_Y  = solucion_completacion(Xa, ya,     0.8, include = "?", row.dom = TRUE),
  pars_nY = solucion_completacion(Xa, 1 - ya, 0.8, include = "?", row.dom = TRUE),
  cons_Y  = solucion_completacion(Xa, ya,     0.8, include = "",  row.dom = TRUE),
  cons_nY = solucion_completacion(Xa, 1 - ya, 0.8, include = "",  row.dom = TRUE))
ok <- mapply(identical, rep, publicadas)
for (k in names(rep)) cat(sprintf("\n%-8s publicada : %s\n         reproducida: %s\n         IDENTICA: %s\n",
                                  k, publicadas[[k]], rep[[k]], ok[[k]]))
if (!all(ok)) { sink(); stop("La solucion publicada NO se reproduce: se detiene el reanalisis.") }

## parametros de ajuste de la solucion parsimoniosa Y (publicado: M1 inclS 0.746, covS 0.538)
m <- minimize(tt, include = "?", row.dom = TRUE); ic <- try(m$IC$sol.incl.cov, silent = TRUE)
cat("\nAjuste pars_Y reproducido: incl", try(round(ic[1, "inclS"], 3), silent = TRUE),
    "cov", try(round(ic[1, "covS"], 3), silent = TRUE), "(publicado 0.746 / 0.538)\n")

## ---- intervalos para TODO el universo: publicados primero, eliminados despues ----
ord <- c(which(an), which(!an))
L <- U <- X[ord, ]
L$DEC[is.na(L$DEC)] <- 0; U$DEC[is.na(U$DEC)] <- 1
saveRDS(list(L = L, U = U, Y = Y[ord], n_pub = sum(an), eliminados = rownames(d)[!an],
             publicadas = publicadas), file.path(OUT, "c02_intervalos.rds"))
cat("\nCeldas faltantes por variable entre los eliminados:\n")
print(colSums(U[-(1:sum(an)), ] - L[-(1:sum(an)), ] > 0))
sink()
