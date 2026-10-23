#!/usr/bin/env Rscript
# Reproducible README figures, drawn from the synthetic example reference.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: make_figures.R DEMO_OUTPUT ASSET_DIRECTORY")
suppressPackageStartupMessages(library(data.table))
library(grid)
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
ink <- "#142D42"; teal <- "#007F87"; gold <- "#DD9A32"; muted <- "#566E81"
txt <- function(label, x, y, size = 18, color = ink, face = "plain", just = "left") {
  grid.text(label, x, y, just = just, gp = gpar(fontsize = size, col = color, fontface = face, fontfamily = "sans"))
}
box <- function(x, y, w, h, fill = "white", border = NA) {
  grid.roundrect(x, y, w, h, r = unit(0.015, "npc"), gp = gpar(fill = fill, col = border))
}
png(file.path(args[2], "methmap-workflow.png"), width = 1800, height = 1100, res = 150, type = "cairo")
grid.newpage()
grid.rect(gp = gpar(fill = "#F1F6FA", col = NA))
grid.rect(.5, .845, 1, .31, gp = gpar(fill = ink, col = NA))
txt("METHMAP", .055, .91, 38, "white", "bold")
txt("De lecturas largas a referencias de metilación por tejido", .055, .84, 21, "#C8E0E8")
txt("hg38  /  5-metilcitosina  /  resolución de una base", .055, .775, 15, "#77D9D1")
xs <- c(.175, .5, .825)
titles <- c("01  ENTRADAS", "02  PROCESAMIENTO", "03  REFERENCIA")
subtitles <- c("BAM con MM/ML", "Modkit + filtro de cobertura", "Distribución entre muestras")
details <- c("FASTA hg38 + índice .fai\nMuestras agrupadas por tejido\nAlternativa: bedMethyl de Modkit",
             "5mC: código m\nCobertura válida ≥ 10\nUn valor por muestra y citosina",
             "Percentiles p5–p95\nNúmero de muestras + IQR\nRDS, TSV, bedGraph y BigWig")
for (i in 1:3) {
  box(xs[i], .505, .29, .285)
  txt(titles[i], xs[i]-.12, .61, 13, teal, "bold")
  txt(subtitles[i], xs[i]-.12, .555, 13, ink, "bold")
  txt(details[i], xs[i]-.12, .46, 12, muted)
}
for (i in 1:2) grid.lines(c(xs[i]+.147, xs[i+1]-.147), c(.51,.51),
  arrow = arrow(length = unit(2.5,"mm")), gp=gpar(col=teal,lwd=2))
box(.5, .225, .89, .19, "#DCEEEF")
txt("UNA REFERENCIA, DOS OPERACIONES", .08, .277, 13, teal, "bold")
txt("build", .08, .217, 21, ink, "bold")
txt("Construye una cohorte desde cero", .19, .217, 13, ink)
txt("add", .55, .217, 21, ink, "bold")
txt("Añade muestras en una nueva versión", .63, .217, 13, ink)
txt("Los datos ausentes se excluyen del cálculo; la metilación cero se conserva.", .08, .165, 12, muted)
txt("MethMap · 2026", .055, .055, 12, muted)
txt("Esquema del método · sin datos de individuos", .945, .055, 12, muted, just="right")
dev.off()

ref <- readRDS(file.path(args[1], "reference", "chr1.rds"))
sites <- as.data.table(ref$sites)[start >= 101 & start < 701]
q <- t(vapply(sites$values, quantile, numeric(5), probs=c(.05,.25,.5,.75,.95), names=FALSE, type=7))
matrix_values <- sapply(ref$sample_ids, function(id) vapply(sites$values,
    function(values) if (id %in% names(values)) unname(values[id]) else NA_real_, numeric(1)))
png(file.path(args[2], "methmap-example.png"), width=1800, height=1400, res=160, type="cairo")
layout(matrix(1:3, ncol=1), heights=c(2.6,2.0,1.15))
par(oma=c(1,1,4,1), mar=c(2,5,2,2), family="sans", col.axis=muted, col.lab=ink, fg=ink)
pos <- sites$start
plot(pos, q[,3], type="n", ylim=c(0,100), xlab="", ylab="Metilación (%)", xaxt="n", bty="l", yaxs="i")
abline(h=seq(0,100,25), col="#E6EDF3")
polygon(c(pos,rev(pos)), c(q[,1],rev(q[,5])), col="#CFE8EC", border=NA)
polygon(c(pos,rev(pos)), c(q[,2],rev(q[,4])), col="#7BC5CB", border=NA)
lines(pos, q[,3], col=teal, lwd=2.5)
points(pos, matrix_values[,1], col=adjustcolor(gold,.85), pch=16, cex=.45)
legend("topright", c("p5–p95","p25–p75","Mediana","Muestra 01"),
       col=c("#CFE8EC","#7BC5CB",teal,gold), lwd=c(8,8,2,NA), pch=c(NA,NA,NA,16),
       bg="white", box.col=NA, horiz=FALSE, cex=.82)
title("A  Perfil poblacional y una muestra individual", adj=0, cex.main=1, col.main=ink)
par(mar=c(2,5,2,2))
image(pos, seq_along(ref$sample_ids), matrix_values, zlim=c(0,100),
      col=colorRampPalette(c("#EDF6F8","#73BEC7","#07576C"))(101),
      xlab="", ylab="", xaxt="n", yaxt="n", useRaster=TRUE)
axis(2, seq_along(ref$sample_ids), labels=paste("Muestra",seq_along(ref$sample_ids)), las=1, tick=FALSE, cex.axis=.8)
title("B  Variación entre muestras", adj=0, cex.main=1, col.main=ink)
mtext("Claro: 0% · oscuro: 100% · blanco: dato ausente", side=3, line=.3, adj=1, cex=.75, col=muted)
par(mar=c(4,5,2,2))
plot(pos, lengths(sites$values), type="h", lwd=2, col=teal, ylim=c(0,length(ref$sample_ids)),
     xlab="Posición en chr1 de la referencia sintética (coordenada BED)", ylab="Muestras", bty="l", yaxs="i")
title("C  Número de muestras con cobertura suficiente por posición", adj=0, cex.main=1, col.main=ink)
mtext("MethMap | Ejemplo de una referencia por tejido", outer=TRUE, line=2, cex=1.3, font=2, col=ink)
mtext("DATOS SINTÉTICOS · seis muestras · ilustración del método, sin interpretación biológica", outer=TRUE, line=.7, cex=.85, col=muted)
dev.off()
