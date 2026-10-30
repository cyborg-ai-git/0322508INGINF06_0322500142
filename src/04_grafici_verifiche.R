# GRAFICI E CONTROLLI - grafica base R, pheatmap e ggplot2.
library(pheatmap)
library(ggplot2)
library(jsonlite)
library(digest)
O <- "risultati"
F <- "figure"
dir.create(F, showWarnings = FALSE)
leggi <- function(nome) read.delim(file.path(O, nome), check.names = FALSE)
meta <- leggi("coorte_pazienti.txt")
res <- leggi("tutti_geni_testati.txt")
deg <- res[res$DEG, ]
G <- as.matrix(read.delim(gzfile(file.path(O,
  "espressione_geni_log2_normalizzata.txt.gz")),
  row.names = 1, check.names = FALSE))
iq <- leggi("filtro_IQR.txt")
riepilogo <- leggi("riepilogo_analisi.txt")
soglia <- riepilogo$value[riepilogo$metric == "IQR_threshold"]
n <- nrow(meta)
N <- "#287C8E"
T <- "#C45C36"
grigio <- "#B7BEC5"
# Ogni figura viene salvata anche in PDF vettoriale.
# Con Rscript non si apre alcuna finestra. Per disabilitare l'anteprima:
# options(elaborato.mostra_plots = FALSE), prima di eseguire il main.
mostra_plots <- interactive() && Sys.getenv("RSTUDIO") == "1" &&
  isTRUE(getOption("elaborato.mostra_plots", TRUE))
stile <- function() par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1),
  family = "sans", las = 1, bty = "l", cex = 0.94,
  col.axis = "#273544", col.lab = "#273544")
salva <- function(nome, larghezza, altezza, disegna) {
  esporta <- function(formato) {
    file <- file.path(F, paste0(nome, ".", formato))
    if (formato == "png") {
      # Quartz su Mac; Cairo o dispositivo Windows sugli altri sistemi.
      tipo <- if (capabilities("aqua")) "quartz" else
        if (capabilities("cairo")) "cairo" else
          if (.Platform$OS.type == "windows") "windows" else
            getOption("bitmapType", "Xlib")
      png(file, width = larghezza, height = altezza, units = "in",
        res = 300, type = tipo)
    } else pdf(file, width = larghezza, height = altezza)
    # Chiudere solo il file appena aperto, anche in caso di errore.
    dispositivo <- dev.cur()
    on.exit(dev.off(dispositivo))
    stile()
    disegna()
  }
  for (formato in c("png", "pdf")) esporta(formato)
  if (mostra_plots) tryCatch({
    dispositivi <- dev.list()
    studio <- dispositivi[names(dispositivi) == "RStudioGD"]
    if (length(studio)) dev.set(studio[1]) else dev.new()
    # Ripristinare il layout: il grafico precedente puo avere piu pannelli.
    stile()
    disegna()
    # Lasciare aperto RStudioGD per conservare grafico e cronologia.
  }, error = function(e) message("Anteprima di ", nome,
    " non disponibile: ", conditionMessage(e),
    ". Allargare Plots e rieseguire. PNG e PDF sono gia salvati."))
}
# 1. Caratteristiche cliniche: un paziente e contato una sola volta.
salva("01_coorte", 8.2, 5.6, function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  hist(as.numeric(meta$age), breaks = seq(35, 95, 10), col = N,
    border = "white", main = "A  Eta", xlab = "Anni", ylab = "Pazienti")
  barre <- function(v, livelli, etichette, titolo) {
    cnt <- table(factor(v, levels = livelli))
    pos <- barplot(cnt, names.arg = etichette, col = N, border = NA,
      ylim = c(0, max(cnt) * 1.2), ylab = "Pazienti", main = titolo)
    text(pos, cnt, labels = cnt, pos = 3, cex = 0.9)
  }
  barre(meta$gender, c("Female", "Male"), c("Donne", "Uomini"),
    "B  Sesso")
  barre(meta$smoking_status, c("Never", "Current"), c("Mai fumatori",
    "Fumatori attuali"), "C  Abitudine al fumo")
  stadi <- c("IA", "IB", "IIA", "IIB", "IIIA", "IV")
  barre(meta$stage, paste("Stage", stadi), stadi, "D  Stadio tumorale")
})
# 2. Quantili salvati prima/dopo la normalizzazione: baffi min e max.
salva("02_normalizzazione", 8.2, 5, function() {
  par(mfrow = c(2, 1), mar = c(3, 4, 2.5, 1))
  for (fase in c("prima", "dopo")) {
    q <- leggi(paste0("quantili_", fase, ".txt"))
    z <- list(stats = t(as.matrix(q[, 3:7])), n = rep(48803, 2 *
      n), conf = matrix(NA_real_, 2, 2 * n), out = numeric(),
      group = numeric(), names = q$sample)
    bxp(z, show.names = FALSE, outline = FALSE, ylim = c(5, 17),
      boxfill = ifelse(q$group == "N", N, T), border = "gray40",
      boxlwd = 0.25, whiskcol = "gray40", staplecol = "gray40",
      whisklty = 1, whisklwd = 0.4, staplelwd = 0.4, medcol = "white",
      medlwd = 0.5, ylab = "Segnale log2")
    axis(1, at = c(1, 29, 57, 86, 114), labels = c("N1", "N29",
      "N57", "T29", "T57"))
    title(if (fase == "prima")
      "A  Dopo log2(x + 1)" else "B  Dopo normalizzazione quantile")
    abline(v = n + 0.5, col = grigio)
  }
})
# 3. IQR senza etichette di gruppo: distribuzione e dettaglio iniziale.
salva("03_filtro_IQR", 8.2, 3.8, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (zoom in c(FALSE, TRUE)) {
    valori <- if (zoom)
      iq$IQR[iq$IQR <= 0.5] else iq$IQR
    hist(valori, breaks = if (zoom)
      seq(0, 0.5, 0.01) else 75, col = N, border = "white",
      xlab = "IQR per gene",
      ylab = "Geni", main = if (zoom)
        "B  Dettaglio degli IQR bassi" else "A  Distribuzione completa")
    abline(v = soglia, col = T, lty = 2, lwd = 2)
    legend("topright", sprintf("10 percentile = %.4f", soglia),
      col = T, lty = 2, bty = "n", cex = 0.8)
  }
})
# 4. Volcano: tutti i geni testati; colore solo se entrambi i criteri valgono.
salva("04_volcano", 8.2, 5.3, function() {
  par(mar = c(4.5, 4.5, 1, 1))
  colore <- ifelse(res$DEG, ifelse(res$logFC >
    0, T, N), grigio)
  plot(res$logFC, -log10(res$p_adj),
    pch = 16, cex = 0.45, col = adjustcolor(colore,
      alpha.f = 0.7), xlim = c(-5.6,
      4.4), ylim = c(0, 42),
    xlab = "log2 fold change: tumorale rispetto a non tumorale",
    ylab = "-log10 p aggiustato BH")
  abline(v = c(-1, 1), h = -log10(0.05),
    col = "gray40", lty = 2)
  geni <- c("GCNT3", "LINC02538",
    "FABP4", "TOP2A", "SPP1")
  righe <- match(geni, res$gene_symbol)
  text(res$logFC[righe], -log10(res$p_adj[righe]),
    geni, pos = c(4, 4, 4, 4, 4),
    cex = 0.8, offset = 0.5)
  legend("topright", c("Non selezionati: 17189",
    "Down: 596", "Up: 323"), col = c(grigio,
    N, T), pch = 16, bty = "n",
    cex = 0.9)
})
# 5. Boxplot dei geni estremi per logFC, con tutte le 57 coppie.
# notch = TRUE aggiunge gli intagli, come nell'esercizio R sui boxplot.
# I limiti sono mediana +/- 1.58 * IQR / sqrt(n): un confronto descrittivo
# delle mediane, che non sostituisce il test t appaiato e la correzione BH.
# Usare TRUE, non T: in questo script T contiene il colore del tumore.
scelti <- deg$gene_symbol[c(which.max(deg$logFC), which.min(deg$logFC))]
salva("05_boxplot_estremi", 8.2, 4.8, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
  set.seed(20260916)  # Sposta solo i punti in orizzontale, non i dati.
  for (gene in scelti) {
    v <- G[gene, ]
    rr <- res[res$gene_symbol == gene, ]
    boxplot(v[1:n], v[(n + 1):(2 * n)], col = adjustcolor(c(N,
      T), 0.35), names = c("Non tumorale", "Tumorale"), notch = TRUE,
      outline = FALSE,
      ylim = range(v), ylab = "Espressione log2 normalizzata")
    segments(rep(1, n), v[1:n], rep(2, n), v[(n + 1):(2 * n)],
      col = adjustcolor(grigio, 0.5), lwd = 0.5)
    points(jitter(rep(1:2, each = n), amount = 0.04), v, pch = 16,
      cex = 0.5, col = rep(c(N, T), each = n))
    title(sprintf("%s\nlogFC = %+.3f   p BH = %.2e", gene, rr$logFC,
      rr$p_adj), cex.main = 0.95)
  }
})
# 6. pheatmap: schema del caso studio, annotazioni su righe e colonne.
# Z-score di riga prima del clustering, correlazione Pearson e average.
# La saturazione +/-2.5 riguarda esclusivamente la scala cromatica.
H <- G[deg$gene_symbol, , drop = FALSE]
Z <- t(scale(t(H)))
albero_righe <- hclust(as.dist(1 - cor(t(Z))), method = "average")
albero_colonne <- hclust(as.dist(1 - cor(Z)), method = "average")
campioni <- leggi("metadati_tutti.txt")
campioni <- campioni[match(colnames(Z), campioni$raw_id), ]
annotation_col <- data.frame(Tessuto = campioni$group,
  Batch = campioni$batch, row.names = colnames(Z))
annotation_row <- data.frame(Direzione = deg$direction,
  row.names = deg$gene_symbol)
colori <- list(Tessuto = c(N = N, T = T), Direzione = c(up = T, down = N),
  Batch = c(`Batch 1` = "#A0A7C9", `Batch 2` = "#D8CCAA"))
salva("06_heatmap_completa", 8.2, 6.7, function() {
  # pheatmap usa un dispositivo temporaneo: ripristinare quello giusto
  # prima di disegnare, anche quando RStudioGD e un file sono aperti.
  dispositivo <- dev.cur()
  mappa <- pheatmap(Z, scale = "none", border_color = NA, silent = TRUE,
    cluster_rows = albero_righe, cluster_cols = albero_colonne,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = colori, show_rownames = FALSE,
    show_colnames = FALSE, color = colorRampPalette(c("#2166AC",
      "#F7F7F7", "#B2182B"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    legend_breaks = c(-2.5, 0, 2.5),
    legend_labels = c("-2.5", "0", "+2.5"),
    fontsize = 10, main = "919 geni x 114 campioni | z-score per gene")
  dev.set(dispositivo)
  grid::grid.newpage()
  grid::grid.draw(mappa$gtable)
})
clusters <- cutree(albero_colonne, k = 2)
write.table(data.frame(sample = colnames(Z), cluster = clusters,
  group = campioni$group), file.path(O, "cluster_campioni_heatmap.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = colnames(Z)[albero_colonne$order]),
  file.path(O, "ordine_campioni_heatmap.txt"), sep = "\t", quote = FALSE,
  row.names = FALSE)
# 7. PCA: prcomp e stato calcolato nello script 02, come nell'esercitazione.
salva("07_PCA", 8.2, 4.3, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (tipo in c("filtrati", "DEG")) {
    pc <- leggi(paste0("PCA_", tipo, "_scores.txt"))
    vr <- leggi(paste0("PCA_", tipo, "_varianza.txt"))$variance_ratio
    plot(pc$PC1, pc$PC2, type = "n", xlab = sprintf("PC1  %.2f%%", vr[1] *
      100), ylab = sprintf("PC2  %.2f%%", vr[2] * 100), main = if (tipo ==
      "filtrati")
      "A  18.108 geni filtrati" else "B  919 DEG")
    abline(h = 0, v = 0, col = "gray90")
    segments(pc$PC1[1:n], pc$PC2[1:n], pc$PC1[(n + 1):(2 * n)], pc$PC2[(n +
      1):(2 * n)], col = grigio, lwd = 0.5)
    points(pc$PC1, pc$PC2, pch = 16, cex = 0.65, col = ifelse(pc$group ==
      "N", N, T))
    if (tipo == "filtrati")
      legend("topleft", c("Non tumorale", "Tumorale"), col = c(N, T),
        pch = 16, bty = "n", cex = 0.75)
  }
})
# 8. Dot plot ggplot2, come Hands on enrichR: x = Gene_count,
# dimensione = Gene_ratio (overlap/totale categoria), colore = p aggiustato.
salva("08_arricchimento", 8.2, 7, function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2,
    1)))
  dbs <- c("GO_Biological_Process_2025", "KEGG_2021_Human")
  for (i in seq_along(dbs)) {
    annotation <- leggi(file.path("Enrichr", paste0("tutti_",
      dbs[i], ".txt")))
    top_term <- 7
    annotation_top <- head(annotation, top_term)
    etichette <- vapply(annotation_top$Term, function(x) paste(strwrap(x,
      width = 45), collapse = "\n"), "")
    annotation_top$Etichetta <- factor(etichette,
      levels = rev(etichette))
    g <- ggplot(annotation_top, aes(x = Gene_count,
      y = Etichetta)) + geom_point(aes(size = Gene_ratio,
      color = Adjusted.P.value)) + scale_colour_gradient(low = "#C45C36",
      high = "#287C8E", labels = scales::label_scientific()) +
      scale_size_continuous(range = c(3, 8)) +
      labs(x = "Geni sovrapposti alla categoria",
        y = NULL, title = if (i == 1)
          "A  GO Biological Process 2025" else "B  KEGG 2021 Human",
        size = "Gene ratio",
        color = "p aggiustato") + theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8.5),
        plot.title = element_text(face = "bold"),
        legend.key.height = grid::unit(0.4, "cm"))
    print(g, vp = grid::viewport(layout.pos.row = i,
      layout.pos.col = 1))
  }
})
# 9. Diagnostica delle differenze: qqnorm/qqline e Shapiro-Wilk.
qa <- list()
salva("09_diagnostica_QQ", 8.2, 3.7, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (gene in scelti) {
    d <- G[gene, (n + 1):(2 * n)] - G[gene, 1:n]
    qqnorm(d, main = gene, pch = 16, cex = 0.6, col = N,
      xlab = "Quantili teorici normali", ylab = "Differenze T - N")
    qqline(d, col = T, lwd = 1.5)
  }
})
for (gene in scelti) {
  d <- G[gene, (n + 1):(2 * n)] - G[gene, 1:n]
  qa[[gene]] <- list(shapiro_p = shapiro.test(d)$p.value,
    positive_pairs = sum(d > 0), negative_pairs = sum(d <
      0))
}
# Controllo alternativo in R: ricostruire t e BH con le formule.
# Non e una seconda coorte ne una conferma sperimentale indipendente.
M <- G[res$gene_symbol, , drop = FALSE]
delta <- M[, (n + 1):(2 * n)] - M[, 1:n]
t_manuale <- rowMeans(delta)/(apply(delta, 1, sd)/sqrt(n))
p_manuale <- 2 * pt(-abs(t_manuale), df = n - 1)
ord <- order(p_manuale)
m <- length(p_manuale)
q_manuale <- numeric(m)
q_manuale[ord] <- pmin(1, rev(cummin(rev(p_manuale[ord] * m/seq_len(m)))))
qa$max_relative_p_error <- max(abs(p_manuale/res$p_value - 1))
qa$max_relative_BH_error <- max(abs(q_manuale/res$p_adj - 1))
stopifnot(qa$max_relative_p_error < 1e-08, qa$max_relative_BH_error < 1e-08,
  max(abs(rowMeans(delta) - res$logFC)) < 1e-10, sum(q_manuale < 0.05 &
    abs(res$logFC) >= 1) == nrow(deg))
qa$tested <- nrow(res)
qa$DEG <- nrow(deg)
qa$pairs <- n
qa$heatmap_k2 <- as.data.frame.matrix(table(clusters, campioni$group))
write_json(qa, file.path(O, "verifica_indipendente.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = 16)
message("Nove figure prodotte in R; controlli numerici superati.")
