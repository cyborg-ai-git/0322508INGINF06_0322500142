# Verifica indipendente dei risultati generati dal progetto.
# Eseguire dopo src/00_mainelaborato.R; salva risultati/esiti_test.txt.
library(jsonlite)
library(digest)
base <- normalizePath(".")
leggi <- function(f) read.delim(file.path(base, "risultati", f),
  check.names = FALSE)
esiti <- list()
ok <- function(nome, condizione, dettaglio) {
  stopifnot(all(condizione))
  esiti[[length(esiti) + 1L]] <<- data.frame(controllo = nome, esito = "OK",
    dettaglio = dettaglio)
  message(nome, ": OK")
}
# Integrita di codice, risultati e figure rispetto al manifesto.
manifest <- read.delim(file.path(base, "manifest_risultati_sha256.txt"))
hash <- vapply(file.path(base, manifest$file), function(f) digest(f,
  file = TRUE, algo = "sha256"), "")
ok("Integrita della consegna", hash == manifest$sha256, paste(length(hash),
  "file"))
res <- leggi("tutti_geni_testati.txt")
deg <- leggi("geni_differenzialmente_espressi.txt")
meta <- leggi("metadati_tutti.txt")
coppie <- leggi("coppie_analizzate.txt")
coorte <- leggi("coorte_pazienti.txt")
G <- as.matrix(read.delim(gzfile(file.path(base, "risultati",
  "espressione_geni_log2_normalizzata.txt.gz")), row.names = 1,
  check.names = FALSE))
ok("Pazienti e gruppi", c(nrow(coppie) == 57,
  !anyDuplicated(coppie$subject), nrow(meta) ==
    116, sum(meta$included) == 114, nrow(coorte) ==
    57), "116 campioni disponibili, 57 coppie complete analizzate")
ok("Assenza di duplicati e dati mancanti",
  c(!anyDuplicated(rownames(G)), all(is.finite(G)),
    ncol(G) == 114), "Matrice di 20.120 geni e 114 campioni")
iq <- apply(G, 1, IQR, type = 7)
keep <- iq >= quantile(iq, 0.1, type = 7)
ok("Filtro IQR", c(sum(keep) == 18108, setequal(rownames(G)[keep],
  res$gene_symbol)), "10 percentile; 18.108 test")
n <- nrow(coppie)
M <- G[res$gene_symbol, , drop = FALSE]
delta <- M[, (n + 1):(2 * n)] - M[, 1:n]
p <- 2 * pt(-abs(rowMeans(delta)/(apply(delta, 1, sd)/sqrt(n))), df = n - 1)
q <- p.adjust(p, "BH")
ok("Test t appaiato e FDR", c(max(abs(p/res$p_value - 1)) < 1e-08,
  max(abs(q/res$p_adj - 1)) < 1e-08, max(abs(rowMeans(delta) - res$logFC)) <
    1e-10), "57 coppie; 56 gradi di liberta; BH sull'intera famiglia")
ok("Selezione e TXT", c(all(c("gene_symbol", "p_adj",
  "logFC", "direction") %in% names(deg)), nrow(deg) ==
  919, sum(deg$direction == "up") == 323, sum(deg$direction ==
  "down") == 596, setequal(deg$gene_symbol,
  res$gene_symbol[q < 0.05 & abs(res$logFC) >=
    1])), "919 DEG: 323 up e 596 down, con tutti i campi richiesti")
# PCA con calcolo indipendente sulla matrice dei prodotti tra campioni.
# eigen e prcomp possono cambiare il segno delle componenti: il controllo
# confronta varianza e correlazione assoluta delle coordinate.
H <- G[deg$gene_symbol, , drop = FALSE]
X <- scale(t(H), center = TRUE, scale = TRUE)
autovalori <- eigen(tcrossprod(X), symmetric = TRUE)
vr <- pmax(autovalori$values, 0)/sum(pmax(autovalori$values, 0))
varianza <- leggi("PCA_DEG_varianza.txt")$variance_ratio
score <- leggi("PCA_DEG_scores.txt")
coordinate <- sweep(autovalori$vectors[, 1:5], 2,
  sqrt(autovalori$values[1:5]), "*")
correlazioni <- abs(diag(cor(coordinate, as.matrix(score[, paste0("PC",
  1:5)]))))
ok("PCA sui soli DEG", c(nrow(X) == 114, ncol(X) == 919,
  max(abs(vr - varianza)) < 1e-10, all(correlazioni >
    1 - 1e-10), identical(score$sample, colnames(H))),
  sprintf("114 campioni x 919 DEG; PC1 %.8f%%; PC2 %.8f%%; somma %.8f%%",
    100 * vr[1], 100 * vr[2], 100 * sum(vr[1:2])))
# Verifica della dimensione e del clustering della heatmap.
Z <- t(scale(t(H)))
albero <- hclust(as.dist(1 - cor(Z)), method = "average")
atteso <- cutree(albero, k = 2)
cluster <- leggi("cluster_campioni_heatmap.txt")
ok("Heatmap completa e clustering", c(nrow(Z) ==
  919, ncol(Z) == 114, identical(cluster$sample,
  colnames(Z)), all(atteso == cluster$cluster)),
  "919 righe e 114 colonne; correlazione Pearson e average linkage")
# I simboli nelle risposte sono resi maiuscoli dal servizio Enrichr.
file_liste <- c(tutti = "geni_DEG_per_Enrichr.txt", up = "geni_up.txt",
  down = "geni_down.txt")
dbs <- c("GO_Biological_Process_2025", "KEGG_2021_Human")
sommario <- leggi("Enrichr/riepilogo.txt")
for (nome in names(file_liste)) {
  geni <- toupper(readLines(file.path(base, "risultati", file_liste[nome])))
  for (db in dbs) {
    a <- leggi(file.path("Enrichr", paste0(nome, "_", db, ".txt")))
    liste <- strsplit(a$Genes, ";", fixed = TRUE)
    alpha <- if (startsWith(db, "GO_"))
      0.01 else 0.05
    ok(paste("Enrichr", nome, db), c(all(unlist(liste) %in%
      geni), all(lengths(liste) == a$Gene_count), max(abs(a$Gene_ratio -
      a$Gene_count/a$Term_size)) < 1e-10, sum(a$Adjusted.P.value <
      alpha) == sommario$significant[sommario$list == nome &
      sommario$library == db]), paste(nrow(a), "termini verificati; soglia",
      alpha))
  }
}
png <- list.files("figure", pattern = "[.]png$")
pdf <- list.files("figure", pattern = "[.]pdf$")
ok("Figure R", c(length(png) == 9, length(pdf) == 9),
  "Nove PNG e nove PDF vettoriali")
write.table(do.call(rbind, esiti), "risultati/esiti_test.txt", sep = "\t",
  quote = FALSE, row.names = FALSE)
message("Tutti i controlli della replica sono superati.")
