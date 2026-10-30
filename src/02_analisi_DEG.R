# Analisi riproducibile di GSE32863, Elementi di Bioinformatica
# Eseguire dalla radice del repository: Rscript src/00_mainelaborato.R
# Usa soltanto R base e gli snapshot dei dati nella cartella data.
# Orientamento delle matrici: righe = sonde o geni; colonne = campioni.
# G contiene i geni aggregati, F quelli ammessi dopo il filtro IQR.
# Le prime 57 colonne sono N e le successive 57 sono T, nello stesso
# ordine di soggetto. res contiene tutti i test, deg la selezione finale.
options(stringsAsFactors = FALSE, digits = 16)
# Soglie principali fissate prima del confronto tra i gruppi.
# Impostazione come nel caso studio: prc_IQR, thr_FC, thr_pval, paired.
# Il fold change 2 corrisponde a logFC=1 (o -1 nella direzione opposta).
prc_IQR <- 0.1
thr_FC <- 2
thr_pval <- 0.05
paired <- TRUE
IQR_PERCENTILE <- prc_IQR
FDR_THRESHOLD <- thr_pval
LOGFC_THRESHOLD <- log2(thr_FC)
root <- normalizePath(".")
D <- file.path(root, "data")
O <- file.path(root, "risultati")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
save_tsv <- function(x, name, row.names = FALSE) write.table(x, file.path(O,
  name), sep = "\t", quote = FALSE, row.names = row.names, col.names = TRUE,
  na = "NA")
# Le liste confrontate con gli snapshot devono avere gli stessi byte
salva_lista <- function(geni, nome) {
  con <- file(file.path(O, nome), open = "wb")
  on.exit(close(con))
  writeLines(enc2utf8(geni), con, sep = "\n", useBytes = TRUE)
}

# Metadati: le righe characteristics non hanno sempre la stessa chiave
# per ciascun campione. Ogni cella viene quindi interpretata separatamente.
con <- gzfile(file.path(D, "GSE32863_series_matrix.txt.gz"), "rt")
header <- character()
repeat {
  z <- readLines(con, n = 1)
  if (!length(z) || z == "!series_matrix_table_begin")
    break
  header <- c(header, z)
}
close(con)
# La pipe |> passa le celle estratte a gsub, che rimuove le virgolette.
cells <- function(z) strsplit(z, "\t", fixed = TRUE)[[1]][-1] |>
  gsub(pattern = "^\"|\"$", replacement = "")
gsm <- cells(header[startsWith(header, "!Sample_geo_accession\t")])
titles <- cells(header[startsWith(header, "!Sample_title\t")])
meta <- data.frame(gsm = gsm, title = titles)
for (z in header[startsWith(header, "!Sample_characteristics_ch1\t")]) {
  vals <- cells(z)
  for (j in seq_along(vals)) if (grepl(": ", vals[j], fixed = TRUE)) {
    key <- tolower(sub(": .*", "", vals[j]))
    key <- gsub(" ", "_", key, fixed = TRUE)
    if (!key %in% names(meta))
      meta[[key]] <- NA_character_
    meta[j, key] <- sub("^[^:]+: ", "", vals[j])
  }
}
meta$group <- ifelse(meta$tissue == "Normal lung", "N", "T")
meta$raw_id <- paste(meta$subject, meta$group, sep = "_")
stopifnot(!anyDuplicated(meta$gsm), !anyDuplicated(meta$raw_id))
pair_tab <- table(meta$subject, meta$group)
subjects <- sort(rownames(pair_tab)[pair_tab[, "N"] == 1 & pair_tab[,
  "T"] == 1])
meta$included <- meta$subject %in% subjects
save_tsv(meta, "metadati_tutti.txt")
save_tsv(meta[!meta$included, ], "campioni_esclusi.txt")
mn <- meta[match(paste0(subjects, "_N"), meta$raw_id), ]
mt <- meta[match(paste0(subjects, "_T"), meta$raw_id), ]
stopifnot(identical(mn$subject, mt$subject), nrow(mn) == 57)
for (v in c("age", "gender", "smoking_status",
  "stage")) stopifnot(identical(mn[[v]], mt[[v]]))
save_tsv(data.frame(subject = subjects, normal = mn$gsm, tumor = mt$gsm),
  "coppie_analizzate.txt")
save_tsv(mt, "coorte_pazienti.txt")

# Intensita non normalizzate e detection p-value alternati nel file GEO.
raw <- read.delim(gzfile(file.path(D, "GSE32863_non-normalized.txt.gz")),
  check.names = FALSE, row.names = 1)
# Le colonne Detection Pval misurano il segnale rispetto al fondo del chip.
# Non sono p-value di espressione differenziale e non entrano nel t-test.
intensity <- as.matrix(raw[, seq(1, ncol(raw), 2)])
stopifnot(all(is.finite(intensity)), all(intensity >= 0),
  !anyDuplicated(rownames(intensity)))
selected <- c(mn$raw_id, mt$raw_id)
stopifnot(all(selected %in% colnames(intensity)))
X <- log2(intensity[, selected, drop = FALSE] + 1)
# Normalizzazione quantile: distribuzione comune e gestione dei pareggi.
quantile_normalize <- function(x) {
  # La distribuzione comune e la media dei valori ordinati per ciascun array.
  # Ogni valore riceve quello del rango corrispondente. Per i pareggi si usa
  # la media dei valori obiettivo ai ranghi interessati.
  target <- rowMeans(apply(x, 2, sort))
  ans <- apply(x, 2, function(v) {
    ord <- order(v)
    out <- numeric(length(v))
    groups <- cumsum(c(TRUE, diff(v[ord]) != 0))
    out[ord] <- ave(target, groups, FUN = mean)
    out
  })
  dimnames(ans) <- dimnames(x)
  ans
}
Q <- quantile_normalize(X)
quantiles <- function(m) t(apply(m, 2, quantile, probs = c(0, 0.25, 0.5,
  0.75, 1), type = 7))
save_tsv(data.frame(sample = selected, group = c(rep("N", 57), rep("T",
  57)), quantiles(X), check.names = FALSE), "quantili_prima.txt")
save_tsv(data.frame(sample = selected, group = c(rep("N", 57), rep("T",
  57)), quantiles(Q), check.names = FALSE), "quantili_dopo.txt")

# Identita delle sonde: GPL6884 -> Entrez Gene ID -> simbolo NCBI corrente.
platform <- readLines(file.path(D, "GPL6884_platform.soft"), warn = FALSE)
i <- match("!platform_table_begin", platform)
j <- match("!platform_table_end", platform)
ann <- read.delim(text = paste(platform[(i + 1):(j - 1)],
  collapse = "\n"), quote = "", comment.char = "", check.names = FALSE,
  colClasses = "character")
gi <- read.delim(gzfile(file.path(D, "Homo_sapiens.gene_info.gz")),
  quote = "", comment.char = "", check.names = FALSE,
  colClasses = "character")
ann <- ann[match(rownames(Q), ann$ID), ]
ann$gene_id <- ann$Entrez_Gene_ID
unambiguous <- !is.na(ann$gene_id) & grepl("^[0-9]+$", ann$gene_id)
ann$gene_symbol <- gi$Symbol[match(ann$gene_id, gi$GeneID)]
valid <- unambiguous & !is.na(ann$gene_symbol) & ann$gene_symbol != "-"
save_tsv(ann[, c("ID", "Symbol", "gene_id", "gene_symbol")],
  "annotazione_sonde.txt")
groups <- ann$gene_symbol[valid]
G <- rowsum(Q[valid, , drop = FALSE], groups, reorder = TRUE)
counts <- table(groups)
# rowsum produce somme: la divisione per il numero di sonde restituisce
# la media delle espressioni log2 normalizzate di ciascun gene.
G <- G/as.numeric(counts[rownames(G)])
stopifnot(all(is.finite(G)), !anyDuplicated(rownames(G)))
# Si esplicita type=7, convenzione predefinita di R per i quantili.
# La slide usa type=5: e una convenzione diversa, non una soglia biologica.
# La scelta qui e fissata e mantenuta anche nei controlli di sensibilita.
iqr <- apply(G, 1, IQR, type = 7)
threshold <- unname(quantile(iqr, IQR_PERCENTILE, type = 7))
keep <- iqr >= threshold
F <- G[keep, , drop = FALSE]
save_tsv(data.frame(gene_symbol = names(iqr), IQR = iqr, retained = keep),
  "filtro_IQR.txt")

# Confronto appaiato bilaterale, 57 soggetti e 56 gradi di liberta.
n <- length(subjects)
test_matrix <- function(M) {
  # dataN e dataC richiamano i nomi del caso studio delle dispense.
  # Le colonne dei due gruppi devono corrispondere agli stessi soggetti.
  dataN <- M[, 1:n, drop = FALSE]
  dataC <- M[, (n + 1):(2 * n), drop = FALSE]
  delta <- dataC - dataN
  tests <- t(vapply(seq_len(nrow(M)), function(i) {
    # Una varianza nulla non consente il normale calcolo del test t.
    # Si interrompe per rendere esplicito il problema, senza inventare p.
    stopifnot(sd(delta[i, ]) > 0)
    tt <- t.test(dataC[i, ], dataN[i, ], paired = paired,
      alternative = "two.sided")
    c(t = unname(tt$statistic), p = tt$p.value, low = tt$conf.int[1],
      high = tt$conf.int[2])
  }, numeric(4)))
  data.frame(gene_symbol = rownames(M), logFC = rowMeans(delta),
    FC = 2^rowMeans(delta), mean_normal = rowMeans(M[, 1:n,
      drop = FALSE]), mean_tumor = rowMeans(M[, (n + 1):(2 *
      n), drop = FALSE]), t = tests[, "t"], p_value = tests[,
      "p"], p_adj = p.adjust(tests[, "p"], method = "fdr"),
    CI95_low = tests[, "low"], CI95_high = tests[, "high"],
    row.names = NULL)
}
res <- test_matrix(F)
res$n_probes <- as.integer(counts[res$gene_symbol])
res$IQR <- iqr[res$gene_symbol]
res$direction <- ifelse(res$logFC > 0, "up", "down")
# Il volcano descrive la famiglia completa dei test.
# In p.adjust, 'fdr' e un alias di 'BH': la correzione precede il logFC.
res$DEG <- res$p_adj < FDR_THRESHOLD & abs(res$logFC) >= LOGFC_THRESHOLD
res <- res[order(res$p_adj, -abs(res$logFC), res$gene_symbol), ]
deg <- res[res$DEG, ]
save_tsv(res, "tutti_geni_testati.txt")
save_tsv(deg, "geni_differenzialmente_espressi.txt")
for (dir in c("up", "down")) salva_lista(deg$gene_symbol[deg$direction ==
  dir], paste0("geni_", dir, ".txt"))
salva_lista(deg$gene_symbol, "geni_DEG_per_Enrichr.txt")
salva_lista(res$gene_symbol, "background_geni_testati.txt")
top <- rbind(head(deg[order(-deg$logFC), ], 10), head(deg[order(deg$logFC),
  ], 10))
save_tsv(top, "top_geni.txt")
con <- gzfile(file.path(O, "espressione_geni_log2_normalizzata.txt.gz"),
  "wt")
write.table(data.frame(gene_symbol = rownames(G), G, check.names = FALSE),
  con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

# Controlli di sensibilita: IQR e Wilcoxon; la correzione include sempre
# tutti i geni della specifica famiglia di test, prima della soglia logFC.
sens <- do.call(rbind, lapply(c(0, 0.05, 0.1, 0.2), function(q) {
  M <- G[iqr >= quantile(iqr, q, type = 7), , drop = FALSE]
  R <- test_matrix(M)
  set <- R$gene_symbol[R$p_adj < 0.05 & abs(R$logFC) >= 1]
  data.frame(percentile = q * 100, tested = nrow(R), DEG = length(set),
    shared_primary = sum(set %in% deg$gene_symbol))
}))
save_tsv(sens, "sensibilita_IQR.txt")
dif <- F[, (n + 1):(2 * n)] - F[, 1:n]
wp <- apply(dif, 1, function(v) wilcox.test(v, exact = FALSE,
  correct = TRUE)$p.value)
wq <- p.adjust(wp, "BH")
save_tsv(data.frame(gene_symbol = names(wp), p_wilcoxon = wp,
  p_adj_wilcoxon = wq), "sensibilita_Wilcoxon.txt")
survived <- sum(wq[deg$gene_symbol] < 0.05)

# PCA su tutti i geni filtrati e sui soli DEG, con standardizzazione.
# t(M) mette i campioni sulle righe, come richiesto da prcomp.
# center centra ogni gene; scale. lo divide per la deviazione standard.
# pc$x contiene le coordinate dei campioni; pc$sdev^2 le varianze.
# La PCA sui DEG descrive geni selezionati negli stessi campioni:
# la separazione dei gruppi non misura una validazione indipendente.
for (label in c("filtrati", "DEG")) {
  M <- if (label == "filtrati")
    F else F[deg$gene_symbol, , drop = FALSE]
  pc <- prcomp(t(M), center = TRUE, scale. = TRUE)
  vr <- pc$sdev^2/sum(pc$sdev^2)
  save_tsv(data.frame(sample = selected, subject = rep(subjects,
    2), group = rep(c("N", "T"), each = n), pc$x[, 1:5]), paste0("PCA_",
    label, "_scores.txt"))
  save_tsv(data.frame(PC = seq_along(vr), variance_ratio = vr),
    paste0("PCA_", label, "_varianza.txt"))
}
qc <- data.frame(metric = c("samples_available", "subjects_available",
  "pairs_analyzed", "samples_analyzed", "raw_probes", "raw_min",
  "raw_max", "missing_intensities", "probes_unambiguous_Entrez",
  "probes_current_symbol", "genes_before_IQR", "IQR_threshold",
  "genes_removed_IQR", "genes_tested", "DEG", "up", "down",
  "BH_significant", "Wilcoxon_confirmed_DEG", "pairs_same_batch"),
  value = c(nrow(meta), length(unique(meta$subject)), n, 2 *
    n, nrow(Q), min(intensity), max(intensity), sum(is.na(intensity)),
    sum(unambiguous), sum(valid), nrow(G), threshold, sum(!keep),
    nrow(F), nrow(deg), sum(deg$direction == "up"), sum(deg$direction ==
      "down"), sum(res$p_adj < 0.05), survived, sum(mn$batch ==
      mt$batch)))
save_tsv(qc, "riepilogo_analisi.txt")
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo_R.txt"))
print(qc, row.names = FALSE)
print(top[, c("gene_symbol", "logFC", "p_adj")], row.names = FALSE)
print(sens, row.names = FALSE)
