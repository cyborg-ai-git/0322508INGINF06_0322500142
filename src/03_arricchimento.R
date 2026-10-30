# ARRICCHIMENTO
# La replica legge le risposte reali del 16/09/2026. Una nuova query,
# opzionale, usa enrichr() con gli stessi geni e background, in R.
library(jsonlite)
library(digest)
O <- "risultati"
S <- file.path(O, "Enrichr")
dir.create(S, showWarnings = FALSE)
dbs <- c("GO_Biological_Process_2025", "KEGG_2021_Human")
file_liste <- c(tutti = "geni_DEG_per_Enrichr.txt", up = "geni_up.txt",
  down = "geni_down.txt")
background <- readLines(file.path(O, "background_geni_testati.txt"))
# FALSE riproduce la relazione anche senza Internet. TRUE richiede
# Internet e salva una nuova analisi separata: non aggiorna la relazione.
NUOVA_QUERY <- isTRUE(getOption("elaborato.nuova_query", FALSE))
if (NUOVA_QUERY) {
  library(enrichR)
  setEnrichrSite("Enrichr")
  nuova <- file.path(O, paste0("Enrichr_", format(Sys.time(),
    "%Y%m%d_%H%M%S")))
  dir.create(nuova)
  for (nome in names(file_liste)) {
    input_list <- readLines(file.path(O, file_liste[nome]))
    df <- enrichr(input_list, dbs, background = background,
      include_overlap = TRUE)
    stopifnot(all(dbs %in% names(df)))
    saveRDS(df, file.path(nuova, paste0(nome, ".rds")))
    for (db in dbs) write.table(df[[db]], file.path(nuova,
      paste0(nome, "_", db, ".txt")), sep = "\t", quote = FALSE,
      row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(nuova,
    "sessionInfo.txt"))
}
# Verifica congiunta di liste, background, librerie GMT e risposte JSON.
# Non e sufficiente avere file con lo stesso nome: i contenuti devono
# corrispondere agli input dell'interrogazione archiviata.
manifest <- fromJSON("config/manifest_enrichr.json")
for (i in seq_len(nrow(manifest$files))) {
  f <- file.path(O, manifest$files$file[i])
  stopifnot(digest(f, algo = "sha256", file = TRUE) ==
    manifest$files$sha256[i])
}
campi <- c("Rank", "Term", "P.value", "Odds.Ratio", "Combined.Score",
  "Genes", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value")
riepilogo <- list()
for (nome in names(file_liste)) {
  input_list <- readLines(file.path(O, file_liste[nome]))
  stopifnot(!anyDuplicated(input_list), all(input_list %in%
    background))
  for (db in dbs) {
    # I nove campi sono quelli dell'API Speedrichr, utilizzata anche
    # da enrichR quando si specifica un background sperimentale.
    risposta <- fromJSON(file.path("data/Enrichr", paste0(nome,
      "_", db, ".json")), simplifyVector = FALSE)[[db]]
    annotation <- as.data.frame(do.call(rbind,
      lapply(risposta, function(x) {
        stopifnot(length(x) == 9)
        x[[6]] <- paste(unlist(x[[6]]), collapse = ";")
        unlist(x)
      })), stringsAsFactors = FALSE)
    names(annotation) <- campi
    for (v in setdiff(campi, c("Term", "Genes"))) {
      annotation[[v]] <- as.numeric(annotation[[v]])
    }
    stopifnot(all(is.finite(annotation$Adjusted.P.value)))
    # Gene_count e Gene_ratio seguono l'esercitazione. Il denominatore
    # di Gene_ratio e il totale dei geni annotati nella categoria GMT.
    gmt <- strsplit(readLines(file.path("data",
      paste0(db, ".gmt")), warn = FALSE), "\t",
      fixed = TRUE)
    dimensioni <- vapply(gmt, function(x) length(unique(x[-c(1,
      2)])), 0L)
    names(dimensioni) <- vapply(gmt, function(x) x[1],
      "")
    annotation$Gene_count <- lengths(strsplit(annotation$Genes,
      ";"))
    annotation$Term_size <- unname(dimensioni[annotation$Term])
    stopifnot(!anyNA(annotation$Term_size), all(annotation$Term_size >
      0))
    annotation$Overlap <- paste(annotation$Gene_count,
      annotation$Term_size, sep = "/")
    annotation$Gene_ratio <- annotation$Gene_count/annotation$Term_size
    annotation$List_size <- length(input_list)
    # Questo secondo rapporto ha un denominatore diverso ed e distinto.
    annotation$Gene_ratio_input <- annotation$Gene_count/length(input_list)
    annotation <- annotation[order(annotation$Adjusted.P.value,
      annotation$P.value, annotation$Term),
      ]
    thr_pval <- if (startsWith(db, "GO_"))
      0.01 else 0.05
    significativi <- annotation[annotation$Adjusted.P.value <
      thr_pval, ]
    write.table(annotation, file.path(S, paste0(nome,
      "_", db, ".txt")), sep = "\t", quote = FALSE,
      row.names = FALSE)
    riepilogo[[length(riepilogo) + 1L]] <- data.frame(list = nome,
      library = db, input_size = length(input_list),
      background_size = length(background),
      terms_returned = nrow(annotation), alpha = thr_pval,
      significant = nrow(significativi))
    message(nome, " / ", db, ": ", nrow(significativi),
      " termini significativi")
  }
}
write.table(do.call(rbind, riepilogo), file.path(S, "riepilogo.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("Risposte Enrichr originali: 16 settembre 2026.",
  "Tabelle rigenerate e verificate con R:",
  as.character(Sys.time()),
  "Gene_ratio = Gene_count / Term_size, come nell'esercitazione.",
  "Le sei famiglie lista x libreria mantengono i p aggiustati del servizio."),
  file.path(S, "rielaborazione_R.txt"))
