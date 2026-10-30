# ELABORATO DI BIOINFORMATICA - avvio completo, esclusivamente in R.
# In RStudio aprire elaborato_bioinformatica.Rproj, quindi eseguire:
# source("src/00_mainelaborato.R", encoding = "UTF-8")
# Dal terminale: Rscript src/00_mainelaborato.R
# GEO fornisce i dati di espressione; tutti gli input finiscono in data/.
esegui_elaborato <- function() {
  if (getRversion() < "4.1.0") stop("Occorre R 4.1 o successivo.")
  # La radice e la cartella del progetto, non la sottocartella src.
  if (!file.exists("config/manifest_fonti.json")) {
    stop("Aprire il file .Rproj o impostare la radice del repository.")
  }
  # Elimina gli oggetti globali, comprese le variabili nascoste.
  base::rm(list = base::ls(envir = .GlobalEnv, all.names = TRUE),
    envir = .GlobalEnv)
  # Chiude tutti i dispositivi e azzera i grafici precedenti in RStudio.
  grDevices::graphics.off()
  options(stringsAsFactors = FALSE, digits = 16, timeout = 600)
  pacchetti <- c("jsonlite", "digest", "curl", "pheatmap",
    "ggplot2", "enrichR")
  dir.create(".R-library", showWarnings = FALSE)
  .libPaths(c(normalizePath(".R-library"), .libPaths()))
  mancanti <- setdiff(pacchetti, rownames(installed.packages()))
  if (length(mancanti)) {
    install.packages(mancanti, lib = ".R-library",
      repos = "https://cloud.r-project.org")
  }
  stopifnot(all(pacchetti %in% rownames(installed.packages())))
  for (cartella in c("data", "risultati", "figure")) {
    dir.create(cartella, showWarnings = FALSE)
  }
  for (script in c("01_acquisizione.R", "02_analisi_DEG.R",
    "03_arricchimento.R", "04_grafici_verifiche.R")) {
    message("Esecuzione: ", script)
    # Gli ambienti separati evitano di riutilizzare variabili residue.
    source(file.path("src", script), local = new.env(),
      encoding = "UTF-8")
  }
  writeLines(capture.output(sessionInfo()), "risultati/sessionInfo_R.txt")
  versioni <- data.frame(pacchetto = pacchetti,
    versione = vapply(pacchetti,
      function(x) as.character(packageVersion(x)), ""))
  write.table(versioni, "risultati/versioni_pacchetti_R.txt",
    sep = "\t", row.names = FALSE, quote = FALSE)
  # Il manifesto lega sorgenti, configurazione e output effettivi.
  files <- sort(list.files(c("src", "config", "tests", "figure", "risultati"),
    recursive = TRUE, full.names = TRUE))
  # Il test scrive il proprio esito dopo aver verificato questo manifesto:
  # escluderlo evita di invalidare il manifesto a ogni nuova verifica.
  files <- setdiff(files, "risultati/esiti_test.txt")
  hash <- vapply(files, digest::digest, "", algo = "sha256", file = TRUE)
  write.table(data.frame(file = files, sha256 = hash),
    "manifest_risultati_sha256.txt", sep = "\t", quote = FALSE,
    row.names = FALSE)
  message("Analisi, figure e verifiche completate in R.")
}
esegui_elaborato()
