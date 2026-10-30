# ACQUISIZIONE - download dei dati di espressione direttamente da GEO.
# Il primo avvio richiede Internet; i file validi sono riusati in seguito.
library(jsonlite)
library(digest)
library(curl)
options(timeout = 600)
dir.create("data", showWarnings = FALSE)
dir.create("risultati", showWarnings = FALSE)
fonti <- fromJSON("config/manifest_fonti.json")$sources
sha <- function(f) digest(f, algo = "sha256", file = TRUE)
verifica <- function(f, atteso) {
  if (!identical(sha(f), atteso)) {
    stop("SHA256 diverso per ", f,
      ". Conservare il file e verificare la fonte prima di proseguire.")
  }
}
# File temporaneo e controllo prima del trasferimento al nome definitivo.
scarica <- function(url, file, handle = new_handle()) {
  handle_setopt(handle, connecttimeout = 30, timeout = 600)
  curl_download(url, file, quiet = TRUE, handle = handle)
  if (!file.exists(file) || file.info(file)$size == 0) {
    stop("Download vuoto: ", basename(file))
  }
}
# La famiglia GPL include moltissimi campioni estranei allo studio.
# Leggere lo stream gzip fino alla fine della sola tabella GPL evita
# di scaricare l'intera famiglia, che occupa diversi gigabyte.
scarica_piattaforma <- function(url, file) {
  h <- new_handle(connecttimeout = 30, timeout = 600)
  ingresso <- gzcon(curl(url, open = "rb", handle = h))
  on.exit(close(ingresso), add = TRUE)
  # Modalità binaria: i terminatori LF restano identici anche su Windows.
  uscita <- base::file(file, "wb")
  on.exit(close(uscita), add = TRUE)
  repeat {
    righe <- readLines(ingresso, n = 1000, warn = FALSE)
    if (!length(righe)) stop("Tabella GPL6884 incompleta.")
    fine <- match("!platform_table_end", righe)
    if (!is.na(fine)) {
      writeLines(righe[seq_len(fine)], uscita, useBytes = TRUE)
      break
    }
    writeLines(righe, uscita, useBytes = TRUE)
  }
}
geo <- fonti[grepl("/geo/", fonti$url, fixed = TRUE), ]
for (i in seq_len(nrow(geo))) {
  f <- file.path("data", geo$file[i])
  if (!file.exists(f)) {
    message("Download diretto GEO: ", geo$file[i])
    tmp <- tempfile(tmpdir = "data", fileext = ".part")
    if (geo$file[i] == "GPL6884_platform.soft") {
      scarica_piattaforma(geo$url[i], tmp)
    } else scarica(geo$url[i], tmp)
    verifica(tmp, geo$sha256[i])
    stopifnot(file.rename(tmp, f))
  }
  verifica(f, geo$sha256[i])
}
# Le espressioni NON provengono da GitHub. La release conserva soltanto
# annotazioni NCBI, librerie GMT e risposte Enrichr del 16 settembre 2026.
# Queste fonti cambiano nel tempo: fissarle riproduce i numeri del PDF.
m <- fromJSON("config/snapshot.json")
presenti <- file.exists(file.path("data", m$files$file))
if (!all(presenti)) {
  archivio <- file.path("data", m$asset)
  if (!file.exists(archivio)) {
    # La release pubblica si scarica direttamente via HTTPS, come GEO.
    # Tag, nome e hash dell'archivio sono fissati nella configurazione.
    url <- paste0("https://github.com/", m$repository,
      "/releases/download/", m$tag, "/", m$asset)
    tmp <- tempfile(tmpdir = "data", fileext = ".part")
    tryCatch(scarica(url, tmp), error = function(e) {
      stop("Download della release non riuscito: ", conditionMessage(e),
        ". Verificare la connessione o scaricare l'archivio dalla ",
        "release nel browser e salvarlo in data/. Vedere il PDF, 7.1.",
        call. = FALSE)
    })
    verifica(tmp, m$archive_sha256)
    stopifnot(file.rename(tmp, archivio))
  }
  verifica(archivio, m$archive_sha256)
  # Il manifesto ammette soltanto nomi interni alla cartella data/.
  nomi <- untar(archivio, list = TRUE)
  if (!setequal(nomi, m$files$file)) stop("Contenuto archivio inatteso.")
  estratti <- tempfile(pattern = ".estrazione-", tmpdir = "data")
  dir.create(estratti)
  untar(archivio, exdir = estratti)
  for (i in seq_len(nrow(m$files))) {
    origine <- file.path(estratti, m$files$file[i])
    verifica(origine, m$files$sha256[i])
    destinazione <- file.path("data", m$files$file[i])
    # Un file modificato non viene sostituito in silenzio.
    if (file.exists(destinazione)) {
      verifica(destinazione, m$files$sha256[i])
    } else {
      dir.create(dirname(destinazione), showWarnings = FALSE,
        recursive = TRUE)
      stopifnot(file.copy(origine, destinazione))
    }
  }
  unlink(estratti, recursive = TRUE)
}
for (i in seq_len(nrow(m$files))) {
  verifica(file.path("data", m$files$file[i]), m$files$sha256[i])
}
# Il rapporto comprende le quattro fonti GEO/NCBI usate nell'analisi.
v <- data.frame(file = fonti$file,
  sha256 = unname(vapply(file.path("data", fonti$file), sha, "")),
  bytes = file.info(file.path("data", fonti$file))$size, esito = "OK")
stopifnot(identical(v$sha256, fonti$sha256))
write.table(v, "risultati/verifica_fonti_R.txt", sep = "\t",
  quote = FALSE, row.names = FALSE)
message("Dati GEO e annotazioni verificati con SHA256.")
