# ELEMENTI DI BIOINFORMATICA 0322508INGINF06

## ELABORATO_DEL_2026_02_23_20_33_36

## 1. Preparare R e RStudio

1. Scaricare e installare **R 4.1 o successivo** dal [sito ufficiale CRAN](https://cloud.r-project.org/).
2. Scaricare e installare **RStudio Desktop** dalla [pagina ufficiale di Posit](https://docs.posit.co/ide/user/download.html), scegliendo la versione compatibile con il proprio sistema operativo.
3. Aprire **RStudio** dopo aver installato R.

> L'analisi è stata verificata con R 4.6.1 su macOS e su un runner Windows
> di GitHub Actions. Su entrambi i sistemi si usano gli stessi comandi.
> Al primo avvio occorre una connessione Internet per scaricare dati e
> pacchetti mancanti.

## 2. Scaricare e aprire il progetto

### Con un archivio ZIP

1. Aprire il [repository GitHub](https://github.com/cyborg-ai-git/0322508INGINF06_0322500142).
2. Scegliere **Code → Download ZIP** ed estrarre l'archivio.
3. Aprire il file **`elaborato_bioinformatica.Rproj`** contenuto nella cartella estratta.

### Con Git da RStudio

In alternativa, con Git installato, scegliere in RStudio
**File → New Project → Version Control → Git** e inserire:

```text
https://github.com/cyborg-ai-git/0322508INGINF06_0322500142.git
```

Scegliere una cartella locale e completare la creazione del progetto.

## 3. Avviare l'analisi

All'avvio vengono eliminati tutti gli oggetti dell'Environment, comprese le
variabili nascoste, e vengono chiusi i dispositivi grafici precedenti.
Salvare prima gli eventuali oggetti che si desidera conservare.

Nella scheda **Console** di RStudio incollare e premere Invio:

```r
source("src/00_mainelaborato.R", encoding = "UTF-8")
```

Il punto di ingresso installa i pacchetti mancanti, scarica i dati necessari
in `./data` ed esegue tutti gli script nell'ordine corretto. I dati di
espressione vengono scaricati direttamente da GEO.

Attendere il messaggio:

```text
Analisi, figure e verifiche completate in R.
```

## 4. Consultare risultati e grafici

Nel riquadro **Files** di RStudio aprire:

| Cartella | Contenuto |
| --- | --- |
| `risultati/` | Tabelle TXT, liste dei geni e risultati della PCA |
| `figure/` | Nove grafici, ciascuno in PNG e PDF |
| `data/` | File scaricati automaticamente |

Durante l'esecuzione interattiva in RStudio, i nove grafici compaiono anche
nella scheda **Plots**. Al termine è visibile l'ultimo: usare le frecce
**indietro/avanti** per rivedere gli altri e **Zoom** per ingrandirli.
Allargare il riquadro Plots se un grafico con più pannelli non entra.
I file PNG e PDF sono sempre salvati in `figure/`.

Per salvare solo i file, senza anteprime, impostare prima dell'avvio:

```r
options(elaborato.mostra_plots = FALSE)
```

Per riattivare le anteprime, impostare la stessa opzione a `TRUE`.
Una nuova esecuzione rigenera e sovrascrive gli output locali.

## 5. Eseguire i controlli

Dopo l'analisi, nella Console eseguire:

```r
source("tests/verifica_risultati.R", encoding = "UTF-8")
```

Attendere **Tutti i controlli della replica sono superati.** Il riepilogo viene
salvato in `risultati/esiti_test.txt`.

## Se R non trova i file

Aprire il file `.Rproj` e controllare la cartella di lavoro:

```r
getwd()
file.exists("config/manifest_fonti.json")  # Deve restituire TRUE.
```

La cartella corretta contiene `src`, `config` e il file `.Rproj`.
Se un controllo segnala un hash diverso, conservare il file e leggere il
messaggio completo prima di proseguire: non disattivare la verifica.

## Nuove interrogazioni Enrichr

Per aggiungere nuove query online alla replica, eseguire:

```r
options(elaborato.nuova_query = TRUE)
source("src/00_mainelaborato.R", encoding = "UTF-8")
options(elaborato.nuova_query = FALSE)
```

Le nuove risposte sono salvate in una cartella `Enrichr` con data e ora dentro
`risultati`.
