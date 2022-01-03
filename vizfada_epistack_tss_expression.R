# paths: change accordingly ---------------
BASEDIR <- "/work/lmorel/data"


# library loading ----------
library(biomaRt)
library(data.table)

# species <- "pig"
# chip <- "ERX3212576"
# rna <- "ERX3212540"

getTssCoordinates <- function(species, basedir) {
    ensembl <- useEnsembl(biomart = "genes")
    pattern <- switch (species,
        "pig" = "Sscrofa"
    )
    dataset <- searchDatasets(ensembl, pattern)
    if (nrow(dataset) != 1L) {
        stop(
            paste("Could not identify Ensembl dataset unambigously.\n
                  Dataset found:\n", dataset)
        )
    }
    dataset <- useDataset(dataset = dataset$dataset[[1]], mart = ensembl)

    transcripts <- getBM(
        attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                       "chromosome_name", "start_position", "end_position",
                       "strand", "transcription_start_site", "gene_biotype"),
        mart = dataset
    ) |> data.table()

    gene_starts <- transcripts[
        ,
        .(median_tss = round(median(transcription_start_site))),
        by = .(ensembl_gene_id, chromosome_name, strand, gene_biotype)
    ]

    gene_starts <- gene_starts[
        ,
        c("chromosome_name", "median_tss", "strand", "ensembl_gene_id", "gene_biotype")
    ]

    fwrite(
        gene_starts,
        file = file.path(basedir, species, "chipseq", "epistack", "median_tss.tsv"),
        sep = "\t"
    )
}

getTssCoordinates("pig", basedir = BASEDIR)
