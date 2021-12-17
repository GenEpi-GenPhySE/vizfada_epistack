# paths: change accordingly ---------------
EPISTACKPATH <- "/home/gdevailly/R/x86_64-pc-linux-gnu-library/4.1/epistack/epistack.R"
BASEDIR <- "/work/lmorel/data"

# EPISTACKPATH <- "~/mnt/inra_p/projets/vizfada/epistack.R"
# PREFIX <- "~/mnt/genotoul_lmoreldata"

# library loading ----------
# library(purrr)
# library(stringr)
# library(jsonlite)
# library(rtracklayer)

library(biomaRt)
library(data.table)

species <- "pig"
chip <- "ERX3212576"
rna <- "ERX3212540"


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
    # TODO
    # fwrite(
    #     gene_starts,
    #     file = file.path(basedir, species, )
    # )
}
