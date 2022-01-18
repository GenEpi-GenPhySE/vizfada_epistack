# paths: change accordingly ---------------
EPISTACKPATH <- "/home/gdevailly/R/x86_64-pc-linux-gnu-library/4.1/epistack/epistack.R"
PREFIX <- "/work/lmorel/data"

# EPISTACKPATH <- "~/mnt/inra_p/projets/vizfada/epistack.R"
# PREFIX <- "~/mnt/genotoul_lmoreldata"

# library loading ----------
library(purrr)
library(data.table)
library(stringr)
library(jsonlite)


# function definitions ----------------
make_commands <- function(species, epistack_path, data_dir, write_meta_table = FALSE) {

    metapath <- list.files(
        file.path(data_dir, species, "chipseq", "metadata"),
        pattern = "chipseq\\.tsv$", full.names = TRUE
    )
    if (length(metapath) != 1)
        stop("more than one .tsv file found in chipseq/metadata")

    metatarget <- fread(metapath, select = c("accession", "cellType", "experiment"))
    metatarget$cellType <- map_chr(metatarget$cellType, function(x) {
        fromJSON(str_replace_all(x, "'", '"'))$text
    })
    metatarget$experiment <- map_chr(metatarget$experiment, function(x) {
        fromJSON(str_replace_all(x, "'", '"'))$target
    })
    metatarget <- unique(metatarget)

    inputs <- list.files(file.path(data_dir, species, "chipseq"), pattern = "^ERX")

    peaks <- map(set_names(inputs), function(x)
        list.files(
            file.path(data_dir, species, "chipseq", x, "macs", "narrowPeak"),
            pattern = "peaks.narrowPeak$"
        ))

    dfc <- map_dfr(peaks, function(x) {
        data.table(
            bound_id = str_extract(x, "^[:alnum:]+"),
            peaks = x
        )
    }, .id = "input_id")
    dfc$input_bw = paste0(dfc$input_id, "_R1.bigWig")
    dfc$bound_bw = paste0(dfc$bound_id, "_R1.bigWig")

    dfc <- merge(dfc, metatarget, by.x = "bound_id", by.y = "accession", all.x = TRUE)
    colnames(dfc) <- c("bound_id", "input_id", "anchors", "input_bw", "bound_bw", "cellType", "experiment")
    dfc$anchors <- "median_tss.tsv"

    expr1 <- fread(file.path(data_dir, species, paste0(species, "_matching_specimen.csv")), select = c("accession_chip", "accession_rna"))
    expr1$notes <- "Same biosample"
    expr2 <- fread(file.path(data_dir, species, paste0(species, "_matching_celltype.csv")))
    expr2 <- expr2[order(accession_rna), ]
    expr2 <- unique(expr2, by = "accession_chip")[, c("accession_chip", "accession_rna")] # always the RNAseq with the earlier ERX id
    expr2$notes <- "Same cell type, distinct biosample"
    expr <- rbind(expr1, expr2)

    dfc <- merge(dfc, expr, by.x = "bound_id", by.y = "accession_chip", all.x = TRUE)
    dfc[is.na(accession_rna), c("accession_rna", "notes")] <- data.frame(
        accession_rna = min(dfc$accession_rna, na.rm=TRUE),
        notes = "Same species, different cell type and biosamples"
    )

    if (write_meta_table) {
        dff <- dfc
        dff$anchor_type <- "TSS"
        setnames(dff, "accession_rna", "scores")
        dff$species <- species
        dff$png <- file.path(species, "chipseq", "epistack", "TSS", paste0(dff$bound_id, ".png"))
        dff <- dff[, c("bound_id", "input_id", "anchors", "input_bw", "bound_bw", "cellType", "experiment", "anchor_type", "scores", "species", "png", "notes")]
        fwrite(
            dff,
            file = file.path(data_dir, species, "chipseq", "epistack", "list_of_plots_tss.tsv"),
            sep = "\t"
        )
    }

    # warnings are often because of pre-existing directories
    suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", "epistack")))
    suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", "epistack", "TSS")))


    commands <- with(dfc, paste0(
        "Rscript --vanilla ", EPISTACKPATH,
        " -a ", file.path(PREFIX, species, "chipseq", "epistack", "median_tss.tsv"),
        " -s ", file.path(PREFIX, species, "rnaseq", "salmon", accession_rna, "quant.genes.sf"),
        " -b ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", bound_bw),
        " -i ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", input_bw),
        " -p ", file.path(PREFIX, species, "chipseq", "epistack", "TSS", paste0(bound_id, ".png")),
        " -t '", paste(experiment, "ChIP-seq in", cellType),
        "' -r start --xlabs='-2.5kb,TSS,+2.5kb' -y 4 -z 2 -c 2 -v -g 5 -m 99999 -f ci95"
    ))

    commands
}

write_commands <- function(commands, species = "", by = 100) {
    from <- seq(1, length(commands), by = by)
    to <- from + by
    to[length(to)] <- length(commands)

    for(i in seq_along(from)) {
        header <- c(
            "#!/bin/bash",
            paste0("#SBATCH -J es_", species, "_", from[i]),
            paste0("#SBATCH -o epistack_", species, "_", from[i], ".out"),
            "#SBATCH --mem=30G",
            "#SBATCH -c 3",
            "#SBATCH --mail-type=END,FAIL",
            "#SBATCH -t 05:00:00",
            "module purge",
            "module load system/R-4.1.1_gcc-9.3.0"
        )

        write.table(
            c(header, commands[seq(from[i], to[i], by = 1L)]),
            file = paste0("epistack_tss_", species, "_", from[i], ".sh"),
            col.names = FALSE, row.names = FALSE, quote = FALSE
        )
    }
}


# running functions ---------------

pigptss <- make_commands("pig", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(pigptss, species = "pig", by = 25)

cowtss <- make_commands("cow", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(cowtss, species = "cow", by = 25)

horsetss <- make_commands("horse", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(horsetss, species = "horse", by = 25)

chicktss <- make_commands("chicken", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(chicktss, species = "chicken", by = 25)
