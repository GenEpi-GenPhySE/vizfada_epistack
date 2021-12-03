# paths: change accordingly ---------------
EPISTACKPATH <- "/home/gdevailly/R/x86_64-pc-linux-gnu-library/4.1/epistack/epistack.R"
PREFIX <- "/work/lmorel/data"

# EPISTACKPATH <- "~/mnt/inra_p/projets/vizfada/epistack.R"
# PREFIX <- "~/mnt/genotoul_lmoreldata"

# library loading ----------
library(purrr)
library(data.table)
library(stringr)

# function definitions ----------------
make_commands <- function(species, epistack_path, data_dir) {
    meta1 <- fread(file.path(PREFIX, species, paste0(species, "_matching_celltype.csv")))
    meta2 <- fread(file.path(PREFIX, species, paste0(species, "_matching_specimen.csv")))

    setnames(meta1, "cellType.text", "cellType")
    setnames(meta2, "cellType.text", "cellType")

    metatissue <- rbind(
        meta1[, c("accession_chip", "cellType")],
        meta2[, c("accession_chip", "cellType")]
    ) |> unique()

    metatarget <- fread(file.path(PREFIX, species, "chipseq", "metadata", "metadata.tsv"), select = c("experiment_accession", "experiment_target"))

    inputs <- list.files(file.path(PREFIX, species, "chipseq"), pattern = "^ERX")

    peaks <- map(set_names(inputs), function(x)
        list.files(
            file.path(PREFIX, species, "chipseq", x, "macs", "narrowPeak"),
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

    dfc <- merge(dfc, metatarget, by.x = "bound_id", by.y = "experiment_accession", all.x = TRUE)
    dfc <- merge(dfc, metatissue, by.x = "bound_id", by.y = "accession_chip"      , all.x = TRUE)

    walk(inputs, function(x) {
        # warnings are often because of pre-existing directories
        suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", x, "epistack")))
        suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", x, "epistack", "peaks")))
    })

    commands <- with(dfc, paste0(
        "Rscript --vanilla ", EPISTACKPATH,
        " -a ", file.path(PREFIX, species, "chipseq", input_id, "macs", "narrowPeak", peaks),
        " -b ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", bound_bw),
        " -i ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", input_bw),
        " -p ", file.path(PREFIX, species, "chipseq", input_id, "epistack", "peaks", paste0(bound_id, ".png")),
        # " -p ", file.path("~/mnt/inra_p/projets/vizfada/plots", paste0(bound_id, ".png")),
        " -t '", paste(experiment_target, "ChIP-seq in", cellType),
        "' -r center -y 10 -z 5 -c 2 -v -g 5"
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
            paste0("#SBATCH -J epistack_", species, "_", from[i]),
            paste0("#SBATCH -o epistack_", species, "_", from[i], ".out"),
            "#SBATCH --mem=16G",
            "#SBATCH -c 3",
            "#SBACTH --mail-type=END,FAIL",
            "module purge",
            "module load system/R-4.1.1_gcc-9.3.0"
        )

        write.table(
            c(header, commands[seq(from[i], to[i], by = 1L)]),
            file = paste0("epistack_", species, "_", from[i], ".sh"),
            col.names = FALSE, row.names = FALSE, quote = FALSE
        )
    }
}

# running functions ---------------

species <- "cow"
cowpeaks <- make_commands("cow", EPISTACKPATH, PREFIX)
write_commands(cowpeaks, species = "cow")

