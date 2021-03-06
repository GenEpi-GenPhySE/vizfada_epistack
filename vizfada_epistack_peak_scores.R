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

    # warnings are often because of pre-existing directories
    suppressWarnings(dir.create(file.path(data_dir, species, "chipseq", "epistack")))
    suppressWarnings(dir.create(file.path(data_dir, species, "chipseq", "epistack", "peaks")))
    
    if (write_meta_table) {
        dff <- dfc
        colnames(dff) <- c("bound_id", "input_id", "anchors", "input_bw", "bound_bw", "cellType", "experiment")
        dff$anchor_type <- "peak centers"
        dff$scores <- dff$anchors
        dff$species <- species
        dff$png <- file.path(species, "chipseq", "epistack", "peaks", paste0(dff$bound_id, ".png"))
        dff$notes <- ""
        if(species == "chicken") {
            # TODO if new annotations / pipeline version: check if still 0 peaks detected
            dff[bound_id %in% c("ERX3783265", "ERX3783279"), "notes"] <- "No peaks detected."
        }
        fwrite(
            dff,
            file = file.path(data_dir, species, "chipseq", "epistack", "list_of_plots_peaks.tsv"),
            sep = "\t"
        )
    }

    commands <- with(dfc, paste0(
        "Rscript --vanilla ", EPISTACKPATH,
        " -a ", file.path(data_dir, species, "chipseq", input_id, "macs", "narrowPeak", peaks),
        " -b ", file.path(data_dir, species, "chipseq", input_id, "bigwig", bound_bw),
        " -i ", file.path(data_dir, species, "chipseq", input_id, "bigwig", input_bw),
        " -p ", file.path(data_dir, species, "chipseq", "epistack", "peaks", paste0(bound_id, ".png")),
        " -t '", paste(experiment, "ChIP-seq in", cellType),
        "' -r center -y 10 -z 5 -c 2 -v -g 5 -m 99999 -f ci95"
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
            "#SBATCH --mem=40G",
            "#SBATCH -c 3",
            "#SBATCH --mail-type=END,FAIL",
            "#SBATCH -t 05:00:00",
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

# TODO rerun when cow metadat are reformated
cowpeaks <- make_commands("cow", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(cowpeaks, species = "cow", by = 25)

pigpeaks <- make_commands("pig", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(pigpeaks, species = "pig", by = 25)

chickpeaks <- make_commands("chicken", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(chickpeaks, species = "chicken", by = 25)

horsepeaks <- make_commands("horse", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(horsepeaks, species = "horse", by = 25)

sheeppeaks <- make_commands("sheep", EPISTACKPATH, PREFIX, write_meta_table = TRUE)
write_commands(sheeppeaks, species = "sheep", by = 25)

