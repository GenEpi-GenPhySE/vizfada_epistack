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
        file.path(PREFIX, species, "chipseq", "metadata"),
        pattern = "\\.tsv$", full.names = TRUE
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
    
    dfc <- merge(dfc, metatarget, by.x = "bound_id", by.y = "accession", all.x = TRUE)
    
    expr1 <- fread(file.path(PREFIX, species, paste0(species, "_matching_specimen.csv")))
    expr2 <- fread(file.path(PREFIX, species, paste0(species, "_matching_celltype.csv")))
    
    
    if (write_meta_table) {
        dff <- dfc
        colnames(dff) <- c("bound_id", "input_id", "anchors", "input_bw", "bound_bw", "cellType", "experiment")
        dff$anchor_type <- "peak centers"
        dff$png <- file.path(species, "chipseq", "epistack", "peaks", paste0(dff$bound_id, ".png"))
        fwrite(
            dff,
            file = file.path(data_dir, species, "chipseq", "epistack", "list_of_plots_peaks.tsv"),
            sep = "\t"
        )
    }
    
    # warnings are often because of pre-existing directories
    suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", "epistack")))
    suppressWarnings(dir.create(file.path(PREFIX, species, "chipseq", "epistack", "peaks")))
    
    
    commands <- with(dfc, paste0(
        "Rscript --vanilla ", EPISTACKPATH,
        " -a ", file.path(PREFIX, species, "chipseq", input_id, "macs", "narrowPeak", peaks),
        " -b ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", bound_bw),
        " -i ", file.path(PREFIX, species, "chipseq", input_id, "bigwig", input_bw),
        " -p ", file.path(PREFIX, species, "chipseq", "epistack", "peaks", paste0(bound_id, ".png")),
        " -t '", paste(experiment, "ChIP-seq in", cellType),
        "' -r center -y 10 -z 5 -c 2 -v -g 5 -m 99999 -f ci95"
    ))
    
    commands
}
