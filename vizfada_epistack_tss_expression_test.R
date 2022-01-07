# interactive tests ---------------
opt <- list()
opt$anchors <- "/work/lmorel/data/pig/chipseq/epistack/median_tss.tsv"
opt$bound <- "/work/lmorel/data/pig/chipseq/ERX3212574/bigwig/ERX3212576_R1.bigWig"
opt$input <- "/work/lmorel/data/pig/chipseq/ERX3212574/bigwig/ERX3212574_R1.bigWig"
opt$scores <- "/work/lmorel/data/pig/rnaseq/salmon/ERX3212540/quant.genes.sf"
opt$cpu <- 2L
opt$reference <- "start"
opt$extend <- 2500L
opt$group <- 5L
opt$title <- "TSS sorted by expr"
opt$xlabs <- c("-2.5kb", "TSS", "+2.5kb")
opt$ylim <- 1
opt$zlim <- 1
opt$errortype <- "ci95"
opt$png <- "~/work/epistack_test.png"

# table(anchors$gene_type) |> sort(decreasing = T)
#
# protein_coding               lncRNA           pseudogene
# 21280                 6790                 1333
# snRNA               snoRNA                miRNA
# 1104                  591                  388
# processed_pseudogene               scaRNA              Mt_tRNA
# 293                   25                   22
# rRNA             misc_RNA            IG_V_gene
# 21                   17                   10
# Y_RNA            TR_V_gene             ribozyme
# 9                    8                    7
# TR_J_gene            vault_RNA              Mt_rRNA
# 4                    3                    2
# IG_C_gene
# 1

# Rscript tests ---------------

Rscript --vanilla ~/work/projects/vizfada/epistack.R \
    -a /work/lmorel/data/pig/chipseq/epistack/median_tss.tsv \
    -s /work/lmorel/data/pig/rnaseq/salmon/ERX3212540/quant.genes.sf \
    -b /work/lmorel/data/pig/chipseq/ERX3212574/bigwig/ERX3212576_R1.bigWig \
    -i /work/lmorel/data/pig/chipseq/ERX3212574/bigwig/ERX3212574_R1.bigWig \
    -p ~/work/epistack_test_cli.png \
    -t 'H3K27me3 I guess' \
    -r start \
    -y 0.5 -z 0.5 -c 2 -v -g 5 -m 99999 -f ci95
