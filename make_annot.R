outputdir <- "/well/myers/rwdavies/chicago_oxford_2019_08_07/"
chr <- "chr19"
ref <- file.path(outputdir, "mm10.fa")

## stitch most recent commit 94271d3d138b96f0c5718b7f1c11c88c12995057
source("~/proj/STITCH/STITCH/R/extract.R")
source("~/proj/STITCH/STITCH/R/writers.R")
source("~/proj/STITCH/STITCH/R/functions.R")
library("STITCH")

## do this process in chunks!
chrlist <- paste0("chr", c(19, 1:18, "X"))
for(chr in chrlist) {
    load(file.path(outputdir, "RData", paste0("EM.all.", chr, ".RData")))
    load(file.path(outputdir, "RData", paste0("sampleNames.", chr, ".RData")))    
    nCores <- 8
    ## break up into chunks
    parts <- getSampleRange(nrow(pos), round(nrow(pos) / 1000))
    out <- mclapply(1:length(parts), mc.cores = nCores, function(ipart) {
        part <- parts[[ipart]]
        print(paste0("Extracting part:", ipart, " out of ", length(parts)))
        pos_local <- pos[part[1]:part[2], ]
        vcf_file <- file.path(outputdir, paste0("stitch.", chr, ".vcf.gz"))        
        cube <- extract_hd_to_cube(
            vcf_file = vcf_file,
            ref = ref,
            bcftools = "/gpfs0/users/flint/rwdavies/personal/bin/bcftools-1.6/bcftools",
            gatk_jar = "/apps/well/gatk/3.7-0/GenomeAnalysisTK.jar",
            samples = NULL,
            pos = pos_local,
            ram = "-Xmx8g",
            verbose = TRUE,
            field = "GP"
        )
        annot <- get_annot(pos_local, cube, sampleNames, nCores = 2)
        return(annot)
    })
    annot <- data.frame(data.table::rbindlist(out))
    file <- paste0("/well/myers/rwdavies/chicago_oxford_2019_08_07/annot.", chr, ".csv")
    write.table(
        annot,
        file = file,
        sep = ",",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    system(paste0("gzip -1 -f ", file))
    ## shockingly similar? and fine?
    print("compare info")
    print(table(round(annot[, "info.oxford"], 1), round(annot[, "info.chicago"], 1)))
    ## pretty close!
    print("passing QC")
    print(table(annot[, "info.both"] > 0.4 & annot[, "hwe.both"] > 1e-6))
    ##
    print("colmeans")
    print(colMeans(annot[, -c(1:4)], na.rm = TRUE))
    gc(reset = TRUE);gc(reset = TRUE);gc(reset = TRUE);
}
