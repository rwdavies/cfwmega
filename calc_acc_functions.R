calculate_accuracy_for_vcfs <- function(
   chrlist,
   nCores,
   outputdir
) {
    mclapply(chrlist, mc.cores = min(1, floor(nCores / 2)), function(chr) {
        mclapply(c("oxford", "chicago"), mc.cores = 2, function(whose_samples) {
            comparison_save_file <- paste0(
                outputdir,
                "acc.data.", chr, ".", whose_samples, ".RData"
            )
        command <- paste0(
            "cd ~/proj/STITCH && ./scripts/compare_vcf_to_truth.R ",
            "--test-file=", outputdir, "/stitch.", chr, ".vcf.gz ",
            "--chr=", chr, " ",
            "--compare-against=megamuga ",
            "--verbose ",
            "--whose-samples=", whose_samples, " ",
            "--cfw-data-dir=/well/myers/rwdavies/cfw_truth/ ",
            "--mega-save-file=", outputdir, "/mega.", chr, ".", whose_samples, ".RData ",
            "--comparison-save-file=", comparison_save_file
        )
        ##if (!file.exists(comparison_save_file)) {
            system(command)
        ##}
        })
    })
    return(NULL)
}


get_specific_average <- function(info, r2, cutoffs) {
    t(sapply(1:(length(cutoffs) - 1), function(i) {
        s <- cutoffs[i]
        e <- cutoffs[i + 1]        
        w <- info >= s & info <= e
        c(sum(w), mean(r2[w], na.rm = TRUE))
    }))
}

get_cumulative_average <- function(info, r2, cutoffs) {
    t(sapply(cutoffs, function(c) {
        c(sum(info >= c), mean(r2[info >= c], na.rm = TRUE))
    }))
}



load_accuracy_results <- function(
   chrlist,
   nCores,
   outputdir
) {
    ##
    ## load data back in, perform bespoke comparisons
    ##
    setwd(outputdir)
    all <- parallel::mclapply(chrlist, mc.cores = nCores, function(chr) {
        ## load oxford data
        load(paste0("acc.data.", chr, ".oxford.RData"))
        outO <- out
        ## load chicago data
        load(paste0("acc.data.", chr, ".chicago.RData"))
        outC <- out
        ## 
        o <- rownames(outO[["callsS"]])
        c <- rownames(outC[["callsS"]])
        both <- intersect(o, c)
        r2.o <- outO[["r2"]][match(both, o)]
        r2.c <- outC[["r2"]][match(both, c)]
        # also load annot
        annot_file <- paste0(outputdir, "/annot.", chr, ".csv.gz")
        if (!file.exists(annot_file)) {
            annot <- NULL
        } else {
            annot <- data.table::fread(cmd = paste0("gunzip -c ", annot_file), data.table = FALSE)
            rownames(annot) <- paste(annot[, 1], annot[, 2], annot[, 3], annot[, 4], sep = "-")
            ## argh - intersect everything
            annot <- annot[match(both, rownames(annot)), ]
            if (nrow(annot) != length(r2.o)) {
                stop(chr)
            }
            if (nrow(annot) != length(r2.c)) {
                stop(chr)
            }
        }
        return(
            list(
                r2.o = r2.o,
                r2.c = r2.c,
                annot = annot
            )
        )
        ## callsO = outO[["callsS"]],
        ## dosagesO = outO[["dosagesS"]],
        ## callsC = outC[["callsS"]],
        ## dosagesC = outC[["dosagesS"]],
    })
    return(all)
}
