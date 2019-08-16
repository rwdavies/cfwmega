## calculate accuracy of CFW oxford and chicago
## stratify by QC metrics

##
## 1 - get data using existing script
## 
chrlist <- paste0("chr", c(1:19, "X"))
parallel::mclapply(chrlist, mc.cores = 8, function(chr) {
    parallel::mclapply(c("oxford", "chicago"), mc.cores = 2, function(whose_samples) {
        comparison_save_file <- paste0("/well/myers/rwdavies/chicago_oxford_2019_08_07/acc.data.", chr, ".", whose_samples, ".RData")
        command <- paste0(
            "cd ~/proj/STITCH && ./scripts/compare_vcf_to_truth.R ",
            "--test-file=/well/myers/rwdavies/chicago_oxford_2019_08_07/stitch.", chr, ".vcf.gz ",
            "--chr=", chr, " ",
            "--compare-against=megamuga ",
            "--verbose ",
            "--whose-samples=", whose_samples, " ",
            "--cfw-data-dir=/well/myers/rwdavies/cfw_truth/ ",
            "--mega-save-file=/well/myers/rwdavies/chicago_oxford_2019_08_07/mega.", chr, ".", whose_samples, ".RData ",
            "--comparison-save-file=", comparison_save_file
        )
        ##if (!file.exists(comparison_save_file)) {
            system(command)
        ##}
    })
})

quit()

##
## load data back in, perform bespoke comparisons
##
setwd("/well/myers/rwdavies/chicago_oxford_2019_08_07/")
chrlist <- paste0("chr", c(1:19, "X"))
all <- parallel::mclapply(chrlist, mc.cores = 8, function(chr) {
    ## load oxford data
    load(paste0("acc.data.", chr, ".oxford.RData"))
    outO <- out
    ## load chicago data
    load(paste0("acc.data.", chr, ".chicago.RData"))
    outC <- out
    ## also load annot
    annot <- data.table::fread(cmd = paste0("gunzip -c /well/myers/rwdavies/chicago_oxford_2019_08_07/annot.", chr, ".csv.gz"), data.table = FALSE)
    rownames(annot) <- paste(annot[, 1], annot[, 2], annot[, 3], annot[, 4], sep = "-")
    ## argh - intersect everything
    o <- rownames(outO[["callsS"]])
    c <- rownames(outC[["callsS"]])
    both <- intersect(o, c)
    r2.o <- outO[["r2"]][match(both, o)]
    r2.c <- outC[["r2"]][match(both, c)]
    annot <- annot[match(both, rownames(annot)), ]
    if (nrow(annot) != length(r2.o)) {
        stop(chr)
    }
    if (nrow(annot) != length(r2.c)) {
        stop(chr)
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

## callsO <- do.call("rbind", lapply(all, function(x) x[["callsO"]]))
## dosagesO <- do.call("rbind", lapply(all, function(x) x[["dosagesO"]]))
## callsC <- do.call("rbind", lapply(all, function(x) x[["callsC"]]))
## dosagesC <- do.call("rbind", lapply(all, function(x) x[["dosagesC"]]))
## merge
r2.o <- unlist(lapply(all, function(x) x[["r2.o"]]))
r2.c <- unlist(lapply(all, function(x) x[["r2.c"]]))
annot <- do.call("rbind", lapply(all, function(x) x[["annot"]]))
## apply blanket call rate
keep <-
    (annot[, "eaf.oxford"] > 0.01) & (annot[, "eaf.oxford"] < 0.99) &
    (annot[, "eaf.chicago"] > 0.01) & (annot[, "eaf.chicago"] < 0.99) &
    (annot[, "eaf.both"] > 0.01) & (annot[, "eaf.both"] < 0.99)
r2.o <- r2.o[keep]
r2.c <- r2.c[keep]
annot <- annot[keep, ]
min_hwe <- 1e-20
annot[annot < min_hwe] <- min_hwe

## now - can make accuracy plots for different groups!

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

    


## accuracy in respective-info score deciles
nRows <- 5
alpha <- 0.25
blue <- rgb(red = 0, green = 0, blue = 1, alpha = alpha)
red <- rgb(red = 1, green = 0, blue = 0, alpha = alpha)
green <- rgb(red = 0, green = 1, blue = 0, alpha = alpha)

for(what_to_plot in c("info", "log_hwe")) {
    if (what_to_plot == "info") {
        cutoffs <- seq(0, 1, length.out = 21)
        cutoffs_mids <- cutoffs[-length(cutoffs)] + diff(cutoffs) /2 
    } else {
        cutoffs <- seq(log10(min_hwe), 0, length.out = 21)
        cutoffs_mids <- cutoffs[-length(cutoffs)] + diff(cutoffs) /2 
    }
    print(paste0("------plot ", what_to_plot, "-------"))
    plotfile <- paste0("/well/myers/rwdavies/chicago_oxford_2019_08_07/", what_to_plot, ".png")
    png(plotfile, height = nRows * 3, width = 9, units = "in", res = 200)
    par(mfcol = c(nRows, 3))
    for(i_which in 1:3) {
        ## give histogram of info
        which_group <- c("oxford", "chicago", "both")[i_which]
        col <- c(blue, red, green)[i_which]
        if (what_to_plot == "info") {
            measure <- annot[, paste0("info.", which_group)]
        } else if (what_to_plot == "log_hwe") {
            measure <- log10(annot[, paste0("hwe.", which_group)])
        }
        hist(measure, col = col, main = paste0("Histogram of ", what_to_plot, "\nfrom ", which_info, " samples"))
        ##
        w <- seq(1, length(measure), length.out = 500)
        plot(sort(measure)[w], ((length(measure):1)/length(measure))[w], col = col, main = "1 - ECDF", type = "l")
        ## accuracy by measure
        plot(measure, r2.o, col = blue, xlab = "Measure", ylab = "r2", pch = 16)
        points(measure, r2.c, col = red, pch = 16)
        legend("bottomleft", c("oxford", "chicago"), col = c(blue, red), lwd = 2)
        ## specific accuracy
        s.o <- get_specific_average(measure, r2.o, cutoffs)
        s.c <- get_specific_average(measure, r2.c, cutoffs)
        ## plot oxford
        plot(cutoffs_mids, s.o[, 2], xlab = what_to_plot, ylab = "Average r2", col = "blue", ylim = c(0, 1), main = "Specific accuracy")
        points(cutoffs_mids, s.c[, 2], col = "red")
        ## no need to plot N - that is basically the histogram
        ## cumulative accuracy of points below that level
        m.o <- get_cumulative_average(measure, r2.o, cutoffs)
        m.c <- get_cumulative_average(measure, r2.c, cutoffs)
        plot(cutoffs, m.o[, 2], xlab = paste0(what_to_plot, " >= cutoff"), ylab = "Average r2", col = "blue", ylim = c(0, 1), main = "Cumulative accuracy")
        points(cutoffs, m.c[, 2], col = "red")
        ## hwe
    }
    dev.off()
}
