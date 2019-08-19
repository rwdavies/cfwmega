source("~/proj/STITCH/STITCH/R/filenames.R")
source("~/proj/STITCH/STITCH/R/sample-reads.R")
source("~/proj/STITCH/STITCH/R/functions.R")
library("STITCH")
library("parallel")
source("~/proj/cfwmega/calc_acc_functions.R")
source("~/proj/cfwmega/check_mega_uniqueness_functions.R")
chrlist <- paste0("chr", 19)
chr <- "chr19"
nCores <- 16
nSNPs <- 152486
## outputdir is normal
## test is with both having unique reads only
## test2 is with mark duplicates
outputdir <- "/well/myers/rwdavies/chicago_oxford_2019_08_07/"
test_outputdir <- "/well/myers/rwdavies/chicago_oxford_test_2019_08_07/"
test_mark_outputdir <- "/well/myers/rwdavies/chicago_oxford_test_mark_2019_08_07/"
dir.create(file.path(test_outputdir, "input"), recursive = TRUE)
dir.create(file.path(test_outputdir, "RData"), recursive = TRUE)
dir.create(file.path(test_mark_outputdir, "input"), recursive = TRUE)
dir.create(file.path(test_mark_outputdir, "RData"), recursive = TRUE)
bamlist <- as.character(read.table("/well/myers/rwdavies/chicago_oxford_2019_08_07/sites/both.bamlist.txt")[, 1])
test_mark_bamlist <- file.path("/well/myers/rwdavies/chicago_oxford_test_mark_2019_08_07/both.bamlist.markduplicates.txt")
write.table(matrix(file.path("/well/myers/rwdavies/chicago_oxford_MarkDuplicatesBams_2019_08_07/", basename(bamlist)), ncol = 1), file = test_mark_bamlist, row.names = FALSE, col.names = FALSE, sep = "", quote = FALSE)


## install.packages("/users/flint/rwdavies/proj/STITCH-private/releases/STITCH_1.5.7.0001.tar.gz")



##
## check how unique they are, in worst case scenario
##
nSamples <- 3234
out <- mclapply(1:nSamples, mc.cores = 16, check_read_uniqueness_and_maybe_remake)
results <- t(sapply(out, function(x) as.numeric(x)))

## take a look
colMeans(results[1:2073, ]) ## Oxford - much longer reads!
colMeans(results[-(1:2073), ]) ## Chicago

##
## prepare hacked input
##
out <- mclapply(1:3234, mc.cores = 16, f, saveNew = TRUE, save_outputdir = test_outputdir)
out <- mclapply(1:3234, mc.cores = 16, f, saveNew = TRUE, save_outputdir = test2_outputdir, downsampleOxford = TRUE)




##
## run STITCH-again on the samples
##

hack_EMall_RData(outputdir, test_outputdir)
hack_EMall_RData(outputdir, test_mark_outputdir)

## run with hacky duplicate read removal
STITCH_again(
    outputdir = test_outputdir,
    originalRegionName = chr, chr = chr, regionStart = NA, regionEnd = NA, buffer = NA, nCores = 16,  posfile = file.path(outputdir, "sites", paste0("pos.", chr, ".txt")), method = "diploid",
    regenerateInput = FALSE,
    regenerateInputWithDefaultValues = TRUE,
    outputSNPBlockSize = 20000 ## juice up for next time
)

calculate_accuracy_for_vcfs(chrlist = chrlist, nCores = nCores, outputdir = test_outputdir)


## run with formal duplicate read marking
STITCH_again(
    outputdir = test_mark_outputdir,
    originalRegionName = chr, chr = chr, regionStart = NA, regionEnd = NA, buffer = NA, nCores = 16,  posfile = file.path(outputdir, "sites", paste0("pos.", chr, ".txt")), method = "diploid",
    regenerateInputWithDefaultValues = TRUE,
    bamlist = test_mark_bamlist,
    outputSNPBlockSize = 20000, ## juice up for next time
    keep_generated_input = TRUE
)

calculate_accuracy_for_vcfs(chrlist = chrlist, nCores = nCores, outputdir = test_mark_outputdir)











##
##
## second part check out the accuracy
##
##

## compare the pileups
pileO <- get_pileups(outputdir, chr, nCores)
pileT <- get_pileups(test_outputdir, chr, nCores)
pileM <- get_pileups(test_mark_outputdir, chr, nCores,input ="input_again")
t(sapply(pileO, colMeans))
t(sapply(pileT, colMeans))
t(sapply(pileM, colMeans))

##
f <- function(x, N) {
    round(100 * sum(x[, 2] > (0.1 * N)) / nrow(x), 2)
}
f(pileO[["ac.o"]], 2073) ## 85% of SNPs have >0.1X
f(pileO[["ac.c"]], 1161) ## 2% of SNPs have >0.1X
f(pileT[["ac.o"]], 2073)
f(pileT[["ac.c"]], 1161)
f(pileM[["ac.o"]], 2073)
f(pileM[["ac.c"]], 1161)




allO <- load_accuracy_results(chrlist = "chr19", nCores = 1, outputdir = outputdir)
allT <- load_accuracy_results(chrlist = "chr19", nCores = 1, outputdir = test_outputdir)
allM <- load_accuracy_results(chrlist = "chr19", nCores = 1, outputdir = test_mark_outputdir)
## basically the same
sapply(allO[[1]], mean, na.rm = TRUE)
sapply(allT[[1]], mean, na.rm = TRUE)
sapply(allM[[1]], mean, na.rm = TRUE)



annot <- read.table("/well/myers/rwdavies/chicago_oxford_2019_08_07/annot.chr19.csv.gz", header = TRUE, sep = ",")

keep <-
    (rowSums(allO[[1]][["annot"]][, c("info.oxford", "info.chicago", "info.both")] > 0.4) == 3) &
    (rowSums(allO[[1]][["annot"]][, c("hwe.oxford", "hwe.chicago", "hwe.both")] > 1e-10) == 3)
f <- function(pile, which, col, N) {
    ## smooth here!
    y <- pile[[which]][, 2] / N
    y2 <- frollmean(y, n = 1)
    y2[y2 > ylim[2]] <- ylim[2]
    points(L2, y2, col = col, pch = 16)
}
L <- allO[[1]][["annot"]][keep, "POS"] / 1e6
L2 <- annot[, "POS"] / 1e6
xlim <- range(L)
ylim <- c(0, 1)
ylim2 <- c(0, 1)
## check locations?

##
## make plot
##
png(paste0(outputdir, "acc.markduplicates.png"), height = 20, width = 20, units = "in", res = 100)
alpha <- 0.25
blue <- rgb(red = 0, green = 0, blue = 1, alpha = alpha)
red <- rgb(red = 1, green = 0, blue = 0, alpha = alpha)
par(mfrow= c(6, 1))
for(i_to_plot in 1:3) {
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "Position", ylab = "r2", main = c("Original", "hacky-mark-duplicates", "Picard-mark-duplicates")[i_to_plot], cex.main = 2)
    all <- list(allO, allT, allM)[[i_to_plot]]
    points(x = L, y = all[[1]][["r2.o"]][keep], col = "blue")
    points(x = L, y = all[[1]][["r2.c"]][keep], col = "red")
    ## pileup
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim2, xlab = "Position", ylab = "Average Depth")
    pile <- list(pileO, pileT, pileM)[[i_to_plot]]
    f(pile, "ac.o", blue, 2073)
    f(pile, "ac.c", red, 1161)
}
dev.off()
