get_annot <- function(pos, cube, sampleNames, nCores = 16) {
    oxford_samples <- grep("CFW", sampleNames)
    oxford_cube <- cube[, oxford_samples, ]
    chicago_cube <- cube[, -oxford_samples, ]
    sampleRanges <- getSampleRange(nrow(pos), nCores)
    out <- parallel::mclapply(sampleRanges, mc.cores = nCores, function(sampleRange) {
        ##
    ## OK, per-variant, calculate metrics!
    annos <- c("info", "hwe", "eaf")
    anno_cols <- paste0(rep(annos, each = 3), ".", c("oxford", "chicago", "both"))
    ##
    annot <- pos[(sampleRange[1]):sampleRange[2], ]
    for(col in anno_cols) {
        annot[, col] <- NA
    }
    for(ii_snp in 1:(sampleRange[2] - sampleRange[1] + 1)) {
        i_snp <- sampleRange[1] + ii_snp - 1
        for(j in 1:3) {
            ##
            if (j == 1) {  gp <- oxford_cube[i_snp, , ]}
            if (j == 2) {  gp <- chicago_cube[i_snp, , ]}
            if (j == 3) {  gp <- cube[i_snp, , ]}
            N <- c(2073, 1161, 3234)[j]
            who <- c("oxford", "chicago", "both")[j]
            ##
            ## Gi \in {0, 1, 2}
            ## pijk = P(Gij = k|H, G)
            ## eij = p_{ij1} + 2 *{p_{ij2} and define fij = p_{ij1} + 4_{pij2}
            eij <- round(gp[, 2] + 2 * gp[, 3], 3)
            fij <- round(gp[, 2] + 4 * gp[, 3], 3)
            thetaHat <- sum(eij) / 2 / N
            info <- 1 - sum(fij - eij ** 2) / (2 * N * thetaHat * (1 - thetaHat))
            info[(round(thetaHat, 2) == 0 ) | (round(thetaHat, 2) == 1)] <- 1
            annot[ii_snp, paste0("info.", who)] <- info
            ##
            max_gen <- get_max_gen_rapid(t(gp))
            x <- c(sum(max_gen[, 2] == 1), sum(max_gen[, 2] == 2), sum(max_gen[, 2] == 3))
            annot[ii_snp, paste0("hwe.", who)] <- rcpp_calculate_hwe_p(x)
            ## freq
            annot[ii_snp, paste0("eaf.", who)] <- sum(0.5 * gp[, 2] + gp[, 3]) / N
        }
    }
    return(annot)
})
annot <- data.frame(data.table::rbindlist(out))
return(annot)
}


get_max_gen_rapid <- function(x) {
    ## assume matrix 3 columns >1 row
    z <- rep(1, ncol(x))
    y <- x[1, ]
    for(i in 2:3) {
        w <- x[i, ] > y
        z[w] <- i
        y[w] <- x[i, w]
    }
    return(cbind(1:ncol(x), z))
}


