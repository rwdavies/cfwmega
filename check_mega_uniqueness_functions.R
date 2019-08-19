hack_EMall_RData <- function(outputdir, test_outputdir) {
    dir.create(paste0(test_outputdir, "/RData/"))
    system(paste0("rsync -av ", outputdir, "/RData/*chr19* ", test_outputdir, "/RData/"))
    hack_file <- file.path(test_outputdir, "/RData/EM.all.chr19.RData")
    load(hack_file)
    gridWindowSize <- NA
    ref_alleleCount <- alleleCount ## not really, but not needed
    ## cat(ls(), sep = ", ")
save(K, L, L_grid, N, alleleCount, alphaMatCurrent_t, chr, eHapsCurrent_t, estimatedAlleleFrequency, gen, gen_imp, grid, grid_distances, hapSumCurrent_t, hwe, hweCount, inRegionL, info, nCores, nGen, nGrids, nSNPs, niterations, passQC, pos, priorCurrent, reference_panel_SNPs, sigmaCurrent, snps_in_grid_1_based, start_and_end_minus_buffer, gridWindowSize, ref_alleleCount, file = hack_file)
}

check_read_uniqueness_and_maybe_remake <- function(
    iSample,
    save_outputdir = NULL,
    saveNew = FALSE,
    downsampleOxford = FALSE
) {
    if ((iSample %% 100) == 0) {
        print(paste0(iSample, ", ", date()))
    }
    ## 
    load(paste0(outputdir, "/input/sample.", iSample, ".input.", chr, ".RData"))
    a <- sapply(sampleReads, function(x) {
        paste0(sapply(x, paste0, collapse = "-"), collapse = "_")
    })
    ## base_positions <- unlist(sapply(sampleReads, function(x) x[[4]]))
    if (saveNew) {
        sampleReads <- sampleReads[match(unique(a), a)]
        if ((downsampleOxford) & (iSample <= 2073)) {
            ## downsample to approximately match
            1
        }
        save(sampleReads, file = paste0(save_outputdir, "/input/sample.", iSample, ".input.", chr, ".RData"))
        return(NULL)
    } else {
        return(
            list(
                nReads = length(sampleReads),
                nUniqueReads = length(unique(a)),
                nUniqueBases = sum(sapply(sampleReads[match(unique(a), a)], function(x) x[[1]]) + 1)
            )
        )
    }
    ## base_positions = base_positions
    ## use only unique reads
}


get_pileups <- function(
    outputdir,
    chr,
    nCores,
    input = "input"
) {
    bundling_info <- get_bundling_position_information(
        N = 2073,
        nCores = nCores,
        blockSize = NA
    )
    alleleCountOx <- buildAlleleCount(
        nSNPs = nSNPs,
        N = 2073,
        nCores = nCores,
        regionName = chr,
        tempdir = paste0(file.path(outputdir, input), "/"),
        bundling_info = bundling_info,
        allSampleReads = NULL
    )
    bundling_info <- get_bundling_position_information(
        N = 3234,
        nCores = nCores,
        blockSize = NA
    )
    alleleCountAll <- buildAlleleCount(
        nSNPs = nSNPs,
        N = 3234,
        nCores = nCores,
        regionName = chr,
        tempdir = paste0(file.path(outputdir, input), "/"),
        bundling_info = bundling_info,
        allSampleReads = NULL
    )
    alleleCountChi <- alleleCountAll - alleleCountOx
    return(
        list(
            ac.o = alleleCountOx,
            ac.c = alleleCountChi,
            ac.a = alleleCountAll
        )
    )
}

