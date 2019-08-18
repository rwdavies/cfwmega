## before hand run
## ~/.conda_activate
## source activate

library("STITCH")
chr <- paste0("chr", commandArgs(trailingOnly = TRUE)[1])
print(paste0("Using chr:", chr))
STITCH(
  outputdir = "/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/",
  nCores = 48,
  K = 4,
  chr = chr,
  posfile = paste0("/Net/dense/data/outbredmice/imputation/ancestral_haps/sites/pos.", chr, ".txt"),
  genfile = paste0("/Net/dense/data/outbredmice/imputation/ancestral_haps/sites/gen.", chr, ".txt"),  
  bamlist = "/Net/dense/data/outbredmice/imputation/ancestral_haps/sites/both.bamlist.txt",
  nGen = 100,
  tempdir = "/dev/shm/rwdavies/"
)


quit()

## for a chr, calculate, per-cohort
## estimatedAlleleFrequency, hwe, info
outputdir <- "/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/"
chr <- "chr19"
extract_hd_to_cube(
    vcf_file,
    ref,
    bcftools = "bcftools",
    gatk_jar = "/data/smew1/rdavies/stitch_richard_paper/bin/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar",
    samples = NULL,
    pos = NULL,
    ram = "-Xmx4g",
    verbose = TRUE,
    field = "GP"
)
