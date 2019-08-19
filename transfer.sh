## script to be run on WCHG for synchronize data between four locations

## make directory on /well
mkdir /well/myers/rwdavies/chicago_oxford_2019_08_07

## reference data
rsync --progress -av sparse:/Net/dense/data/outbredmice/refs/mm10/mm10.dict /well/myers/rwdavies/chicago_oxford_2019_08_07/

## imputed data
rsync --exclude=input --progress -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/ /well/myers/rwdavies/chicago_oxford_2019_08_07/

## get posfile, genfile
rsync -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/sites /well/myers/rwdavies/chicago_oxford_2019_08_07/

## now transfer from /well to ucsc
rsync -av --progress /well/myers/rwdavies/chicago_oxford_2019_08_07/* ucsc:/projects/ps-palmer/robbie/chicago_oxford_2019_08_07/

## get truth data from smew
mkdir /well/myers/rwdavies/cfw_truth
rsync --progress -av smew:/data/smew1/rdavies/stitch_development/truth/cfw/megamuga_2018_05_18.tgz /well/myers/rwdavies/cfw_truth/
tar -xzvf megamuga_2018_05_18.tgz

##
mkdir /well/myers/rwdavies/chicago_oxford_MarkDuplicatesBams_2019_08_07/
rsync -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/mark_duplicate_bams/*  /well/myers/rwdavies/chicago_oxford_MarkDuplicatesBams_2019_08_07/

##
## scratch
## 

## copy chromosome 19 input files
## rsync --progress -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/input/*chr19* /well/myers/rwdavies/chicago_oxford_2019_08_07/input/



## rsync --progress -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/input/*chr19* /well/myers/rwdavies/chicago_oxford_2019_08_07/input/
