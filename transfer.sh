## script to be run on WCHG for synchronize data between three locations

## make directory on /well
mkdir /well/myers/rwdavies/chicago_oxford_2019_08_07

## reference data
rsync --progress -av sparse:/Net/dense/data/outbredmice/refs/mm10/mm10.dict /well/myers/rwdavies/chicago_oxford_2019_08_07/

## imputed data
rsync --progress -av sparse:/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/ /well/myers/rwdavies/chicago_oxford_2019_08_07/

## now transfer from /well to ucsc
rsync -av --progress /well/myers/rwdavies/chicago_oxford_2019_08_07/* ucsc:/projects/ps-palmer/robbie/chicago_oxford_2019_08_07/

