. ~/.conda_activate
source activate

chrlist="16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 X"
chrlist="12 1 2"
for chr in `echo ${chrlist}`
do
    echo --- ${chr} --- 
    date
    log_file=/Net/dense/data/outbredmice/imputation/ancestral_haps/chicago_oxford_2019_08_07/log.${chr}.txt
    R -f ~/proj/outbred/chicago_oxford.R --args ${chr} &> ${log_file}
    echo --- done ${chr}
done