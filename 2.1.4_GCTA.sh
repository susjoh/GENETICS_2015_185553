#!/bin/sh

# Run GCTA to extract genomic related matrices and inbreeding coefficients.
# Author: Susan Johnston


#### Get inbreeding coefficients

./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --autosome --ibc --out 20160205_autosomal_IBD

#### Create a script for all

echo '#!/bin/sh'                   >  all.sh
echo ''                            >> all.sh
echo '#$ -cwd'                     >> all.sh
echo '#$ -l h_rt=2:00:00'          >> all.sh
echo ''                            >> all.sh
echo '. /etc/profile.d/modules.sh' >> all.sh
echo ''                            >> all.sh
  
echo './gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --autosome --make-grm-gz --out 160205_autoGRM' >> all.sh
echo './gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --autosome ---grm-cutoff 0.025 --make-grm-gz --out 160205_autoGRM_cut' >> all.sh
echo './gcta64 --grm-gz 160205_autoGRM --grm-adj 0 --make-grm-gz --out 160205_autoGRM_adj'  >> all.sh

qsub all.sh

### Create a script for all + PAR

echo '#!/bin/sh'                   >  allwpar.sh
echo ''                            >> allwpar.sh
echo '#$ -cwd'                     >> allwpar.sh
echo '#$ -l h_rt=2:00:00'          >> allwpar.sh
echo ''                            >> allwpar.sh
echo '. /etc/profile.d/modules.sh' >> allwpar.sh
echo ''                            >> allwpar.sh
  
echo './gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 27 --exclude chr27nonparsnplist.txt --make-grm-gz --out 160205_autowparGRM' >> allwpar.sh
echo './gcta64 --grm-gz 160205_autowparGRM --grm-adj 0 --make-grm-gz --out 160205_autowparGRM_adj'  >> allwpar.sh

qsub allwpar.sh

### Run each individual chromosome

for i in {1..26}
do

  echo '#!/bin/sh' > chr${i}.sh
  echo '' >> chr${i}.sh
  echo '#$ -cwd' >> chr${i}.sh
  echo '#$ -l h_rt=2:00:00' >> chr${i}.sh
  echo '' >> chr${i}.sh
  echo '. /etc/profile.d/modules.sh' >> chr${i}.sh
  echo '' >> chr${i}.sh
    
  echo "./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --chr ${i} --make-grm-gz --out 160205_chr${i}_GRM" >> chr${i}.sh
  echo "./gcta64 --grm-gz 160205_chr${i}_GRM --grm-adj 0 --make-grm-gz --out 160205_chr${i}_GRM_adj" >> chr${i}.sh

  echo "./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --autosome --exclude chr${i}snplist.txt --make-grm-gz --out 160205_WOchr${i}_GRM" >> chr${i}.sh
  echo "./gcta64 --grm-gz 160205_WOchr${i}_GRM --grm-adj 0 --make-grm-gz --out 160205_WOchr${i}_GRM_adj" >> chr${i}.sh

#  echo "./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 27 --exclude chr${i}snplist.txt --exclude chr27nonparsnplist.txt --make-grm --out 160205_WOchr${i}_incPAR_GRM" >> chr${i}.sh
#  echo "./gcta64 --grm 160205_WOchr${i}_incPAR_GRM --grm-adj 0 --make-grm --out 160205_chr${i}_GRM_adj" >> chr${i}.sh

  qsub chr${i}.sh
  
done

### Sex chromosome

./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --make-grm-xchr-gz --out 160205_chr27nonpar_GRM

### PseudoAutosomalRegion

./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 27 --chr 27 --exclude chr27nonparsnplist.txt --make-grm-gz --out 160205_chr27par_GRM
./gcta64 --grm-gz 160205_chr27par_GRM --grm-adj 0 --make-grm-gz --out 160205_chr27par_GRM_adj



