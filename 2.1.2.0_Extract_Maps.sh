for i in {1..27}
do

  echo '#!/bin/sh' > crim${i}.sh
  echo '' >> crim${i}.sh
  echo '#$ -cwd' >> crim${i}.sh
  echo '#$ -l h_rt=4:00:00' >> crim${i}.sh
  echo '' >> crim${i}.sh
  echo '. /etc/profile.d/modules.sh' >> crim${i}.sh
  echo '' >> crim${i}.sh
    
  echo "./gcta64 --bfile 20150129merged1_66nodups.QC2 --autosome-num 26 --chr ${i} --make-grm-gz --out 150129_chr${i}_GRM" >> crim${i}.sh
  echo "./crimap ${i}g build > chr${i}g.map"

  qsub crim${i}.sh
  
done