# Determine multiple testing threshold based on LD
# Author: Susan Johnston

setwd("data")

# system("cmd", input = "plink --file 20150129merged1_66nodups.QC2 --sheep --recode12 --transpose --out t20150129merged1_66nodups.QC2")
# 
# readLines("t20150129merged1_66nodups.QC2.tped", n = 1)
# 
# system("cmd", input = "cut -d \" \" -f 2,5-  t20150129merged1_66nodups.QC2.tped > test.txt")
# beepr::beep()

x <- readLines("test.txt", n = 1)
x <- strsplit(x, split = " ")[[1]]
(length(x) - 1)/2

library(R.utils)
y <- countLines("t20150129merged1_66nodups.QC2.tped")

#~~ Run through Keff


# Keffective.exe FileName Nmar Nind SigLevel WindowSize

system("cmd", input = paste("Keffective.exe test.txt", y, (length(x) - 1)/2, 0.05, 50))
