library(devtools)
source('c:/RForge/Nuweb.R')
mc <- as.package("c:/RForge/multiCA")


nuweb(mc)
document(mc)
run_examples(mc)  # or dev_example("ran.CMData")
load_all(mc)

check(mc, cleanup = FALSE, check_dir = "c:/RForge", check_version = TRUE, cran = TRUE)


#create data set
strk <- data.matrix(read.delim("z:/EOGeorge/MultiTrend/StrokeData.txt", row.names=1))
colnames(strk) <- gsub("X", "", colnames(strk))
stroke <- as.data.frame.table(strk)
names(stroke) <- c("Type", "Year", "Freq")
stroke$Year <- as.numeric(as.character(stroke$Year))

use_data(stroke, pkg=mc)
promptData(stroke)