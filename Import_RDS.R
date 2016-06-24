CSH_PRC2RBRID <- data.frame(readRDS("~/RBR ID/160419.RBRID.table.PRC2.rds"))

cols <- c(3,4,5,6,7,37:41)

PRC2RBRID <-PRC2RBRID[,cols]
View(PRC2RBRID)
