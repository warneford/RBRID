CSH_PRC2RBRID <- data.frame(readRDS("~/RBR ID/160419.RBRID.table.PRC2.rds"))

cols <- c(3,4,5,6,7,37:41)

PRC2RBRID <-PRC2RBRID[,cols]
View(PRC2RBRID)

#Import EZH2 IP/RBR-ID with full proteome search
library(gdata)
installXLSXsupport()
Prot_table <- read.xls(xls = "Robert test 1 with mouse proteome.xlsx", sheet = 1)

#Import EZH2 IP/RBR-ID with custom PRC2 database search
PRC2_table <- read.xls(xls = "Robert test 1 with custom database.xlsx", sheet = 1)
str(PRC2_table)

#Sanitize Sample names + Add description
colnames(Prot_table)[8:15] <- c("-4su.Quant_01","-4su.Quant_02","-4su_1:10_Quant_01", "-4su_1:10_Quant_02","+4su.Quant_01","+4su.Quant_02","+4su_1:10_Quant_01","+4su_1:10_Quant_02")

# normalize peptide intensity relative to total




View(Prot_table)

