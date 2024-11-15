###############
# Enrichment RNA
###############

# This script is used to perform enrichment analyses from the results from the TWAS analyses
# The complete TWAS analyses can be found in the script ...

#Packages
library(enrichR)
library(ggplot2)
require(gridExtra)
library(ggpubr)

#############
# Enrichr checks
#############
#Prepare enrichr and check if the connection is live
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
# This should give you a list with all available databases.
# If not, check you connection to internet

# Select the databases you are interested in, in this project these are the GO 2023 and KEGG 2021 databases. 
dsb1 <- c("GO_Molecular_Function_2023", "KEGG_2021_Human")

#############
# Data
#############

# This is a file that should contain a list of all the genes you want to investigate in the enrichment analyses
# For example, a list of genes that were significantly differently expressed between two groups or a list of the top 5% most stronlgy differenltially expressed genes
# This file should not contain a header. 

sig <- read.table("INSERT_FILE_NAMES.txt", header = F)

############
# Enrichment analyses
############

# This command will perform the enrichment in the earlier defined databases (line 28)
if (websiteLive) {
  enriched <- enrichr(c(sig$V1
  ), dsb1)
}

############
# Plots
############

# This command will make a plot of the top 10 most strongly, based on p-value, enriched pathways in the indicated database
bp1 <- if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 10, numChar = 85, y = "Count", orderBy = "P.value", title = "GO Molecular Function 2023")+theme(text=element_text(size=15))
}
# This command will make a plot of the top 10 most strongly, based on p-value, enriched pathways in the indicated database
bp2 <- if (websiteLive) {
  plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG Human 2021")+theme(text=element_text(size=15))
}

#combine the figures
tiff(filename = paste0(as.character(Sys.Date()),"INSERT_FILE_NAME.tiff"),
     res = 300, width = 550, height = 180, units = "mm")
ggarrange(bp1,bp2,# list of plots
          common.legend = T, # COMMON LEGEND
          legend = "bottom", # legend position
          align = "hv", # Align them both, horizontal and vertical
          nrow = 1, ncol = 2)  # number of rows & columns
dev.off()

