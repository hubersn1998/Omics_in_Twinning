#####
# MA Metabolomics
####

# This script can be used to perform a meta-analyses on metabolite levels and zygosity.
# This script requires input that can be made with the script ...
# This script is based on two cohorts, but more can be added if necessary.

#Libraries
library(ggplot2)
library(grid)
library(metafor)
library(readxl)
library(xlsx)
library(data.table)

##############
# Load the data
##############
# For each cohort, the data set should at least contain the columns; study, metabolite, N, BETA, SE, P

#dataset 1
df1 <- read_excel('INSERT_FILE_NAME.xlsx',
                     sheet='LinearResults')
# Rename the columns if needed
df1[,
       'Zyg_BETA'] <- df1$zyg_BETA
df1[,
       'Zyg_SE'] <- df1$zyg_SE
df1[,
       'Zyg_PVALUE'] <- df1$zyg_PVALUE
df1 [,
        'Study'] <- 'df1'

# Keep only the needed columns
df1 <- df1[,
                 c('Study',
                   'Metabolite',
                   'N',
                   'Zyg_BETA',
                   'Zyg_SE',
                   'Zyg_PVALUE')]
############################################

####################################################################
#dataset 2
df2 <- read_excel('INSERT_FILE_NAME.xlsx',
                     sheet='LinearResults')
# Adapt the columnnames
df2[,
       'Metabolite'] <- df_nl$METABOLITE
df2[,
       'Study'] <- 'df2'

# Keep only the needed columns
df2 <- df2[,
                 c('Study',
                   'Metabolite',
                   'N',
                   'Zyg_BETA',
                   'Zyg_SE',
                   'Zyg_PVALUE')]
#################################################################

#####################
# Merge the data
#####################

# Merge NTR and old
MDF <- merge(x=df1[,
                      c('Metabolite',
                        'N')],
             y=df2[,
                      c('Metabolite',
                        'N')],
             by='Metabolite', 
             suffixes=c('_df1','_df2'),
             all=TRUE)

# Create a list with all shared metabolites. 
SharedTraits <- MDF[complete.cases(MDF),]

# Create an empty data frame where you can copy all shared metabolites
MADF <- data.frame(row.names=1:length(SharedTraits))
MADF[,
     'Metabolite'] <- SharedTraits 
# Select the shared metabolites from the two data sets
df1_selected <- df1[df1$Metabolite %in% SharedTraits,]
df2_selected <- df2[df2$Metabolite %in% SharedTraits,]

#####################
# Run an REMA and FEMA
###################

# this loop takes one metabolite and then runs a meta analyses with a random effect model (REMA) and a fixed effect model (FEMA)
# afterwards the loop takes the next metabolites.
# The loop also immediately binds the outcomes of the two models for each metabolite to a data frame called MADF

for(i in 1:length(SharedTraits)){
    Metabolite2 <- SharedTraits[i]
    MAInpuDF <- as.data.frame(c(df1_selected[df1_selected$Metabolite==Metabolite2,
                                       ]))
    MAInpuDF <- rbind(MAInpuDF,                          
                      as.data.frame(c(df2_selected[df2_selected$Metabolite==Metabolite2,
                                             ])))
    remove(Metabolite2)
    REMAModel <- rma(yi=Zyg_BETA,
                     sei=Zyg_SE,
                     weights=N,
                     data=MAInpuDF)
    FEMAModel <- rma(yi=Zyg_BETA,
                     sei=Zyg_SE,
                     weights=N,
                     method='FE',
                     data=MAInpuDF)
    MADF[i,
         'N_MA'] <- sum(MAInpuDF$N)
    #remove(MAInpuDF)
    MADF[i,
         'Zyg_HETISQ_MA'] <- REMAModel$I2
    MADF[i,
         'Zyg_HETPVALUE_MA'] <- REMAModel$QEp
    MADF[i,
         'Zyg_BETA_REMA'] <- REMAModel$beta
    MADF[i,
         'Zyg_CILB_REMA'] <- REMAModel$ci.lb
    MADF[i,
         'Zyg_CIUB_REMA'] <- REMAModel$ci.ub
    MADF[i,
         'Zyg_SE_REMA'] <- REMAModel$se
    MADF[i,
         'Zyg_PVALUE_REMA'] <- REMAModel$pval
    #remove(REMAModel)
    MADF[i,
         'Zyg_BETA_FEMA'] <- FEMAModel$beta
    MADF[i,
         'Zyg_CILB_FEMA'] <- FEMAModel$ci.lb
    MADF[i,
         'Zyg_CIUB_FEMA'] <- FEMAModel$ci.ub
    MADF[i,
         'Zyg_SE_FEMA'] <- FEMAModel$se
    MADF[i,
         'Zyg_PVALUE_FEMA'] <- FEMAModel$pval
    #remove(FEMAModel)
}
remove(i)
remove(SharedTraits_oldNTR)
#########################################################################

###############
# Create output file and save
##############

# Bind the complete data sets
MDF <- merge(x=df1,
             y=df2,
             by='Metabolite',
             suffixes=c('_df1',
                        '_df2'),
             all=TRUE)
# Bind the complete data sets with the outcomes from the REMA and FEMA model
MDF <- merge(x=MDF,
             y=MADF,
             by='Metabolite',
             all=TRUE)

# Save the complete output data. 
write.xlsx(x=MDF,
           file='INSERT_FILE_NAME.xlsx',
           row.names=FALSE,
           showNA=FALSE)
