#
# README
# ------
#
# This script should be utilized for analyzing how blood metabolite levels depend on zygosuty.
# The analysis plan by René Pool (r.pool#vu.nl), Nikki Hubers (n.hubers@vu.nl) and Leonie Bogl (leonie-helen.bogl@helsinki.fi)
# on which this script is based, is included in directory "../Documents".
#
# This script needs editing enabling it to run correctly for your cohort.
# Lines that need editing below are preceded by "### <-- EDIT" and followed
# by "### EDIT -->".
#
# Make sure your phenotype data file and metabolomics data file are located in
# "../Phenotypes".
#
# The analyses results are written to  "../Results/<CURRENT_DATE_IN_YYYYMMDD>".
# The results that need to be provided are written to 
# "../Uploaded/<CURRENT_DATE_IN_YYYYMMDD>". 
# the generated file "../Uploaded/<CURRENT_DATE_IN_YYYYMMDD>.zip".
#
# This script also requires the analyses from Metabolomics_Twinning_AnalysesFunctions.R
# Please direct any questions or remarks on this script to Nikki Hubers (n.hubers@vu.nl).
#
# Happy computing :-)
#
# René Pool & Nikki Hubers
#

######################
# + Import libraries #
######################
library('summarytools')
library('xlsx')
######################
# - Import libraries #
######################

################################################
# + Specify analysis variables for your cohort #
################################################
## <-- EDIT
Cohort <- 'cohort' # Your cohort abbreviation
## EDIT -->
Date <- strftime(Sys.time(),
                 format="%Y%m%d") # Set automatically
## <-- EDIT
GLMFunctionName <- 'geeglm' # Currently the only possibilities are 'glm' or
# 'geeglm'

IDString <- NULL # Do not change this line!
if(GLMFunctionName == 'geeglm'){
  library('geepack')
  IDString <- 'FamilyNumber' # Set the appropriate identitystring
  # IMPORTANT: Sort data by this variable otherwise
  # clusters will become ill defined!
}
InDepVar <- 'Zyg' # Dependent variable

# Inspect the list of covariates for your cohort below and change accordingly
CohortSpecificCovariates <- c('covariats')
# CohortSpecificCovariates <- c()
AnalyisCovariates = c('Age',
                      'Sex',
                      'Smoke',
                      'BMI',
                      'Fasting',
                      'LipMed',
                      CohortSpecificCovariates)

remove(CohortSpecificCovariates)
################################################
# - Specify analysis variables for your cohort #
################################################

###############
# + Set paths #
###############
ParentWd <- strsplit(getwd(),
                     split='/')[[1]]
ParentWd <- paste(ParentWd[1:(length(ParentWd)-1)],
                  collapse = '/')
ScriptsPath <- paste(ParentWd,
                     'Scripts',
                     sep='/')
PhenotypesPath <- paste(ParentWd,
                        'Phenotypes',
                        sep='/')
ResultsPath <- paste(ParentWd,
                     'Results',
                     paste(Date,
                           sep='.'),
                     sep='/')
dir.create(ResultsPath,
           showWarnings=FALSE)
UploadedPath <- paste(ParentWd,
                      'Uploaded',
                      Date,
                      sep='/')
dir.create(UploadedPath,
           showWarnings=FALSE)
setwd(ScriptsPath)
###############
# - Set paths #
###############

########################
# + Source R functions #
########################
FunctionsFName <- 'AnalysisFunctions.R'
FunctionsFName <- paste(ScriptsPath,
                        FunctionsFName,
                        sep='/')
source(FunctionsFName)
remove(FunctionsFName)
########################
# - Source R functions #
########################

##########################
# + Parse phenotype data #
##########################
setwd("")
PhenotypeDF <- read.table('',
                          header=T,
                          sep='\t',
                          na.strings='NA',
                          as.is=T) # Depending on the file type containing 
                                   # your phenotype data, the type and/or 
                                   # argument(s) of the parsing function 
remove(PhenotypeFName)

MetaboliteList <- colnames(PhenotypeDF)[StartMetabolitesColumn:length(colnames(PhenotypeDF))]

PhenotypeDF[,
            'Zyg'] <- 0
PhenotypeDF[!is.na(PhenotypeDF$DZTwinCode),
            'Zyg'] <- 1
PhenotypeDF[!is.na(PhenotypeDF$MZTwinCode),
            'Zyg'] <- 2
PhenotypeDF <- PhenotypeDF[PhenotypeDF$Zyg!=0,
                           ]
row.names(PhenotypeDF) <- NULL
##########################
# - Parse phenotype data #
##########################

#############################################
# + Determine the number of functional test #
#############################################
Metabolites <- PhenotypeDF[StartMetabolitesColumn:length(colnames(PhenotypeDF))]

DetermineMEffLi <- function(Data,Columns){
  MeanCenteredAutoScaledData <- Data
  for(i in 1:length(Columns)){
    Entry <- Columns[i]
    Mean <- mean(Data[,Entry])
    Std <- sd(Data[,Entry])
    MeanCenteredAutoScaledData[,Entry] <- (MeanCenteredAutoScaledData[,Entry] - Mean)/Std
  }
  
  
  SVD <- svd(x=as.matrix(MeanCenteredAutoScaledData[,Columns]))
  EigenValues <- SVD$d^2 / dim(SVD$u)[1]
  M <- length(EigenValues)
  L <- M - 1
  Var <- var(EigenValues)
  MEff <- M*(1.0-(L*Var/M^2))
  IntEVals <- as.numeric(EigenValues>=1.0)
  NonIntEVals <- EigenValues - floor(EigenValues)
  MEffLi <- ceiling(sum(IntEVals) + sum(NonIntEVals))
  MEffLi
}

#running MSD
meta_list <- colnames(Metabolites)
MEffLiC <- DetermineMEffLi(Data=Metabolites, Columns=meta_list) 
print (MEffLiC)

#bonferonni correction 
print(0.05/MEffLiC)

#############################################
# - Determine the number of functional test #
#############################################

###########################################################
# + Generate distributions of metabolomics data           #
###########################################################
pdf(file=paste(ResultsPath,
               '',
               sep='/'))
for(i in 1:length(MetaboliteList)){
  M <- MetaboliteList[i]
  PlotHistogram(x=PhenotypeDF[,
                              M],
                Metabolite=M,
                Legend=TRUE)
  remove(M)
}
dev.off()

###########################################################
# - Generate distributions of metabolomics data           #
###########################################################

#############################################################
# + Get NaN decriptives for al metabolite columns           #
#############################################################
Result <- GetNaN(PhenotypeDF[,
                             MetaboliteList])
save(Result,
     file=paste(ResultsPath,
                '.RData',
                sep='/'))
write.table(Result,
            file=paste(ResultsPath,
                       '.tsv',
                       sep='/'),
            sep='\t',
            col.names=NA)
remove(Result)
#############################################################
# - Get NaN decriptives for al metabolite columns           #
#############################################################

################################
# + Perform analyses           #
################################
RegressionDF <- PhenotypeDF
if(GLMFunctionName == 'geeglm'){
  RegressionDF <- RegressionDF[order(RegressionDF[,
                                                  IDString],
                                     na.last=TRUE),]
  rownames(RegressionDF) <- NULL
}

# Single metabolite linear regression Model:
LinearResults <- GLMAnalysisMetabolitesDependent(RegressionData=RegressionDF,
                                                 DependentVariables=MetaboliteList,
                                                 IndependentVariable=InDepVar,
                                                 Covariates=AnalyisCovariates,
                                                 GlmFunction=GLMFunctionName,
                                                 Id=IDString)
FName <- paste(Cohort,
               'Zyg',
               'cohort',
               Date,
               'RData',
               sep='.')
save(LinearResults,
     file=paste(ResultsPath,
                FName,
                sep='/'))
remove(FName)
FName <- paste(Cohort,
               'Zyg',
               'cohort',
               Date,
               'tsv',
               sep='.')
write.table(LinearResults,
            file=paste(ResultsPath,
                       FName,
                       sep='/'),
            sep='\t',
            row.names=FALSE)
remove(FName)
FName <- paste(Cohort,
               'Zyg',
               'cohort',
               Date,
               'xlsx',
               sep='.')
write.xlsx(LinearResults,
           file=paste(ResultsPath,
                      FName,
                      sep='/'),
           sheetName='LinearResults',
           row.names=FALSE,
           col.names=TRUE,
           showNA=FALSE)
remove(FName)
remove(LinearResults)

################################
# - Perform analyses           #
################################

