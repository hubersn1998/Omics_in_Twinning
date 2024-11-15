#########################
# First analysis GEE
# Gene Expression data
########################

#############
# Preparations
#############

#load R packages
library(haven)
library(gee)
library(geepack)
library(foreach)
library(doMC)
library(ggplot2)
registerDoMC(15)  #number of cpu's to use

# Load the data, these should be one file containing all the RNA probes as columns and IDs as rows. 
# The second file should contain the phenotype and covariate information as column and the IDs as rows. 
# Ensure you change the command to the appropriate input data
RNA <- load("RNA_FILE.RData")
df <- load("PHENOTYPE_FILE.RData")

############
# Select all the covariates and phenotypes
############

# Make sure you add as.numeric, as.factor, etc. depending on your data. 

#ID and Fam ID
ID <- df$"SELECT_COLUMN"
fam <- df$"SELECT_COLUMN"
# Outcome of zygosity
zyg <- df$"SELECT_COLUMN"
# Covariates
sex <- df$"SELECT_COLUMN"
age <- df$"SELECT_COLUMN"
plate <- df$"SELECT_COLUMN"
bmi <-df$"SELECT_COLUMN"
smoking <- df$"SELECT_COLUMN"
well <- df$"SELECT_COLUMN"
# Cell Counts
wbc <- df$"SELECT_COLUMN"
neutro <- df$"SELECT_COLUMN"
lympho <- df$"SELECT_COLUMN"
mono <- df$"SELECT_COLUMN"
eos <- df$"SELECT_COLUMN"
baso <- df$"SELECT_COLUMN"

#############
# Combine the data
#############

# Make a new list 'probes', which will be used to fill in a specific probe in the analyses step
probes <- rep(NA,times=nrow(RNA))
probes

# Combine all the information
data <- data.frame(ID,fam, zyg, probes, sex,age,plate,bmi,smoking,well,wbc,
                   neutro, lympho, mono, eos)
head(data)

# make a list of all rows that have complete data
select <- which(!is.na(data$zyg) &!is.na(data$sex) & !is.na(data$age) & !is.na(data$plate) & !is.na(data$bmi) &!is.na(data$smoking) &!is.na(data$well)
                &!is.na(data$wbc)&!is.na(data$neutro)&!is.na(data$lympho)&!is.na(data$mono)&!is.na(data$eos)&!is.na(data$baso)&!is.na(data$fam))
select
length(select)

# Keep only the rows with complete data
mydata <- data[select,]
dim(mydata)
# Make sure the row names of mydata match the IDs
row.names(mydata) <- mydata$ID
# select these samples from the RNA data, which should also have the ID as rownames. 
RNA <- RNA[rownames(mydata),]
dim(RNA)  

##########
# Run linear model for all probes
#########

# Order the data by the relatedness structure e.g. family/pedigree number
mydata <- mydata[order(mydata$fam),]
RNA <- RNA[match(rownames(mydata), rownames(RNA)), ]
# Check if the row names are the same
table(rownames(RNA)==rownames(mydata))
# Needs to be TRUE before continuing. 

# linear model for a continuous outcome (RNA level) 
# in case of a related sample, the 'id=fam' statement needs to include the relatedness structure of your data
# Make sure that not all cell counts are included in the model as these will add up to 100% and therefore create multi-colinearity problems. 
# We apply this model to all RNA probes 

# Define the model
lmmodel <- function (dat,probes) 
{
  r1=gee(as.formula(paste0(probes,"~zyg+age+sex+bmi+smoking+wbc+neutro+lympho+mono+eos+plate+well")), 
      id=fam, family=gaussian,data=dat)
  coeff1 <- summary(r1)
  coeff <- as.data.frame(coeff1[["coefficients"]])
  beta = coeff[,1]
  naive.se = coeff[,2]
  naive.z = coeff[,3]
  robust.se = coeff[,4]
  robust.z = coeff[,5]
  pvalue= 2 * (1 - pnorm(abs(coeff[,5])))
  out = cbind(beta,naive.se,naive.z, robust.se, robust.z, pvalue)
  return(out)
}

# test Run 15 probes
TWAS <- foreach(i=1:15) %dopar%
  {
    lmmodel(dat=mydata, probes = RNA[i])
  }

# Make table of the results
# Adapted the order/names of the rows to match your covariates 
TWASresults <- as.data.frame(matrix(NA,nrow = 15,ncol = 0))
row.names(TWASresults) <- c(colnames(RNA[1:15]))

# Betas
TWASresults$beta.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,1] })))
TWASresults$beta.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,1] })))
TWASresults$beta.age <- c(unlist(lapply(TWAS, function(x) { x[3,1] })))
TWASresults$beta.sex <- c(unlist(lapply(TWAS, function(x) { x[4,1] })))
TWASresults$beta.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,1] })))
TWASresults$beta.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,1] })))
TWASresults$beta.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,1] })))
TWASresults$beta.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,1] })))
TWASresults$beta.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,1] })))
TWASresults$beta.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,1] })))
TWASresults$beta.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,1] })))
TWASresults$beta.mono <- c(unlist(lapply(TWAS, function(x) { x[12,1] })))
TWASresults$beta.eos <- c(unlist(lapply(TWAS, function(x) { x[13,1] })))

# SE
TWASresults$naive.se.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,2] })))
TWASresults$naive.se.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,2] })))
TWASresults$naive.se.age <- c(unlist(lapply(TWAS, function(x) { x[3,2] })))
TWASresults$naive.se.sex <- c(unlist(lapply(TWAS, function(x) { x[4,2] })))
TWASresults$naive.se.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,2] })))
TWASresults$naive.se.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,2] })))
TWASresults$naive.se.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,2] })))
TWASresults$naive.se.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,2] })))
TWASresults$naive.se.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,2] })))
TWASresults$naive.se.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,2] })))
TWASresults$naive.se.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,2] })))
TWASresults$naive.se.mono <- c(unlist(lapply(TWAS, function(x) { x[12,2] })))
TWASresults$naive.se.eos <- c(unlist(lapply(TWAS, function(x) { x[13,2] })))

# naive.z
TWASresults$naive.z.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,3] })))
TWASresults$naive.z.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,3] })))
TWASresults$naive.z.age <- c(unlist(lapply(TWAS, function(x) { x[3,3] })))
TWASresults$naive.z.sex <- c(unlist(lapply(TWAS, function(x) { x[4,3] })))
TWASresults$naive.z.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,3] })))
TWASresults$naive.z.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,3] })))
TWASresults$naive.z.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,3] })))
TWASresults$naive.z.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,3] })))
TWASresults$naive.z.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,3] })))
TWASresults$naive.z.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,3] })))
TWASresults$naive.z.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,3] })))
TWASresults$naive.z.mono <- c(unlist(lapply(TWAS, function(x) { x[12,3] })))
TWASresults$naive.z.eos <- c(unlist(lapply(TWAS, function(x) { x[13,3] })))

# robust.se
TWASresults$robust.se.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,4] })))
TWASresults$robust.se.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,4] })))
TWASresults$robust.se.age <- c(unlist(lapply(TWAS, function(x) { x[3,4] })))
TWASresults$robust.se.sex <- c(unlist(lapply(TWAS, function(x) { x[4,4] })))
TWASresults$robust.se.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,4] })))
TWASresults$robust.se.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,4] })))
TWASresults$robust.se.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,4] })))
TWASresults$robust.se.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,4] })))
TWASresults$robust.se.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,4] })))
TWASresults$robust.se.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,4] })))
TWASresults$robust.se.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,4] })))
TWASresults$robust.se.mono <- c(unlist(lapply(TWAS, function(x) { x[12,4] })))
TWASresults$robust.se.eos <- c(unlist(lapply(TWAS, function(x) { x[13,4] })))

# robust.z
TWASresults$robust.z.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,5] })))
TWASresults$robust.z.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,5] })))
TWASresults$robust.z.age <- c(unlist(lapply(TWAS, function(x) { x[3,5] })))
TWASresults$robust.z.sex <- c(unlist(lapply(TWAS, function(x) { x[4,5] })))
TWASresults$robust.z.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,5] })))
TWASresults$robust.z.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,5] })))
TWASresults$robust.z.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,5] })))
TWASresults$robust.z.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,5] })))
TWASresults$robust.z.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,5] })))
TWASresults$robust.z.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,5] })))
TWASresults$robust.z.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,5] })))
TWASresults$robust.z.mono <- c(unlist(lapply(TWAS, function(x) { x[12,5] })))
TWASresults$robust.z.eos <- c(unlist(lapply(TWAS, function(x) { x[13,5] })))

# pvalue
TWASresults$pvalue.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,6] })))
TWASresults$pvalue.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,6] })))
TWASresults$pvalue.age <- c(unlist(lapply(TWAS, function(x) { x[3,6] })))
TWASresults$pvalue.sex <- c(unlist(lapply(TWAS, function(x) { x[4,6] })))
TWASresults$pvalue.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,6] })))
TWASresults$pvalue.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,6] })))
TWASresults$pvalue.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,6] })))
TWASresults$pvalue.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,6] })))
TWASresults$pvalue.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,6] })))
TWASresults$pvalue.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,6] })))
TWASresults$pvalue.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,6] })))
TWASresults$pvalue.mono <- c(unlist(lapply(TWAS, function(x) { x[12,6] })))
TWASresults$pvalue.eos <- c(unlist(lapply(TWAS, function(x) { x[13,6] })))

head(TWASresults)
# If this works well you can go ahead an run the analyses for all probes

# Run for all of the probes
TWAS <- foreach(i=1:ncol(RNA)) %dopar%
  {
    lmmodel(dat=mydata, probes = RNA[i])
  }

# Make a table with the results
# Make table of the results
TWAS_All <- as.data.frame(matrix(NA,nrow = ncol(RNA),ncol = 0))
row.names(TWAS_All) <- c(colnames(RNA))

# Betas
TWAS_All$beta.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,1] })))
TWAS_All$beta.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,1] })))
TWAS_All$beta.age <- c(unlist(lapply(TWAS, function(x) { x[3,1] })))
TWAS_All$beta.sex <- c(unlist(lapply(TWAS, function(x) { x[4,1] })))
TWAS_All$beta.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,1] })))
TWAS_All$beta.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,1] })))
TWAS_All$beta.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,1] })))
TWAS_All$beta.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,1] })))
TWAS_All$beta.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,1] })))
TWAS_All$beta.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,1] })))
TWAS_All$beta.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,1] })))
TWAS_All$beta.mono <- c(unlist(lapply(TWAS, function(x) { x[12,1] })))
TWAS_All$beta.eos <- c(unlist(lapply(TWAS, function(x) { x[13,1] })))

# SE
TWAS_All$naive.se.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,2] })))
TWAS_All$naive.se.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,2] })))
TWAS_All$naive.se.age <- c(unlist(lapply(TWAS, function(x) { x[3,2] })))
TWAS_All$naive.se.sex <- c(unlist(lapply(TWAS, function(x) { x[4,2] })))
TWAS_All$naive.se.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,2] })))
TWAS_All$naive.se.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,2] })))
TWAS_All$naive.se.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,2] })))
TWAS_All$naive.se.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,2] })))
TWAS_All$naive.se.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,2] })))
TWAS_All$naive.se.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,2] })))
TWAS_All$naive.se.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,2] })))
TWAS_All$naive.se.mono <- c(unlist(lapply(TWAS, function(x) { x[12,2] })))
TWAS_All$naive.se.eos <- c(unlist(lapply(TWAS, function(x) { x[13,2] })))

# naive.z
TWAS_All$naive.z.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,3] })))
TWAS_All$naive.z.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,3] })))
TWAS_All$naive.z.age <- c(unlist(lapply(TWAS, function(x) { x[3,3] })))
TWAS_All$naive.z.sex <- c(unlist(lapply(TWAS, function(x) { x[4,3] })))
TWAS_All$naive.z.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,3] })))
TWAS_All$naive.z.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,3] })))
TWAS_All$naive.z.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,3] })))
TWAS_All$naive.z.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,3] })))
TWAS_All$naive.z.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,3] })))
TWAS_All$naive.z.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,3] })))
TWAS_All$naive.z.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,3] })))
TWAS_All$naive.z.mono <- c(unlist(lapply(TWAS, function(x) { x[12,3] })))
TWAS_All$naive.z.eos <- c(unlist(lapply(TWAS, function(x) { x[13,3] })))

# robust.se
TWAS_All$robust.se.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,4] })))
TWAS_All$robust.se.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,4] })))
TWAS_All$robust.se.age <- c(unlist(lapply(TWAS, function(x) { x[3,4] })))
TWAS_All$robust.se.sex <- c(unlist(lapply(TWAS, function(x) { x[4,4] })))
TWAS_All$robust.se.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,4] })))
TWAS_All$robust.se.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,4] })))
TWAS_All$robust.se.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,4] })))
TWAS_All$robust.se.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,4] })))
TWAS_All$robust.se.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,4] })))
TWAS_All$robust.se.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,4] })))
TWAS_All$robust.se.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,4] })))
TWAS_All$robust.se.mono <- c(unlist(lapply(TWAS, function(x) { x[12,4] })))
TWAS_All$robust.se.eos <- c(unlist(lapply(TWAS, function(x) { x[13,4] })))

# robust.z
TWAS_All$robust.z.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,5] })))
TWAS_All$robust.z.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,5] })))
TWAS_All$robust.z.age <- c(unlist(lapply(TWAS, function(x) { x[3,5] })))
TWAS_All$robust.z.sex <- c(unlist(lapply(TWAS, function(x) { x[4,5] })))
TWAS_All$robust.z.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,5] })))
TWAS_All$robust.z.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,5] })))
TWAS_All$robust.z.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,5] })))
TWAS_All$robust.z.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,5] })))
TWAS_All$robust.z.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,5] })))
TWAS_All$robust.z.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,5] })))
TWAS_All$robust.z.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,5] })))
TWAS_All$robust.z.mono <- c(unlist(lapply(TWAS, function(x) { x[12,5] })))
TWAS_All$robust.z.eos <- c(unlist(lapply(TWAS, function(x) { x[13,5] })))

# pvalue
TWAS_All$pvalue.intercept <- c(unlist(lapply(TWAS, function(x) { x[1,6] })))
TWAS_All$pvalue.zyg <- c(unlist(lapply(TWAS, function(x) { x[2,6] })))
TWAS_All$pvalue.age <- c(unlist(lapply(TWAS, function(x) { x[3,6] })))
TWAS_All$pvalue.sex <- c(unlist(lapply(TWAS, function(x) { x[4,6] })))
TWAS_All$pvalue.bmi <- c(unlist(lapply(TWAS, function(x) { x[5,6] })))
TWAS_All$pvalue.smokingcurrent <- c(unlist(lapply(TWAS, function(x) { x[6,6] })))
TWAS_All$pvalue.smokingex <- c(unlist(lapply(TWAS, function(x) { x[7,6] })))
TWAS_All$pvalue.smokingnooit <- c(unlist(lapply(TWAS, function(x) { x[8,6] })))
TWAS_All$pvalue.wbc <- c(unlist(lapply(TWAS, function(x) { x[9,6] })))
TWAS_All$pvalue.neutro <- c(unlist(lapply(TWAS, function(x) { x[10,6] })))
TWAS_All$pvalue.lympho <- c(unlist(lapply(TWAS, function(x) { x[11,6] })))
TWAS_All$pvalue.mono <- c(unlist(lapply(TWAS, function(x) { x[12,6] })))
TWAS_All$pvalue.eos <- c(unlist(lapply(TWAS, function(x) { x[13,6] })))

head(TWAS_All)

# Add your own multiple correction, here we added the bonfferoni correction as example 
alpha_bonf <- 0.05 / ncol(RNA)
sig <- TWAS_All[which(TWAS_All$pvalue.zyg <= alpha_bonf),]

save(TWAS_All, file = "INSERT_FILE_NAME.RData")
save(sig, file = "INSERT_FILE_NAME.RData")     

