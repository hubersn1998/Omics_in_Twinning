######################
# Histogram function #
######################
PlotHistogram <- function(x,
                          Metabolite,
                          Legend=FALSE){
  h <- hist(x,
            plot=FALSE)
  xfit <- seq(min(x,
                  na.rm=TRUE),
              max(x,
                  na.rm=TRUE),
              length=100)
  yfit <- dnorm(xfit,
                mean=mean(x,
                          na.rm=TRUE),
                sd=sd(x,
                      na.rm=TRUE))
  d <- density(x[complete.cases(x)],
               from=min(x[complete.cases(x)]),
               to=max(x[complete.cases(x)]))
  XLim = c(min(x[complete.cases(x)]),
           max(x[complete.cases(x)]))
  YLim = c(0.0,
           max(max(d$y),
               max(yfit),
               max(h$density)))
  h <- hist(x,
            xlab=Metabolite,
            main=paste0('Histogram of ',
                        Metabolite),
            freq=FALSE,
            xlim=XLim,
            ylim=YLim)
  lines(xfit,
        yfit,
        col='blue',
        lwd=2)
  lines(d$x,
        d$y,
        col='red',
        lwd=2)
  if(Legend){
    legend('topright',
           xjust=1,
           yjust=1,
           lty=c(1,
                 1,
                 1),
           lwd=c(1,
                 2,
                 2),
           col=c('black',
                 'blue',
                 'red'),
           legend=c('histogram',
                    'normal fit',
                    'kde fit'),
           bty='n')
  }
}

##############################
# Calculate NaN desrcivtives #
##############################
GetNaN <- function(MetabolomicsData) {
  PercentageNA <- matrix()
  NumberOfNA <- matrix()
  for (i in 1:ncol(MetabolomicsData)) {
    Percentage <- (length(which(is.na(MetabolomicsData[,
                                                       i])))/length(MetabolomicsData[,
                                                                                     i]))*100
    NumberOfNA[[i]] <- length(which(is.na(MetabolomicsData[,
                                                           i])))
    PercentageNA[[i]]<- matrix(Percentage)
  }
  NAInformation <- cbind(NumberOfNA,
                         PercentageNA)
  NAInformation <- as.data.frame(NAInformation)
  rownames(NAInformation) <- colnames(MetabolomicsData)
  return(NAInformation)
}

######################################
# Single metabolite linear regression #
#######################################
GLMAnalysisMetabolitesDependent <- function(RegressionData,
                                            DependentVariables,
                                            IndependentVariable,
                                            Covariates,
                                            GlmFunction,
                                            Id){
  # The metabolites M are the DependentVariables
  ResultsDF <- data.frame()
  ResultsDF[,
            'METABOLITE'] <- character(0)
  ResultsDF[,
            'N'] <- integer(0)
  ResultsDF[,
            paste(IndependentVariable,
                  'BETA',
                  sep='_')] <- numeric(0)
  ResultsDF[,
            paste(IndependentVariable,
                  'SE',
                  sep='_')] <- numeric(0)
  ResultsDF[,
            paste(IndependentVariable,
                  'PVALUE',
                  sep='_')] <- numeric(0)
  for(C in Covariates){
    ResultsDF[,
              paste(C,
                    'BETA',
                    sep='_')] <- numeric(0)
    ResultsDF[,
              paste(C,
                    'SE',
                    sep='_')] <- numeric(0)
    ResultsDF[,
              paste(C,
                    'PVALUE',
                    sep='_')] <- numeric(0)
  }
  for(M in DependentVariables){
    Formula <- as.formula(paste(M,
                                paste(c(IndependentVariable,
                                        Covariates),
                                      collapse=' + '),
                                sep=' ~ '))
    if(GlmFunction=='glm'){
      Data <- RegressionData[,
                             c(M,
                               IndependentVariable,
                               Covariates)]
      Data <- Data[complete.cases(Data),
                   ] # use complete cases only!
      DropList <- c()
      for(C in Covariates){
        if(length(unique(Data[,
                              C]))==1)
          DropList <- c(DropList,
                        C)
        if(min(table(Data[,
                          C])) / nrow(Data) < 0.05)
          DropList <- c(DropList,
                        C)
      }
      DropList <- unique(DropList)
      if(length(DropList) > 0){
        for(C in DropList){
          Covariates <- Covariates[Covariates != C]
          ResultsDF[,
                    paste(C,
                          'BETA',
                          sep='_')] <- NA
          ResultsDF[,
                    paste(C,
                          'SE',
                          sep='_')] <- NA
          ResultsDF[,
                    paste(C,
                          'PVALUE',
                          sep='_')] <- NA
        }
        Data <- RegressionData[,
                               c(M,
                                 IndependentVariable,
                                 Covariates,
                                 Id)]
        Data <- Data[complete.cases(Data),
                     ] # use complete cases only!
        Formula <- as.formula(paste(M,
                                    paste(c(IndependentVariable,
                                            Covariates),
                                          collapse=' + '),
                                    sep=' ~ '))
      }
      Model <- glm(formula=Formula,
                   data=Data,
                   family=gaussian(link='identity'))
    }
    else if(GlmFunction=='geeglm'){
      Data <- RegressionData[,
                             c(M,
                               IndependentVariable,
                               Covariates,
                               Id)]
      Data <- Data[complete.cases(Data),
                   ] # use complete cases only!
      DropList <- c()
      for(C in Covariates){
        if(length(unique(Data[,
                              C]))==1)
          DropList <- c(DropList,
                        C)
        Table <- (table(Data[,
                             C])) / nrow(Data)
        if((C == 'BATCH') &&
           (min(table(Data[,
                           C])) / nrow(Data) < 0.05))
          DropList <- c(DropList,
                        C)
      }
      DropList <- unique(DropList)
      if(length(DropList) > 0){
        for(C in DropList){
          Covariates <- Covariates[Covariates != C]
          ResultsDF[M,
                    paste(C,
                          'BETA',
                          sep='_')] <- NA
          ResultsDF[M,
                    paste(C,
                          'SE',
                          sep='_')] <- NA
          ResultsDF[M,
                    paste(C,
                          'PVALUE',
                          sep='_')] <- NA
        }
        Data <- RegressionData[,
                               c(M,
                                 IndependentVariable,
                                 Covariates,
                                 Id)]
        Data <- Data[complete.cases(Data),
                     ] # use complete cases only!
        Formula <- as.formula(paste(M,
                                    paste(c(IndependentVariable,
                                            Covariates),
                                          collapse=' + '),
                                    sep=' ~ '))
      }
      Model <- geeglm(formula=Formula,
                      data=Data,
                      family=gaussian(link='identity'),
                      id=Data[,
                              Id],
                      corst='exchangeable')
    }
    Output <- coef(summary(Model))
    ResultsDF[M,
              'METABOLITE'] <- M
    ResultsDF[M,                
              'N'] <- nrow(Data)
    ResultsDF[M,
              c(paste(IndependentVariable,
                      'BETA',
                      sep='_'),
                paste(IndependentVariable,
                      'SE',
                      sep='_'),
                paste(IndependentVariable,
                      'PVALUE',
                      sep='_'))] <- as.numeric(Output[IndependentVariable,
                                                      c(1,
                                                        2,
                                                        ncol(Output))])
    for(C in Covariates){
      ResultsDF[M,
                c(paste(C,
                        'BETA',
                        sep='_'),
                  paste(C,
                        'SE',
                        sep='_'),
                  paste(C,
                        'PVALUE',
                        sep='_'))] <- as.numeric(Output[C,
                                                        c(1,
                                                          2,
                                                          ncol(Output))])
    }
  }
  return(ResultsDF)
}
