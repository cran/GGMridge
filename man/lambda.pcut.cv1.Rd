% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda.pcut.cv1.R
\name{lambda.pcut.cv1}
\alias{lambda.pcut.cv1}
\title{Choose the Tuning Parameter of the Ridge Inverse and 
 Thresholding Level of the Empirical p-Values.
 
Calculate total prediction error for test data after fitting partial 
  correlations from train data for all values of lambda and pcut.}
\usage{
lambda.pcut.cv1(train, test, lambda, pcut)
}
\arguments{
\item{train}{An n x p data matrix from which the model is fitted.}

\item{test}{An m x p data matrix from which the model is evaluated.}

\item{lambda}{A vector of candidate tuning parameters.}

\item{pcut}{A vector of candidate cutoffs of pvalues.}
}
\value{
Total prediction error for all the candidate lambda and pvalue 
  cutoff values.
}
\description{
Choose the Tuning Parameter of the Ridge Inverse and 
 Thresholding Level of the Empirical p-Values.
 
Calculate total prediction error for test data after fitting partial 
  correlations from train data for all values of lambda and pcut.
}
\examples{
 p <- 100 # number of variables
 n <- 50 # sample size
 
 ###############################
 # Simulate data
 ###############################
 simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
 data <- simulation$data[[1L]]
 
 ###############################
 # Split into train/test sets
 ###############################
 testindex <- sample(1L:n, 10L)
 
 train <- data[-testindex,,drop = FALSE]
 stdTrain <- scale(x = train, center = TRUE, scale = TRUE)
 
 test <- data[testindex,,drop = FALSE]
 stdTest <- scale(x = test, center = TRUE, scale = TRUE)
 
 ###############################
 # Calculate total prediction 
 # errors for all candidate 
 # lambda and p-value cutoffs
 ###############################
 lambda.array <- seq(from = 0.1, to = 5, length = 10) * (n - 1.0)
 pcut.array <- seq(from = 0.01, to = 0.05, by = 0.01)
 tpe <- lambda.pcut.cv1(train = stdTrain,
                        test = stdTest,
                        lambda = lambda.array,
                        pcut = pcut.array)
  
}
\references{
Ha, M. J. and Sun, W. 
 (2014).
 Partial correlation matrix estimation using ridge penalty followed 
 by thresholding and re-estimation.
 Biometrics, 70, 762--770.
}
\author{
Min Jin Ha
}
