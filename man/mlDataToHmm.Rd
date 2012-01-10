\name{mlDataToHmm}
\alias{mlDataToHmm}
\title{Postprocessing of summarized data with VanillaICE}
\usage{
  mlDataToHmm(mlData, runmed = 5, tau = 1e+08,
    copynumberstates = c(0.1, 1, 2, 3, 4),
    states = c("homDel", "hemDel", "normal", "3copy", "4copy"),
    normalindex = 3,
    initiallogprobs = log(rep(1/length(states), length(states))))
}
\arguments{
  \item{mlData}{data from cn.farms}

  \item{runmed}{Integer. Size of window for running median.
  A running median of the total copy number is used to
  estimate the probability that a copy number estimate is
  an outlier.}

  \item{tau}{Numeric. Factor for scaling the distance
  weighted transition probability.}

  \item{copynumberstates}{The mean values for each of the
  hidden states.}

  \item{states}{Character vector indicating the names of
  the hidden states}

  \item{normalindex}{Integer. Indicates which element of
  the states vector corresponds to the 'normal' state.}

  \item{initiallogprobs}{Numeric vector of the same length
  as the number of states. Specifies the initial state
  probabilities on the log scale.}
}
\value{
  Some data
}
\description{
  Postprocessing of summarized data with VanillaICE
}
\author{
  Andreas Mitterecker
}

