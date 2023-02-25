#' @title Average Relative Interval Length (ARIL)
#' @description
#' This function calculates the average relative interval length (ARIL)
#' given a set of predictions 'X', the target 'Y', and a confidence level
#' of 1 - alpha.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @param alpha specifies 1 - alpha confidence level \[0 < alpha < 1\]
#' @return average relative interval length \[0 < scalar < 1\]
#' @references
#' Dogulu, N., P. Lopez Lopez, D.P. Solomatine, A.H. Weerts and
#' D.L. Shreshta (2015), Estimation of predictive hydrologic uncertainty
#' using the quantile regression and UNEEC methods and their comparison
#' on constrasting catchments, Hydrol. Earth. Syst. Sci., 19, 3181-3201.
#' @rdname aril
#' @export
aril <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate average relative interval length (ARIL)
  res = mean( (Qlu[2,] - Qlu[1,]) / Y)

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' @title S index of average asymmetry degree
#' @description
#' This function calculates the S index of average asymmetry degree given a
#' set of predictions 'X', the target 'Y', and a confidence level of
#' 1 - alpha.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @param alpha specifies 1 - alpha confidence level \[0 < alpha < 1\], Default: 0.5
#' @return S index of average asymmetry degree \[scalar\]
#' @references
#' Xiong, L., N. Wan, X. Wei, and K. O'Connor (2009), Indices
#' for assessing the prediction bounds of hydrological models
#' and application by generalised likelihood uncertainty
#' estimation, Hydrological Sciences Journal, 54:5, 852-871.
#' @rdname avg.asym.deg.S
#' @export
avg.asym.deg.S <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate S index of average asymmetry degree
  res = mean((Qlu[2,]-Y)/(Qlu[2,] - Qlu[1,]))

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' @title T index of average asymmetry degree
#' @description
#' This function calculates the T index of average asymmetry degree given a
#' set of predictions 'X', the target 'Y', and a confidence level of
#' 1 - alpha.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @param alpha specifies 1 - alpha confidence level \[0 < alpha < 1\], Default: 0.5
#' @return T index of average asymmetry degree \[scalar\]
#' @references
#' Xiong, L., N. Wan, X. Wei, and K. O'Connor (2009), Indices
#' for assessing the prediction bounds of hydrological models
#' and application by generalised likelihood uncertainty
#' estimation, Hydrological Sciences Journal, 54:5, 852-871.
#' @rdname avg.asym.deg.T
#' @export
avg.asym.deg.T <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate T index of average asymmetry degree
  res = mean((abs( (Qlu[2,]-Y)^3 + (Qlu[1,]-Y)^3 ) / (Qlu[2,] - Qlu[1,])^3)^1/3)

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' @title Average Deviation Amplitude
#' @description
#' This function calculates the average deviation amplitude given a
#' set of predictions 'X', the target 'Y', and a confidence level of
#' 1 - alpha.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @param alpha specifies 1 - alpha confidence level \[0 < alpha < 1\], Default: 0.5
#' @return average deviation amplitude \[scalar\]
#' @references
#' Xiong, L., N. Wan, X. Wei, and K. O'Connor (2009), Indices
#' for assessing the prediction bounds of hydrological models
#' and application by generalised likelihood uncertainty
#' estimation, Hydrological Sciences Journal, 54:5, 852-871.
#' @rdname avg.dev.amp
#' @export
avg.dev.amp <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate average deviation amplitude (D)
  res = mean(abs(0.5*(Qlu[2,] - Qlu[1,])-Y))

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' @title Average Relative Deviation Amplitude
#' @description
#' This function calculates the average relative deviation amplitude given
#' a set of predictions 'X', the target 'Y', and a confidence level of
#' 1 - alpha.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @param alpha specifies 1 - alpha confidence level \[0 < alpha < 1\], Default: 0.5
#' @return average relative deviation amplitude \[scalar\]
#' @references
#' Xiong, L., N. Wan, X. Wei, and K. O'Connor (2009), Indices
#' for assessing the prediction bounds of hydrological models
#' and application by generalised likelihood uncertainty
#' estimation, Hydrological Sciences Journal, 54:5, 852-871.
#' @rdname avg.rel.dev.amp
#' @export
avg.rel.dev.amp <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate average relative deviation amplitude (RD)
  res = mean(abs(((0.5*(Qlu[2,] - Qlu[1,]))/Y)-1))

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' @title Coverage Probability Plot (CPP)
#' @description
#' This function calculates the coverage probability plot (CPP)
#' given a set of predictions 'X' and the target 'Y'.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X predictions matrix \[N x D\]
#' @return list containing:
#' \itemize{
#'  \item{"fq"}{forecast quantiles(0,0.1,...,1)}
#'  \item{"tq"}{theoretical quantiles(0,0.1,...,1)}
#'  \item{"cpp_mse"}{mean square error between forecast and theoretical quantiles}
#' }
#' @references
#' F. Laio and S. Tamea (2007), Verification tools for probabilistic
#' forecasts of continuous hydrological variables,
#' Hydrol. Earth Syst. Sci., 11, pp. 1267-1277.
#'
#' B. Renard, D. Kavetski, G. Kuczera, M. Thyer, and S.W. Franks (2009), Understanding
#' predictive uncertainty in hydrologic modeling: the challenge of identifying input
#' and structural errors. Water. Resour. Res., 46, p. W05521.
#' @rdname cpp
#' @export
cpp <- function(Y, X){


  homemade_ecdf = function(dat){

    n = length(dat)

    v_sort_dat = sort(dat)
    v_uni_dat = unique(v_sort_dat)
    nu = length(v_uni_dat)

    v_dat_ecdf = matrix(0,nrow=nu,ncol=1)

    for(i in 1:nu){

      current_dat = v_uni_dat[i]
      v_dat_ecdf[i] = sum(v_sort_dat <= current_dat)/n

    }

    v_x = rbind(v_uni_dat[1], as.matrix(v_uni_dat))
    v_f = rbind(0, v_dat_ecdf)


    list (v_f = v_f, v_x=v_x)

  }

  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  z = matrix(n, nrow=n, ncol=1)

  for(i in 1:n){

    #' evaluate the cumulative distribution function of the predictions
    #' in correspondence to each observed value (i.e., calculate ecdf of
    #' forecasts for each time point (1:n) and evaluate at each
    #' observation)

    #' calculate ecdf
    tmp = homemade_ecdf(X[i,])
    F = tmp$v_f
    V = tmp$v_x


    # calculate probability of observation
    probs = (V <= Y[i])

    # evaluate ecdf for probability of observation

    # catch 0 probabilities
    if( sum(probs) == 0){

      z[i] = 0

    }

    else{

      z[i] = max(F[probs])

    }

  }

  #' sort z and find rank of each z[i] then calculate the empirical
  #' cumulative distribution of z, i.e., rank[i]/n for i in 1:n

  # S = sort(z)
  # R = match(z, S)

  # ecdf of z
  # Ri_n = R/n
  Ri_n = rank(z, ties.method = "first") / n

  tq = Ri_n # theoretical quantiles
  fq = z # forecast quantiles

  # tq = seq(from = 0, to = 1, length.out = n) # theoretical quantiles
  # fq = sort(z) # forecast quantiles

  cpp_mse = mean((tq - fq)^2) # coverage probability plot mse
  alpha_index = 1 - 2*mean(abs(tq-fq)) # alpha index (Renard et al., 2009)
  eta_index = 1 - (sum(fq==0 | fq==1)/n) # eta index (Renard et al., 2009)

  plot(fq,tq,xlim=c(0,1),ylim=c(0,1))
  abline(0,1)


  list( fq = fq, tq = tq,
        cpp_mse = cpp_mse, alpha_index = alpha_index, eta_index = eta_index)


} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the Liu-Mean Efficiency (LME), a modified version of
#' the Nash-Sutcliffe (NSE) and Kling-Gupta Efficiency (KGE), given a vector of
#' predictions 'X' and targets 'Y'.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - prediction vector  [N x 1]
#'
#' Output:
#' @param res - Liu-Mean Efficiency (LME) [-inf < scalar <= 1]

#' References:
#'
#' Liu, D. (2020), A rational performance criterion for hydrological model, J. Hydrol.,
#' 580, 125488.
#'
#' Created on: Feb. 5, 2021 by JMQ
#' Updated on: Feb. 5, 2021 by JMQ
#'
#' Usage:
#'
#'

lme <- function(Y, X){

  X = as.numeric(X)
  Y = as.numeric(Y)

  k1 = cor(Y,X)*(var(X)/var(Y))
  beta = mean(X)/mean(Y)

  #' Liu-Mean Efficiency (LME)
  res = 1 - sqrt( (k1-1)^2 + (beta-1)^2 )

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the mean interval score (meanIS)
#' given a set of predictions 'X', the target 'Y', and a
#' confidence level of 1 - alpha.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - predictions matrix [N x D]
#' @param alpha - specifies 1 - alpha confidence level [0 < alpha < 1]
#'
#' Output:
#' @param res - mean interval score [scalar]

#' References:
#'
#' T. Gneiting and  A.E. Raferty (2007), Strictly proper score rules,
#' prediction, and estimation, Journal of the American Statistical
#' Association, 102, pp. 359-378.

#' F. Bourgin et al. (2015), Transferring global uncertainty estimates
#' from gauged to ungauged catchments, Hydrol. Earth Syst. Sci., 19, pp.
#' 2535-2546.

#' C. Wan et al. (2014), Probabiistic forecasting of wind power
#' generation using extreme learning machine, IEEE Transactions on Power
#' Systems, 29, pp. 1033-1044.
#'
#' Created on: Aug. 24, 2019 by JMQ
#' Updated on: Aug. 24, 2019 by JMQ
#'
#' Usage:
#'
#'

meanIS <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)

  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  WPI = Qlu[2,] - Qlu[1,]

  S = matrix(n, nrow=n, ncol=1)

  for(i in 1:n){


    if(Y[i] < Qlu[1,i]){ # beneath lower pi

      S[i] = WPI[i] + (2/alpha)*(Qlu[1,i] - Y[i])

    }

    else if(Y[i] >= Qlu[1,i] & Y[i] <= Qlu[2,i]){ # within pi

      S[i] = WPI[i]

    }

    else if(Y[i] > Qlu[2,i]){ # above upper pi

      S[i] = WPI[i] + (2/alpha)*(Y[i] - Qlu[2,i])

    }


  }

  res = mean(S) # mean interval score

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the mean prediction interval (MPI) given a
#' set of predictions 'X', the target 'Y', and a confidence level of
#' 1 - alpha.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param X - predictions matrix [N x D]
#' @param alpha - specifies 1 - alpha confidence level [0 < alpha < 1]
#'
#' Output:
#' @param res - mean prediction interval [scalar]

#' References:
#'
#' Shrestha, D.L., N. Kayastha, D. Solomatine, and R. Price (2014),
#' Encapsulation of parametric uncertainty statistics by various
#' predictive machine leanring models: MLUE method, J. Hydroinformatics,
#' 16, 95-113.
#'
#' Created on: Oct. 16, 2018 by JMQ
#' Updated on: Oct. 16, 2018 by JMQ
#'
#' Usage:
#'
#'

mpi <- function(X, alpha=0.5){

  X = as.matrix(X)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate mean prediction interval (MPI)
  res = mean(Qlu[2,] - Qlu[1,])

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the normalized uncertainty efficiency (NUE)
#' given a set of predictions 'X', the target 'Y', a weighting factor,
#' and a confidence level of 1 - alpha.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - predictions matrix [N x D]
#' @param w - weighting factor [0 < w < 1]
#' @param alpha - specifies 1 - alpha confidence level [0 < alpha < 1]
#'
#' Output:
#' @param res - normalized uncertainty efficiency [0 < scalar < 1]

#' References:
#'
#' Dogulu, N., P. Lopez Lopez, D.P. Solomatine, A.H. Weerts and
#' D.L. Shreshta (2015), Estimation of predictive hydrologic uncertainty
#' using the quantile regression and UNEEC methods and their comparison
#' on constrasting catchments, Hydrol. Earth. Syst. Sci., 19, 3181-3201.
#'
#' Created on: Oct. 16, 2018 by JMQ
#' Updated on: Oct. 16, 2018 by JMQ
#'
#' Usage:
#'
#'

nue <- function(Y, X, w=1, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate prediction interval coverage probability (PICP)
  pcp = sum(Y >= Qlu[1,] & Y <= Qlu[2,])/length(Y)

  #' calculate average relative interval length (ARIL)
  arl = mean( (Qlu[2,] - Qlu[1,]) / Y)

  #' calculate normalized uncertainty efficiency (NUE)
  res = pcp / (w * arl)

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the prediction interval coverage probability
#' (PICP) given a set of predictions 'X', the target 'Y', and a
#' confidence level of 1 - alpha.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - predictions matrix [N x D]
#' @param alpha - specifies 1 - alpha confidence level [0 < alpha < 1]
#'
#' Output:
#' @param res - prediction interval coverage probability [0 < scalar < 1]

#' References:
#'
#' Shrestha, D.L., N. Kayastha, D. Solomatine, and R. Price (2014),
#' Encapsulation of parametric uncertainty statistics by various
#' predictive machine leanring models: MLUE method, J. Hydroinformatics,
#' 16, 95-113.
#'
#' Created on: Oct. 16, 2018 by JMQ
#' Updated on: Oct. 16, 2018 by JMQ
#'
#' Usage:
#'
#'

picp <- function(Y, X, alpha=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)

  #' get lower and upper qualtiles at confidence level
  Qlu = apply(X, 1, quantile, probs = c(alpha/2, 1 - alpha/2),
              na.rm = TRUE)

  #' calculate prediction interval coverage probability (PICP)
  res = sum(Y >= Qlu[1,] & Y <= Qlu[2,])/length(Y)

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the p-quantile loss (QLp) or p-risk given a set of predictions
#' 'X', the target 'Y', and a quantile 'p' (0, 1).
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - predictions matrix [N x D]
#' @param p - quantile [0 < p < 1]
#'
#' Output:
#' @param res - p-quantile loss [scalar]

#' References:
#'
#' Seeger, M.W., Salinas, D., Flunkert, V., 2016. Bayesian Intermittent Demand Forecasting
#' for Large Inventories, in: Lee, D.D., Sugiyama, M., Luxburg, U. V, Guyon, I., Garnett, R.
#' (Eds.), Advances in Neural Information Processing Systems 29. Curran Associates, Inc.,
#' pp. 4646-4654.
#'
#' #' Li, S., Jin, X., Xuan, Y., Zhou, X., Chen, W., Wang, Y.-X., Yan, X., 2019. Enhancing
#' the Locality and Breaking the Memory Bottleneck of Transformer on Time Series Forecasting.
#' arXiv.
#'
#' D. Salinas, V. Flunkert, J. Gasthaus et al., 2019, DeepAR: Probabilistic forecasting
#' with autoregressive recurrent networks. International Journal of Forecasting,
#' https://doi.org/10.1016/j.ijforecast.2019.07.001.
#'
#'
#' Created on: Jan. 18, 2020 by JMQ
#' Updated on: Jan. 18, 2020 by JMQ
#'
#' Usage:
#'
#'

qlp <- function(Y, X, p=0.5){

  X = as.matrix(X)
  Y = as.matrix(Y)
  pvec = matrix(p, nrow = length(Y),1)
  onevec = matrix(1, nrow = length(Y),1)

  #' get quantile
  Qp = apply(X, 1, quantile, probs = p,
             na.rm = TRUE)

  #' calculate p-quantile loss (QLp) (Li et al., 2019) also known as p-risk in Seeger et al. (2016)
  #' and Salinas et al. (2019)
  #'
  #' Note that Eqns from Seeger et al. (2016) and Li et al. (2019) result in the same answer while
  #' Salinas et al. (2019) differs

  res = 2 * sum( (pvec - (Y <= Qp)) * (Y - Qp) ) / sum(Y) # Li et al. (2019)

  # res = 2 * sum( (Y - Qp) * ( pvec*(Y > Qp) - (onevec - pvec)*(Y <= Qp)) ) / sum(Y) # Seeger et al. (2016) (same result as Li et al. (2019))
  # res = 2 * sum( (Qp - Y) * ( pvec*(Qp > Y) - (onevec - pvec)*(Qp <= Y)) ) / sum(Y) # Salinas et al. (2019) (gives different results than the above)

  return(res)

} # EOF

# ------------------------------------------------------------------------------

#' This function calculates the non-parametric Kling-Gupta Efficiency (RNP), a
#' modified version of the Kling-Gupta Efficiency (KGE), given a vector of
#' predictions 'X' and targets 'Y'.
#'
#'
#' Inputs (must be a matrix (i.e., ?as.matrix):
#' @param Y - target vector (response variable) [N x 1]
#' @param X - prediction vector  [N x 1]
#'
#' Output:
#' @param res - non-parametric Kling-Gupta Efficiency (RNP) [-inf < scalar <= 1]

#' References:
#'
#' Pool, S., Vis, M., Seibert, J. (2018), Evaluating model performance: towards a
#' non-parametric variant of the Kling-Gupta effficiency, Hydrological Sci. J.,
#' 63 (13-14), 1941-1953.
#'
#' Created on: Feb. 5, 2021 by JMQ
#' Updated on: Feb. 5, 2021 by JMQ
#'
#' Usage:
#'
#'

rnp <- function(Y, X){

  X = as.numeric(X)
  Y = as.numeric(Y)

  # mean of predictions and target
  mean.X = mean(X)
  mean.Y = mean(Y)

  # normalized flow duration curves
  fdc.X = sort(X / (mean.X * length(X)))
  fdc.Y = sort(Y / (mean.Y * length(Y)))

  # alpha component
  rnp.alpha = 1 - 0.5 * sum(abs(fdc.X - fdc.Y))

  # beta component
  rnp.beta = mean.X / mean.Y

  # r component
  rnp.r = cor(X, Y, method="spearman")

  #' non-parametric Kling-Gupta Efficiency (RNP)
  res = 1 - sqrt( (rnp.alpha-1)^2 + (rnp.beta-1)^2 + (rnp.r - 1)^2)

  return(res)

} # EOF
