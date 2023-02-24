#' @title Calibrate (train) a data-driven model (DDM)
#' @description
#' This function calibrates a data-driven model (DDM) for a given input-
#' output dataset.
#' @param yc calibration target \[N x 1\]
#' @param xc calibration inputs \[N x D\]
#' @param ddm data-driven model \[string: 'spov','knnr','grnn','rrf'\]
#' @param ddm_param data-driven model hyper-parameter(s)
#' \itemize{
#'  \item{"spov"}{ddm_param\[1\] = model order \[1 < integer <= 3\]}
#'  \item{"knnr"}{ddm_param\[1\] = no. of nearest neighbours \[1 < integer < sample size-1\]}
#'  \item{"grnn"}{no ddm_param as kernel bandwidth determined automatically}
#'  \item{"rrf"}{ddm_param\[1\] = ntrees \[(1,500)\]}
#' }
#' @return integer, list, etc. that contains calibrated model structure
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @references
#' Quilty, J., J. Adamowski, B. Khalil, and M. Rathinasamy (2016), Bootstrap rank-
#' ordered conditional mutual information (broCMI): A nonlinear input variable
#' selection method for water resources modeling, Water Resour. Res., 52,
#' doi:10.1002/2015WR016959.
#'
#' Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
#' based hydrological and water resources forecasting models for real-world
#' applications with best practices and a new forecasting framework, J. Hydrol.,
#' doi:10.1016/j.jhydrol.2018.05.003.
#' @rdname calibrateDDM
#' @export
calibrateDDM <- function(yc,xc,ddm,ddm_param){


  if(ddm != 'rrf'){

    inputs_c = as.matrix(xc)

  }

  else{

    inputs_c = as.data.frame(xc)

  }


  switch(ddm,

         # sparse p-th order Volterra Series model
         spov={

           mdl_params = spoV_estimate(yc,inputs_c,
                                      ddm_param) # estimate parameters

         },

         # k nearest neighbours regression (KNNR)

         knnr={

           mdl_params = ddm_param # use specified k


         },

         # general regression neural network (GRNN)

         grnn={

           #' perform k-fold cross-validation to select best kernel
           #' bandwith parameter

           k = 5
           n = length(yc)
           inds = seq(1,n,by=1)
           rinds = sample(n)
           folds = cut(inds, breaks = k, labels=FALSE)
           tmp_params = vector(mode="numeric",k)
           msqerr = vector(mode="numeric",k)

           for(i in 1:k){

             tstinds = rinds[which(folds==i,arr.ind=TRUE)]

             tmp_params[i] = grnn_estimate(yc[-tstinds],
                                           inputs_c[-tstinds,],
                                           yc[tstinds],
                                           inputs_c[tstinds,])

             p = grnn_predict(yc, inputs_c, inputs_c, tmp_params[i])

             msqerr[i] = mse(obs = as.matrix(yc), sim = p)

           }

           # select grnn kernel bandwidth that provides lowest MSE over
           # calibration set

           mdl_params =  tmp_params[which(msqerr == min(msqerr))]


         },

         # regularized random forests regression (RRF)

         rrf={

           # estimate parameters
           mdl_params = RRF(y=yc, x=inputs_c,
                            ntree=ddm_param, # usually set between 128 and 1000, 500 is typical
                            mtry = max(floor(ncol(inputs_c)/3), 1), # usually set to max(floor(ncol(x)/3), 1)
                            flagReg=1, # set to 0 for standard RF
                            coefReg=1) # set (0,1], smaller values force sparsity

         }

  )

  return(mdl_params)

}

# ------------------------------------------------------------------------------

#' @title Predict using a calibrated data-driven model (DDM)
#' @description
#' This function makes predictions using a calibrated data-driven model
#' (DDM) for a given input vector/matrix.
#' @param x inputs \[M x D\]
#' @param xc calibration inputs \[N x D\]
#' @param yc calibration target \[N x 1\]
#' @param ddm data-driven model \[string: 'spov','knnr','grnn','rrf'\]
#' @param ddm_param data-driven model hyper-parameter(s)
#' \itemize{
#'  \item{"spov"}{ddm_param\[1\] = model order \[1 < integer <= 3\]}
#'  \item{"knnr"}{ddm_param\[1\] = no. of nearest neighbours \[1 < integer < sample size-1\]}
#'  \item{"grnn"}{no ddm_param as kernel bandwidth determined automatically}
#'  \item{"rrf"}{ddm_param\[1\] = ntrees \[(1,500)\]}
#' }
#' @param mdl_params calibrated model structure \[integer, list, etc.\]
#' @return predictions \[N x 1\]
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @references
#' Quilty, J., J. Adamowski, B. Khalil, and M. Rathinasamy (2016), Bootstrap rank-
#' ordered conditional mutual information (broCMI): A nonlinear input variable
#' selection method for water resources modeling, Water Resour. Res., 52,
#' doi:10.1002/2015WR016959.
#'
#' Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
#' based hydrological and water resources forecasting models for real-world
#' applications with best practices and a new forecasting framework, J. Hydrol.,
#' doi:10.1016/j.jhydrol.2018.05.003.
#' @rdname predictDDM
#' @export
predictDDM <- function(x,xc=NULL,yc=NULL,ddm,ddm_param,mdl_params){


  isxcnull = is.null(xc) # check if calibration inputs were provided
  if (isxcnull){ # calibration inputs NOT provided

    xc = x
    yc = yc

  }

  if(ddm != 'rrf'){

    inputs_c = as.matrix(xc)
    inputs = as.matrix(x)

  }

  else{

    inputs_c = as.data.frame(xc)
    inputs = as.data.frame(x)

  }


  switch(ddm,

         # sparse p-th order Volterra Series model
         spov={

           pred = spoV_predict(inputs,mdl_params,ddm_param) # make predictions

         },

         # k nearest neighbours regression (KNNR)

         knnr={

           pred = knnr_predict(yc,inputs_c,
                               inputs,ddm_param,"kd_tree") # make predictions


         },

         # general regression neural network (GRNN)

         grnn={

           pred = grnn_predict(yc,inputs_c,
                               inputs,mdl_params) # make predictions


         },

         # regularized random forests regression (RRF)

         rrf={

           pred = predict(mdl_params,inputs) # make predictions

         }

  )

  return(pred)

}

# ------------------------------------------------------------------------------

# require(pracma) # required for Euclidean distance calculations and for optimizing kernel bandwidth (via optimize)
# require(DEoptim) # required for optimizing kernel bandwidth (via DEoptim)
# require(hydroGOF) # required for objective function for optimizing kernel bandwidth

#' @title Calibrate a Generalized Regression Neural Newtork (GRNN)
#' @description
#' This function calibrates a Generalized Regression Neural Network
#' model's kernel bandwidth 'sigma' based on the radial basis
#' function (RBF, also known as Gaussian) kernel using the 'optimize'
#' function from the 'stats' package for a set of calibration data
#' pairs (Y,X) and a set of validation data pairs (Yv,Xv).
#' @param Y target vector (response variable) \[N x 1\]
#' @param X input matrix (explanatory variables) \[N x D\]
#' @param Yv validation target vector \[Nv x 1\]
#' @param Xv validation input matrix for generating validation predictions \[Nv x D\]
#' @return optimized kernel bandwith \[scalar\]
#' @details
#' Note that the number of distance calculations required for estimating the RBF
#' kernels (parameters) in the GRNN model grows O(N^2).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @references
#'  D.F. Specht (1991), A general regression neural network, IEEE Transactions
#'  on Neural Networks, 2, 6, pp. 568-576, DOI: 10.1109/72.97934.
#'
#'  T.I. Harrold, A. Sharma, and S. Sheather (2001), Selection of a kernel
#'  bandwidth for measuring dependence in hydrologic time series using the
#'  mutual information criterion, Stoch. Environ. Res. Risk Assess., 15, pp.
#'  310-324.
#' @rdname grnn_estimate
#' @export
grnn_estimate <- function(Y,X,Yv,Xv){

  Y = as.matrix(Y)
  X = as.matrix(X)
  Yv = as.matrix(Yv)
  Xv = as.matrix(Xv)

  # set kernel bandwidth interval sufficiently wide (adjust if necessary)
  bw_interval = c(10^-5, 10^5)


  # N = nrow(X) # number of points in calibration set
  # Nd = ncol(X) # number of inputs

  # calculate kernel matrix
  D = pdist2(X,Xv) # Euclidean distance (see Eq. 4 in Specht, 1991)

  #' sub-function 1: calculate GRNN predictions for given target (Y), distance
  #' matrix (D), and kernel bandwidth (bw)

  get_grnn_prediction <- function(Y1,D1,bw1){

    N1 = ncol(D1)
    K1 = exp(-(D1^2) / (2 * (bw1^2))) # RBF kernel function (see Eq. 5 in Specht, 1991)

    P = matrix(0,N1,1)
    for (i in 1:N1){

      num = sum(Y1 * K1[,i]) # (see numerator in Eq. 5 in Specht, 1991)
      den = sum(K1[,i]) # (see denominator in Eq. 5 in Specht, 1991)
      den = max(den, 1e-6)
      P[i,] = num/den # GRNN prediction

    }

    get_grnn_prediction <- P # return predictions

  }

  #' sub-function 2: calculate performance of GRNN model predictions (Pv) against
  #' validation target (Yv)

  eval_grnn_perf <- function(bw2,Y2,D2,Yv2){

    Pv = get_grnn_prediction(Y2,D2,bw2) # get GRNN prediction

    eval_grnn_perf <- mse(as.matrix(Yv2),Pv) # evaluate objective function

  }

  # optimize kernel bandwith using 'optimize' or 'DEoptim' uncomment one or the other
  # bw_opt = optimize(eval_grnn_perf,bw_interval,Y,D,Yv) # Line Search optimizer
  #
  # grnn_estimate <- bw_opt$minimum # return optimized kernel bandwidth

  #' uncomment lines below (and 'bw_opt' and 'grnn_estimate' lines above) here
  #' to use Differential Evolution (DE) -based optimizer 'DEoptim'

  bw_grr = ((4/(nrow(X)+2))^(1/(nrow(X)+4)))*
    nrow(X)^(-1/(ncol(X)+4)) # GRR bandwidth
  bw_opt = DEoptim(eval_grnn_perf,
                   bw_grr,Y,D,Yv,
                   lower=bw_interval[1],
                   upper=bw_interval[2]) # Differential Evolution optimzer

  grnn_estimate <- bw_opt$optim$bestmem



} # EOF

# ------------------------------------------------------------------------------

# require(pracma) # required for Euclidean distance calculations

#' @title Predict using a calibrated Generalized Regression Neural Network (GRNN) model
#' @description
#' This function calculates predictions for a set of model inputs 'Xtest'
#' using a Generalized Regression Neural Network model for a given kernel
#' bandwidth 'sigma' and using the radial basis function (RBF, also known
#' as Gaussian) kernel.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X input matrix (explanatory variables) \[N x D\]
#' @param Xtest test input matrix for generating predictions \[Ntest x D\]
#' @param bw RBF kernel bandwidth \[scalar\], Default: ((4/(nrow(X) + 2))^(1/(nrow(X) + 4))) * N^(-1/(ncol(X) + 4))
#' @return prediction for input matrix Xtest \[Ntest x 1\]
#' @details Default value of the `bw` input argument set to the Gaussian Reference Rule (Harrold et al., 2001),
#' see below in the references.
#' Note that the number of distance calculations required for estimating the RBF
#' kernels (parameters) in the GRNN model grows O(N^2).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @references
#'  D.F. Specht (1991), A general regression neural network, IEEE Transactions
#'  on Neural Networks, 2, 6, pp. 568-576, DOI: 10.1109/72.97934.
#'
#'  T.I. Harrold, A. Sharma, and S. Sheather (2001), Selection of a kernel
#'  bandwidth for measuring dependence in hydrologic time series using the
#'  mutual information criterion, Stoch. Environ. Res. Risk Assess., 15, pp.
#'  310-324.
#' @rdname grnn_predict
#' @export
grnn_predict <- function(Y,X,Xtest,
                         bw =
                           ((4/(nrow(X)+2))^(1/(nrow(X)+4)))*
                           N^(-1/(ncol(X)+4))){

  Y = as.matrix(Y)
  X = as.matrix(X)
  Xtest = as.matrix(Xtest)

  N = nrow(Xtest) # number of points in Xtest to predict
  D = pdist2(X,Xtest) # Euclidean distance (see Eq. 4 in Specht, 1991)
  K = exp(-(D^2) / (2 * (bw^2))) # RBF kernel function (see Eq. 5 in Specht, 1991)

  P = matrix(0,N,1)
  for (i in 1:N){

    num = sum(Y * K[,i]) # (see numerator in Eq. 5 in Specht, 1991)
    den = sum(K[,i]) # (see denominator in Eq. 5 in Specht, 1991)
    den = max(den, 1e-6)
    P[i,] = num/den # GRNN prediction

  }

  grnn_predict <- P # return predictions

} # EOF

# ------------------------------------------------------------------------------

# require(FNN) # required for nearest neighbour estimation

#' @title Predict using a K-Nearest Neighbors (KNN) regression model
#' @description
#' This function calculates predictions for a set of model inputs 'Xtest'
#' using a k nearest neighbour model whose neighbour searching algorithm can
#' be selected according to different methods. The prediction is calculated
#' by assuming equally weighted predictions based on the k nearest neighbors.
#' @param Y target vector (response variable) \[N x 1\]
#' @param X input matrix (explanatory variables) \[N x D\]
#' @param Xtest test input matrix for generating predictions \[Ntest x D\]
#' @param k number of nearest neighbours to use in prediction \[scalar < N\]
#' @param method `c("kd_tree", "cover_tree", "CR", "brute")`, Default: `"kd_tree"`
#' @return prediction for input matrix Xtest \[Ntest x 1\]
#' @details
#' The cover tree is O(n) space data structure which allows us to
#' answer queries in the same O(log(n)) time as kd tree given a
#' fixed intrinsic dimensionality. Templated code from
#' http://hunch.net/~jl/projects/cover_tree/cover_tree.html is used.
#' The kd tree algorithm is implemented in the Approximate Near
#' Neighbor (ANN) C++ library (see http://www.cs.umd.edu/~mount/ANN/).
#' The exact nearest neighbors are searched in this package.
#' The CR algorithm is the VR using distance 1-x'y assuming x and y
#' are unit vectors. The brute algorithm searches linearly. It is a
#' naive method.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @references
#'  Lall, U., and A. Sharma (1996), A nearest neighbor bootstrap for time series
#'  resampling, Water. Resour. Res., 32, 679-693.
#' @rdname knnr_predict
#' @export
knnr_predict <- function(Y,X,Xtest,k=3,method="kd_tree"){

  # get k nearest neighbours to Xtest

  nn_inds = get.knnx(X,Xtest,k,"kd_tree")$nn.index

  # identify target related to each of the nearest neighbours

  Yr = matrix(Y[nn_inds],nrow(as.matrix(Xtest)),k)

  # assign probabilities to k nearest neighbours

  probs = (1 / k) * rep(1,k)

  # get predictions according to probabilities

  P = rowSums(probs*Yr)

  knnr_predict <- P # return predictions

} # EOF

# ------------------------------------------------------------------------------

# require(pracma) # required for Moore-Penrose generalized inverse function pinv()

#' @title Calibrate (train) a sparse p-th order Volterra series model
#' @description
#' This function calculates a  sparse p-th order multiple input single output
#' (MISO) Volterra series model for given input-output data pairs X and Y,
#' using the Moore-Penrose generalized inverse to estimate the different
#' triangular kernels (see Eq. 8 in Wu and Kareem \[2014\]).
#' @param Y target vector (response variable) \[N x 1\]
#' @param X input matrix (explanatory variables) \[N x D\]
#' @param p order (polynomial degree) of Volterra model \[scalar\]
#' @return p-th order kernels of Volterra model \[Q x 1\]
#' @details
#' NOTE: the number of kernels (parameters) in the Volterra series model grows exponentially
#' based on p (the model order).
#' @examples
#' \dontrun{
#' if(interactive()){
#' # Usage:
#' p=2
#' X=matrix(runif(1000),500,20)
#' Y=matrix(runif(500),500,1)
#' H = spoV_estimate(Y,X,2) # sparse second-order Volterra series kernels
#'  }
#' }
#' @references
#'  T. Wu and A. Kareem (2014), Simulation of nonlinear bridge aerodynamics: A sparse
#'  third-order Volterra model, Journal of Sound and Vibration, 333, 1, pp. 178-188
#'  https://doi.org/10.1016/j.jsv.2013.09.003.
#'
#'  https://stackoverflow.com/questions/49538911/r-given-a-matrix-and-a-power-produce-multiple-matrices-containing-all-unique
#' @rdname spoV_estimate
#' @export
spoV_estimate <- function(Y,X,p){

  N = nrow(as.matrix(X)) # number of samples

  #' create p-th degree design matrices transformed from original inputs X
  #' (see Eqs. 5a - 5d in Wu and Kareem, 2014)

  mat = as.data.frame(X) # convert to data frame for calculating
  const = matrix(rep(1,N),N,1) # constant term (not included in Wu and Kareem, 2014)

  # preallocate space for p = 1:p
  res_mat = list()
  res_mat[[1]] = mat # first order kernels

  if (p > 1) { # triangular kernels higher than first order (p > 2)

    for(i in 2:p) {

      combns = combn(ncol(mat) + i - 1, i) - 0:(i - 1) # "unique" permutations only (see Eq. 8 in Wu and Kareem [2014])
      #colnames(combs) <- apply(combs, 2, paste, collapse = "-") # Just for display of output, we keep info of combinations in colnames

      # p-th order triangular kernels
      res_mat[[i]] = apply(combns, 2, function(x) Reduce(`*`, mat[ ,x])) # multiply the relevant columns

    }

  }

  Xt = do.call(cbind,res_mat) # unroll list into matrix
  Xtc = as.matrix(cbind(const,Xt)) # include constant term

  # calculate spoV kernels using Moore-Penrose generalized inverse pinv()
  H = pinv(Xtc) %*% Y
  # prediction = Xtc %*% H

  spoV_estimate <- H # return parameters

} # EOF

# ------------------------------------------------------------------------------

#' This function generates predictions for a sparse p-th order multiple input single
#' output (MISO) Volterra series model for given input X and parameters H.
#'
#' Inputs:
#' @param X - input matrix (explanatory variables) [N x D] (must be a matrix (i.e., ?as.matrix)
#' @param H - p-th order kernels of Volterra model [Q x 1]
#' @param p - order (polynomial degree) of Volterra model [scalar]
#'
#' Output:
#' @param P - prediction [N x 1]

#' Reference:
#'
#'  T. Wu and A. Kareem (2014), Simulation of nonlinear bridge aerodynamics: A sparse
#'  third-order Volterra model, Journal of Sound and Vibration, 333, 1, pp. 178-188
#'  https://doi.org/10.1016/j.jsv.2013.09.003.
#'
#' Created on: Apr. 4, 2018 by JMQ
#' Updated on: Apr. 25, 2018 by JMQ
#'
#' Usage:
#' p=2
#' X=matrix(runif(1000),500,20)
#' Y=matrix(runif(500),500,1)
#' H = spoV_estimate(Y,X,2); # sparse second-order Volterra series kernels
#' P = spoV_predict(X,H,p); # predictions
#'

spoV_predict <- function(X,H,p){

  N = nrow(as.matrix(X)) # number of samples

  #' create p-th degree design matrices transformed from original inputs X
  #' (see Eqs. 5a - 5d in Wu and Kareem, 2014)

  mat = as.data.frame(X) # convert to data frame for calculating
  const = matrix(rep(1,N),N,1) # constant term (not included in Wu and Kareem, 2014)

  # preallocate space for p = 1:p
  res_mat = list()
  res_mat[[1]] = mat # first order kernels

  if (p > 1) { # triangular kernels higher than first order (p > 2)

    for(i in 2:p) {

      combns = combn(ncol(mat) + i - 1, i) - 0:(i - 1)  # "unique" permutations only (see Eq. 8 in Wu and Kareem [2014])
      #colnames(combs) <- apply(combs, 2, paste, collapse = "-") # Just for display of output, we keep info of combinations in colnames

      # p-th order triangular kernels
      res_mat[[i]] = apply(combns, 2, function(x) Reduce(`*`, mat[ ,x])) # multiply the relevant columns

      if(N==1){ res_mat[[i]] = matrix(res_mat[[i]], ncol=length(res_mat[[i]]))} # catch case where only single observation

    }

  }

  Xt = do.call(cbind,res_mat) # unroll list into matrix
  Xtc = as.matrix(cbind(const,Xt)) # include constant term

  # calculate spoV prediction
  Yhat = Xtc %*% H

  spoV_predict <- Yhat # return prediction

} # EOF
