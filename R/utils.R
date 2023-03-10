#' @title Column Bind Target and Input Variables
#' @description
#' This utility function column binds target and input variables together in a data frame
#' and removes 'NA' values.
#' @param y target time series \[N x 1\]
#' @param x input time series \[N x D\]
#' @return resulting matrix \[N x D+1\]
#' @rdname bindYX
#' @export
bindYX <- function(y,x){

  yx = cbind(y,x)

  return(yx[stats::complete.cases(yx), ]) # return only rows w/o 'NA' value

} # EOF...

# ------------------------------------------------------------------------------

#' @title Column Bind Target, Input and Auxiliary Variables
#' @description
#' This function column binds target, input variables, and auxiliary variables together
#' in a data frame and removes 'NA' values.
#' @param y target time series \[N x 1\]
#' @param x input time series \[N x D\]
#' @param z auxiliary variables \[N x M\]
#' @return resulting matrix \[N x D+M+1\]
#' @rdname bindYXZ
#' @export
bindYXZ <- function(y,x,z){

  yxz = cbind(y,x,z)

  return(yxz[stats::complete.cases(yxz), ]) # return only rows w/o 'NA' value

} # EOF...

# ------------------------------------------------------------------------------

#' @title Scale Columns of a Matrix Based on Another Matrix
#' @description
#' This function scales an input matrix 'X2' between values 'a' and 'b' with the
#' option to exclude NA values from the scaling (see example below) according to
#' the content (max and min values) of scaling matrix 'X1'.
#' @param X1 input matrix whose properties are used in scaling \[N x D\]
#' @param X2 input matrix to scaled based on properties of X1 \[M x D\]
#' @param a lower limit value for scaled data \[scalar\]
#' @param b upper limit value for scaled data \[scalar\]
#' @param ... further arguments to the `min()` and `max()` functions
#' @return scaled input matrix X2 \[M x D\]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X1 = matrix(runif(16),4,4) # square matrix
#'  X1[2:3,3:4] = NA           # inject some NA's into X1
#'  X2 = matrix(rnorm(16),4,4) # another square matrix
#'
#'  # Compute statistics with X1, apply scaling to X2
#'  X3 = scale_ab(X1,X2,-1,1,na.rm=TRUE)
#'
#'  # NOTE: Typically, X1 will be the training set input features,
#'  #       X2 will be the validation/test set.
#'  }
#' }
#' @references
#' https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points
#' @rdname scale_ab
#' @export
scale_ab <- function(X1,X2,a,b,...){

  X1 = as.matrix(X1)
  X2 = as.matrix(X2)

  X1min = apply(X1, MARGIN = 2,FUN = function(X1) (min(X1,...)))
  X1max = apply(X1, MARGIN = 2,FUN = function(X1) (max(X1,...)))

  Xs = t(t(sweep(X2, 2, X1min)) / (X1max - X1min)) * (b - a) + a

  # check for any columns that are completely NA and replace with Os

  test.nona = colSums(is.na(Xs)) == 0

  if (any(test.nona == FALSE)){

    Xs[,as.integer(which(test.nona == FALSE))] = 0

  }


  scale_ab <- Xs

} # EOF

# ------------------------------------------------------------------------------

#' @title De-scale Columns of a Matrix Based on Another Matrix
#' @description
#' This function de-scales an input matrix 'X2' from between values 'a' and 'b' with
#' the option to disclude NA values from the scaling (see example below) according to
#' the content (max and min values) of scaling matrix 'X1'. In other words, this function
#' undoes what `scale_ab` does.
#' @param X1 input matrix whose properties are used in de-scaling \[N x D\]
#' @param X2 input matrix to de-scale based on properties of X1 \[M x D\]
#' @param a lower limit value for scaled data \[scalar\]
#' @param b upper limit value for scaled data \[scalar\]
#' @param ... further arguments to the `min()` and `max()` functions
#' @return de-scaled input matrix X2 \[M x D\]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X1 = matrix(runif(16),4,4)
#'  X2 = matrix(rnorm(16),4,4)
#'  X1[2:3,3:4] = NA
#'  X3 = de_scale_ab(X1,X2,-1,1,na.rm=TRUE)
#'  }
#' }
#' @rdname de_scale_ab
#' @references
#' https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points
#' @export
de_scale_ab <- function(X1,X2,a,b,...){

  X1 = as.matrix(X1)
  X2 = as.matrix(X2)

  X1min = apply(X1, MARGIN = 2,FUN = function(X1) (min(X1,...)))
  X1max = apply(X1, MARGIN = 2,FUN = function(X1) (max(X1,...)))

  Xds = t(t(sweep(X2 - a,2,(X1max - X1min) / (b - a),'*')) + X1min)

  #Xs = t(t(sweep(X2, 2, X1min)) / (X1max - X1min)) * (b - a) + a

  # check for any columns that are completely NA and replace with Os

  test.nona = colSums(is.na(Xds)) == 0

  if (any(test.nona == FALSE)){

    Xds[,as.integer(which(test.nona == FALSE))] = 0

  }


  de_scale_ab <- Xds

} # EOF

# ------------------------------------------------------------------------------

#' @title Standardize the Columns of a Matrix Based on Another Matrix
#' @description
#' This function standardizes an input matrix 'X2' with the option to exclude
#' NA values from the scaling (see example below) according to the content
#' (mean and standard deviation values) of standardizing matrix 'X1'.
#' @param X1 input matrix whose properties are used in standardization \[N x D\]
#' @param X2 input matrix to standardize based on properties of X1 \[M x D\]
#' @param ... further arguments to `mean()` and `std()`
#' @return standardized input matrix X2 \[M x D\]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X1 = matrix(runif(16),4,4)
#'  X2 = matrix(rnorm(16),4,4)
#'  X1[2:3,3:4] = NA
#'  X3 = standardize(X1,X2,na.rm=TRUE)
#'  }
#' }
#' @references
#' https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points
#' @rdname standardize
#' @export
standardize <- function(X1,X2,...){

  X1 = as.matrix(X1)
  X2 = as.matrix(X2)

  X1mean = apply(X1, MARGIN = 2,FUN = function(X1) (mean(X1,...)))
  X1sd = apply(X1, MARGIN = 2,FUN = function(X1) (stats::sd(X1,...)))

  Xs = t(t(sweep(X2, 2, X1mean)) / X1sd)

  # check for any columns that are completely NA and replace with Os

  test.nona = colSums(is.na(Xs)) == 0

  if (any(test.nona == FALSE)){

    Xs[,as.integer(which(test.nona == FALSE))] = 0

  }


  standardize <- Xs

} # EOF

# ------------------------------------------------------------------------------

#' @title De-standardize the Columns of a Matrix Based on Another Matrix
#' @description
#' This function de-standardizes an input matrix 'X2' with the option to exclude
#' NA values from the de-scaling (see example below) according to the content
#' (mean and standard deviation values) of standardizing matrix 'X1'.
#' @param X1 input matrix whose properties are used in de-standardization \[N x D\]
#' @param X2 input matrix to de-standardize based on properties of X1 \[M x D\]
#' @param ... further arguments to `mean()` and `std()`
#' @return de-standardized input matrix X2 \[M x D\]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X1 = matrix(runif(16),4,4)
#'  X2 = rnorm(rnorm(16),4,4)
#'  X1[2:3,3:4] = NA
#'  X3 = de_standardize(X1,X2,na.rm=TRUE)
#'  }
#' }
#' @references
#' https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points
#' @rdname de_standardize
#' @export
de_standardize <- function(X1,X2,...){

  X1 = as.matrix(X1)
  X2 = as.matrix(X2)

  X1mean = apply(X1, MARGIN = 2,FUN = function(X1) (mean(X1,...)))
  X1sd = apply(X1, MARGIN = 2,FUN = function(X1) (stats::sd(X1,...)))

  Xds = t(t(sweep(X2,2,X1sd,'*')) + X1mean)
  #Xds = t(t(sweep(X2, 2, X1mean)) / X1sd)

  # check for any columns that are completely NA and replace with Os

  test.nona = colSums(is.na(Xds)) == 0

  if (any(test.nona == FALSE)){

    Xds[,as.integer(which(test.nona == FALSE))] = 0

  }


  de_standardize <- Xds

} # EOF

# ------------------------------------------------------------------------------

# This function time lags the target variable based on the lead time and
# required time delay.
#' @keywords internal
lagmatrix <- function(x,max.lag){stats::embed(c(rep(NA,max.lag),x),max.lag+1)}

# Example
# lagmatrix(1:10,2)
#      [,1] [,2] [,3]
# [1,]    1   NA   NA
# [2,]    2    1   NA
# [3,]    3    2    1
# [4,]    4    3    2
# [5,]    5    4    3
# [6,]    6    5    4
# [7,]    7    6    5
# [8,]    8    7    6
# [9,]    9    8    7
# [10,]  10    9    8

# ------------------------------------------------------------------------------

#' @title Get train, validation and test set indices
#' @description
#' This function provides train, validation and test set indices in a sequential fashion.
#' The user must provide the desired validation set size `nval` and test set size `ntst`.
#' @param N number of time series records \[scalar\]
#' @param nval number of validation records \[scalar\]
#' @param ntst number of test records \[scalar\]
#' @return list that holds training, validation, and test indices
#' @rdname partitionInds
#' @export
partitionInds <- function(N,nval,ntst){

  ind_t = seq(N - ntst + 1, N, by=1) # test indices
  ind_v = seq(min(ind_t) - nval, min(ind_t) - 1, by=1) # validation indices
  ind_c = seq(1, min(ind_v) - 1, by=1) # calibration indices

  return(list(indc = ind_c,
              indv = ind_v,
              indt = ind_t))


}

# ------------------------------------------------------------------------------

# This function creates a periodic sequence at the yearly scale using
# either monthly, daily, or hourly dates.  If desired, a phase component
# can also be included.
#' @keywords internal
time2sine = function(dte=seq(as.Date("1910/1/1"), as.Date("1999/12/30"), "days"),
                     phase=0,timescale = "daily"){

  switch(timescale,

         hourly = {

           dy = as.numeric(strftime(
             dte,format = "%j"))

           hr = as.numeric(strftime(
             dte,format = "%H"))/24

           res = sin(2*pi*(phase + dy + hr)/366)

         },

         daily = {

           res = sin(2*pi*(phase + as.numeric(strftime(
             dte,format = "%j")))/366)

         },

         monthly = {

           res = sin(2*pi*(phase + as.numeric(strftime(
             dte,format = "%m")))/12)

         }


  )

  return(res)

} # EOF

# ------------------------------------------------------------------------------

# copy and paste data from R Studio into Excel
#' @keywords internal
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  utils::write.table(x,"clipboard-16384",sep="\t",row.names=row.names,col.names=col.names,...)
}

# ------------------------------------------------------------------------------

# This function time lags the input variables based on the lead time and
# required time delay.

# Input(s):
# @param x: input time series [N x D]
# @param leadtime: forecast lead time [scalar >= 0]
# @param timedelays: number of time delays [vector (of scalar >= 1)]

# Output:
# @param y: input matrix based forecast lead time and time delays [N x D +sum(timedelays)]
#' @keywords internal
xLeadTimeDelay <- function(x,leadtime=1,timedelays=NULL){

  isvec = is.vector(x)
  istdnull = is.null(timedelays)

  # check if timedelays is null and assign time delays based on x being a vector or not
  if(istdnull){
    if(isvec){timedelays = 1}
    else{timedelays = rep(1,ncol(x))}
  }

  # check to ensure there is a time delay for each column of x
  if(!isvec){
    if(length(timedelays) != ncol(x)){
      timedelays = rep(max(timedelays),ncol(x))}
  }

  y = list()

  if (!isvec){

    for(i in 1:ncol(x)){

      tmp = stats::embed(c(rep(NA,leadtime + timedelays[i]),x[,i]),
                  leadtime + timedelays[i])
      tmp = as.data.frame(tmp) # convert to data frame (embed only operates on vectors/matrices)

      if(!is.null(colnames(x))){

        colnames(tmp) <- paste0(colnames(x)[i],'_', 0:(ncol(tmp)-1)) # 0 represents target forecast

      }

      else{

        colnames(tmp) <- paste0('x',i,'_', 0:(ncol(tmp)-1)) # 0 represents target forecast

      }


      y[[i]] = tmp[2:nrow(tmp),(leadtime + 1):(leadtime + timedelays[i])]

    }

    y = do.call(cbind,y)

  }

  else{

    tmp = stats::embed(c(rep(NA,leadtime + timedelays),x),leadtime + timedelays)
    tmp = as.data.frame(tmp) # convert to data frame (embed only operates on vectors/matrices)

    if(!is.null(colnames(x))){

      colnames(tmp) <- paste0(colnames(x)[i],'_', 0:(ncol(tmp)-1)) # 0 represents target forecast

    }

    else{

      colnames(tmp) <- paste0('x',i,'_', 0:(ncol(tmp)-1)) # 0 represents target forecast

    }


    y = tmp[2:nrow(tmp),(leadtime + 1):(leadtime + timedelays)]

  }

  return(y)

} # EOF...

# ------------------------------------------------------------------------------

# This function time lags the target variable based on the lead time and
# required time delay.

# Input(s):
# @param y: time series [N x 1]
# @param leadtime: forecast lead time [scalar >= 0]
# @param timedelay: number of time delays [scalar >= 1]

# Output:
# @param y: forecast lead time and time delay matrix [N x timedelay+1]
#' @keywords internal
yLeadTimeDelay <- function(x,leadtime=1,timedelay=1){


  tmp = stats::embed(c(rep(NA,leadtime + timedelay),x),leadtime + timedelay)

  if(!is.null(colnames(x))){

    colnames(tmp) = paste0(colnames(x),'_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

  }

  else{

    colnames(tmp) = paste0('y','_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

  }

  y = as.data.frame(tmp) # convert to data frame (embed only operates on vectors/matrices)

  return(y[2:nrow(y),c(1,(leadtime + 1):(leadtime + timedelay))])

} # EOF...

# ------------------------------------------------------------------------------

#' @title Lag Input Variables
#' @description
#' This function time lags the input variables based on the lead time and
#' required time delay.
#' @param x input time series \[N x D\]
#' @param leadtime forecast lead time \[scalar >= 0\], Default: 1
#' @param timedelays number of time delays \[vector (of scalar >= 0)\], Default: NULL
#' @return input matrix based forecast lead time and time delays \[N x D + sum(timedelays)\]
#' @details
#' Based on the internal 'xLeadTimeDelay' function (used only for forecasting); See source code. This
#' function can be used for both modeling and forecasting.
#' @rdname xLTD
#' @export
xLTD <- function(x,leadtime=1,timedelays=NULL){

  x = as.matrix(x)
  istdnull = is.null(timedelays)

  # check if timedelays is null then assign time delays based on x being a vector or not
  if(istdnull){

    timedelays = rep(1,ncol(x)) # set time delays to 1 per input

  }

  # check to ensure there is a time delay for each column of x
  if(length(timedelays) != ncol(x)){

    timedelays = rep(max(timedelays),ncol(x)) # set time delays to max specified time delay

  }

  y = list()

  for(i in 1:ncol(x)){

    # case 1 #############################################################

    if( (leadtime == 0) & (timedelays[i] == 0) ){

      tmp = as.data.frame(x[,i])

      if(!is.null(colnames(x))){

        colnames(tmp) <- paste0(colnames(x)[i],'_0') # t=0 represents target forecast

      }

      else{

        colnames(tmp) <- paste0('x',i,'_0') # t=0 represents target forecast

      }

      y[[i]] = tmp # time-delayed matrix for current variable i

    }

    # case 2 #############################################################

    else if( (leadtime == 0) & (timedelays[i] > 0) ){

      # tmp = cbind(x[,i],stats::embed(c(rep(NA,timedelays[i]),x[,i]),
      #                         timedelays[i]))
      tmp = cbind(x[,i],stats::embed(c(rep(NA,timedelays[i]),x[,i]),
                              timedelays[i])[1:(length(x[,i])),])

      if(!is.null(colnames(x))){

        colnames(tmp) <- paste0(colnames(x)[i],'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      else{

        colnames(tmp) = paste0('x',i,'_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      #df = as.data.frame(tmp[1:(nrow(tmp)-1),1:(1 + timedelays[i])]) # convert to data frame (embed only operates on vectors/matrices)
      df = as.data.frame(tmp[,1:(1 + timedelays[i])]) # convert to data frame (embed only operates on vectors/matrices)
      colnames(df) = colnames(tmp)[1:(1 + timedelays[i])]
      y[[i]] = df # time-delayed matrix for current variable i
    }

    # case 3 #############################################################

    else if( (leadtime > 0) & (timedelays[i] == 0) ){

      tmp = stats::embed(c(rep(NA,leadtime + 1),x[,i]),leadtime + 1)

      if(!is.null(colnames(x))){

        colnames(tmp) <- paste0(colnames(x)[i],'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      else{

        colnames(tmp) <- paste0('x',i,'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      df = as.data.frame(tmp[2:nrow(tmp),leadtime + 1]) # convert to data frame (embed only operates on vectors/matrices)
      colnames(df) = colnames(tmp)[leadtime + 1]
      y[[i]] = df # time-delayed matrix for current variable i

    }

    # case 4 #############################################################

    else{ # (leadtime > 0) & (timedelay > 0) )

      tmp = stats::embed(c(rep(NA,leadtime + timedelays[i]),x[,i]),
                  leadtime + timedelays[i])

      if(!is.null(colnames(x))){

        colnames(tmp) <- paste0(colnames(x)[i],'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      else{

        colnames(tmp) <- paste0('x',i,'_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

      }

      df = as.data.frame(tmp[2:nrow(tmp),
                             (leadtime + 1):(leadtime + timedelays[i])]) # convert to data frame (embed only operates on vectors/matrices)
      colnames(df) = colnames(tmp)[(leadtime + 1):(leadtime + timedelays[i])]
      y[[i]] = df # time-delayed matrix for current variable i

    }

  }

  y = do.call(cbind,y)

  return(y)

} # EOF...

# ------------------------------------------------------------------------------

#' @title Lag Target Variables
#' @description
#' This function time lags the target variable based on the lead time and
#' required time delay.
#' @param x time series \[N x 1\]
#' @param leadtime forecast lead time \[scalar >= 0\], Default: 1
#' @param timedelay number of time delays \[scalar >= 0\], Default: 1
#' @return forecast lead time and time delay matrix \[N x timedelay+1\]
#' @details
#' Based on the internal 'yLeadTimeDelay' function (used only for forecasting);See source code. This
#' function can be used for both modeling and forecasting.
#' @rdname yLTD
#' @export
yLTD <- function(x,leadtime=1,timedelay=1){

  # case 1 #############################################################

  if( (leadtime >= 0) & (timedelay == 0) ){

    tmp = as.data.frame(x)

    if(!is.null(colnames(x))){

      colnames(tmp) = paste0(colnames(x),'_0') # t=0 represents target forecast

    }

    else{

      colnames(tmp) = 'y_0' # t=0 represents target forecast

    }

    y = tmp

    return(y)

  }

  # case 2 #############################################################

  else if( (leadtime == 0) & (timedelay > 0) ){

    #tmp = cbind(x,stats::embed(c(rep(NA,timedelay),x),timedelay))
    tmp = cbind(x,
                stats::embed(c(rep(NA,timedelay),x),timedelay)[1:length(x),])

    if(!is.null(colnames(x))){

      colnames(tmp) = paste0(colnames(x),'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

    }

    else{

      colnames(tmp) = paste0('y','_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

    }

    # y = as.data.frame(tmp) # convert to data frame (embed only operates on vectors/matrices)
    # return(y[1:(nrow(y)-1),c(1,2:(1 + timedelay))])

    #df = as.data.frame(tmp[1:(nrow(tmp)-1),1:(1 + timedelay)]) # convert to data frame (embed only operates on vectors/matrices)
    df = as.data.frame(tmp[,1:(1 + timedelay)]) # convert to data frame (embed only operates on vectors/matrices)
    colnames(df) = colnames(tmp)[1:(1 + timedelay)]
    y = df # time-delayed matrix for current variable i

    return(y)

  }

  # case 3 #############################################################

  else{ # (leadtime > 0) & (timedelay > 0) )

    tmp = stats::embed(c(rep(NA,leadtime + timedelay),x),leadtime + timedelay)

    if(!is.null(colnames(x))){

      colnames(tmp) = paste0(colnames(x),'_',0:(ncol(tmp)-1)) # t=0 represents target forecast

    }


    else{

      colnames(tmp) = paste0('y','_', 0:(ncol(tmp)-1)) # t=0 represents target forecast

    }

    df = as.data.frame(tmp[2:nrow(tmp),
                           c(1,(leadtime + 1):(leadtime + timedelay))]) # convert to data frame (embed only operates on vectors/matrices)
    colnames(df) = colnames(tmp)[c(1,(leadtime + 1):(leadtime + timedelay))]
    y = df # time-delayed matrix for current variable i

    return(y)

  }


} # EOF...
