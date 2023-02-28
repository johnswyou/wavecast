#' @title Wavelet Data-Driven Forecasting Framework (WDDFF)
#' @description
#' This function calibrates a model and generates predictions using the
#' Wavelet Data-Driven Forecasting Framework (WDDFF) based on the work in
#' Quilty and Adamowski (2018).
#' @param y target time series \[N x 1\]
#' @param x input time series \[N x D\]
#' @param z auxilliary input time series \[N x C\]
#' @param leadtime forecast lead time \[scalar >= 0\]
#' @param lag_Y number of time delays for target \[scalar >= 1\]
#' @param lag_X number of time delays \[vector (or scalar >= 1)\]
#' @param nval number of validation records \[scalar\]
#' @param ntst number of test records \[scalar\]
#' @param wfm wavelet-based forecasting method \[string: 'none', 'single', 'within', 'across','single-hybrid','within-hybrid','across-hybrid'\]
#' @param wt type of wavelet transform \[string; 'at','modwt'\]
#' @param wavelet scaling filter name \[string; 'haar', 'db2', 'db4', 'coif1', etc.\]
#' @param decomp_level decomposition level \[1 < integer << N/2\]
#' @param max_decomp_level maximum decomposition level \[scalar >= decomp_level\]
#' @param max_wavelet_length maxmimum wavelet filter length \[scalar >= length(wavelet)\]
#' @param ivsm input variable selection (IVS) method
#'   \[string:'none','boruta','rrf','ea_cmi_htc','ea_cmi_tol','knn_cmi_tol','knn_cmi_bi_tol',pmis_bic','pcis_bic'\]
#' @param ivs_param parameters for input variable selection method
#' \itemize{
#' \item `rrf`:
#' \itemize{
#' \item ivs_param\[1\] = ntrees \[integer > 1)\] usually between 128 and 1000
#' \item ivs_param\[2\] = regularization parameter \[0 < scalar <= 1\]}
#' \item `ea_cmi`:
#' \itemize{
#' \item ivs_param\[1\] =  Hampel test criterion threshold \[scalar > 0\]}
#' \item `ea_cmi_tol`:
#' \itemize{
#' \item ivs_param\[1\] = CMI/MI threshold \[0 < scalar <= 1\]
#' \item ivs_param\[2\] = no. of nearest neighbours \[1 < integer < sample size-1\]}
#' \item `knn_cmi_tol`:
#' \itemize{
#' \item ivs_param\[1\] = CMI/MI threshold \[0 < scalar <= 1\]
#' \item ivs_param\[2\] = no. of nearest neighbours \[1 < integer < sample size-1\]}
#' }
#' @param ddm data-driven model \[string: 'spov','knnr','grnn','rrf'\]
#' @param ddm_param data-driven model hyper-parameter(s)
#' - 'spov': ddm_param\[1\] = model order \[1 < integer <= 3\]
#' - 'knnr': ddm_param\[1\] = no. of nearest neighbours \[1 < integer < sample size-1\]
#' - 'grnn': no ddm_param as kernel bandwidth determined automatically
#' - 'rrf': ddm_param\[1\] = ntrees \[(1,500)\]
#' @param scale_inputs scale inputs (predictors) \[0,1\] \[binary\]
#' @param scale_target scale target (predictand) \[0,1\] \[binary\]
#' @param cutoff0 ensure that all predictions > 0 \[binary\]
#' @param light store all function inputs or neglect 'y', 'x', and 'z' \[binary\]
#' @param savefile save file using function name and current Sys.time() \[binary\]
#' @return
#' mdl_rslt \[list\] that contains:
#' \itemize{
#' \item mdl \[list\] that includes all function inputs
#' (this list does not include 'y', 'x', and 'z' if 'light' = `TRUE`)
#' \item rslt \[list\],
#' \itemize{
#' \item predictions,
#' \item target,
#' \item residual,
#' \item n.inds,
#' \item indc,
#' \item indv,
#' \item indt,
#' \item ivs,
#' \item ddm,
#' \item perf_c,
#' \item perf_v,
#' \item perf_t}
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
#' @rdname wddff
#' @export
wddff <- function(y, x=NULL, z=NULL, leadtime=1, lag_Y=1, lag_X=1,
                  nval=1, ntst=1,
                  wfm='single_hybrid', wt='modwt', wavelet='haar',
                  decomp_level=1, max_decomp_level=decomp_level, max_wavelet_length=length(fastWavelets::wavelet_filter(wavelet)),
                  ivsm = 'none', ivs_param = NULL,
                  ddm='spov', ddm_param=1,
                  scale_inputs=FALSE, scale_target=FALSE,
                  cutoff0=FALSE, light=TRUE, savefile=FALSE){

  ###########################################################################
  # HIGH LEVEL SETTINGS
  ###########################################################################

  # set.seed(123) # make reproducible

  y = as.matrix(y)
  if( !is.null(x) ){ x = as.matrix(x)}
  if( !is.null(z) ) {z = as.matrix(z)}

  # store model settings

  mdl = c(as.list(environment()), list())

  if( light ){

    mdl$y = NULL
    mdl$x = NULL
    mdl$z = NULL

  }


  ###########################################################################
  # CREATE INPUT-OUTPUT (IO) DATASET
  ###########################################################################

  # create IO dataset (lags input-output dataset and peforms wavelet-
  # decomposition, if needed)

  IO.pd = prepareData(y=y,x=x,z=z,leadtime,lag_Y,lag_X,
                      wfm,wt,wavelet,decomp_level)

  # adjust IO dataset for boundary-condition (very important for wavelet-
  # based models)

  IO.od = organizeIOData(IO.pd,max_decomp_level,max_wavelet_length)
  nbc = ((2^max_decomp_level) - 1) * (max_wavelet_length - 1) + 1
  mdl$num_boundary_coefs = nbc

  ###########################################################################
  # PARTITION INPUT-OUTPUT (IO) DATASET
  ###########################################################################

  # partition IO dataset

  N = nrow(IO.od[[1]])
  IO.inds = partitionInds(N,nval,ntst)

  ###########################################################################
  ###########################################################################

  # PERFORM:
  #          1) DATA-SCALING (OPTIONAL)
  #          2) INPUT VARIBALE SELECTION (IVS);
  #          3) DATA-DRIVEN MODEL (DDM) CALIBRATION;
  #          4) DDM PREDICTION (FOR CALIBRATION, VALIDATION, TEST SETS)

  # NOTE!
  #
  # IO.od contains a set of lists (at least one) that depends on:
  #
  # 1) whether a wavelet-based forecasting model is used; and
  # 2) the type of wavelet-based forecast model ('within','across',
  #    'within-hybrid', and 'across-hybrid') contain J+1 lists,
  #    each related to a different scale for the wavelet and scaling
  #    coefficients
  #
  # Therefore, data-scaling (optional), IVS, DDM calibration and
  # prediction is done in a for loop with bounds between 1 and J + 1 times.

  ###########################################################################
  ###########################################################################

  # for each item in the list perform:
  #
  #   1) DATA-SCALING (OPTIONAL);
  #   2) IVS;
  #   3) DDM CALIBRATION; AND
  #   4) DDM PREDICTION (FOR CALIBRATION, VALIDATION, TEST SETS)

  # preallocate space for the model predictions (this vector is aggregated)
  # through the loop below in the case of particular wavelet-based
  # forecasting models such as 'within' or 'across')

  pred = matrix(0,N,1)
  targ = matrix(0,N,1)
  IO.ivs = vector(mode="list", length = length(IO.od))
  IO.ddm = vector(mode="list", length = length(IO.od))

  for ( i in 1 : length (IO.od)){

    ###########################################################################
    # DATA-SCALING (OPTIONAL)
    ###########################################################################

    IO.s = scale_ab(IO.od[[i]][IO.inds$indc,],IO.od[[i]],0,1)

    # input scaling

    if( !scale_inputs ){

      # IO.s[,-1] = IO.od[[i]][,-1]
      IO.s[,-1] = as.matrix(IO.od[[i]][,-1])

    }

    # target scaling

    if( !scale_target ){

      IO.s[,1] = IO.od[[i]][,1]

    }


    ###########################################################################
    # PERFORM INPUT VARIABLE SELECTION (IVS)
    ###########################################################################

    IO.ivs[[i]] = hydroIVS::ivsIOData(y=IO.s[IO.inds$indc,1],
                                      x=IO.s[IO.inds$indc,-1],
                                      ivsm,ivs_param)

    ###########################################################################
    # CALIBRATE DATA-DRIVEN MODEL (DDM)
    ###########################################################################

    IO.ddm[[i]] = calibrateDDM(yc=IO.s[IO.inds$indc,1],
                               xc=IO.s[IO.inds$indc,
                                       IO.ivs[[i]]$names_sel_inputs],
                               ddm,ddm_param)

    ###########################################################################
    # MAKE PREDICTIONS USING CALIBRATED DATA-DRIVEN MODEL (DDM)
    ###########################################################################

    # get target (if certain types of wavelet-based forecasting models are
    # used (e.g., 'within' or 'across'), then target is aggregated across
    # the different scales)

    targ = targ + IO.od[[i]][,1]

    # make predictions using the calibrated model

    ptemp = predictDDM(x=IO.s[,IO.ivs[[i]]$names_sel_inputs],
                       xc=IO.s[IO.inds$indc,
                               IO.ivs[[i]]$names_sel_inputs],
                       yc=IO.s[IO.inds$indc,1],
                       ddm,ddm_param,IO.ddm[[i]])

    # if required, rescale predictions to target scale

    if(scale_target){

      ptemp = de_scale_ab(as.matrix(IO.od[[i]][IO.inds$indc,1]),
                          as.matrix(ptemp),0,1)

    }

    # final model predictions (if certain types of wavelet-based
    # forecasting models are used (e.g., 'within' or 'across' ), then
    # 'pred' will aggregate predictions across scales)

    pred = pred + ptemp

  }

  # if required, ensure only non-negative values are used

  if(cutoff0){ pred[pred < 0] = 0}

  # get wddff performance

  perf_c = hydroGOF::gof(obs=targ[IO.inds$indc], sim=pred[IO.inds$indc])
  perf_v = hydroGOF::gof(obs=targ[IO.inds$indv], sim=pred[IO.inds$indv])
  perf_t = hydroGOF::gof(obs=targ[IO.inds$indt], sim=pred[IO.inds$indt])

  ###########################################################################
  # STORE (AND, IF REQUIRED, SAVE) RESULTS
  ###########################################################################

  # results of wddff model

  rslt = list()
  rslt$predictions = pred
  rslt$target = targ
  rslt$n.inds = N
  rslt$indc = IO.inds$indc
  rslt$indv = IO.inds$indv
  rslt$indt = IO.inds$indt
  rslt$ivs = IO.ivs
  rslt$ddm = IO.ddm
  rslt$perf_c = perf_c
  rslt$perf_v = perf_v
  rslt$perf_t = perf_t

  # if required, save results to .rda file

  if( savefile ){

    # filename = paste0('wddff_',trunc(
    # as.numeric(Sys.time())*1000, digits=13),'.rda')
    st = format(Sys.time(), "%Y-%m-%d-%H-%M-%OS3")
    st = gsub("[.]", "-", st)
    filename = paste0('wddff_',st,'.rda')
    mdl$filename = filename
    mdl_rslt = list(mdl = mdl, rslt = rslt)
    save(mdl_rslt, file=filename) # save file as '.rda'

  }

  # return results to user

  return( list(

    mdl = mdl,
    rslt = rslt

  ))

}

# ------------------------------------------------------------------------------

#' @title Update WDDFF With New Data
#' @references
#' This function takes in new time series measurements for the target and
#' input time series (standard and/or auxilliary), performs time-lagging,
#' (wavelet-decomposition,) and generates predictions given an exisiting
#' Wavelet Data-Driven Forecasting Framework (WDDFF) model.
#' @param mdl_rslt list that contains all essential inputs from 'wddff'
#' @param y new target time series \[N x 1\]
#' @param x new input time series \[N x D\]
#' @param z new auxilliary input time series \[N x C\]
#' @param light store all function inputs or neglect 'y', 'x', and 'z' \[binary\]
#' @param savefile save file using function name and current Sys.time() \[binary\]
#' @return
#' list that contains:
#' - mdl
#' - predictions
#' - target
#' - residual
#' - n.inds
#' - indc
#' - indv
#' - indt
#' - ivs
#' - ddm
#' - perf_c
#' - perf_v
#' - perf_t
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
#' @rdname upd.wddff
#' @export
upd.wddff <- function(mdl_rslt,
                      y, x=NULL, z=NULL,
                      light=TRUE, savefile=FALSE){

  ###########################################################################
  # HIGH LEVEL SETTINGS
  ###########################################################################

  # set.seed(123) # make reproducible

  ###########################################################################
  # SET NEW MODEL STRUCTURE AND UPDATE VARIABLES
  ###########################################################################

  # set model structure

  mdl = mdl_rslt$mdl # transfer old model settings to new mdl
  rslt = mdl_rslt$rslt # transfer old model results to new mdl
  rm(mdl_rslt) # remove old model settings and results

  # update target and inputs (if any)

  y = as.matrix(y)
  y = rbind(mdl$y, y)
  mdl$y = y

  if( !is.null(x) ){

    x = as.matrix(x)
    x = rbind(mdl$x, x)
    mdl$x = x

  }

  if( !is.null(z) ) {

    z = as.matrix(z)
    z = rbind(mdl$z, z)
    mdl$z = z

  }

  if( light ){

    mdl$y = NULL
    mdl$x = NULL
    mdl$z = NULL

  }


  ###########################################################################
  # CREATE INPUT-OUTPUT (IO) DATASET
  ###########################################################################

  # create IO dataset (lags input-output dataset and peforms wavelet-
  # decomposition, if needed)

  IO.pd = prepareData(y=y,x=x,z=z,mdl$leadtime,mdl$lag_Y,mdl$lag_X,
                      mdl$wfm,mdl$wt,mdl$wavelet,mdl$decomp_level)

  # adjust IO dataset for boundary-condition (very important for wavelet-
  # based models)

  IO.od = organizeIOData(IO.pd,mdl$max_decomp_level,mdl$max_wavelet_length)

  ###########################################################################
  # IDENTIFY INDICES OF NEW DATA AND CREATE NEW TEST PARTITION
  ###########################################################################

  # create a new 'forecast' parition based on new data and update 'test'
  # partition

  old.N = rslt$n.inds
  N = nrow(IO.od[[1]])
  indf = seq(old.N + 1, N, by = 1)
  rslt$indv = c(rslt$indv,rslt$indt)
  rslt$indt = indf
  mdl$nval = length(rslt$indv)
  mdl$ntst = length(rslt$indt)
  rslt$n.inds = N
  rslt$n.inds.old = old.N

  ###########################################################################
  ###########################################################################

  # PERFORM:
  #          1) DATA-SCALING (OPTIONAL)
  #          2) DATA-DRIVEN MODEL (DDM) PREDICTION (FOR CALIBRATION,
  #             VALIDATION, TEST, FORECAST SETS)

  # NOTE!
  #
  # IO.od contains a set of lists (at least one) that depends on:
  #
  # 1) whether a wavelet-based forecasting model is used; and
  # 2) the type of wavelet-based forecast model ('within','across',
  #    'within-hybrid', and 'across-hybrid') contain J+1 lists,
  #    each related to a different scale for the wavelet and scaling
  #    coefficients
  #
  # Therefore, data-scaling (optional) and DDM prediction is done in a
  # for loop with bounds between 1 and J + 1 times.

  ###########################################################################
  ###########################################################################

  # for each item in the list perform:
  #
  #   1) DATA-SCALING (OPTIONAL);
  #   2) DDM PREDICTION (FOR CALIBRATION, VALIDATION, TEST, FORECAST SETS)

  # preallocate space for the model predictions (this vector is aggregated)
  # through the loop below in the case of particular wavelet-based
  # forecasting models such as 'within' or 'across')

  pred = matrix(0,length(indf),1)
  targ = matrix(0,length(indf),1)

  for ( i in 1 : length (IO.od)){

    ###########################################################################
    # DATA-SCALING (OPTIONAL)
    ###########################################################################

    IO.s = scale_ab(IO.od[[i]][rslt$indc,],IO.od[[i]],0,1)

    # input scaling
    if( !mdl$scale_inputs ){

      IO.s[,-1] = IO.od[[i]][,-1]

    }

    # target scaling
    if( !mdl$scale_target ){

      IO.s[,1] = IO.od[[i]][,1]

    }

    ###########################################################################
    # MAKE PREDICTIONS USING PRE-CALIBRATED DATA-DRIVEN MODEL (DDM)
    ###########################################################################

    # get target (if certain types of wavelet-based forecasting models are
    # used (e.g., 'within' or 'across'), then target is aggregated across
    # the different scales)

    targ = targ + IO.od[[i]][indf,1]

    # make predictions using the calibrated model

    ptemp = predictDDM(x=IO.s[indf,rslt$ivs[[i]]$names_sel_inputs],
                       xc=IO.s[rslt$indc,rslt$ivs[[i]]$names_sel_inputs],
                       yc=IO.s[rslt$indc,1],
                       mdl$ddm,mdl$ddm_param,rslt$ddm[[i]])

    # if required, rescale predictions to target scale

    if(mdl$scale_target){

      ptemp = de_scale_ab(as.matrix(IO.od[[i]][rslt$indc,1]),
                          as.matrix(ptemp),0,1)

    }

    # final model predictions (if certain types of wavelet-based
    # forecasting models are used (e.g., 'within' or 'across' ), then
    # 'pred' will aggregate predictions across scales)

    pred = pred + ptemp

  }

  # if required, ensure only non-negative values are used

  if(mdl$cutoff0){ pred[pred < 0] = 0}

  # get wddff performance for new test set

  target = rbind(rslt$target,targ)
  predictions = rbind(rslt$predictions,pred)
  perf_v = hydroGOF::gof(obs=target[rslt$indv], sim=predictions[rslt$indv])
  perf_t = hydroGOF::gof(obs = target[rslt$indt], sim = predictions[rslt$indt])

  ###########################################################################
  # STORE (AND, IF REQUIRED, SAVE) RESULTS
  ###########################################################################

  # results of updated wddff model

  rslt$predictions = predictions
  rslt$target = target
  rslt$perf_v = perf_v
  rslt$perf_t = perf_t

  # if required, save results to .rda file

  if( savefile ){

    # append update number to original WDDFF model

    filename = paste0('upd_',length(mdl$filename),'_',mdl$filename[1])
    mdl$filename = c(mdl$filename,filename)
    mdl_rslt = list(mdl = mdl, rslt = rslt)
    save(mdl_rslt, file=filename) # save file as '.rda'

  }

  # return results to user

  return( list(

    mdl = mdl,
    rslt = rslt

  ))

}

# ------------------------------------------------------------------------------

#' @title Model Agnostic WDDFF
#' @description
#' This function creates input-output (IO) datasets for a given wavelet-based
#' forecasting method, the number of time series lags to consider as predictors,
#' and the number of boundary-effected coefficients to remove.
#' @param ts time series/matrix/vector \[N x 1\]
#' @param W_ts wavelet decomposed time series \[N x J+1\]
#' @param wfm wavelet-based forecasting method \[string; 'single', 'within', 'across','single_hybrid','within_hybrid','across_hybrid'\]
#' @param lags lag length to include as predictors \[integer >= 0\]
#' @param maxlag maximum lag length \[1 < integer <= lags\]
#' @param nbc number of boundary-effected coefficients \[integer >= 0\]
#' @return list of input-output dataset(s) according to wfm \[N x D x J+1\]
#' @references
#' Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
#' based hydrological and water resources forecasting models for real-world
#' applications with best practices and a new forecasting framework, J. Hydrol.,
#' doi:10.1016/j.jhydrol.2018.05.003.
#' @rdname wddff.wfm
#' @export
wddff.wfm <- function(ts,W_ts,wfm,lags,maxlag,nbc){

  J = ncol(W_ts)

  if(nbc > 0){ # boundary-correction req'd

    W_ts_bc = as.matrix(W_ts[-(1:nbc),]) # boundary-corrected coefficients
    ts_bc = ts[-(1:nbc),] # boundary-corrected target

  }

  else{ # not req'd

    W_ts_bc = W_ts
    ts_bc = ts

  }

  switch(wfm,

         # none

         none={

           IO = list()
           IO[[1]] = lagmatrix(W_ts_bc,lags)[-(1:maxlag),] # boundary-corrected & lag-corrected target and inputs

         },

         # single - applicable with AT and MODWT

         single={

           W_ts_bc_lc_inputs = list()
           IO = list()

           for(j in 1:J){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients at each scale
             W_ts_bc_lc_inputs[[j]] = W_ts_bc_lc[,-1] # boundary-corrected & lag-corrected coefficients for each scale (inputs only!)

           }
           W_ts_bc_lc_inputs = do.call(cbind,W_ts_bc_lc_inputs) # unroll list to matrix

           ts_bc_lc = as.matrix(ts_bc)[-(1:maxlag),] # boundary-corrected, lag-corrected target
           IO[[1]] = cbind(ts_bc_lc,W_ts_bc_lc_inputs) # input-output dataset (first column in target, remaining columns are predictors)

         },

         # within - applicable with AT only

         within={

           IO = list()

           for(j in 1:J){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients for each scale
             IO[[j]] = W_ts_bc_lc # input-output dataset (first column in target, remaining columns are predictors)

           }

         },

         # across - applicable with AT onnly

         across={

           W_ts_bc_lc_inputs = list()
           W_ts_bc_lc_output = list()
           IO = list()

           for(j in 1:J){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients for each scale
             W_ts_bc_lc_inputs[[j]] = W_ts_bc_lc[,-1] # boundary-corrected & lag-corrected coefficients for each scale (inputs only!)
             W_ts_bc_lc_output[[j]] = W_ts_bc_lc[,1] # boundary-corrected & lag-corrected coefficients for each scale (output only!)

           }

           for(j in 1:J){

             W_ts_bc_lc_inputs_tmp = do.call(cbind,W_ts_bc_lc_inputs) # unroll list to matrix
             IO[[j]] = cbind(W_ts_bc_lc_output[[j]],W_ts_bc_lc_inputs_tmp) # input-output dataset (first column in target, remaining columns are predictors)


           }

         },

         # single_hybrid - applicable with AT and MODWT

         single_hybrid={

           W_ts_bc_lc_inputs = list()
           IO = list()

           ts_lag_bc_lc = lagmatrix(ts_bc,lags)[-(1:maxlag),-1] # boundary-corrected & lag-corrected lagged target

           for(j in 1:(J)){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients at each scale
             W_ts_bc_lc_inputs[[j]] = W_ts_bc_lc[,-1] # boundary-corrected & lag-corrected coefficients for each scale (inputs only!)

           }
           W_ts_bc_lc_inputs = do.call(cbind,W_ts_bc_lc_inputs) # unroll list to matrix

           ts_bc_lc = as.matrix(ts_bc)[-(1:maxlag),] # boundary-corrected, lag-corrected target
           IO[[1]] = cbind(ts_bc_lc,ts_lag_bc_lc,W_ts_bc_lc_inputs) # input-output dataset (first column in target, remaining columns are predictors)


         },

         # within_hybrid - applicable with AT only

         within_hybrid={

           IO = list()

           ts_lag_bc_lc = lagmatrix(ts_bc,lags)[-(1:maxlag),-1] # boundary-corrected & lag-corrected lagged target

           for(j in 1:J){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients for each scale
             IO[[j]] = cbind(W_ts_bc_lc[,1],ts_lag_bc_lc,W_ts_bc_lc[,-1]) # input-output dataset (first column in target, remaining columns are predictors)

           }

         },

         # across_hybrid - applicable with AT onnly

         across_hybrid={

           W_ts_bc_lc_inputs = list()
           W_ts_bc_lc_output = list()
           IO = list()

           ts_lag_bc_lc = lagmatrix(ts_bc,lags)[-(1:maxlag),-1] # boundary-corrected & lag-corrected lagged target

           for(j in 1:J){

             W_ts_bc_lc = lagmatrix(W_ts_bc[,j],lags)[-(1:maxlag),] # boundary-corrected & lag-corrected coefficients for each scale
             W_ts_bc_lc_inputs[[j]] = W_ts_bc_lc[,-1] # boundary-corrected & lag-corrected coefficients for each scale (inputs only!)
             W_ts_bc_lc_output[[j]] = W_ts_bc_lc[,1] # boundary-corrected & lag-corrected coefficients for each scale (output only!)

           }

           for(j in 1:J){

             W_ts_bc_lc_inputs_tmp = do.call(cbind,W_ts_bc_lc_inputs) # unroll list to matrix
             IO[[j]] = cbind(W_ts_bc_lc_output[[j]],ts_lag_bc_lc,W_ts_bc_lc_inputs_tmp) # input-output dataset (first column in target, remaining columns are predictors)


           }

         }

  )

  return(IO)

}

# ------------------------------------------------------------------------------

# This function performs wavelet decomposition using the à trous algorithm (AT) or
# the maximal overlap discrete wavelet transform (MODWT).
#
# # Inputs:
# @param ts - time series [N x 1]
# @param wt - type of wavelet transform [string; at','modwt']
# @param wavelet - scaling filter name [string]
# @param decomp_level - decomposition level [1 < integer < N/2]
#
# Output:
# @param W - wavelet and scaling coefficients [N x D x J+1] (wavelet coefs in first J columns)
#
# Date created: May 15, 2018
# Date updated: May 15, 2018
#
# Reference(s):
#
# Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
# based hydrological and water resources forecasting models for real-world
# applications with best practices and a new forecasting framework, J. Hydrol.,
# doi:10.1016/j.jhydrol.2018.05.003.
#

# wddff.wt <- function(ts=matrix(runif(100),100,1),
#                      wt='at',wavelet='haar',decomp_level=1){
#
#   switch(wt,
#
#          # AT
#
#          at={
#
#            W = atrous_dwt(ts,wavelet,decomp_level) # use AT algoritm
#
#          },
#
#          # MODWT
#
#          modwt={
#
#            W = mo_dwt(ts,wavelet,decomp_level) # use MODWT algoritm
#
#          }
#
#   )
#
#   return(W)
#
# }

