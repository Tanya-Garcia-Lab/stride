
####################################
####################################
##
##
## Functions for the main method
##
####################################
##
####################################

stride.estimator.wrapper <- function(
  n,m,p,qvs,q,
  x,delta,ww,zz,
  run.NPMLEs,
  run.NPNA,
  run.NPNA_avg,
  run.NPNA_wrong,
  run.OLS,
  run.WLS,
  run.EFF,
  tval,tval0,
  z.use,w.use,
  update.qs=FALSE,
  know.true.groups=FALSE,
  true.group.identifier=NULL,
  run.prediction.accuracy=FALSE,
  do_cross_validation_AUC_BS=FALSE){

  ###############################
  ## check for errors in input ##
  ###############################
  common_error_messages(n,m,p,qvs,q,
                        x,delta,ww,zz,
                        run.NPMLEs,
                        run.NPNA,
                        run.NPNA_avg,
                        run.NPNA_wrong,
                        run.OLS,
                        run.WLS,
                        run.EFF,
                        tval,tval0,
                        z.use,w.use,
                        update.qs,
                        know.true.groups,
                        true.group.identifier,
                        run.prediction.accuracy)




  ##############################
  ## setup output for storage ##
  ##############################
  est.names <- get_method_label(run.NPMLEs,
                                run.NPNA,
                                run.NPNA_avg,
                                run.NPNA_wrong,
                                run.OLS,
                                run.WLS,
                                run.EFF
  )

  method.label <- est.names


  ###############
  ## Compute r ##
  ###############
  r <- compute.r(n,m,p,qvs,q)



  ###################
  ## Make data set ##
  ###################
  data <- make.data.set(
    n,m,p,qvs,q,
    x,delta,ww,zz,
    know.true.groups,true.group.identifier)



  #######################
  ## main computations ##
  #######################
  problem <- NULL


  ####################
  ## run estimators ##
  ####################
  data.out <- estimator.main(data,
                             n,p,m,r,qvs,
                             tval,tval0,
                             method.label,
                             z.use,w.use,
                             update.qs,
                             run.prediction.accuracy,
                             do_cross_validation_AUC_BS
  )

  return(data.out)
}

#' Common error messages
#'
#' Function to check for common errors in the input.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#' @param x a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.
#' @param delta a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).
#' @param ww a numeric vector of length \code{n} containing the values of the continuous
#' covariate for each person in the sample.
#' @param zz a numeric vector of length \code{n} containing the values of the discrete
#' covariate for each person in the sample.
#' @param run.NPMLEs a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data based on the type-I and type II nonparametric maximum likelihood
#' estimators. The type I nonparametric maximum likelihood estimator is referred
#' to as the "Kaplan-Meier" estimator in Garcia and Parast (2020). Neither the type I nor type II
#' include covariates or dynamic landmarking.
#' @param run.NPNA a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that accounts for covariates and dynamic
#' landmarking. This estimator is called "NPNA" in Garcia and Parast (2020).
#' @param run.NPNA_avg a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that averages out over the observed covariates.
#' This is referred to as NPNA_marg in Garcia and Parast (2020).
#' @param run.NPNA_wrong a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that ignores landmarking. This is referred to as
#' NPNA_{t_0=0} in the paper.
#' @param run.OLS a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using an ordinary least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.WLS  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using a weighted least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.EFF  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using the efficient influence function based on Hilbert space projection theory results.
#' The estimator adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param tval numeric vector of time points at which the distribution function is evaluated, all values must
#' be non-negative.
#' @param tval0 numeric vector of time points representing the landmark times. All values must be non-negative
#' and smaller than the maximum of \code{tval}.
#' @param z.use numeric vector at which to evaluate the discrete covariate \eqn{Z} at in the estimated distribution function.
#' The values of \code{z.use} must be in the range of the observed \code{zz}.
#' @param w.use numeric vector at which to evaluate the continuous covariate \eqn{W} at in the estimated distribution function.
#' The values of \code{w.use} must be in the range of the observed \code{ww}.
#' @param update.qs logical indicator. If TRUE, the mixture proportions \code{q} will be updated. This is currently not implemented.
#' @param know.true.groups logical indicator. If TRUE, then we know the population identifier for each person in the sample.
#' This option is only used for simulation studies. Default is FALSE.
#' @param true.group.identifier numeric vector of length \code{n} denoting the population identifier for each person in the sample.
#' Default is NULL.
#' @param run.prediction.accuracy logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS). Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#'
#' @references
#' Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#' Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
#' In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.
#'
#' Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
#' Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
#' Annals of Applied Statistics, 8(2), 1182-1208.
#'
#' Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
#' application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association,
#' 107, 1324-1338.
#'
#' @return Error or warning messages when input is not appropriate for the methods.
#' @export
common_error_messages <- function(n,m,p,qvs,q,
                                  x,delta,ww,zz,
                                  run.NPMLEs,
                                  run.NPNA,
                                  run.NPNA_avg,
                                  run.NPNA_wrong,
                                  run.OLS,
                                  run.WLS,
                                  run.EFF,
                                  tval,tval0,
                                  z.use,w.use,
                                  update.qs,
                                  know.true.groups,
                                  true.group.identifier,
                                  run.prediction.accuracy)
{
  ## check if n>=1
  error_positive_value(n,get_variable_name(n))

  ## check m>=1
  error_positive_value(m,get_variable_name(m))

  ## check p>=1
  error_positive_value(p,get_variable_name(p))

  ## check if qvs of size p x m
  stop_length_message(qvs,c(p,m))
  error_nonpositive_value(qvs,get_variable_name(qvs))

  ## check if q is of size p x n
  stop_length_message(q,c(p,n))
  error_nonpositive_value(q,get_variable_name(q))

  ## check if x is of length n
  stop_length_message(x,n)

  ## check if delta is of length n
  stop_length_message(delta,n)

  ## check if ww is of length n
  if(!is.null(ww)){
    stop_length_message(ww,n)
  }

  ## check if zz is of length n
  if(!is.null(zz)){
    stop_length_message(zz,n)
  }

  ## check if tval >=0
  error_nonpositive_value(tval,get_variable_name(tval))

  ## check if tval0 >=0
  error_nonpositive_value(tval0,get_variable_name(tval0))

  ## check if max(tval) > tval0
  if(any(tval0>max(tval))){
    stop("All values of tval0 must be smaller than max(tval).")
  }

  ## check if z.use in the range of observed zz
  if(!is.null(z.use)){
    if(is.null(zz)){
      stop("To evaluate the method(s) at z.use, you need to provide z-covariate values for each subject.")
    }

    if(any(!z.use%in% range(zz) ) ){
      warning("Not all values in z.use are in the range of observed Z.")
    }
  }

  ## check if w.use in the range of observed ww
  if(!is.null(w.use)){

    if(is.null(ww)){
      stop("To evaluate the method(s) at w.use, you need to provide w-covariate values for each subject.")
    }

    if(any(w.use < min(ww) | w.use > max(ww) ) ){
      warning("Not all values in w.use are in the range of observed W.")
    }

  }

  ## update.qs must be FALSE
  if(update.qs==TRUE){
    stop("Method not implemented yet.")
  }

  ## data needed when know.true.groups=TRUE
  if(know.true.groups==TRUE & is.null(true.group.identifier)){
    stop("User must provide true.group.identifier.")
  }

  if(run.prediction.accuracy==TRUE){
    if(know.true.groups==FALSE | is.null(true.group.identifier)){
      stop("To run prediction accuracy, user must set
           know.true.groups=TRUE and provide true.group.identifier.")
    }

  }

  ## add warning message about NPNA methods and need for covariates
  if(run.NPNA==TRUE | run.NPNA_avg==TRUE | run.NPNA_wrong==TRUE){
    if(is.null(zz) | is.null(ww)){
      stop("For estimators of NPNA, NPNA_avg, NPNA_wrong, the user must provide covariate values z and w
           for each subject.")
    }

    if(is.null(z.use) | is.null(w.use)){
      stop("For estimators of NPNA, NPNA_avg, NPNA_wrong, the user must provide covariate z.use and w.use.")
    }
  }

  if(!any(c(run.NPMLEs,run.NPNA,run.NPNA_avg,run.NPNA_wrong,
            run.OLS,run.WLS,run.EFF))){
    stop("No methods are set to run.")
  }
}

#' Dynamic landmark prediction estimator for mixture data with covariates
#'
#' Estimates the distribution function for mixture data where
#' the population identifiers are unknown, but the probability of belonging
#' to a population is known. The distribution functions are evaluated at
#' time points \code{tval} and adjust for dynamic landmark prediction and one
#' discrete covariate (\code{zz}) and one continuous covariate (\code{ww}).
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#' @param x a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.
#' @param delta a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).
#' @param ww a numeric vector of length \code{n} containing the values of the continuous
#' covariate for each person in the sample.
#' @param zz a numeric vector of length \code{n} containing the values of the discrete
#' covariate for each person in the sample.
#' @param run.NPMLEs a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data based on the type-I and type II nonparametric maximum likelihood
#' estimators. The type I nonparametric maximum likelihood estimator is referred
#' to as the "Kaplan-Meier" estimator in Garcia and Parast (2020). Neither the type I nor type II
#' include covariates or dynamic landmarking.
#' @param run.NPNA a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that accounts for covariates and dynamic
#' landmarking. This estimator is called "NPNA" in Garcia and Parast (2020).
#' @param run.NPNA_avg a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that averages out over the observed covariates.
#' This is referred to as NPNA_marg in Garcia and Parast (2020)..
#' @param run.NPNA_wrong a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that ignores landmarking. This is referred to as
#' NPNA_{t_0=0} in the paper.
#' @param run.OLS a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using an ordinary least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.WLS  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using a weighted least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.EFF  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using the efficient influence function based on Hilbert space projection theory results.
#' The estimator adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param tval numeric vector of time points at which the distribution function is evaluated, all values must
#' be non-negative.
#' @param tval0 numeric vector of time points representing the landmark times. All values must be non-negative
#' and smaller than the maximum of \code{tval}.
#' @param z.use numeric vector at which to evaluate the discrete covariate \eqn{Z} at in the estimated distribution function.
#' The values of \code{z.use} must be in the range of the observed \code{zz}.
#' @param w.use numeric vector at which to evaluate the continuous covariate \eqn{W} at in the estimated distribution function.
#' The values of \code{w.use} must be in the range of the observed \code{ww}.
#' @param update.qs logical indicator. If TRUE, the mixture proportions \code{q} will be updated. This is currently not implemented.
#' @param know.true.groups logical indicator. If TRUE, then we know the population identifier for each person in the sample.
#' This option is only used for simulation studies to check prediction accuracy. Default is FALSE.
#' @param true.group.identifier numeric vector of length \code{n} denoting the population identifier for each person in the sample.
#' Default is NULL.
#' @param run.prediction.accuracy logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS). Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#' @param do_cross_validation_AUC_BS logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS) using cross-validation. Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#'
#' @section Details:
#' We estimate the distribution function for mixture data  where
#' the population identifiers are unknown, but the probability of belonging
#' to a population is known. The distribution functions are evaluated at
#' time points \code{tval} and adjust for dynamic landmark prediction and one
#' discrete covariate (\code{zz}) and one continuous covariate (\code{ww}).
#' Dynamic landmark prediction means that the distribution function is computed knowing
#' that the survival time, \eqn{T}, satisfies \eqn{T >t_0}
#' where \eqn{t_0} are the time points in \code{tval0}.
#'
#' @references
#' Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#' Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
#' In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.
#'
#' Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
#' Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
#' Annals of Applied Statistics, 8(2), 1182-1208.
#'
#' Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
#' application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association,
#' 107, 1324-1338.
#'
#'
#' @return \code{stride.estimator} returns a list containing
#' \itemize{
#'    \item{problem: }{a numeric indicator of errors in the NPNA estimator. If NULL, no error is reported.
#'    Otherwise, there is an error in the computation of the NPNA estimator.}
#'
#'    \item{Ft.estimate: }{a numeric array containing the estimated distribution functions for all methods for all
#'    \code{p} populations. The distribution function is evaluated at each \code{tval},
#'    \code{tval0}, \code{z.use}, \code{w.use}, and for all \code{p} populations.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{p}.  The distribution function is only valid for \eqn{t\geq t_0}, so
#'    \code{Ft.estimate} shows NA for any combination for which \eqn{t<t_0}.
#'    }
#'
#'    \item {St.estimate: }{a numeric array containing the estimated distribution functions for all methods
#'    for all \code{m} mixture proportion subgroups. The distribution function is evaluated
#'    at each \code{tval}, \code{tval0}, \code{z.use}, \code{w.use}, and for all \code{m} mixture
#'    proportion subgroups.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{m}.  The distribution function is only valid for \eqn{t\geq t_0}, so
#'    \code{St.estimate} shows NA for any combination for which \eqn{t<t_0}.
#'    }
#'
#'    \item{Ft.AUC.BS: }{a numeric array containing the
#'    area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{p} populations. The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#'
#'    Results for both the estimated distributon functions and prediction
#'    accuracy measures (AUC, BS) are only valid when \eqn{t\geq t_0}, so
#'    arrays show NA for any combination for which \eqn{t<t_0}.
#'
#'    }
#'
#'    \item {St.AUC.BS: }{a numeric array containing the results are the
#'    area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{m} mixture proportion groups.
#'    The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#'
#'    Results for both the estimated distributon functions and prediction
#'    accuracy measures (AUC, BS) are only valid when \eqn{t\geq t_0}, so
#'    arrays show NA for any combination for which \eqn{t<t_0}.
#'
#'    }
#'
#' }
#'
#' @example man/examples/stride_ex.R
#'
#' @export
stride.estimator <- function(n,m,p,qvs,q,
                              x,delta,ww,zz,
                              run.NPMLEs,
                              run.NPNA,
                              run.NPNA_avg,
                              run.NPNA_wrong,
                              run.OLS,
                              run.WLS,
                              run.EFF,
                              tval,tval0,
                              z.use,w.use,
                              update.qs,
                              know.true.groups=FALSE,
                              true.group.identifier=NULL,
                              run.prediction.accuracy,
                              do_cross_validation_AUC_BS){

  ############################
  ## run the main procedure ##
  ############################
  estimators.out <- stride.estimator.wrapper(
    n,m,p,qvs,q,
    x,delta,ww,zz,
    run.NPMLEs,
    run.NPNA,
    run.NPNA_avg,
    run.NPNA_wrong,
    run.OLS,
    run.WLS,
    run.EFF,
    tval,tval0,
    z.use,w.use,
    update.qs,
    know.true.groups,
    true.group.identifier,
    run.prediction.accuracy=FALSE,
    do_cross_validation_AUC_BS=FALSE
  )

  ## check if we have a problem with NPNA estimator.
  Ft.estimate <- estimators.out$Ft.store
  St.estimate <- estimators.out$Sout.store
  Ft.AUC.BS <- NULL 	## place holder for now
  St.AUC.BS <- NULL ## place holder for now

  problem <- estimators.out$problem
  if(is.null(problem)){
    #############################
    ## get prediction accuracy ##
    #############################
    ## only run prediction computation on simulated data
    if(run.prediction.accuracy==TRUE & know.true.groups==TRUE){
      prediction.out <- stride.estimator.wrapper(
        n,m,p,qvs,q,
        x,delta,ww,zz,
        run.NPMLEs,
        run.NPNA,
        run.NPNA_avg,
        run.NPNA_wrong,
        run.OLS,
        run.WLS,
        run.EFF,
        tval,tval0,
        z.use,w.use,
        update.qs,
        know.true.groups,
        true.group.identifier,
        run.prediction.accuracy,
        do_cross_validation_AUC_BS)

      Ft.AUC.BS <- prediction.out$Ft.store
      St.AUC.BS <- prediction.out$Sout.store
    }
  }


  return(list(problem=problem,
              Ft.estimate=Ft.estimate,
              St.estimate=St.estimate,
              Ft.AUC.BS=Ft.AUC.BS,
              St.AUC.BS=St.AUC.BS))
}



######################################
## Functions for variance estimates ##
######################################
#' Bootstrap variance estimator
#'
#' Computes the bootstrap variance for all methods.
#'
#' @param nboot number of bootstrap replicates, must be at least 1.
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#' @param x a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.
#' @param delta a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).
#' @param ww a numeric vector of length \code{n} containing the values of the continuous
#' covariate for each person in the sample.
#' @param zz a numeric vector of length \code{n} containing the values of the discrete
#' covariate for each person in the sample.
#' @param run.NPMLEs a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data based on the type-I and type II nonparametric maximum likelihood
#' estimators. The type I nonparametric maximum likelihood estimator is referred
#' to as the "Kaplan-Meier" estimator in Garcia and Parast (2020). Neither the type I nor type II
#' include covariates or dynamic landmarking.
#' @param run.NPNA a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that accounts for covariates and dynamic
#' landmarking. This estimator is called "NPNA" in Garcia and Parast (2020).
#' @param run.NPNA_avg a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that averages out over the observed covariates.
#' This is referred to as NPNA_marg in Garcia and Parast (2020).
#' @param run.NPNA_wrong a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that ignores landmarking. This is referred to as
#' NPNA_{t_0=0} in the paper.
#' @param run.OLS a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using an ordinary least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.WLS  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using a weighted least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.EFF  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using the efficient influence function based on Hilbert space projection theory results.
#' The estimator adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param tval numeric vector of time points at which the distribution function is evaluated, all values must
#' be non-negative.
#' @param tval0 numeric vector of time points representing the landmark times. All values must be non-negative
#' and smaller than the maximum of \code{tval}.
#' @param z.use numeric vector at which to evaluate the discrete covariate \eqn{Z} at in the estimated distribution function.
#' The values of \code{z.use} must be in the range of the observed \code{zz}.
#' @param w.use numeric vector at which to evaluate the continuous covariate \eqn{W} at in the estimated distribution function.
#' The values of \code{w.use} must be in the range of the observed \code{ww}.
#' @param update.qs logical indicator. If TRUE, the mixture proportions \code{q} will be updated. This is currently not implemented.
#' @param know.true.groups logical indicator. If TRUE, then we know the population identifier for each person in the sample.
#' This option is only used for simulation studies. Default is FALSE.
#' @param estimator_Ft a numeric array of the estimated distribution functions for all \code{p} populations.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{p}.
#' 	Results are only valid when \eqn{t\geq t_0}, so arrays show NA for any combination for which \eqn{t<t_0}.
#' @param estimator_St a numeric array of the estimated distribution functions for all \code{m} mixture proportion groups.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{m}.
#' 	Results are only valid when \eqn{t\geq t_0}, so arrays show NA for any combination for which \eqn{t<t_0}.
#' @param AUC_BS_Ft a numeric array of the area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{p} populations. The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#' 	Results are only valid when \eqn{t\geq t_0}, so arrays show NA for any combination for which \eqn{t<t_0}.
#' @param AUC_BS_St a numeric array of the area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{m} mixture proportion groups.
#'    The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#'    Results are only valid when \eqn{t\geq t_0}, so
#'    arrays show NA for any combination for which \eqn{t<t_0}.
#'
#' @param true.group.identifier numeric vector of length \code{n} denoting the population identifier for each person in the sample.
#' Default is NULL.
#' @param run.prediction.accuracy logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS). Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#' @param do_cross_validation_AUC_BS logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS) using cross-validation. Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#'
#' @references
#' Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#' Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
#' In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.
#'
#' Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
#' Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
#' Annals of Applied Statistics, 8(2), 1182-1208.
#'
#' Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
#' application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association,
#' 107, 1324-1338.
#'
#' @return \code{stride.bootstrap.variance} returns a list containing
#' \itemize{
#'    \item{Ft.estimate.boot: }{a numeric array containing the estimated bootstrap variances,
#'    the 2.5% bootstrap quantile, and the 97.5% bootstrap quantile for all methods
#'    of the distribution function at all \code{p} populations. The results are shown at each \code{tval},
#'    \code{tval0}, \code{z.use}, \code{w.use}, and for all \code{p} populations.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{p} by 3. The last dimension
#'    corresponds to the bootstrap variance estimate, 2.5% bootstrap quantile,
#'    and the 97.5% bootstrap quantile.  The results are only valid for \eqn{t\geq t_0}, so
#'    \code{Ft.estimate} shows NA for any combination for which \eqn{t<t_0}.
#'    }
#'
#'    \item {St.estimate.boot: }{a numeric array containing the estimated bootstrap variances,
#'    the 2.5% bootstrap quantile, and the 97.5% bootstrap quantile for all methods
#'    of the distribution function at all \code{m} \code{m} mixture proportion subgroups. The
#'    results are shown at each \code{tval}, \code{tval0}, \code{z.use}, \code{w.use}, and for all \code{m} mixture
#'    proportion subgroups.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{m} by 3.
#'    The last dimension
#'    corresponds to the bootstrap variance estimate, 2.5% bootstrap quantile,
#'    and the 97.5% bootstrap quantile.  The results are only valid for \eqn{t\geq t_0}, so
#'    \code{Sout.data.out} shows NA for any combination for which \eqn{t<t_0}.
#'    }
#'
#'
#'    \item{Ft.AUC.BS.boot: }{a numeric array containing the estimated bootstrap variances,
#'    the 2.5% bootstrap quantile, and the 97.5% bootstrap quantile for all methods
#'    of the area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{p} populations. The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2 by 3. The second to last stores the AUC and BS results.
#'	  The last dimension corresponds to the bootstrap variance estimate, 2.5% bootstrap quantile,
#'    and the 97.5% bootstrap quantile.
#'    }
#'
#'    \item {St.AUC.BS.boot: }{a numeric array containing the estimated bootstrap variances,
#'    the 2.5% bootstrap quantile, and the 97.5% bootstrap quantile for all methods
#'	  of the area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{m} mixture proportion subgroups. The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2 by 3. The second to last stores the AUC and BS results.
#'	  The last dimension corresponds to the bootstrap variance estimate, 2.5% bootstrap quantile,
#'    and the 97.5% bootstrap quantile.
#'    }
#'
#' }
#'
#' @example man/examples/stride_ex.R
#'
#' @export
#' @import stats
stride.bootstrap.variance <- function(nboot,n,m,p,qvs,q,
                                       x,delta,ww,zz,
                                       run.NPMLEs,
                                       run.NPNA,
                                       run.NPNA_avg,
                                       run.NPNA_wrong,
                                       run.OLS,
                                       run.WLS,
                                       run.EFF,
                                       tval,tval0,
                                       z.use,w.use,
                                       update.qs,
                                       know.true.groups,
                                       true.group.identifier,
                                       estimator_Ft,
                                       estimator_St,
                                       AUC_BS_Ft,
                                       AUC_BS_St,
                                       run.prediction.accuracy,
                                       do_cross_validation_AUC_BS){


  ###############################
  ## check for errors in input ##
  ###############################
  common_error_messages(n,m,p,qvs,q,
                        x,delta,ww,zz,
                        run.NPMLEs,
                        run.NPNA,
                        run.NPNA_avg,
                        run.NPNA_wrong,
                        run.OLS,
                        run.WLS,
                        run.EFF,
                        tval,tval0,
                        z.use,w.use,
                        update.qs,
                        know.true.groups,
                        true.group.identifier,
                        run.prediction.accuracy)

  ##############################
  ## setup output for storage ##
  ##############################
  est.names <- get_method_label(run.NPMLEs,
                                run.NPNA,
                                run.NPNA_avg,
                                run.NPNA_wrong,
                                run.OLS,
                                run.WLS,
                                run.EFF)



  method.label <- est.names

  ####################
  ## output storage ##
  ####################
  ## for F_p(t) data
  Ft.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                  first.label.name=method.label,
                                  tval,tval0,z.use,w.use,Ft.name="Ft",p,
                                  label.dim.simus=nboot,
                                  label.name.simus=paste("boot",1:nboot,sep=""))

  Ft.boot.out <- Ft.null.theta$null.theta.simus$estimator
  Ft.estimate.boot <- Ft.null.theta$null.theta.ci$estimator

  ## for S_m(t) data
  Sout.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                    first.label.name=method.label,
                                    tval,tval0,z.use,w.use,Ft.name="St",p=m,
                                    label.dim.simus=nboot,
                                    label.name.simus=paste("boot",1:nboot,sep=""))
  Sout.boot.out <- Sout.null.theta$null.theta.simus$estimator
  St.estimate.boot <- Sout.null.theta$null.theta.ci$estimator



  ################################################
  ## setup output for storage: AUC, BS analysis ##
  ################################################
  method.label.AUCBS <- c(method.label,"NPNA-NPMLE1","NPNA-NPNA_wrong","NPMLE1-NPNA_avg")


  AUC.BS.Ft.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                         first.label.name=method.label.AUCBS,
                                         tval,tval0,z.use,w.use,Ft.name="Ft",p,
                                         label.dim.simus=nboot,
                                         label.name.simus=paste("boot",1:nboot,sep=""))

  predict.boot.out <- AUC.BS.Ft.null.theta$null.theta.simus$prediction
  Ft.AUC.BS.boot <- AUC.BS.Ft.null.theta$null.theta.ci$prediction

  AUC.BS.Sout.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                           first.label.name=method.label.AUCBS,
                                           tval,tval0,z.use,w.use,Ft.name="St",p=m,
                                           label.dim.simus=nboot,
                                           label.name.simus=paste("boot",1:nboot,sep=""))
  Sout.predict.boot.out <- AUC.BS.Sout.null.theta$null.theta.simus$prediction
  St.AUC.BS.boot <- AUC.BS.Sout.null.theta$null.theta.ci$prediction

  ## function to compute AUC_BS differences between estimators
  get_AUC_BS_differences <- function(bb,name1,name2,boot.out){
    boot.out[bb,paste(name1,name2,sep="-"),,,] <-
      boot.out[bb,name1,,,] - boot.out[bb,name2,,,]
    return(boot.out)
  }





  ######################################################
  ## create array to store estimators and differences ##
  ######################################################
  AUC.BS.Ft.null.theta.diff <- all.null.theta(theta.names=
                                                c("estimator","prediction"),
                                              first.label.name=method.label.AUCBS,
                                              tval,tval0,z.use,w.use,Ft.name="Ft",p,
                                              label.dim.simus=1,
                                              label.name.simus=paste("boot",1:1,sep=""))

  predict.differences <- AUC.BS.Ft.null.theta.diff$null.theta.simus$prediction

  AUC.BS.Sout.null.theta.diff <- all.null.theta(theta.names=
                                                  c("estimator","prediction"),
                                                first.label.name=method.label.AUCBS,
                                                tval,tval0,z.use,w.use,Ft.name="St",p=m,
                                                label.dim.simus=1,
                                                label.name.simus=paste("boot",1:1,sep=""))
  Sout.predict.differences <- AUC.BS.Sout.null.theta$null.theta.simus$prediction


  ## store in the results
  if(run.prediction.accuracy==TRUE){
    predict.differences[1,method.label,,,] <- AUC_BS_Ft[method.label,,,]
    Sout.predict.differences[1,method.label,,,] <- AUC_BS_St[method.label,,,]


    if(run.NPMLEs==TRUE & run.NPNA==TRUE){
      name1 <- "NPNA"
      name2 <- "NPMLE1"
      predict.differences <- get_AUC_BS_differences(1,name1,name2,predict.differences)
      Sout.predict.differences <- get_AUC_BS_differences(1,name1,name2,Sout.predict.differences)
    }

    if(run.NPMLEs==TRUE & run.NPNA_avg==TRUE){
      name1 <- "NPMLE1"
      name2 <- "NPNA_avg"
      predict.differences <- get_AUC_BS_differences(1,name1,name2,predict.differences)
      Sout.predict.differences <- get_AUC_BS_differences(1,name1,name2,Sout.predict.differences)

    }

    if(run.NPNA==TRUE & run.NPNA_wrong==TRUE){
      name1 <- "NPNA"
      name2 <- "NPNA_wrong"
      predict.differences <- get_AUC_BS_differences(1,name1,name2,predict.differences)
      Sout.predict.differences <- get_AUC_BS_differences(1,name1,name2,Sout.predict.differences)

    }

    predict.differences <- 	predict.differences[1,,,,]
    Sout.predict.differences <- Sout.predict.differences[1,,,,]
  }
  ###################
  ## run bootstrap ##
  ###################
  ## for testing
  ##for(bb in 1:36){

  bb <- 1
  while(bb <= nboot){
    cat("\n\n bb=",bb,"\n\n")
    mysamp <- sample(1:n,replace=TRUE)
    x.boot <- x[mysamp]
    delta.boot <- delta[mysamp]
    ww.boot <- ww[mysamp]
    zz.boot <- zz[mysamp]
    true.group.identifier.boot <- true.group.identifier[mysamp]
    q.boot <- q[,mysamp]

    ## get new qvs
    qvs.boot <- t(as.matrix(unique(t(q.boot))))

    ## re-name qvs.boot with matching names as in qvs
    qvs.boot.names <- match.names(qvs.boot,qvs)
    colnames(qvs.boot) <- paste("m",qvs.boot.names,sep="")
    qvs.boot <- qvs.boot[,order(qvs.boot.names)]

    ## get new m
    m.boot <- ncol(qvs.boot)

    ##################
    ## extract data ##
    ##################
    estimators.out <- stride.estimator(
      n,m=m.boot,
      p,
      qvs=qvs.boot,
      q=q.boot,
      x=x.boot,
      delta=delta.boot,
      ww=ww.boot,
      zz=zz.boot,
      run.NPMLEs,
      run.NPNA,
      run.NPNA_avg,
      run.NPNA_wrong,
      run.OLS,
      run.WLS,
      run.EFF,
      tval,tval0,
      z.use,w.use,
      update.qs,
      know.true.groups,
      true.group.identifier.boot,
      run.prediction.accuracy,
      do_cross_validation_AUC_BS)


    ## check if we have a problem
    problem <- estimators.out$problem
    if(is.null(problem)){
      ## there is no problem, so get the results.
      if(!is.null(z.use) & !is.null(w.use)){
        Ft.boot.out[bb,,,,,,] <- estimators.out$Ft.estimate
        Sout.boot.out[bb,,,,,,] <- estimators.out$St.estimate
      } else if(is.null(z.use) & is.null(w.use)) {
        Ft.boot.out[bb,,,,] <- estimators.out$Ft.estimate
        Sout.boot.out[bb,,,,] <- estimators.out$St.estimate
      }

      if(run.prediction.accuracy==TRUE){
        predict.boot.out[bb,method.label,,,] <- estimators.out$Ft.AUC.BS
        Sout.predict.boot.out[bb,method.label,,,] <- estimators.out$St.AUC.BS


        if(run.NPMLEs==TRUE & run.NPNA==TRUE){
          name1 <- "NPNA"
          name2 <- "NPMLE1"
          predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,predict.boot.out)
          Sout.predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,Sout.predict.boot.out)
        }

        if(run.NPMLEs==TRUE & run.NPNA_avg==TRUE){
          name1 <- "NPMLE1"
          name2 <- "NPNA_avg"

          predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,predict.boot.out)
          Sout.predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,Sout.predict.boot.out)

        }

        if(run.NPNA==TRUE & run.NPNA_wrong==TRUE){
          name1 <- "NPNA"
          name2 <- "NPNA_wrong"
          predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,predict.boot.out)
          Sout.predict.boot.out <- get_AUC_BS_differences(bb,name1,name2,Sout.predict.boot.out)

        }

      }

      ## update bb
      bb <- bb + 1
    }
  }

  get.my.variance <- function(theta.boot.out, theta.est){
    varest <- apply.index(theta.boot.out,"iters",var,na.rm=TRUE)

    ##varlo <-  apply.index(theta.boot.out,"iters",myquantiles.lo)
    ##varhi <-  apply.index(theta.boot.out,"iters",myquantiles.hi)
    varlo <- theta.est + qnorm(0.025)*sqrt(varest)
    varhi <- theta.est + qnorm(0.975)*sqrt(varest)
    list(varest=varest,varlo=varlo,varhi=varhi)
  }

  Ft.variance <- get.my.variance(Ft.boot.out,estimator_Ft)
  Sout.variance <- get.my.variance(Sout.boot.out,estimator_St)

  if(!is.null(z.use) & !is.null(w.use)){
    Ft.estimate.boot[,,,,,,"varest"]  <- Ft.variance$varest
    Ft.estimate.boot[,,,,,,"varlo"]  <- Ft.variance$varlo
    Ft.estimate.boot[,,,,,,"varhi"]  <- Ft.variance$varhi

    St.estimate.boot[,,,,,,"varest"]  <- Sout.variance$varest
    St.estimate.boot[,,,,,,"varlo"]  <- Sout.variance$varlo
    St.estimate.boot[,,,,,,"varhi"]  <- Sout.variance$varhi

  } else if(is.null(z.use) & is.null(w.use)){
    Ft.estimate.boot[,,,,"varest"]  <- Ft.variance$varest
    Ft.estimate.boot[,,,,"varlo"]  <- Ft.variance$varlo
    Ft.estimate.boot[,,,,"varhi"]  <- Ft.variance$varhi

    St.estimate.boot[,,,,"varest"]  <- Sout.variance$varest
    St.estimate.boot[,,,,"varlo"]  <- Sout.variance$varlo
    St.estimate.boot[,,,,"varhi"]  <- Sout.variance$varhi
  }



  if(run.prediction.accuracy==TRUE){
    predict.variance <- get.my.variance(predict.boot.out,predict.differences)
    Ft.AUC.BS.boot[,,,,"varest"]  <- predict.variance$varest
    Ft.AUC.BS.boot[,,,,"varlo"]  <- predict.variance$varlo
    Ft.AUC.BS.boot[,,,,"varhi"]  <- predict.variance$varhi

    Sout.predict.variance <- get.my.variance(Sout.predict.boot.out,Sout.predict.differences)
    St.AUC.BS.boot[,,,,"varest"]  <- Sout.predict.variance$varest
    St.AUC.BS.boot[,,,,"varlo"]  <- Sout.predict.variance$varlo
    St.AUC.BS.boot[,,,,"varhi"]  <- Sout.predict.variance$varhi
  }
  list(Ft.estimate.boot=Ft.estimate.boot,
       St.estimate.boot=St.estimate.boot,
       Ft.AUC.BS.boot=Ft.AUC.BS.boot,
       St.AUC.BS.boot=St.AUC.BS.boot)
}




#' Data matrix for dynamic landmark prediction of mixture data
#'
#' Create data set based on user input.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#' @param x a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.
#' @param delta a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).
#' @param ww a numeric vector of length \code{n} containing the values of the continuous
#' covariate for each person in the sample.
#' @param zz a numeric vector of length \code{n} containing the values of the discrete
#' covariate for each person in the sample.
#' @param know.true.groups logical indicator. If TRUE, then we know the population identifier for each person in the sample.
#' This option is only used for simulation studies. Default is FALSE.
#' @param true.group.identifier numeric vector of length \code{n} denoting the population identifier for each person in the sample.
#' Default is NULL.
#'
#' @return a matrix of the concatenated data. The columns are "x", "delta", "q1" to "qp" (correspond to the mixture proportions),
#' "w" (for \code{ww}), "z" (for \code{zz}), "uset" (for "u-set" meaning the uth mixture proportion group), and
#' "group" (for the true population group as defined by \code{true.group.identifier}). If \code{know.true.groups} is
#' FALSE, then the column "group" is a vector of 0's.

make.data.set <- function(
  n,m,p,qvs,q,
  x,delta,ww,zz,
  know.true.groups,true.group.identifier){

  #####################
  ## Compute u-group ##
  #####################
  uset <- compute.uset(n,m,p,qvs,q)

  if(know.true.groups==FALSE){
    true.group.identifier <- rep(0,length(delta))
  }

  #data <- data.frame(x=x,delta=delta,t(q),ww,zz,uset,true.group.identifier)
  #colnames(data) <- c("x","delta",paste("q",1:p,sep=""),"w","z","uset","group")

  dataset <-data.frame(x=x,delta=delta,t(q),uset,true.group.identifier)
  colnames_dataset <- c("x","delta",paste("q",1:p,sep=""),"uset","group")

  ## add z covariates
  if(!is.null(zz)){
    dataset <- data.frame(dataset,zz)
    colnames_dataset <- c(colnames_dataset,"z")
  }

  ## add w covariates
  if(!is.null(ww)){
    dataset <- data.frame(dataset,ww)
    colnames_dataset <- c(colnames_dataset,"w")
  }

  ## add column names
  colnames(dataset) <- colnames_dataset

  ## sort in increasing order by x
  dataset <- dataset[order(dataset$x),]
  return(dataset)
}

## function to sum NAs from array
sum_array_na <- function(m1,m2){
  #return(Reduce(`+`, lapply(list(m1,m2),function(x) {x[is.na(x)] <-0;x})))
  return(ifelse(is.na(m1), ifelse(is.na(m2), NA, m2), ifelse(is.na(m2), m1, m1 + m2)))
}

#' Estimator implementation
#'
#' Main function to estimate the distribution function for mixture data where
#' the population identifiers are unknown, but the probability of belonging
#' to a population is known. The distribution functions are evaluated at
#' time points \code{tval} and adjust for dynamic landmark prediction and one
#' discrete covariate (\code{zz}) and one continuous covariate (\code{ww}).
#'
#' @param data data matrix obtained from \code{make.data.set}
#' @param n sample size, must be at least 1.
#' @param p number of populations, must be at least 2.
#' @param m number of different mixture proportions, must be at least 2.
#' @param r numeric vector including the number of individuals in each mixture proportion
#' group.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param tval numeric vector of time points at which the distribution function is evaluated, all values must
#' be non-negative.
#' @param tval0 numeric vector of time points representing the landmark times. All values must be non-negative
#' and smaller than the maximum of \code{tval}.
#' @param method.label character vector of methods implemented. This is the result from \code{get_method_label()}/
#' @param z.use numeric vector at which to evaluate the discrete covariate \eqn{Z} at in the estimated distribution function.
#' The values of \code{z.use} must be in the range of the observed \code{zz}.
#' @param w.use numeric vector at which to evaluate the continuous covariate \eqn{W} at in the estimated distribution function.
#' The values of \code{w.use} must be in the range of the observed \code{ww}.
#' @param update.qs logical indicator. If TRUE, the mixture proportions \code{q} will be updated. This is currently not implemented.
#' @param run.prediction.accuracy logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS). Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#' @param do_cross_validation_AUC_BS logical indicator. If TRUE, then we compute the prediction accuracy measures, including the
#' area under the receiver operating characteristic curve (AUC) and the Brier Score (BS) using cross-validation. Prediction accuracy is only valid
#' in simulation studies where \code{know.true.groups}=TRUE and \code{true.group.identifier} is available.
#'
#' @section Details:
#' We estimate the distribution function for mixture data  where
#' the population identifiers are unknown, but the probability of belonging
#' to a population is known. The distribution functions are evaluated at
#' time points \code{tval} and adjust for dynamic landmark prediction and one
#' discrete covariate (\code{zz}) and one continuous covariate (\code{ww}).
#' Dynamic landmark prediction means that the distribution function is computed knowing
#' that the survival time, \eqn{T}, satisfies \eqn{T >t_0}
#' where \eqn{t_0} are the time points in \code{tval0}.
#'
#' @references
#' Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#' Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
#' In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.
#'
#' Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
#' Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
#' Annals of Applied Statistics, 8(2), 1182-1208.
#'
#' Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
#' application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association,
#' 107, 1324-1338.
#'
#'
#' @return \code{estimator.main} returns a list containing
#' \itemize{
#'    \item{Ft.store: }{a numeric array. When
#'    \code{run.prediction.accuracy} is FALSE, then the results are the
#'    the estimated distribution functions for all \code{p} populations.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{p}.
#'
#'    When \code{run.prediction.accuracy} is TRUE, then the results are the
#'    area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{p} populations. The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#'
#'    Results for both the estimated distributon functions and prediction
#'    accuracy measures (AUC, BS) are only valid when \eqn{t\geq t_0}, so
#'    arrays show NA for any combination for which \eqn{t<t_0}.
#'
#'    }
#'
#'    \item {St.store: }{a numeric array. When
#'    \code{run.prediction.accuracy} is FALSE, then the results are the
#'    the estimated distribution functions for all \code{m} mixture proportion
#'    groups.
#'    The dimension of the array is \# of methods by \code{length(tval)} by \code{lenth(tval0)} by
#'    \code{length(z.use)} by \code{length(w.use)} by \code{m}.
#'
#'    When \code{run.prediction.accuracy} is TRUE, then the results are the
#'    area under the receiver operating characteristic curve (AUC) and
#'    Brier Score (BS) for the \code{m} mixture proportion groups.
#'    The dimension of the array is \# of methods by
#'    \code{length(tval)} by \code{length(tval0)} by 2,
#'    where the last dimension stores the AUC and BS results.
#'
#'    Results for both the estimated distributon functions and prediction
#'    accuracy measures (AUC, BS) are only valid when \eqn{t\geq t_0}, so
#'    arrays show NA for any combination for which \eqn{t<t_0}.
#'
#'    }
#'
#'    \item{problem: }{a numeric indicator of errors in the NPNA estimator. If NULL, no error is reported.
#'    Otherwise, there is an error in the computation of the NPNA estimator.}
#'
#' }
#'
#' @export
estimator.main <- function(data,
                           n,p,m,r,qvs,
                           tval,tval0,
                           method.label,
                           z.use,w.use,
                           update.qs,
                           run.prediction.accuracy,
                           do_cross_validation_AUC_BS){

  #############################
  ## dummy variables needed  ##
  #############################
  problem <- NULL	## for NPNA estimator
  bootvar <- FALSE  	## for kincohort estimator

  ##################
  ## store output ##
  ##################

  ## for F_p(t) data
  nsimu <- 1
  Ft.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                  first.label.name=method.label,
                                  tval,tval0,z.use,w.use,Ft.name="Ft",p,
                                  label.dim.simus=nsimu,
                                  label.name.simus=paste("iters",1:nsimu,sep=""))

  ## for S_m(t) data
  Sout.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                    first.label.name=method.label,
                                    tval,tval0,z.use,w.use,Ft.name="St",p=m,
                                    label.dim.simus=nsimu,
                                    label.name.simus=paste("iters",1:nsimu,sep=""))

  if(run.prediction.accuracy==FALSE){
    ## Run estimator evaluated at specified (z,w) pairs

    ## storage for output
    Ft.store <- Ft.null.theta$null.theta$estimator
    Sout.store <- Sout.null.theta$null.theta$estimator

  } else {
    ## Run predictions

    ## storage for output
    Ft.store <- Ft.null.theta$null.theta$prediction
    Sout.store <- Sout.null.theta$null.theta$prediction
  }


  if(do_cross_validation_AUC_BS==TRUE){
    cv_folds <- 100
  } else {
    cv_folds <- 1
  }

  ## BEGIN OLD CODE FOR CV
  ## create data folds
  #  if(cv_folds > 1){
  #	  random_sample <- sample(1:3,size=nrow(data),
  #			replace=TRUE,prob=rep(1,cv_folds)/cv_folds)
  #	  data_folds <- setNames(split(data,random_sample), paste0("samp",1:cv_folds))
  # }
  ## END OLD CODE FOR CV



  for(cvv in 1:cv_folds){

    ###########################
    ## split data into folds ##
    ###########################
    if(cv_folds>1){
      ## BEGIN OLD CODE FOR CV
      #	data_train <- dplyr::bind_rows(data_folds[-cvv])
      #	  data_test <- data_folds[[cvv]]
      ## END OLD CODE FOR CV

      num_train <- round(0.6*n)
      train_ind <- sample(1:n, num_train, replace = FALSE)
      test_ind <- setdiff(1:n, train_ind)

      data_train <- data[train_ind,]
      data_test <- data[test_ind,]
    } else {
      data_train <- data
      data_test <- data
    }

    ## sort in increasing order by x
    data_train <- data_train[order(data_train$x),]
    data_test <- data_test[order(data_test$x),]

    n_use <- nrow(data_train)
    p_use <- p
    q_use <- t(data_train[,c("q1","q2")])
    qvs_use <- t(as.matrix(unique(t(q_use))))
    tmp_info <- get_new_qvs(qvs_use,qvs)
    qvs_use <- tmp_info$qvs_new
    m_use<- tmp_info$m_new
    r_use <- compute.r(n_use,m_use,p_use,qvs_use,q_use)

    #################
    ## update q's? ##
    #################
    get.info.use <- get.new.qs(data_train,n_use,p_use,m_use,r_use,qvs_use,tval0,
                               method.label,update.qs)


    ##########################
    ## estimation procedure ##
    ##########################
    for(tt0 in 1:length(tval0)){
      ##cat("tt0=",tt0)


      for(kk in 1:length(method.label)){
        ##cat("kk=",kk)


        ########################
        ## Index used for tt0 ##
        ########################
        if(method.label[kk]=="NPNA_wrong"){
          tt0.use <- which(tval0==0)
        } else {
          tt0.use <- tt0
        }


        ####################################################
        ## subset of data where subjects survived past t0 ##
        ####################################################
        #index_train_use <- which(data$x>tval0[tt0.use] & data$delta>0)## 8/24/2016: wrong
        index_train_use <- which(data_train$x>tval0[tt0.use]) ## 8/24/2016: this is the correct one.
        data_train_subset <- data_train[index_train_use,]

        index_test_use <- which(data_test$x>tval0[tt0.use])
        data_test_subset <- data_test[index_test_use,]

        ############################
        ## information used at t0 ##
        ############################
        info.tmp <- get.info.use[[tt0.use]]

        ## qvs.use depends on method of Ft estimator
        qvs.tmp <- info.tmp$qvs[[method.label[kk]]]

        ## r info for weights
        r.tmp <- info.tmp$r[[method.label[kk]]]

        #####################
        ## set up the data ##
        #####################
        tmp.out <- get.tmp.storage(run.prediction.accuracy,
                                   method=method.label[kk],
                                   z.use,w.use,
                                   data_test_subset,p,m)

        data_test_zw <- tmp.out$data.zw

        Ft_test_out <- tmp.out$Ft.tmp.out
        Sout_test_out <- tmp.out$Sout.tmp.out

        ## information to run kin-cohort estimators. From data_train
        q.tmp <- t(info.tmp$q[[method.label[kk]]])
        n.tmp <- nrow(data_train)
        m.tmp <- length(table(data_train$q1))
        x.tmp <- data_train[,"x"]
        delta.tmp <- data_train[,"delta"]


        ##########################
        ## NPMLE estimator at t0 ##
        ##########################
        if(!grepl("NPNA*",method.label[kk])){
          set.run <- set.kin.run(method.label[kk])

          ## run kin-cohort estimator
          kin.out0 <- kincohort.estimators(n.tmp,
                                           q=q.tmp,
                                           x=x.tmp,
                                           delta=delta.tmp,
                                           tval0[tt0],
                                           qvs=qvs.tmp,
                                           p,m.tmp,
                                           r=r.tmp,
                                           boot=0,bootvar,
                                           useOLS=set.run$run.OLS,
                                           useWLS=set.run$run.WLS,
                                           useEFF=set.run$run.EFF,
                                           useNPMLEs=set.run$run.NPMLEs)
        }

        for(tt in 1:length(tval)){
          ##cat("tt=",tt)

          ###########################################
          ## Valid to compute F(t|t0) when t> t0. ##
          ###########################################
          if(tval[tt] > tval0[tt0]){

            ##################################
            ## NPMLE,OLS,WLS, EFF estimator ##
            ##################################

            if(!grepl("NPNA*",method.label[kk])){
              ## run kin-cohort estimator
              kin.out <- kincohort.estimators(n.tmp,
                                              q=q.tmp,
                                              x=x.tmp,
                                              delta=delta.tmp,
                                              tval[tt],
                                              qvs=qvs.tmp,
                                              p,m.tmp,
                                              r=r.tmp,
                                              boot=0,bootvar,
                                              useOLS=set.run$run.OLS,
                                              useWLS=set.run$run.WLS,
                                              useEFF=set.run$run.EFF,
                                              useNPMLEs=set.run$run.NPMLEs)

              if(method.label[kk]=="NPMLE1"){
                ## 1-S(t|t0,z,w)
                myHout <- (kin.out$hts0-kin.out0$hts0)/
                  (1-kin.out0$hts0)

                mySout <- 1- myHout

                ## F(t|t0,z,w)
                u.tmp <- t(qvs.tmp)

                myFest <- solve(t(u.tmp) %*% diag(r.tmp) %*% u.tmp) %*%
                  t(u.tmp) %*% diag(r.tmp) %*% myHout

              } else{
                myFest <- (kin.out$Fest[,method.label[kk]]-kin.out0$Fest[,method.label[kk]])/
                            (1-kin.out0$Fest[,method.label[kk]])

                ## We put 0's for Sout_test_out. We do  not compute this for OLS, WLS, EFF, NPMLE2 estimators
                mySout <- rep(0,m.tmp)
              }
              ## Repeat the estimates for Ft and St since these are the same for all covariate values.
              print(method.label[kk])
              print(myFest)
              Ft_test_out[,paste("Ft",1:p,sep="")] <- rep.row(myFest,nrow(Ft_test_out))


              Sout_test_out[,paste("St",1:m,sep="")] <- rep.row(mySout,nrow(Sout_test_out))



            } else {

              ####################
              ## NPNA estimator ##
              ####################

              #####################################
              ## what(z,w,u) data to evaluate at ##
              #####################################
              data_test_zw_evaluate <- get.zw.evaluate(data_test_zw,m)

              ######################################
              ## get data for landmark estimation ##
              ######################################
              data.land <- data_train_subset[,c("x","delta","z","w","uset")]

              ## run NPNA estimator at tt, and report values at data.evaluate
              S.NPNA.zw.out <- S.NPNA.zw(t=tval[tt],
                                         data=data.land,
                                         newdata=data_test_zw_evaluate)

              ## report any problems from S.NPNA.zw
              problem.tmp <- S.NPNA.zw.out$problem
              if(!is.null(problem.tmp)){
                problem <- problem.tmp
              }

              ## report output
              NPNA.out <-  as.data.frame(S.NPNA.zw.out$data.out)

              ## Form F(t),St.estimate
              for(ll in 1:nrow(Ft_test_out)){

                ## get appropriate rows for (z,w) estimates
                z.tmp <- Ft_test_out[ll,"z"]
                w.tmp <- Ft_test_out[ll,"w"]

                ## get appropriate subset of NPNA.out
                NPNA.out.tmp <- NPNA.out[which(NPNA.out$w==w.tmp & NPNA.out$z==z.tmp),]
                Sout.tmp <- NPNA.out.tmp$survival.v.new
                Sout_test_out[ll,paste("St",1:m,sep="")] <- Sout.tmp

                u.index <- NPNA.out.tmp$u
                u.tmp <- t(qvs.tmp[,u.index])

                ## FOR TESTING. NEED TO REMOVE!!!
                ##r.tmp <- rep(1,length(r.tmp))

                ## form (U^TU)^{-1} U^T (1-S(t|z,w))
                Ft_test_out[ll,paste("Ft",1:p,sep="")] <-
                  solve(t(u.tmp) %*% diag(r.tmp) %*% u.tmp) %*%
                  t(u.tmp) %*%  diag(r.tmp) %*% (1-Sout.tmp)


              }

            }

            if(run.prediction.accuracy==TRUE){
              ###################################
              ## output results for prediction ##
              ###################################


              ########################
              ## NPNA-avg estimator ##
              ########################
              if(method.label[kk]=="NPNA_avg"){
                Ft.tmp.avg <- apply(Ft_test_out[,paste("Ft",1:p,sep="")],2,mean)
                Ft_test_out[,paste("Ft",1:p,sep="")] <- rep.row(Ft.tmp.avg,nrow(Ft_test_out))

                Sout.tmp.avg <- apply(Sout_test_out[,paste("St",1:m,sep="")],2,mean)
                Sout_test_out[,paste("St",1:m,sep="")] <-
                  rep.row(Sout.tmp.avg,nrow(Sout_test_out))
              }

              ## need to feed in data for all samples
              Ft.predict.tmp <- matrix(NA,nrow(data_test))
              Sout.predict.tmp <- matrix(NA, nrow(data_test))

              for(ii in 1:length(index_test_use)){
                u <- data_test_subset[ii,"uset"]
                Sout.predict.tmp[index_test_use[ii]] <- Sout_test_out[ii,paste("St",u,sep="")]
                Ft.predict.tmp[index_test_use[ii]] <-
                  Ft_test_out[ii,paste("Ft",
                                       data_test[index_test_use[ii],"group"],sep="")]
              }

              Sout.predict.tmp.use <- Sout.predict.tmp[index_test_use]
              Ft.predict.tmp.use <- Ft.predict.tmp[index_test_use]

              ## feed in original data
              data.land.predict <- data_test_subset[,c("x","delta","z","w","uset")]
              Sout.prediction <- get.predictions(data.land.predict,
                                                 Ft=1-Sout.predict.tmp.use,
                                                 t0=tval0[tt0],tau=tval[tt]-tval0[tt0])

              Sout.store[kk,tt,tt0,"AUC"] <- sum_array_na(Sout.store[kk,tt,tt0,"AUC"],
                                                          Sout.prediction$AUC$AUC.est/cv_folds)

              Sout.store[kk,tt,tt0,"BS"] <- sum_array_na(Sout.store[kk,tt,tt0,"BS"],
                                                         Sout.prediction$BS$Brier.score/cv_folds)

              Ft.prediction <- get.predictions(data.land.predict,Ft=Ft.predict.tmp.use,
                                               t0=tval0[tt0],tau=tval[tt]-tval0[tt0])

              Ft.store[kk,tt,tt0,"AUC"] <- sum_array_na(Ft.store[kk,tt,tt0,"AUC"],
                                                        Ft.prediction$AUC$AUC.est/cv_folds)

              Ft.store[kk,tt,tt0,"BS"] <- sum_array_na(Ft.store[kk,tt,tt0,"BS"],
                                                       Ft.prediction$BS$Brier.score/cv_folds)

            } else {

              ######################################
              ## output results at specific (z,w) ##
              ######################################
              if(method.label[kk]!="NPNA_avg"){
                Ft.store[kk,tt,tt0,,,] <- unflatten.array(Ft.store[kk,tt,tt0,,,],
                                                          dim.order=c("zz","ww","Ft"),
                                                          Ft_test_out[,paste("Ft",1:p,sep="")],
                                                          flatten.name="Ft")

                Sout.store[kk,tt,tt0,,,] <- unflatten.array(Sout.store[kk,tt,tt0,,,],
                                                            dim.order=c("zz","ww","Ft"),
                                                            Sout_test_out[,paste("St",1:m,sep="")],
                                                            flatten.name="Ft")
              } else{
                ## take average over all (z,w) combinations
                Ft.tmp.avg <- apply(Ft_test_out[,paste("Ft",1:p,sep="")],2,mean)
                Sout.tmp.avg <- apply(Sout_test_out[,paste("St",1:m,sep="")],2,mean)

                if(!is.null(z.use) & !is.null(w.use)){
                  ## repeat the estimates for all (z,w) combinations
                  Ft.store[kk,tt,tt0,,,] <- repeat.zw(Ft.tmp.avg,z.use,w.use,p)
                  Sout.store[kk,tt,tt0,,,] <- repeat.zw(Sout.tmp.avg,z.use,w.use,m)
                } else if(is.null(z.use) & is.null(w.use)){
                  ## report the estimates
                  Ft.store[kk,tt,tt0,] <- Ft.tmp.avg
                  Sout.store[kk,tt,tt0,] <- Sout.tmp.avg
                }
              }
            }
          }
        }
      }
    }
  }
  list(Ft.store=Ft.store,Sout.store=Sout.store,problem=problem)
}



#' Method names
#'
#' Creates character vector of method names.
#'
#' @param run.NPMLEs a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data based on the type-I and type II nonparametric maximum likelihood
#' estimators. The type I nonparametric maximum likelihood estimator is referred
#' to as the "Kaplan-Meier" estimator in Garcia and Parast (2020). Neither the type I nor type II
#' include covariates or dynamic landmarking.
#' @param run.NPNA a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that accounts for covariates and dynamic
#' landmarking. This estimator is called "NPNA" in Garcia and Parast (2020).
#' @param run.NPNA_avg a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that averages out over the observed covariates.
#' This is referred to as NPNA_marg in Garcia and Parast (2020).
#' @param run.NPNA_wrong a logical indicator. If TRUE, then the output includes the
#' estimated distribution function for mixture data that ignores landmarking. This is referred to as
#' NPNA_{t_0=0} in the paper.
#' @param run.OLS a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using an ordinary least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.WLS  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using a weighted least squares influence function that adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#' @param run.EFF  a logical indicator. If TRUE, then the output includes the estimated distribution function
#' computed using the efficient influence function based on Hilbert space projection theory results.
#' The estimator adjusts for censoring using
#' inverse probability weighting (IPW), augmented inverse probability weighting (AIPW),
#' and imputation (IMP). See details in Wang et al (2012). These estimators do not include covariates nor dynamic landmarking.
#'
#' @references
#' Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#' Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
#' In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.
#'
#' Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
#' Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
#' Annals of Applied Statistics, 8(2), 1182-1208.
#'
#' Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
#' application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association,
#' 107, 1324-1338.
#'
#' @return A character vector of method names. If \code{run.NPMLEs} is TRUE,
#' the character vector includes "NPMLE1" and "NPMLE2". If \code{run.NPNA} is TRUE,
#' the character vector includes "NPNA". If \code{run.NPNA_avg} is TRUE,
#' the character vector includes "NPNA_avg". If \code{run.NPNA_wrong} is TRUE,
#' the character vector includes "NPNA_wrong". If \code{run.OLS} is TRUE,
#' the character vector includes "OLSIPW","OLSAIPW","OLSIMP". If \code{run.WLS} is TRUE,
#' the character vector includes "WLSIPW","WLSAIPW","WLSIMP". If \code{run.EFF} is TRUE,
#' the character vector includes "EFFIPW","EFFAIPW","EFFIMP".
#'
#' @export
get_method_label <- function(run.NPMLEs,
                             run.NPNA,
                             run.NPNA_avg,
                             run.NPNA_wrong,
                             run.OLS,
                             run.WLS,
                             run.EFF){

  ## Names for methods run
  est.names <- NULL
  mylabel <- c("IPW","AIPW","IMP")

  if(run.NPMLEs==TRUE){
    est.names <- c(est.names,paste("NPMLE",1:2,sep=""))
  }

  if(run.OLS==TRUE){
    est.names <- c(est.names,paste("OLS",mylabel,sep=""))
  }

  if(run.WLS==TRUE){
    est.names <- c(est.names,paste("WLS",mylabel,sep=""))
  }

  if(run.EFF==TRUE){
    est.names <- c(est.names,paste("EFF",mylabel,sep=""))
  }

  if(run.NPNA==TRUE){
    est.names <- c(est.names,"NPNA")
  }

  if(run.NPNA_avg==TRUE){
    est.names <- c(est.names,"NPNA_avg")
  }

  if(run.NPNA_wrong==TRUE){
    est.names <- c(est.names,"NPNA_wrong")
  }

  return(est.names)
}

## another function to re-order qvs
get_new_qvs <- function(qvs_new,qvs){
  ## re-name qvs_new with matching names as in qvs

  qvs_new_names <- match.names(qvs_new,qvs)
  colnames(qvs_new) <- paste("m",qvs_new_names,sep="")
  qvs_new <- qvs_new[,order(qvs_new_names)]

  ## get new m
  m_new <- ncol(qvs_new)
  return(list(qvs_new=qvs_new,m_new=m_new))
}



## function to have qvs.boot and qvs have the same column names
## and in the correct order.
match.names <- function(qvs.boot,qvs){
  index.bucket <- NULL
  for(kk in 1:ncol(qvs.boot)){
    for(ii in 1:ncol(qvs)){
      if(all(qvs.boot[,kk]==qvs[,ii])){
        index.bucket <- c(index.bucket,ii)
      }
    }
  }
  return(index.bucket)
}



#############################################
## functions to run kincohort estimators   ##
#############################################
#' Nonparametric kincohort estimators
#'
#' Wrapper function to run kincohort estimators that
#' ignore covariates and landmarking.
#'
#' @param n sample size, must be at least 1.
#' @param q a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#' @param x a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.
#' @param delta a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).
#' @param t numeric value at which the distribution function is evaluated.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param p number of populations, must be at least 2.
#' @param m number of different mixture proportions, must be at least 2.
#' @param r numeric vector including the number of individuals in each mixture proportion
#' group.
#' @param boot number of bootstrap replicates.
#' @param bootvar logical indicator. If TRUE, we compute the bootstrap variance estimates of
#' the estimators.
#' @param useOLS logical indicator. If TRUE, we compute the distribution function
#' for the mixture data where the influence function is the ordinary least squares estimator.
#' @param useWLS logical indicator. If TRUE, we compute the distribution function
#' for the mixture data where the influence function is the weighted least squares estimnator.
#' @param useEFF logical indicator. If TRUE, we compute the distribution function
#' for the mixture data where the influence function is the efficient estimator.
#' @param useNPMLEs logical indicator. If TRUE, we compute the distribution function
#' for the mixture data based on the type I and type II
#' nonparametric maximum likelihood esimators (NPMLEs).
#'
#' @return a list containing
#' \itemize{
#'    \item{hts0: }{numeric vector of length \code{m} containing the
#'    Kaplan-Meier estimates of the survival curves for the \code{m} mixture
#'    proportion groups.}
#'    \item{Fest: }{numeric array containing the estimated distribution functions for
#'    all influence functions of interest.}
#'    \item{var_est: }{numeric array containing the estimated variances for
#'    all influence functions of interest.}
#'    \item{eflag: }{numeric value indicating if errors occurred in computing  the estimated distribution functions.
#'    \code{eflag}<0 indicates an error.}
#'
#' }
#'
#' @export
#' @useDynLib stride
kincohort.estimators <- function(n,q,x,delta,
                                 t,qvs,p,m,r,
                                 boot,bootvar,
                                 useOLS,useWLS,useEFF,
                                 useNPMLEs){
  ## set up for f90
  storage.mode(n) <- "integer"
  storage.mode(q) <- "double"
  storage.mode(x) <- "double"
  storage.mode(delta) <- "double"
  storage.mode(t) <- "double"
  storage.mode(qvs) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(r) <- "integer"
  storage.mode(boot) <- "integer"
  storage.mode(bootvar) <- "logical"
  storage.mode(useOLS) <- "logical"
  storage.mode(useWLS) <- "logical"
  storage.mode(useEFF) <- "logical"
  storage.mode(useNPMLEs) <- "logical"

  num_estimators <- 11 ## fixed
  storage.mode(num_estimators) <- "integer"

  ## set as standard values since we don't make use of these objects
  lim <- 0; storage.mode(lim) <- "double"
  rat <- 0; storage.mode(rat) <- "double"
  setting <- "regcoxcase"; storage.mode(setting) <- "character"
  d <- 0; storage.mode(d) <- "double"
  H0 <- TRUE; storage.mode(H0) <- "logical"
  usetruth <- FALSE; storage.mode(usetruth) <- "logical"

  ## prepare output
  mylabel <- c("IPW","AIPW","IMP")
  Fest <- array(0,dim=c(p,num_estimators),dimnames=list(paste("p",1:p,sep=""),
                                                        c(paste("OLS",mylabel,sep=""),
                                                          paste("WLS",mylabel,sep=""),
                                                          paste("EFF",mylabel,sep=""),
                                                          paste("NPMLE",1:2,sep=""))))
  var_est <- array(0,dim=c(p,p,num_estimators),dimnames=list(paste("p",1:p,sep=""),
                                                             paste("p",1:p,sep=""),
                                                             c(paste("OLS",mylabel,sep=""),
                                                               paste("WLS",mylabel,sep=""),
                                                               paste("EFF",mylabel,sep=""),
                                                               paste("NPMLE",1:2,sep=""))))
  hts0 <- rep(0,m)
  eflag <- 0
  storage.mode(Fest) <- "double"
  storage.mode(var_est) <- "double"
  storage.mode(eflag) <- "integer"
  storage.mode(hts0) <- "double"

  ##if(t>0){
  ## COMMENT THIS OUT FOR THE R PACKAGE
  out <- .Fortran("kincohort_estimators",n,q,x,delta,t,qvs,p,m,r,lim,rat,setting,d,H0,
                  num_estimators,boot,bootvar,usetruth,useOLS,useWLS,useEFF,useNPMLEs,
                  hts0=hts0,
                  Fest=Fest,var_est=var_est,eflag=eflag,
                  PACKAGE="stride")


  ## UNCOMMENT FOR RUNNING CODE IN UNIX
  #  out <- .Fortran("kincohort_estimators",n,q,x,delta,t,qvs,p,m,r,lim,rat,setting,d,H0,
  #                  num_estimators,boot,bootvar,usetruth,useOLS,useWLS,useEFF,useNPMLEs,
  #                  hts0=hts0,
  #                  Fest=Fest,var_est=var_est,eflag=eflag)

  ## } else {
  ##   out <- list(hts0=hts0,Fest=Fest,var_est=var_est,eflag=eflag)
  ## }
  list(hts0=out$hts0,Fest=out$Fest,var_est=out$var_est,eflag=out$eflag)
}


#' Alpha quantile
#'
#' Computes the alpha-quantile of a vector \code{x}.
#'
#' @param x numeric vector of values.
#'
#' @return \code{myquantiles.lo} returns the \eqn{\alpha/2}-quantile.
#'
#' @export
myquantiles.lo <-function(x){
  alpha <- 0.05
  B <- length(x)
  lo <- sort(x)[B*alpha/2]
  return(lo)
}

#' One minus alpha quantile
#'
#' Computes the \eqn{1-\alpha} quantile of a vector \code{x}.
#'
#' @param x numeric vector of values.
#'
#' @return \code{myquantiles.hi} returns the \eqn{1-\alpha/2}-quantile.
#'
#' @export
myquantiles.hi <-function(x){
  alpha <- 0.05
  B <- length(x)
  up <- sort(x)[B*(1-alpha/2)]
  return(up)
}

#' Sample sizes of mixture proportion groups.
#'
#' Computes the number of individuals in each mixture proportion group.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q.use a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#'
#' @export
compute.r <- function(n,m,p,qvs,q.use){
  r <- rep(0,m)		    ## number in each qvs subgroup
  for(jj in 1:n){
    for(ii in 1:m){
      q.tmp <- all(q.use[1:p,jj]==qvs[,ii])
      if(q.tmp==TRUE){
        r[ii] <- r[ii] + 1
      }
    }
  }
  return(r)
}


## how to set run.OLS, run.WLS, run.EFF, run.NPMLEs based on \code{method.label.unit}
#' @importFrom utils glob2rx
set.kin.run <- function(method.label.unit){
  run.OLS <- FALSE
  run.WLS <- FALSE
  run.EFF <- FALSE
  run.NPMLEs <- FALSE

  if(length(grep(utils::glob2rx("OLS*"),method.label.unit))>0){
    run.OLS <- TRUE
  }

  if(length(grep(utils::glob2rx("WLS*"),method.label.unit))>0){
    run.WLS <- TRUE
  }

  if(length(grep(utils::glob2rx("EFF*"),method.label.unit))>0){
    run.EFF <- TRUE
  }

  if(length(grep(utils::glob2rx("NPMLE*"),method.label.unit))>0){
    run.NPMLEs <- TRUE
  }

  list(run.OLS=run.OLS,run.WLS=run.WLS,run.EFF=run.EFF,run.NPMLEs=run.NPMLEs)
}

#############################################################################
## Store data set that will be used to estimate the distribution functions ##
## and compute the prediction accuracy measures.                           ##
#############################################################################
get.tmp.storage <- function(run.prediction.accuracy,method,
                            z.use,w.use,data,p,m){

  if(!is.null(z.use) & !is.null(w.use)){
    if(run.prediction.accuracy==TRUE | method=="NPNA_avg" ){
      ## storage for calculation
      data.zw <- data[,c("z","w")]
    } else{
      ## storage for calculation
      data.zw <- get.zw.combinations(z.use,w.use)
    }

    ## Ft and Sout have n rows
    Ft.tmp.out <- tmp.storage(data.tmp=data.zw,dim.col=p,dim.names="Ft")


    Sout.tmp.out <- tmp.storage(data.tmp=data.zw,dim.col=m,dim.names="St")

    ## remove duplicates
    data.zw <- data.zw[!duplicated(data.zw),]
  } else if(is.null(z.use) & is.null(w.use)){
    data.zw <- NULL
    Ft.tmp.out <- array(NA,dim=c(1,p),dimnames=list(1,paste0("Ft",1:p)))
    Sout.tmp.out <- array(NA,dim=c(1,m),dimnames=list(1,paste0("St",1:m)))
  }
  list(data.zw=data.zw,Ft.tmp.out=Ft.tmp.out,Sout.tmp.out=Sout.tmp.out)
}


## Repeats a vector into rows of a matrix
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


## function to get z,w,u combinations to evaluate S_j(t|z,w) = S(t|z,w,u_j)
get.zw.evaluate <- function(data.evaluate.orig,m){
  data.evaluate <- data.evaluate.orig[rep(seq_len(nrow(data.evaluate.orig)),each=m),]
  data.evaluate <- cbind(data.evaluate,rep(1:m,nrow(data.evaluate)/m))

  colnames(data.evaluate) <- c("z","w","u")
  rownames(data.evaluate) <- 1:nrow(data.evaluate)

  return(data.evaluate)
}




## function to repeat Ft entries
repeat.zw <- function(x,z.use,w.use,p){
  ## entries go in as p, w, z
  ## entries go out as z by w by p
  x.tmp <- array(x,dim=c(p,length(w.use),length(z.use)))
  x.tmp <- aperm(x.tmp,perm=c(3,2,1))
  return(x.tmp)
}

## Function to update the mixture proportions. Not implemented yet.
get.new.qs <- function(data,n,p,m,r,qvs,tval0,method.label,update.qs){

  ############################
  ## separate results by t0 ##
  ############################
  new.output <- vector("list",length(tval0))
  names(new.output) <- paste("t0",tval0,sep="")

  #####################################
  ## initialize values for r, qvs, q ##
  #####################################
  r.use <- vector("list",length(method.label))
  names(r.use) <- method.label
  qvs.use <- r.use
  q.use <- r.use

  for(kk in 1:length(method.label)){
    r.use[[kk]] <- r
    qvs.use[[kk]] <- qvs
    q.use[[kk]] <- data[,paste("q",1:p,sep="")]
  }

  ######################
  ## update r, qvs, q ##
  ######################
  if(update.qs==TRUE){
    ## update q's

    for(tt0 in 1:length(tval0)){
      out.tmp <- kin.updateq(data,n,p,m,r.orig=r.use,qvs.orig=qvs.use,q.orig=q.use,
                             t0=tval0[tt0],
                             method.label)

      new.output[[tt0]] <- list(qvs=out.tmp$qvs.new.out,
                                q=out.tmp$q.new.out,
                                r=out.tmp$r.new.out)

      ## update r.use,qvs.use,q.use
      r.use <- new.output[[tt0]]$r
      qvs.use <- new.output[[tt0]]$qvs
      q.use <- new.output[[tt0]]$q
    }
  } else {
    ## no update to the q's, so at each t0, we use the original values.
    for(tt0 in 1:length(tval0)){
      new.output[[tt0]] <- list(qvs=qvs.use,q=q.use,r=r.use)
    }
  }
  return(new.output)
}





#' Mixture porportion group identifier
#'
#' Computes the mixture proportion group identifier.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q.use a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#'
#' @section Details:
#' The matrix \code{qvs} contains all mixture proportions. The function
#' \code{compute.uset} assigns individual \eqn{i} to mixture proportion \eqn{u},
#' for \eqn{u=1,\ldots,m}, if \eqn{q_i=qvs_{u}}.
#'
#' @return A numeric vector of length \code{n} containing the mixture proportion
#' group to which each person belongs.
#'
#' @export
compute.uset <- function(n,m,p,qvs,q.use){
  uset <- rep(0,n)          ## indicator of which subgroup qvs
  for(jj in 1:n){
    for(ii in 1:m){
      q.tmp <- all(q.use[1:p,jj]==qvs[,ii])
      if(q.tmp==TRUE){
        uset[jj] <- ii
      }
    }
  }
  return(uset)
}




## function to make all (z,w) combinations
get.zw.combinations <- function(z.use,w.use){
  data.evaluate.list <- list(z=z.use,w=w.use)
  data.evaluate.orig <- expand.grid(data.evaluate.list)
  return(data.evaluate.orig)
}



## get storage matrices for Ft, Sout
# dim.col=p
tmp.storage <- function(data.tmp,dim.col,dim.names="Ft"){
  out <- matrix(NA,nrow=nrow(data.tmp),ncol=dim.col)
  out.name <- paste(dim.names,1:dim.col,sep="")
  out <- cbind(data.tmp,out)
  colnames(out) <- c(colnames(data.tmp),out.name)
  return(out)
}






## update q,qvs, r. Function not implemented yet.
kin.updateq <- function(data,n,p,m,r.orig,qvs.orig,q.orig,t0,method.label){

  bootvar <- FALSE  ## no need to run variance

  ####################
  ## output storage ##
  ####################
  q.new.out <- vector("list",length(method.label))
  names(q.new.out) <- method.label

  qvs.new.out <- q.new.out
  r.new.out <- q.new.out


  for(kk in 1:length(method.label)){
    set.run <- set.kin.run(method.label[kk])

    ################################################
    ## set up q.orig, qvs.orig, r.orig for method ##
    ################################################
    q.orig.tmp <- q.orig[[method.label[kk]]]
    qvs.orig.tmp <- qvs.orig[[method.label[kk]]]
    r.orig.tmp <- r.orig[[method.label[kk]]]

    ####################
    ## get estimators ##
    ####################
    kin.out <- kincohort.estimators(n,
                                    q=t(q.orig.tmp),
                                    x=data$x,
                                    delta=data$delta,
                                    t0,
                                    qvs=qvs.orig.tmp,
                                    p,m,
                                    r=r.orig.tmp,
                                    boot=0,bootvar,
                                    useOLS=set.run$run.OLS,
                                    useWLS=set.run$run.WLS,
                                    useEFF=set.run$run.EFF,
                                    useNPMLEs=set.run$run.NPMLEs)

    Fest.tmp <- t(kin.out$Fest[,method.label[kk]])

    ## update q's
    q.new <- (q.orig.tmp[,1]*(1-Fest.tmp[1]))/
      (q.orig.tmp[,1]*(1-Fest.tmp[1])+q.orig.tmp[,2]*(1-Fest.tmp[2]))
    q.new.out[[method.label[kk]]] <- cbind(q.new,1-q.new)
    colnames(q.new.out[[method.label[kk]]]) <- paste("q",1:p,sep="")

    ## update qvs
    qvs.new <- (qvs.orig.tmp[1,]*(1-Fest.tmp[1]))/
      (qvs.orig.tmp[1,]*(1-Fest.tmp[1])+qvs.orig.tmp[2,]*(1-Fest.tmp[2]))
    qvs.new.out[[method.label[kk]]] <- rbind(qvs.new,1-qvs.new)

    ## update r
    r.new.out[[method.label[kk]]] <- t(as.numeric(table(q.new)))
    colnames(r.new.out[[method.label[kk]]]) <- paste("m",1:m,sep="")
  }

  list(q.new.out=q.new.out,qvs.new.out=qvs.new.out,r.new.out=r.new.out)
}

####################################
####################################
##
##
## Modified functions from landpred
##
##
####################################
####################################

## Computes the modified nonparametric Nelson-Aalen estimator.
#' @import stats
pred.smooth.surv <- function(w.vector=NULL, t, data.use, covariate.value, group.value, weight)
{ Xi.use = data.use[,1]
Di.use = data.use[,2]
Zi.use = data.use[,3]
Wi.use = data.use[,4]
Ui.use = data.use[,5]
if(is.null(w.vector)) {w.vector = Wi.use}
h = bw.nrd(Wi.use[Zi.use == covariate.value & Ui.use == group.value])
bandwidth = h*length(Xi.use[Zi.use == covariate.value & Ui.use == group.value])^(-.10)
K = Kern.FUN(Wi.use[Zi.use == covariate.value & Ui.use == group.value],w.vector,bandwidth)
Xi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,1]
Di.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,2]
Wi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,4]
tmpind = (Xi.temp<=t)&(Di.temp==1)
tj = Xi.temp[tmpind];
kerni.1 = t(weight[Zi.use == covariate.value & Ui.use == group.value]*t(K))
pihamyt0.tj.ss = sum.I(tj, "<=", Xi.temp, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##
dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss;
#dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
ret = apply(dLamhat.tj.ss,2,sum)
S.return  =exp(-ret)
return(S.return)
}



###################################################
##SURVIVAL USING Z and W COVARIATE INFORMATION######
###################################################
#data structure: first column should be observed event time (or censoring time), second column should be indicator of whether an event was observed (Delta), third column should be discrete covariate Z, fourth column should be continuous covariate W, fifth column should indicate the U group that each person is in (for example if there are 3 "u" groups, u_1, u_2, u_3, then the fifth column would either have "1","2", or "3"); returns the data matrix with an extra column, the extra column is the survival probability at that Z and W
#newdata is an optional n by 3 matrix where the first column is the discrete covariate Z, the second column is the continuous covariate W and the third column is the U group. Predicted survival probabilities are estimated for these data;returns the data matrix with an extra column, the extra column is the survival probability at that Z and W

S.NPNA.zw = function(t, data, newdata = NULL, weight = NULL){
  problem <- NULL
  Xi = data[,1]
  Di = data[,2]
  Zi = data[,3]
  Wi = data[,4]
  Ui = data[,5]

  if(sum(Xi > t) == 0) {print(paste("No observations past time t=",t))}
  if(is.null(weight)) {weight = rep(1,length(Xi))}


  zi.cat = unique(Zi)
  ui.cat = unique(Ui)

  ## 3/16/2017: Below will only run if newdata is NULL
  if(is.null(newdata)){
    survival.v <- rep(NA, length = dim(data)[1]) ## Changed 1/30/2017
    ##survival.v = vector(length = dim(data)[1])
    for(j in 1:length(ui.cat)) {
      ##print(j)
      for(k in 1:length(zi.cat)) {
        ##print(k)
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]

        if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
          print(paste("Warning: Very few individuals with covariate value = ",
                      zi.value, ",
                      and in U group = ",ui.value))}

        if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
          print(paste("Warning: Very few individuals
                      observed to survive past t=",t," with covariate value = ",
                      zi.value, ", and in U group = ",ui.value))
        }

        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){  ## Changed 1/30/2017
          problem <- TRUE
          print(paste("Warning-Error: Less than 2 individuals with W covariate for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0.  W is",
                      Wi[Zi == zi.value & Ui == ui.value ], sep=""))
          P.return <- 0
        } else {
          P.return <- pred.smooth.surv(w.vector = Wi[Zi == zi.value & Ui == ui.value], t=t,
                                       data.use = data, covariate.value = 	zi.value,
                                       group.value = ui.value, weight = weight)
        }
        survival.v[Zi == zi.value & Ui == ui.value] <- P.return
        }
      }
    data = cbind(data, survival.v)
}

  if(!is.null(newdata)) {
    Zi.new = newdata[,1]
    Wi.new = newdata[,2]
    Ui.new = newdata[,3]

    survival.v.new = rep(NA,length = dim(newdata)[1]) ## Changed 1/30/2017
    ##survival.v.new = vector(length = dim(newdata)[1])
    for(j in 1:length(ui.cat)) {
      for(k in 1:length(zi.cat)) {
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]

        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
        #  print(paste("Warning: Very few individuals with covariate value = ",
        #  			zi.value, ",
        #    and in U group = ",ui.value))}

        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
        #  print(paste("Warning: Very few individuals
        #                  observed to survive past t=",t," with covariate value = ",
        #                  zi.value, ", and in U group = ",ui.value))
        #}


        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){ ## Changed 1/30/2017
          problem <- TRUE
          print(paste("Warning-Error: Less than 2 individuals with W covariate
                      for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0. W is",
                      Wi[Zi == zi.value & Ui == ui.value],
                      sep=""))
          P.return <- 0
        } else {
          P.return = pred.smooth.surv(w.vector = Wi.new[Zi.new == zi.value &
                                                          Ui.new == ui.value],
                                      t=t,data.use = data, covariate.value =  zi.value, group.value = ui.value,
                                      weight = weight)
        }
        survival.v.new[Zi.new == zi.value & Ui.new == ui.value] = P.return
      }
    }
    newdata.matrix = cbind(newdata, survival.v.new)
  }
  if(is.null(newdata)) {return(list(data.out=data,problem=problem))}
  if(!is.null(newdata)) {return(list(data.out=newdata.matrix,problem=problem))}
}


## function to get AUC, and BS calculations
#' @importFrom landpred AUC.landmark
#' @importFrom landpred BS.landmark
get.predictions <- function(data.land,Ft,t0,tau){
  data.for.prediction <- cbind(data.land$x,data.land$delta,
                               rep(NA,length(data.land$delta)),Ft)
  AUC <- landpred::AUC.landmark(t0=t0,tau=tau,short=FALSE,data=data.for.prediction)
  BS <- landpred::BS.landmark(t0=t0,tau=tau,short=FALSE,data=data.for.prediction)
  list(AUC=AUC,BS=BS)
}


cumsum2 <- function(mydat)     #cumsum by row, col remains the same
{
  if(is.null(dim(mydat))) return(cumsum(mydat))
  else{
    out <- matrix(cumsum(mydat), nrow=nrow(mydat))
    out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
    return(out)
  }
}

########################
##HELPER FUNCTION######
########################

sum.I <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{
  if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
  if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
  pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
  if(is.null(Vi)){return(pos)}else{
    Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
    out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
    out[pos!=0,] <- Vi[pos,]
    if(is.null(dim(Vi))) out <- c(out)
    return(out) ## n.y x p
  }
}




########################
##HELPER FUNCTION######
########################
Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix
{
  out = (VTM(zz,length(zi))- zi)/bw
  norm.k = dnorm(out)/bw
  norm.k
}


########################
##HELPER FUNCTION######
########################
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

####################################################
####################################################
##                                                ##
##                                                ##
## Functions to organize parameter storage        ##
##                                                ##
##                                                ##
####################################################
####################################################


## Produces arrays with NA components.
## The arrays will be used to store the results for the estimated
## distribution function and prediction accuracy measures (e.g., AUC, BS).
get.null.theta <- function(theta.names=c("estimator","prediction"),
                           first.label.name,
                           tval,tval0,z.use,w.use,Ft.name="Ft",p){


  null.theta <- get.empty.list(theta.names)

  ## modify length of null.theta depending on inclusion of w.use, z.use
  if(!is.null(z.use) & !is.null(w.use)){
    null.theta[["estimator"]] <- array(NA,dim=c(length(first.label.name),
                                                length(tval),
                                                length(tval0),
                                                length(z.use),
                                                length(w.use),p),
                                       dimnames=list(
                                         method=first.label.name,
                                         time=paste("tt",tval,sep=""),
                                         time0=paste("tt0",tval0,sep=""),
                                         zz=z.use,
                                         ww=w.use,
                                         Ft=paste(Ft.name,1:p,sep="")))
  } else if(is.null(z.use) & is.null(w.use)){
    null.theta[["estimator"]] <- array(NA,dim=c(length(first.label.name),
                                                length(tval),
                                                length(tval0),p),
                                       dimnames=list(
                                         method=first.label.name,
                                         time=paste("tt",tval,sep=""),
                                         time0=paste("tt0",tval0,sep=""),
                                         Ft=paste(Ft.name,1:p,sep="")))

  }

  null.theta[["prediction"]] <- array(NA,dim=c(length(first.label.name),
                                               length(tval),
                                               length(tval0),2),
                                  dimnames=list(
                                    method=first.label.name,
                                    time=paste("tt",tval,sep=""),
                                    time0=paste("tt0",tval0,sep=""),
                                    measure=c("AUC","BS")
                                  ))

  return(null.theta)
}

## add dimension to each array in a list of arrays.
add.dimension.null <- function(null.theta,location,label.dim,label.name){
  null.theta0 <- null.theta
  for(u in 1:length(null.theta)){
    if(!is.null(null.theta0[[u]])){
      if(location=="first"){
        new.names <- list(iters=list(label.name),dimnames(null.theta0[[u]]))

        null.theta[[u]] <- array(NA,dim=c(label.dim,dim(null.theta0[[u]])),
                                 dimnames=unlist(new.names,recursive=FALSE))
      } else if(location=="last"){
        new.names <- list(dimnames(null.theta0[[u]]),val=list(label.name))
        null.theta[[u]] <- array(NA,dim=c(dim(null.theta[[u]]),label.dim),
                                 dimnames=unlist(new.names,recursive=FALSE))
      }
    }
  }
  return(null.theta)
}

## get all list of arrays used in the analysis.
all.null.theta <- function(theta.names=c("estimator","prediction"),
                           first.label.name,
                           tval,tval0,z.use,w.use,Ft.name,p,
                           label.dim.simus,label.name.simus){

  null.theta <- get.null.theta(theta.names,
                               first.label.name,
                               tval,tval0,z.use,w.use,Ft.name,p)
  ## remove method
  method.names <- which(names(dimnames(null.theta$estimator))=="method")
  null.theta.nomethod <- array(NA,dim=dim(null.theta$estimator)[-method.names],
                               dimnames=dimnames(null.theta$estimator)[-method.names])


  ## add simus layer
  null.theta.simus <- add.dimension.null(null.theta,location="first",
                                         label.dim=label.dim.simus,label.name=label.name.simus)

  ## ci results layer
  null.theta.ci <-  add.dimension.null(null.theta,location="last",
                                       label.dim=3,label.name=c("varest","varlo","varhi"))

  ## est, ci results
  ##null.theta.est.ci <-  add.dimension.null(null.theta,location="last",
  ##              label.dim=4,label.name=c("est","varest","varlo","varhi"))

  ## simus and ci results layer to null.theta
  ##null.theta.simus.ci <-  add.dimension.null(null.theta.simus,location="last",
  ##              label.dim=3,label.name=c("varest","varlo","varhi"))

  ## simus, est, ci results
  null.theta.simus.est.ci <- add.dimension.null(null.theta.simus,location="last",
                                                label.dim=4,label.name=c("est","varest","varlo","varhi"))

  list(null.theta=null.theta,null.theta.simus=null.theta.simus,
       ##null.theta.simus.ci=null.theta.simus.ci, ## not used
       null.theta.ci=null.theta.ci,
       null.theta.simus.est.ci=null.theta.simus.est.ci,
       ##null.theta.est.ci =  null.theta.est.ci, ## not used
       null.theta.nomethod=null.theta.nomethod )
}

## Change a matrix to an array.
unflatten.array <- function(x,dim.order,data.change,
                            flatten.name="time"){

  ## get correct permutation order of array
  xorig <- x
  if(!is.null(dim.order)){
    x <- aperm(x,dim.order)
  }

  ## get index of flatten variable
  dim.flatten <- which(names(dimnames(x))==flatten.name)

  ## find differences
  dim.use <- 1:length(dimnames(x))
  dim.use <- setdiff(dim.use,dim.flatten)

  out  <- as.matrix(data.change[, unlist(dimnames(x)[dim.flatten])])
  dim(out) <- dim(x)[c(dim.use,dim.flatten)]
  dimnames(out) <- dimnames(x)[c(dim.use,dim.flatten)]

  ## put order as original
  array.order <-  match(names(dimnames(xorig)),names(dimnames(out)))
  out <- aperm(out,array.order)

  return(out)
}


find.apply.index <- function(x,names.use){
  out <- which(names(dimnames(x)) %in% names.use)
  out <- setdiff(1:length(dimnames(x)),out)
  return(out)
}

apply.index <- function(x,names.use,fun,...){
  out <- apply(x,find.apply.index(x,names.use),fun,...)
  return(out)
}

#' Create empty list
#'
#' Returns an empty list where list names are determined by \code{names}.
#'
#' @param names character vector of names for the list
#'
#' @examples
#' list_names <- c("kitty","dog","bird")
#' get.empty.list(list_names)
#'
#' @export
get.empty.list <- function(names){
  out <-  vector("list",length(names))
  names(out) <- names
  return(out)
}




#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
##
## Functions to produce Error Messages
##
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+

#' Check number of predictors are valid
#'
#' Reports an error if the number_of_functional_predictors =0 and
#' number_of_nonfunctional_predictors=0. Our model runs only if we have at least one
#' functional_predictor or at least one nonfunctional_predictor.
#'
#' @param number_of_functional_predictors Number of functional predictors.
#' @param number_of_nonfunctional_predictors Number of nonfunctional predictors.
#'
#' @return Error if  number_of_functional_predictors =0 and
#' number_of_nonfunctional_predictors=0.
#'
#' @examples
#' ## Returns error (not run)
#' #number_of_functional_predictors <- 0
#' #number_of_nonfunctional_predictors <- 0
#' #error_number_of_predictors(number_of_functional_predictors,number_of_nonfunctional_predictors)
#'
#' @export
error_number_of_predictors <- function(number_of_functional_predictors,number_of_nonfunctional_predictors){
  return(invisible((
    if(number_of_functional_predictors==0 & number_of_nonfunctional_predictors==0){
      stop("No parameters to estimate.
           Change number_of_functional_predictors>=1 and/or number_of_nonfunctional_predictors>=1")
    }


  )))
}


#' Check positivity
#'
#' Stops code if  \code{x} is not >0.
#'
#' @param x numeric vector.
#' @param variable_name name of object x.
#'
#' @return Error message if any component of \code{x} is <= 0.
#'
#' @examples
#' ## Returns error (not run)
#' #test <- c(-1,2)
#' #error_positive_value(test,get_variable_name(test))
#'
#' @export
error_positive_value <- function(x,variable_name){
  return(invisible(if(any(x <=0)){
    stop(variable_name," must be >0.")
  }
  ))
}





#' Incorrect type
#'
#' Reports an error if \code{x} is not of type \code{y}.
#'
#' @param x character vector or array that must be of type \code{y}.
#' @param y character vector of allowed types.
#'
#' @return Error message if any component of \code{x} is not of type \code{y}.
#'
#' @examples
#' ## Reports an error (not run)
#' #x <- c("cat","dog")
#' #y <- c("cat","bird","fish")
#'
#' #stop_options_message(x,y)
#'
#' @export
stop_options_message <- function(x,y){
  if(is.vector(x)){
    return(invisible(
      if(!all(x %in% y)){
        stop("Options of ",deparse(substitute(x))," must be one of ", paste(y,collapse=" or "),".")
      }
    ))
  } else {
    ## arrays/matrices
    return(invisible(
      if(any(!all(x%in% y))==TRUE){
        stop("Options of ",deparse(substitute(x))," must be one of ", paste(y,collapse=" or "),".")
      }
    ))

  }
}

#' Incorrect size or type
#'
#' Reports error if \code{type} is not the correct size specified by \code{number}
#' or if \code{type} is not one of the \code{type_options}.
#'
#' @param type character vector or array.
#' @param number required length of \code{type} if \code{type} is a vector, or
#' required dimension of \code{type} if \code{type} is an array.
#' @param type_options character of what \code{type} is allowed to be.
#'
#'
#' @return Error message if any component of \code{type} is not the correct size or
#' allowed type.
#'
#' @examples
#' ## Reports an error because of wrong type. (not run)
#' #number_of_functional_predictors <- 3
#' #functional_coefficients_basis_type <- rep("cubic-spline",number_of_functional_predictors)
#' #functional_basis_options <- c("fpca","quantlet")
#' #error_type_number(functional_coefficients_basis_type,number_of_functional_predictors,
#' #functional_basis_options)
#'
#' @export
error_type_number <- function(type,number,type_options){
  return(invisible((
    if(number>=0){
      stop_length_message(type,number)
      stop_options_message(type,type_options)
    }
  )))
}



#' Error from exceeding maximum
#'
#' Reports error if x > max_x
#'
#' @param x numeric vector or array object.
#' @param max_x maximum allowed value in the vector or array.
#' @param variable_name name of object x.
#'
#' @return Error message if any component of \code{x} is larger than \code{max_x}.
#'
#' @examples
#' ## Reports an error (not run)
#' ##love <- c(1,2)
#' ##max_love <- 1.5
#' ##error_max_value(love,max_love,get_variable_name(love))
#'
#' @export
error_max_value <- function(x,max_x,variable_name){
  return(invisible(if(any(x >max_x)){
    stop(variable_name," must be <=", max_x)
  }
  ))
}



#' Variable name
#'
#' Reports variable name.
#'
#' @param variable any object type.
#'
#' @return Character name of variable.
#'
#' @examples
#' love <- c(1,2)
#' get_variable_name(love)
#'
#' dog <- matrix(1:4,2,2)
#' get_variable_name(dog)
#'
#' @export
get_variable_name <- function(variable) {
  return(deparse(substitute(variable)))
}


#' Check non-negativity
#'
#' Stops code if  \code{x} is not >=0.
#'
#' @param x numeric vector.
#' @param variable_name name of object x.
#'
#' @return Error message if any component of \code{x} < 0.
#'
#' @examples
#' ## Returns error (not run)
#' ##test <- c(-1,2)
#' ##error_nonpositive_value(test,get_variable_name(test))
#'
#' @export
error_nonpositive_value <- function(x,variable_name){
  return(invisible(if(any(x <0)){
    stop(variable_name," must be >=0.")
  }
  ))
}



#' Check dimension requirement
#'
#' Stops code if dimension of \code{x} is not as specified by \code{y}.
#'
#' @param x vector that should be length \code{y}, or an array that should have
#' dimensions specified by \code{y}.
#' @param y numeric scalar indicating required length of x,
#' or vector of scalars indicating required dimensions of x.
#'
#' @return Error message if \code{x} does not have dimensions specified by \code{y},
#' otherwise, nothing is returned.
#'
#' @examples
#' ## Returns error (not run)
#' #x <- c(1,2,3)
#' #y <- 4
#' #stop_length_message(x,y)
#'
#'
#' ##Returns error (not run)
#' #x <- matrix(1:4,2,2)
#' #y <- c(1,2)
#' #stop_length_message(x,y)
#'
#' @export
stop_length_message <- function(x,y){
  if(is.vector(x)){
    return(invisible(
      if(length(x)!=y){
        stop("Length of ",deparse(substitute(x))," must be ",y,".")
      }
    ))
  } else {
    ## arrays/matrices
    return(invisible(
      if(any(dim(x)!=y)==TRUE){
        stop("Dimensions of ",deparse(substitute(x))," must be ",paste0(y,collapse=" by "),".")
      }
    ))

  }
}


