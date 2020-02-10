#' Generate data
#'
#' Produces data from different populations with the probability of
#' belonging to a population. Also produces one discrete covariate and one continuous covariate.
#'
#' @param n sample size, must be at least 1.
#' @param p number of populations, must be at least 2.
#' @param m number of different mixture proportions, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param censoring.rate a scalar indicating the censoring proportion. Options are 0 or 40.
#' @param simu.setting Character indicating simulation setting.
#' Options are "Log-Normal-No-Covariates", "Log-Normal-With-Covariates", "HD-No-Covariates","HD-With-Covariates".
#' Setting "Log-Normal-No-Covariates" and "Log-Normal-With-Covariates" refer to simulation setting 1 in Garcia and Parast (2020).
#' "Log-Normal-No-Covariates" means the
#' survival outcomes do NOT depend on the covariates, and "Log-Normal-With-Covariates" means the
#' survival outcomes do depend on the covariates.
#' Setting "HD-No-Covariates" and "HD-With-Covariates" refer to Simulation setting 2 in Garcia and Parast (2020), "HD-No-Covariates" means the
#' survival outcomes do NOT depend on the covariates, and "HD-With-Covariates" means the
#' survival outcomes do depend on the covariates.
#'
#' @return Returns a list containing
#' \itemize{
#'    \item{x: }{a numeric vector of length \code{n} containing the observed event times
#' for each person in the sample.}
#'    \item{delta: }{a numeric vector of length \code{n} that denotes
#' censoring (1 denotes event is observed, 0 denotes event is censored).}
#'    \item{q: }{a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.}
#'    \item{ww: }{a numeric vector of length \code{n} containing the values of the continuous
#' covariate for each person in the sample.}
#'    \item{zz: }{ a numeric vector of length \code{n} containing the values of the discrete
#' covariate for each person in the sample.}
#'    \item{true.group.identifier: }{numeric vector of length \code{n} denoting the population identifier for each person in the sample.
#'    This is only used simulation study to evaluate the prediction accuracy of our methods.}
#' }
#'
#'@references
#'Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.
#'
#'
#' @export
GenerateData <- function(n,p,m,qvs,censoring.rate,simu.setting){


  q <- array(0,dim=c(p,n))  ## mixture proportions
  x <- rep(0,n)		    ## observed event time: min(T,C)
  delta <- rep(0,n)	    ## censoring indicator
  uset <- rep(0,n)	    ## indicator of which subgroup qvs
  r <- rep(0,m)		    ## number in each qvs subgroup
  true.group.identifier <- rep(0,n)   ## to which group each person belongs

  ###################
  ## Simulate data ##
  ###################

  ## set evenly spaced subgroup sizes
  r0 <- seq(0,1,length=m+1)

	for(i in 1:n){
      a <- runif(1)
      for(j in 1:m){
        if( (a >= r0[j]) && (a < r0[j+1]) ){
          q[,i] <- qvs[,j]
	  r[j]  <- r[j] + 1
	  uset[i] <- j
        }
      }
    }

    ## set covariates
    zw.data <- gendata.zw(n,simu.setting)
    zz <- zw.data$zz
    ww <- zw.data$ww

    ## set dependence on covariates
    if(simu.setting=="Log-Normal-No-Covariates" | simu.setting=="HD-No-Covariates"){
      covariate.dependent <- FALSE
    } else {
      covariate.dependent <- TRUE
    }


    for(i in 1:n){
      q1<- q[,i]

	  negative_value <- TRUE
	  while(negative_value==TRUE){
		t1 <- trueinvFt(p,ww[i],zz[i],simu.setting,covariate.dependent)
		if(all(t1>0)){
			negative_value <- FALSE
		}
	  }

      a <- runif(1)
      tmp <- q1[1]
      j <- 1
      while(j <= p){
        if(a <= tmp){
          s <- t1[j]
	  true.group.identifier[i] <- j
	  j <- p+1
        } else {
          tmp <- tmp + q1[j+1]
	  j <- j+1
        }
      }

      cens <- genc(s,censoring.rate,simu.setting)
      x[i] <- min(s,cens)
      delta[i] <- which.min(c(cens,s)) - 1
    }

  list(x=x,delta=delta,q=q,ww=ww,zz=zz,true.group.identifier=true.group.identifier)
}


################################
## Functions to generate F(t) ##
################################

## function to get sd's used in F(t)
getsd <- function(p){
  sd.use <- seq(1,2,length.out=p)
  return(sd.use)
}

mult.w <- function(covariate.dependent){
  if(covariate.dependent==TRUE){
    out <- 1
  } else {
    out <- 0
  }
  return(out)
}

mult.z <- function(covariate.dependent){
  if(covariate.dependent==TRUE){
    out <- 0.5
  } else {
    out <- 0
  }
  return(out)
}

#' @import stats
Ft.form <- function(tt,tt0,ww,zz,p,simu.setting){
	out <- rep(0,p)

	## set dependence on covariates
	if(simu.setting=="Log-Normal-No-Covariates" | simu.setting=="HD-No-Covariates"){
	  covariate.dependent <- FALSE
	} else {
	  covariate.dependent <- TRUE
	}


	if(simu.setting == "Log-Normal-No-Covariates" | simu.setting=="Log-Normal-With-Covariates"){
		sd.use <- getsd(p)
      	constant.z <- mult.z(covariate.dependent)
      	constant.w <- mult.w(covariate.dependent)

		for(jj in 1:p){
			out[jj] <-
			(pnorm(log(tt)- constant.w * ww -
				constant.z * zz,mean=0,sd=sd.use[jj])-

			 pnorm(log(tt0)- constant.w * ww -
				   constant.z * zz,mean=0,sd=sd.use[jj]))/
    	   (1-pnorm(log(tt0)- constant.w * ww -
				   constant.z * zz, mean=0,sd=sd.use[jj]))
		}
	}else if(simu.setting=="HD-No-Covariates" | simu.setting=="HD-With-Covariates"){
      	           out[1] <-
		( F1(tt,ww,zz,covariate.dependent)- F1(tt0,ww,zz,covariate.dependent) ) /
      		      	 ( 1 - F1(tt0,ww,zz,covariate.dependent) )

                   out[2] <-
		( F2(tt,ww,zz,covariate.dependent)- F2(tt0,ww,zz,covariate.dependent) ) /
      		      	( 1 - F2(tt0,ww,zz,covariate.dependent) )
	}
	return(out)
}


## Function to produce true F(t)
trueFt <- function(z.use,w.use,tval,tval0,
			p,real_data,simu.setting,
			Ft.null.theta, method.label){

   ## storage
   ##out <- Ft.null.theta$null.theta.nomethod

   out <- Ft.null.theta$null.theta$estimator

   ## set dependence on covariates
   if(simu.setting=="Log-Normal-No-Covariates" | simu.setting=="HD-No-Covariates"){
     covariate.dependent <- FALSE
   } else {
     covariate.dependent <- TRUE
   }

  ####################
  ## simu setting 1 ##
  ####################
      	## function to generate T from true F(t)
      	## Since log(T)=W+0.5Z+N, where N~Normal(0,\sigma_u). Then, T=exp(W+0.5Z+N).
      	## S(t|z,w,sigma_u) = Pr(T > t |z,w,sigma_u) = Pr{ exp(W+0.5Z+N) > t | z, w, sigma_u}
      ##  = Pr(W+0.5Z+N > log(t)|z,w,sigma_u)
      ##  = Pr(N > log(t) - W - 0.5Z|sigma_u)
      ##  = 1- pnorm(log(t)-W-0.5Z,mean=0,sd=sigma_u)
      ## So F(t|z,w,u)=pnorm(log(t)-W-0.5Z,mean=0,sd=sigma_u)
      ## note: pnorm() is the distribution function.

      ## In the general case,
      ## pr(T <= t | T> t_0,z,w,sigma_u)= (pr(T<=t|z,w,sigma_u) - pr(T<= t_0|z,w,sigma_u))/pr(T>t_0|z,w,sigma_u)

  ####################
  ## simu setting 2 ##
  ####################
  ## real data

  if(real_data==FALSE){
	for(ee in 1:length(method.label)){
      for(tt in 1:length(tval)){
        for(tt0 in 1:length(tval0)){

			if(tval[tt]> tval0[tt0]){
   			    if(method.label[ee]!="NPNA_avg"){

				  for(zz in 1:length(z.use)){
					  for(ww in 1:length(w.use)){
						out[ee,tt,tt0,zz,ww,] <-
							Ft.form(tval[tt],tval0[tt0],
							w.use[ww],z.use[zz],p,simu.setting)
				  	 }
				  }
			    } else {
				  ## NPNA-avg

				   ## generate z,w data
				  zw.data <- gendata.zw(n=100000,simu.setting)
				  zz.tmp <- zw.data$zz
				  ww.tmp <- zw.data$ww

				  ## compute average over (z,w) pairs
				  avg.tmp <- rep(0,p)
				  for(ii in 1:length(zz.tmp)){
					 avg.tmp <- avg.tmp + Ft.form(tval[tt],tval0[tt0],
							ww.tmp[ii],zz.tmp[ii],
							p,simu.setting)
				  }

				  ## report the avg
				  for(zz in 1:length(z.use)){
					for(ww in 1:length(w.use)){
						out[ee,tt,tt0,zz,ww,] <- avg.tmp/100000
					}
				  }

			  }
			}
		}
	  }
    }
  }
  return(out)
}

## randomly choose subgroup sizes
#' @import stats
random.subgroup.sizes <- function(m,n){
  rc <- rep(0,m+1)
  rc[2:m] <- round(runif((m-1))*(n-m))
  rc[m+1] <- n

  rc <- sort(rc)
  rc[2:m] <- sort(rc[2:m])+2:m-1

  return(rc)
}


## function to generate mixture proportions
set.qvs.values <- function(p,m,qvs.setting,mydata=NULL){
  qvs <- matrix(0,nrow=p,ncol=m)
  colnames(qvs) <- paste("m",1:m,sep="")
  rownames(qvs) <- paste("p",1:p,sep="")

  if(qvs.setting=="original"){
	qvs <- qvs.values(p,m)
  } else if(qvs.setting=="HD"){
    ##qvs[1,1:m] <- unique(mydata$q1)  ## probability of being carrier
    qvs[,1:m] <- t(as.matrix(unique(mydata[,paste("q",1:p,sep="")])))
  }

}

#' Generate finite set of mixture proportions
#'
#' Produces the finite set of mixture proportions for simulated data.
#'
#' @param p number of populations, must be at least 2.
#' @param m number of different mixture proportions, must be at least 2.
#'
#' @return Returns a \code{p} by \code{m} matrix of mixture proportions.
#'
#' @export
qvs.values <- function(p,m){
  qvs <- matrix(0,nrow=p,ncol=m)
  colnames(qvs) <- paste("m",1:m,sep="")
  rownames(qvs) <- paste("p",1:p,sep="")

  tmp <- 0.2
  qvs[1,1] <- 1
  qvs[1,2] <- (1+tmp)/2
  qvs[1,3] <- tmp
  qvs[1,m] <- 1+tmp*tmp/4-tmp/2-3/4

  qvs[1,] <- sort(qvs[1,])
  qvs[2,] <- 1-qvs[1,]

  return(qvs)
}

## function for censoring
#' @import stats
genc <- function(s,censoring.rate,simu.setting){
  if(censoring.rate==0){
    out <- s + 1
  } else if(censoring.rate==20){
    if(simu.setting=="Log-Normal-No-Covariates"){
      out <- runif(1,0,25)
    } else if(simu.setting=="Log-Normal-With-Covariates"){
      out <- runif(1,0,22)
    } else if(simu.setting=="HD-No-Covariates"){
      #out <- runif(1,50,90)
      ##out <- runif(1,0,80)
      out <- runif(1,0,210)
    } else if(simu.setting=="HD-With-Covariates"){
      out <- runif(1,0,220)
    }
  } else if(censoring.rate==40){
    if(simu.setting=="Log-Normal-No-Covariates"){
      out <- runif(1,0,8)
    } else if(simu.setting=="Log-Normal-With-Covariates"){
      out <- runif(1,0,4)
    } else if(simu.setting=="HD-No-Covariates"){
      #out <- runif(1,50,90)
      ##out <- runif(1,0,80)
      out <- runif(1,10,100)
    } else if(simu.setting=="HD-With-Covariates"){
      out <- runif(1,15,100)
    }
  }
  return(out)
}





############################################
## Functions used to replicate HD setting ##
############################################
F.exp <- function(t,w,z,mu0,mu1,mu2,sigma,r,w.constant){
  out <- ( 1+exp( - ((t-mu0)+w.constant * (w-mu1)+mu2*z)/sigma)  )^r
  return(out)
}

normalizing.constant <- function(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,w.constant){
  out <- a*(t.low) + F.exp(t.high,w,z,mu0,mu1,mu2,sigma,r,w.constant)-F.exp(t.low,w,z,mu0,mu1,mu2,sigma,r,w.constant)
  return(out)
}

intersection.constant <- function(t.low,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  if(p==2){
    out <- a*t.low - F.exp(t.low,w,z,mu0,mu1,mu2,sigma,r,w.constant)
  } else {
    out <- 0
  }
  return(out)
}

time.generate <- function(u,t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  out <- mu0 - w.constant*(w-mu1) - mu2*z -
       sigma * log ((u *
       	     normalizing.constant(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,w.constant)-
       	       	       intersection.constant(t.low,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant) )^(1/r)-1)
  return(out)
}

log.limit <- function(u,t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  out <- (u * normalizing.constant(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,
      	    a,w.constant)-
                       intersection.constant(t.low,w,z,mu0,mu1,mu2,sigma,r,a,
		       p,w.constant) )^(1/r)-1
  return(out)
}

F.constants <- function(p,covariate.dependent=TRUE){
  if(p==1){
    mu1 <- 40
    mu0 <- 43
    if(covariate.dependent==TRUE){
      ## covariates impact onset ages
      mu2 <- -0.5
      w.constant <- 1
    } else {
      mu2 <- 0
      w.constant <- 0
    }
    sigma <- 7
    a <- 0
    t.low <- 0
    t.high <- 100
    r <- -0.9
  } else {
    mu1 <- 42
    mu0 <- 48
    if(covariate.dependent==TRUE){
      ## covariates impact onset ages
      w.constant <- 1
      mu2 <- -0.5
    } else {
      mu2 <- 0
      w.constant <- 0
    }
    sigma <- 10.5
    a <- 0.0007
    t.low <- 13
    t.high <- 100
    r <- -2
  }
  list(mu0=mu0,mu1=mu1,mu2=mu2,sigma=sigma,a=a,t.low=t.low,t.high=t.high,r=r,
	w.constant=w.constant)
}

## function to set w
get.w.set <-function(w,mu1,covariate.dependent){
  if(covariate.dependent==TRUE){
    ## survival times depend on covariates
    w.set <- w
  } else {
    ## survival times DO NOT depend on covariates
    w.set <- mu1
  }
 return(w.set)
}

F1 <- function(t,w,z,covariate.dependent){
  constants <- F.constants(1,covariate.dependent)
  mu0 <- constants$mu0
  mu1 <- constants$mu1
  mu2 <- constants$mu2
  w.constant <- constants$w.constant
  sigma <- constants$sigma
  a <- constants$a
  t.low <- constants$t.low
  t.high <- constants$t.high
  r <- constants$r
  ## set w-value
  w.set <- get.w.set(w,mu1,covariate.dependent)
  out <- F.exp(t,w=w.set,z,mu0,mu1,mu2,sigma,r,w.constant)/
      	 normalizing.constant(t.low,t.high,w=w.set,z,mu0,mu1,mu2,sigma,r,a,w.constant)
  return(out)
}

F2 <- function(t,w,z,covariate.dependent){
  constants <- F.constants(2,covariate.dependent)
  mu0 <- constants$mu0
  mu1 <- constants$mu1
  mu2 <- constants$mu2
  w.constant <- constants$w.constant
  sigma <- constants$sigma
  a <- constants$a
  t.low <- constants$t.low
  t.high <- constants$t.high
  r <- constants$r

  ## set w-value
  w.set <- get.w.set(w,mu1,covariate.dependent)
  if(t<=t.low){
    out <- a * t
  } else {
    out <- intersection.constant(t.low,w=w.set,z,mu0,mu1,mu2,sigma,r,a,p=2,w.constant) + F.exp(t,w=w.set,z,mu0,mu1,mu2,sigma,r,w.constant)
  }
  out <- out/normalizing.constant(t.low,t.high,w=w.set,z,mu0,mu1,mu2,sigma,r,a,w.constant)
  return(out)
}



## function to generate T from true F(t)
#' @import stats
trueinvFt <- function(p,w,z,simu.setting,covariate.dependent){
  if(simu.setting =="Log-Normal-No-Covariates" | simu.setting=="Log-Normal-With-Covariates"){
    out <- rep(0,p)
    constant.z <- mult.z(covariate.dependent)
    constant.w <- mult.w(covariate.dependent)
    sd.use <- getsd(p)
    for(jj in 1:p){
      ## Model: log(T)=W+0.5Z+N, so T=exp(W+0.5Z+N)
      out[jj] <- exp(constant.w * w + constant.z * z + rnorm(1,0,sd.use[jj]))
    }
  } else if(simu.setting=="HD-No-Covariates" | simu.setting=="HD-With-Covariates"){
    out <- rep(0,p) ## p=2 (fixed)

    ## get constants needed
    constants.F1 <- F.constants(1,covariate.dependent)
    mu0.F1 <- constants.F1$mu0
    mu1.F1 <- constants.F1$mu1
    mu2.F1 <- constants.F1$mu2
    w.constant.F1 <- constants.F1$w.constant
    sigma.F1 <- constants.F1$sigma
    a.F1 <- constants.F1$a
    t.low.F1 <- constants.F1$t.low
    t.high.F1 <- constants.F1$t.high
    r.F1 <- constants.F1$r

    ## set w-value
    w.set.F1 <- get.w.set(w,mu1.F1,covariate.dependent)

    use.value <- FALSE
    while(use.value==FALSE){
      u <- runif(1)
      if( log.limit(u,t.low.F1,t.high.F1,w=w.set.F1,z,mu0.F1,
		mu1.F1,mu2.F1,sigma.F1,r.F1,a.F1,p=1,w.constant.F1) > 0 ){
        use.value <- TRUE
      }
    }
    out[1] <- time.generate(u,t.low.F1,t.high.F1,w=w.set.F1,z,mu0.F1,
  	    mu1.F1,mu2.F1,sigma.F1,r.F1,a.F1,p=1,w.constant.F1)

    ## get constants needed
    constants.F2 <- F.constants(2,covariate.dependent)
    mu0.F2 <- constants.F2$mu0
    mu1.F2 <- constants.F2$mu1
    mu2.F2 <- constants.F2$mu2
    w.constant.F2 <- constants.F2$w.constant
    sigma.F2 <- constants.F2$sigma
    a.F2 <- constants.F2$a
    t.low.F2 <- constants.F2$t.low
    t.high.F2 <- constants.F2$t.high
    r.F2 <- constants.F2$r

    ## set w-value
    w.set.F2 <- get.w.set(w,mu1.F2,covariate.dependent)

    use.value <- FALSE
    while(use.value==FALSE){
      u <- runif(1)
      if(u < a.F2 * t.low.F2 / normalizing.constant(t.low.F2,
    	       t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,r.F2,
	       a.F2,w.constant.F2) |
	       log.limit(u,t.low.F2,t.high.F2,w=w.set.F2,z,mu0.F2,
                mu1.F2,mu2.F2,sigma.F2,r.F2,a.F2,p=2,w.constant.F2) > 0 ) {
        use.value <- TRUE
      }
    }

    if(u <  a.F2 * t.low.F2 / normalizing.constant(t.low.F2,
               t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,r.F2,
	       a.F2,w.constant.F2)){
      out[2] <- u * normalizing.constant(t.low.F2,
               t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,
	       r.F2,a.F2,w.constant.F2)/ a.F2
    } else {
      out[2] <- time.generate(u,t.low.F2,t.high.F2,w=w.set.F2,z,mu0.F2,
  	    mu1.F2,mu2.F2,sigma.F2,r.F2,a.F2,p=2,w.constant.F2)
    }
  }
  return(out)
}


################################
## Function to generate (z,w) ##
################################

## function to get w-distribution
w.sample <- function(n){
  cag.values <- seq(30,60,by=1)
  out <- sample(x=cag.values,n,replace=TRUE,prob=rep(1,length(cag.values))/length(cag.values))
  return(out)
}

#' @import stats
gendata.zw <- function(n,simu.setting){
  zz <- rbinom(n,1,0.5)

  if(simu.setting=="Log-Normal-No-Covariates" | simu.setting=="Log-Normal-With-Covariates"){
    ww <- runif(n,0,1)
  } else if(simu.setting=="HD-No-Covariates" | simu.setting=="HD-With-Covariates"){
    ww <- w.sample(n)
    ##ww <- runif(n,30,60)
  }


  ## setting without covariates: set z=w=0
  ##zz <- rep(0,n)
  ##ww <- rep(0,n)

  list(zz=zz,ww=ww)
}

