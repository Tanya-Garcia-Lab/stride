# Setup parameters to generate the data
set.seed(1)
censoring.rate <- 40
p <- 2
n <- 2000
m <- 4
tval <- c(20,40,60,80)
tval0 <- c(0)
z.use <- c(0,1)
w.use <- seq(35,55,by=1)
run.prediction.accuracy <- FALSE
update.qs <- FALSE
simu.setting <- "HD-No-Covariates"
covariate.dependent <- TRUE
run.NPMLEs <- FALSE
run.NPNA <- TRUE
run.OLS <- FALSE
run.WLS <- FALSE
run.EFF <- FALSE
run.NPNA_avg <- FALSE
run.NPNA_wrong <- FALSE
do_cross_validation_AUC_BS <- FALSE
know.true.groups <- TRUE


## compute the finite set of mixture proportions
qvs <- qvs.values(p,m)

## generate the data

data.gen <- GenerateData(n,p,m,qvs,censoring.rate,simu.setting)

x <- data.gen$x
delta <- data.gen$delta
q <- data.gen$q
ww <- data.gen$ww
zz <- data.gen$zz


## true group membership (needed to compute the AUC/BS for simulated data)
true.group.identifier <- data.gen$true.group.identifier

## Perform the estimation
estimators.out <- stride.estimator(n,m,p,qvs,q,
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

## Show results for the estimates
## estimators.out$Ft.estimate
## estimators.out$St.estimate

## Show results for prediction accuracy AUC and BS measures (only valid for simulated data)
## estimators.out$Ft.AUC.BS
## estimators.out$St.AUC.BS


## NOT RUN
## Do bootstrap variance
#nboot <- 100
#variance.estimation <- TRUE

#varboot <- stride.bootstrap.variance(
#						nboot,n,m,p,qvs,q,
#						x,delta,ww,zz,
#						run.NPMLEs,
#						run.NPNA,
#						run.NPNA_avg,
#						run.NPNA_wrong,
#           run.OLS,
#           run.WLS,
#           run.EFF,
#						tval,tval0,
#						z.use,w.use,
#						update.qs,
#						know.true.groups,
#						true.group.identifer,
#						estimator_Ft=estimators.out$Ft.estimate,
#						estimator_St=estimators.out$St.estimate,
#						AUC_BS_Ft=estimators.out$Ft.AUC.BS,
#						AUC_BS_St=estimators.out$St.AUC.BS,
#						run.prediction.accuracy,
#						do_cross_validation_AUC_BS=FALSE)

## Show results for the bootstrap variances of the estimates
## varboot$Ft.estimate.boot
## varboot$St.estimate.boot


## Show results for the bootstrap variances of the prediction accuracy measures, AUC and BS
## varboot$Ft.AUC.BS.boot
## varboot$St.AUC.BS.boot



