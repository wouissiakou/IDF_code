# This file contains:
# -IDF-package description
# -IDF.agg for preparing the data
# -IDF.plot for plotting of IDF curves at a chosen station


#### IDF-package ####

#' Introduction
#' @name IDF-package
#' @docType package
#' @description This package provides functions to estimate IDF relations for given
#' precipitation time series on the basis of a duration-dependent
#' generalized extreme value distribution (d-GEV). 
#' The central function is \code{\link{gev.d.fit}}, which uses the method 
#' of maximum-likelihood estimation for the d-GEV parameters, whereby it is 
#' possible to include generalized linear modeling 
#' for each parameter. This function was implemented on the basis of \code{\link[ismev]{gev.fit}}.
#' For more detailed information on the methods and the application of the package for estimating 
#' IDF curves with spatial covariates, see Ulrich et. al (2020). 
#' @details 
#' * The __d-GEV__ is defined following Koutsoyiannis et al. (1998): 
#' \deqn{G(x)= \exp[-( 1+\xi(x/\sigma(d)- \tilde{\mu}) ) ^{-1/\xi}] } 
#' defined on \eqn{ \{ x: 1+\xi(x/\sigma(d)- \tilde{\mu} > 0) \} },
#' with the duration dependent scale parameter \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta > 0},
#' modified location parameter \eqn{\tilde{\mu}=\mu/\sigma(d)\in R}
#' and shape parameter \eqn{\xi\in R}, \eqn{\xi\neq 0}.
#' The parameters \eqn{\theta \leq 0} and \eqn{0<\eta<1} are duration offset and duration exponent
#' and describe the slope and curvature in the resulting IDF curves, respectively.
#' * The dependence of scale and location parameter on duration, \eqn{\sigma(d)} and \eqn{\mu(d)}, can be extended by multiscaling
#' and flattening, if requested. Multiscaling introduces a second duration exponent \eqn{\eta_2}, enabling the model to change slope
#' linearly with return period. Flattening adds a parameter \eqn{\tau}, that flattens the IDF curve for long durations:
#' \deqn{\sigma(x)=\sigma_0(d+\theta)^{-(\eta+\eta_2)}+\tau }
#' \deqn{\mu(x)=\tilde{\mu}(\sigma_0(d+\theta)^{-\eta_1}+\tau)}
#' * A useful introduction to __Maximum Likelihood Estimation__ for fitting for the 
#' generalized extreme value distribution (GEV) is provided by Coles (2001). It should be noted, however, that this method uses
#' the assumption that block maxima (of different durations or stations) are independent of each other. 
#' @references 
#' * Ulrich, J.; Jurado, O.E.; Peter, M.; Scheibel, M.; 
#' Rust, H.W. Estimating IDF Curves Consistently over Durations with Spatial Covariates. Water 2020, 12, 3119,
#' https://doi.org/10.3390/w12113119
#' * Demetris Koutsoyiannis, Demosthenes Kozonis, Alexandros Manetas,
#' A mathematical framework for studying rainfall intensity-duration-frequency relationships,
#' Journal of Hydrology,
#' Volume 206, Issues 1–2,1998,Pages 118-135,ISSN 0022-1694, https://doi.org/10.1016/S0022-1694(98)00097-3
#' * Coles, S.An Introduction to Statistical Modeling of Extreme Values;   Springer:  New York, NY, USA, 2001,
#' https://doi.org/10.1198/tech.2002.s73
#' @md
#' 
#' @examples 
#' ## Here are a few examples to illustrate the order in which the functions are intended to be used.
#' 
#' ## Step 0: sample 20 years of example hourly 'precipitation' data
# dates <- seq(as.POSIXct("2000-01-01 00:00:00"),as.POSIXct("2019-12-31 23:00:00"),by = 'hour')
# sample.precip <- rgamma(n = length(dates), shape = 0.05, rate = 0.4)
# precip.df <- data.frame(date=dates,RR=sample.precip)
# 
# ## Step 1: get annual maxima
# durations <- 2^(0:6) # accumulation durations [h] 
# ann.max <- IDF.agg(list(precip.df),ds=durations,na.accept = 0.1)
# # plotting the annual maxima in log-log representation
# plot(ann.max$ds,ann.max$xdat,log='xy',xlab = 'Duration [h]',ylab='Intensity [mm/h]')
# 
# ## Step 2: fit d-GEV to annual maxima
# fit <- gev.d.fit(xdat = ann.max$xdat,ds = ann.max$ds,sigma0link = make.link('log'))
# # checking the fit 
# gev.d.diag(fit,pch=1,legend = FALSE)
# # parameter estimates 
# params <- gev.d.params(fit)
# print(params)
# # plotting the probability density for a single duration 
# q.min <- floor(min(ann.max$xdat[ann.max$ds%in%1:2]))
# q.max <- ceiling(max(ann.max$xdat[ann.max$ds%in%1:2]))
# q <- seq(q.min,q.max,0.2)
# plot(range(q),c(0,0.55),type = 'n',xlab = 'Intensity [mm/h]',ylab = 'Density')
# for(d in 1:2){ # d=1h and d=2h
#   hist(ann.max$xdat[ann.max$ds==d],main = paste('d=',d),q.min:q.max
#        ,freq = FALSE,add=TRUE,border = d)   # sampled data
#   lines(q,dgev.d(q,params$mut,params$sigma0,params$xi,params$theta,params$eta,d = d),col=d) # etimated prob. density
# }
# legend('topright',col=1:2,lwd=1,legend = paste('d=',1:2,'h'),title = 'Duration')
# 
# ## Step 3: adding the IDF-curves to the data
# plot(ann.max$ds,ann.max$xdat,log='xy',xlab = 'Duration [h]',ylab='Intensity [mm/h]')
# IDF.plot(durations,params,add=TRUE)
NULL

#### IDF.agg ####

#' Aggregation and annual maxima for chosen durations
#' @description Aggregates several time series for chosen durations and finds annual maxima 
#' (either for the whole year or chosen months). Returns data.frame that can be used for
#' the function \code{\link{gev.d.fit}}.
#'
#' @param data list of data.frames containing time series for every station. 
#' The data.frame must have the columns 'date' and 'RR' unless other names 
#' are specified in the parameter `names`. The column 'date' must contain strings with 
#' standard date format.
#' @param ds numeric vector of aggregation durations in hours. 
#' (Must be multiples of time resolution at all stations.)
#' @param na.accept numeric in [0,1) giving maximum percentage of missing values 
#' for which block max. should still be calculated.
#' @param which.stations optional, subset of stations. Either numeric vector or character vector 
#' containing names of elements in data. If not given, all elements in `data` will be used.
#' @param which.mon optional, subset of months (as list containing values from 0 to 11) of which to calculate the annual maxima from. 
#' @param names optional, character vector of length 2, containing the names of the columns to be used. 
#' @param cl optional, number of cores to be used from \code{\link[parallel]{parLapply}} for parallel computing.
#'
#' @details If data contains stations with different time resolutions that need to be aggregated at
#' different durations, IDF.agg needs to be run separately for the different groups of stations. 
#' Afterwards the results can be joint together using `rbind`.
#'
#' @return data.frame containing the annual intensity maxima [mm/h] in `$xdat`, the corresponding duration in `$ds`,
#' the `$year` and month (`$mon`) in which the maxima occurred 
#' and the station id or name in `$station`.
#' 
#' @seealso \code{\link{pgev.d}}
#' 
#' @export
#' @importFrom parallel parLapply makeCluster stopCluster
#' @importFrom pbapply pblapply 
#' @importFrom RcppRoll roll_sum
#' @importFrom fastmatch ctapply
#'
#' @examples
#' dates <- as.Date("2019-01-01")+0:729
#' x <- rgamma(n = 730, shape = 0.4, rate = 0.5)
#' df <- data.frame(date=dates,RR=x)
#' 
#' # get annual maxima
#' IDF.agg(list('Sample'= df),ds=c(24,48),na.accept = 0.01)
#' 
#' ##      xdat ds year  mon station
#' ## 0.2853811 24 2019 0:11  Sample
#' ## 0.5673122 24 2020 0:11  Sample
#' ## 0.1598448 48 2019 0:11  Sample
#' ## 0.3112713 48 2020 0:11  Sample
#' 
#' # get monthly maxima for each month of june, july and august
#' IDF.agg(list('Sample'=df),ds=c(24,48),na.accept = 0.01,which.mon = list(5,6,7))
#' 
#' # get maxima for time range from june to august
#' IDF.agg(list('Sample'=df),ds=c(24,48),na.accept = 0.01,which.mon = list(5:7))
#' 
    IDF.agg <- function(data,ds,na.accept = 0,
                        which.stations = NULL,which.mon = list(0:11),names = c('date','RR'),cl = 1){
      
      if(!inherits(data, "list"))stop("Argument 'data' must be a list, instead it is a: ", class(data))
      
      # function 2: aggregate station data over durations and find annual maxima:                                
      agg.station <- function(station){
        data.s <- data[[station]]
        if(!is.data.frame(data.s)){stop("Elements of 'data' must be data.frames. But element "
                                        ,station," contains: ", class(data.s))}
        if(sum(is.element(names[1:2],names(data.s)))!=2){stop('Dataframe of station ', station 
                                                              ,' does not contain $', names[1]
                                                              ,' or $', names[2], '.')}
        dtime<-as.numeric((data.s[,names[1]][2]-data.s[,names[1]][1]),units="hours")
        
        if(any((ds/dtime)%%1 > 10e-8)){
          stop('At least one of the given aggregation durations is not multiple of the time resolution = '
               ,dtime,'hours at station ',station,'.')}
        
        # function 1: aggregate over single durations and find annual maxima:
        agg.ts <- function(ds){
          runsum = RcppRoll::roll_sum(data.s[,names[2]],round(ds/dtime),fill=NA,align='right')
          #runmean <- rollapplyr(as.zoo(data.s[,names[2]]),ds/dtime,FUN=sum,fill =NA,align='right')
          runsum <- runsum/ds #intensity per hour
          max.subset <- lapply(1:length(which.mon),function(m.i){
            subset <- is.element(as.POSIXlt(data.s[,names[1]])$mon,which.mon[[m.i]])
            max <- fastmatch::ctapply(runsum[subset],(as.POSIXlt(data.s[,names[1]])$year+1900)[subset],
                          function(vec){
                            n.na <- sum(is.na(vec))
                            max <- ifelse(n.na <= na.accept*length(vec),max(vec,na.rm = TRUE),NA)
                            return(max)})
            df <- data.frame(xdat=max,ds=ds,year=as.numeric(names(max)),mon=deparse(which.mon[[m.i]]),
                             station= station,
                             stringsAsFactors = FALSE)
            return(df)})
          df <- do.call(rbind,max.subset)  
          return(df) # maxima for single durations
        }
        # call function 1 in lapply to aggregate over all durations at single station
        if(cl>1){
          clust <- parallel::makeCluster(cl,type='FORK')
          data.agg <- parallel::parLapply(cl = clust,ds,agg.ts)  
          parallel::stopCluster(clust)
        }else{data.agg <-lapply(ds,agg.ts)}

        df <- do.call(rbind,data.agg)
        return(df) # maxima for all durations at one station
      }
      # which stations should be used?
      if(is.null(which.stations))which.stations <- if(is.null(names(data))){1:length(data)}else{names(data)}
      # call function 2 in lapply to aggregate over all durations at all stations
      station.list <- pbapply::pblapply(which.stations,agg.station)
      
      return(do.call('rbind',station.list))
    }

#### IDF.plot ####

#' Plotting of IDF curves at a chosen station
#'
#' @param durations vector of durations for which to calculate the quantiles. 
#' @param fitparams vector containing parameters mut, sigma0, xi, theta, eta
#' (modified location, scale offset, shape, duration offset, duration exponent) for chosen station
#' as obtained from \code{\link{gev.d.fit}}
#' (or \code{\link{gev.d.params}} for model with covariates).
#' @param probs vector of non-exceedance probabilities for which to plot IDF curves (p = 1-1/(Return Period))
#' @param cols vector of colors for IDF curves. Should have same length as \code{probs}
#' @param add logical indicating if plot should be added to existing plot, default is FALSE
#' @param legend logical indicating if legend should be plotted (TRUE, the default)
#' @param ... additional parameters passed on to the \code{plot} function
#'
#' @export
#' @importFrom grDevices rgb
#' @importFrom graphics axis box lines plot points 
#' @examples
#' data('example',package = 'IDF')
#' # fit d-gev
#' fit <- gev.d.fit(example$dat,example$d,ydat = as.matrix(example[,c("cov1","cov2")])
#'                  ,mutl = c(1,2),sigma0l = 1)
#' # get parameters for cov1 = 1, cov2 = 1
#' par <- gev.d.params(fit = fit, ydat = matrix(1,1,2))
#' # plot quantiles
#' IDF.plot(durations = seq(0.5,35,0.2),fitparams = par)
#' # add data points
#' points(example[example$cov1==1,]$d,example[example$cov1==1,]$dat)
IDF.plot <- function(durations,fitparams,probs=c(0.5,0.9,0.99),
                     cols=4:2,add=FALSE,
                     legend=TRUE,...){
  # if cols is too short, make longer    
  if(length(cols)!=length(probs))cols <- rep_len(cols,length.out=length(probs))
  ## calculate IDF values for given probability and durations
  qs <- lapply(durations,qgev.d,p=probs,mut=fitparams[1],sigma0=fitparams[2],xi=fitparams[3],
         theta=fitparams[4],eta=fitparams[5], eta2=fitparams[6], tau=fitparams[7])
  idf.array <- matrix(unlist(qs),length(probs),length(durations)) # array[probs,durs]
  if(!add){ #new plot
    ## initialize plot window
    # check if limits were passed
    if(is.element('ylim',names(list(...)))){
      ylim <- list(...)[['ylim']]
    }else{ylim <- range(idf.array,na.rm=TRUE)}
    if(is.element('xlim',names(list(...)))){
      xlim <- list(...)[['xlim']]
    }else{xlim <- range(durations,na.rm=TRUE)}
    if(is.element('main',names(list(...)))){
      main <- list(...)[['main']]
    }else{main <- ''}
    
    # empty plot
    plot(NA,xlim=xlim,ylim=ylim,xlab="Duration [h]",ylab="Intensity [mm/h]",log="xy",main=main)
  }
  ## plot IDF curves
  for(i in 1:length(probs)){
    lines(durations,idf.array[i,],col=cols[i],...)
  }
  
  if(legend){## plot legend
    # check if lwd, lty were passed
    if(is.element('lwd',names(list(...)))){
      lwd <- list(...)[['lwd']]
    }else{lwd <- 1}
    if(is.element('lty',names(list(...)))){
      lty <- list(...)[['lty']]
    }else{lty <- 1}
    
    legend(x="topright",title = 'p-quantile',legend=probs,
           col=cols,lty=lty,lwd=lwd)
  }
}

#################################################################
## gev.d.fit function code
########################################################################
# This file contains the functions:
# - gev.d.fit, gev.d.init for fitting
# - gev.d.diag for diagnostic plots
# - gev.d.params for calculation of parameters
# and the documentation of the example data

#### gev.d.fit ####

#' @title Maximum-likelihood Fitting of the duration-dependent GEV Distribution
#' @description Modified \code{\link[ismev]{gev.fit}} function for Maximum-likelihood fitting
#' for the duration-dependent generalized extreme
#' value distribution, following Koutsoyiannis et al. (1998), including generalized linear
#' modeling of each parameter.
#' @param xdat A vector containing maxima for different durations.
#' This can be obtained from \code{\link{IDF.agg}}.
#' @param ds A vector of aggregation levels corresponding to the maxima in xdat.
#' 1/60 corresponds to 1 minute, 1 corresponds to 1 hour.
#' @param ydat A matrix of covariates for generalized linear modeling of the parameters
#' (or NULL (the default) for stationary fitting). The number of rows should be the same as the
#' length of xdat.
#' @param  mutl,sigma0l,xil,thetal,etal,taul,eta2l Numeric vectors of integers, giving the columns of ydat that contain
#'  covariates for generalized linear modeling of the parameters (or NULL (the default)
#'  if the corresponding parameter is stationary).
#'  Parameters are: modified location, scale offset, shape, duration offset, duration exponent, respectively.
#' @param mutlink,sigma0link,xilink,thetalink,etalink,eta2link,taulink Link functions for generalized linear
#' modeling of the parameters, created with \code{\link{make.link}}. The default is \code{make.link("identity")}.
#' @param init.vals list, giving initial values for all or some parameters
#' (order: mut, sigma0, xi, theta, eta, eta2, tau). If one of those parameters shall not be used (see theta_zero, eta2_zero, tau_zero),
#' the number of parameters has to be reduced accordingly. If some or all given values in init.vals are NA or 
#' no init.vals at all is declared (the default), initial parameters are obtained
#' internally by fitting the GEV separately for each duration and applying a linear model to obtain the
#' duration dependency of the location and shape parameter.
#' Initial values for covariate parameters are assumed as 0 if not given.
#' @param theta_zero Logical value, indicating whether theta should be estimated (FALSE, the default) or
#' should stay zero.
#' @param tau_zero,eta2_zero Logical values, indicating whether tau,eta2 should be estimated (TRUE, the default).
#' @param show Logical; if TRUE (the default), print details of the fit.
#' @param method The optimization method used in \code{\link{optim}}.
#' @param maxit The maximum number of iterations.
#' @param ... Other control parameters for the optimization.
#' @return A list containing the following components.
#' A subset of these components are printed after the fit.
#' If \code{show} is TRUE, then assuming that successful convergence is indicated,
#' the components nllh, mle and se are always printed.
#' \item{nllh}{single numeric giving the negative log-likelihood value}
#' \item{mle}{numeric vector giving the MLE's for the modified location, scale_0, shape,
#' duration offset and duration exponent, resp. If requested, contains also second duration exponent and intensity-offset}
#' \item{se}{numeric vector giving the standard errors for the MLE's (in the same order)}
#' \item{trans}{A logical indicator for a non-stationary fit.}
#' \item{model}{A list with components mutl, sigma0l, xil, thetal and etal. If requested, contains also eta2l and taul}
#' \item{link}{A character vector giving link functions.}
#' \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}.
#' A zero indicates successful convergence.}
#' \item{data}{data is standardized to standard Gumbel.}
#' \item{cov}{The covariance matrix.}
#' \item{vals}{Parameter values for every data point.}
#' \item{init.vals}{Initial values that were used.}
#' \item{ds}{Durations for every data point.}
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' @seealso \code{\link{IDF-package}}, \code{\link{IDF.agg}}, \code{\link{gev.fit}}, \code{\link{optim}}
#' @export
#' @importFrom stats optim
#' @importFrom stats make.link
#'
#' @examples
#' # sampled random data from d-gev with covariates
#' # GEV parameters:
#' # mut = 4 + 0.2*cov1 +0.5*cov2
#' # sigma0 = 2+0.5*cov1
#' # xi = 0.5
#' # theta = 0
#' # eta = 0.5
#' # eta2 = 0
#' # tau = 0
#'
#' data('example',package ='IDF')
#'
#' gev.d.fit(xdat=example$dat,ds = example$d,ydat=as.matrix(example[,c('cov1','cov2')])
#' ,mutl=c(1,2),sigma0l=1)

gev.d.fit <- function(xdat, ds, ydat = NULL, mutl = NULL, sigma0l = NULL, xil = NULL, thetal = NULL, etal = NULL, taul = NULL, eta2l = NULL,
                      mutlink = make.link("identity"), sigma0link = make.link("identity"), xilink = make.link("identity"),
                      thetalink = make.link("identity"), etalink = make.link("identity"), taulink = make.link("identity"), eta2link = make.link("identity"),
                      init.vals = NULL, theta_zero = FALSE, tau_zero=TRUE, eta2_zero=TRUE,
                      show = TRUE, method = "Nelder-Mead", maxit = 10000,...)
{
  if (length(xdat) != length(ds)) {
    stop(paste0('The length of xdat is ',length(xdat),', but the length of ds is ',length(ds),'.'))
  }
  z <- list()
  # number of parameters (betas) to estimate for each parameter:
  npmu <- length(mutl) + 1
  npsc <- length(sigma0l) + 1
  npsh <- length(xil) + 1
  npth <- ifelse(!theta_zero,length(thetal) + 1,0)
  npet <- length(etal) + 1
  npta <- ifelse(!tau_zero,  length(taul)   + 1,0)
  npe2 <- ifelse(!eta2_zero,  length(eta2l)   + 1,0)
  z$trans <- FALSE  # indicates if fit is non-stationary
  z$model <- list(mutl, sigma0l, xil, thetal, etal, eta2l, taul)
  z$link <- list(mutlink=mutlink, sigma0link=sigma0link, xilink=xilink, thetalink=thetalink, etalink=etalink, eta2link=eta2link, taulink=taulink)
  
  # test for NA values:
  if(any(is.na(xdat))) stop('xdat contains NA values. NA values need to be removed first.')
  # test for finite values:
  if(any(is.infinite(xdat))) stop('xdat contains non finite values. Inf and -Inf need to be removed first.')
  
  # test if covariates matrix is given correctly
  npar <- max(sapply(z$model,function(x){return(ifelse(is.null(x),0,max(x)))}))
  if(any(npar>ncol(ydat),npar>0 & is.null(ydat)))stop("Not enough columns in covariates matrix 'ydat'.")
  
  # initial values
  init.necessary.length = 4 + as.numeric(!theta_zero) + as.numeric(!eta2_zero) + as.numeric(!tau_zero)  # mut, sigma0, xi, theta, eta, eta2, tau
  if (is.null(init.vals)) {init.vals = as.list(rep(NA,init.necessary.length))
  }else{init.vals = as.list(init.vals)}
  
  if(length(init.vals)!=init.necessary.length | !is.list(init.vals)) {
    print(paste0('Parameter init.vals is not used, because it is no list of length ',init.necessary.length,'.'))
    init.vals <- gev.d.init(xdat,ds,z$link)
    
  }else{ # length of given values is correct
    
    # name given initial values
    names1=c('mu','sigma','xi')                 # standard set of parameters
    if (!theta_zero){names1=c(names1,'theta')}  # add theta (in case)
    names1=c(names1,'eta')                      # add eta   (always)
    if (!eta2_zero){names1=c(names1,'eta2')}    # add eta2  (in case)
    if (!tau_zero){names1=c(names1,'tau')}      # add tau   (in case)
    names(init.vals) <- names1
    
    # add missing initial value (fixed internal number of parameters: 7)
    if (theta_zero) init.vals$theta = 0
    if (eta2_zero) init.vals$eta2 = 0#init.vals$eta
    if (tau_zero) init.vals$tau = 0
    init.vals = init.vals[c("mu","sigma","xi","theta","eta","eta2","tau")]
    iv = init.vals
    init.vals = list(mu=iv$mu,sigma=iv$sigma,xi=iv$xi,theta=iv$theta,eta=iv$eta,eta2=iv$eta2,tau=iv$tau)
    
    if(!any(is.na(init.vals))){ #all initial values are given
      # do nothing
    }else if(any(!is.na(init.vals))) { #some initial values are given
      #if (!eta2_zero) print("autmoatic inital value setting not implemented yet for multiscaling (eta2_zero=FALSE)")
      init.vals.user <- init.vals
      init.vals <- gev.d.init(xdat,ds,z$link) #calculate init.vals using gev.d.init
      for (i in 1:length(init.vals)){ #overwrite the calculated initial values with the ones given by the user
        if(!is.na(init.vals.user[[names(init.vals.user)[i]]])) {
          init.vals[[names(init.vals.user)[i]]]<-init.vals.user[[names(init.vals.user)[i]]]
        } 
      }
    }else{ #no initial values are given
      #if (!eta2_zero) print("autmoatic inital value setting not implemented yet for multiscaling (eta2_zero=FALSE)")
      init.vals <- gev.d.init(xdat,ds,z$link)
    }
  } 
  
  # transform eta2 to eta2~~eta oriented. the optim function does not need eta2~~0, but eta2~~eta
  init.vals$eta2=init.vals$eta + init.vals$eta2
  
  # generate covariates matrices:
  if (is.null(mutl)) { #stationary
    mumat <- as.matrix(rep(1, length(xdat)))
    muinit <- init.vals$mu
  }else { #non stationary
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mutl])
    muinit <- c(init.vals$mu, rep(0, length(mutl)))[1:npmu] #fill with 0s to length npmu
  }
  if (is.null(sigma0l)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    siginit <- init.vals$sigma
  }else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigma0l])
    siginit <- c(init.vals$sigma, rep(0, length(sigma0l)))[1:npsc]
  }
  if (is.null(xil)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    shinit <- init.vals$xi
  }else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, xil])
    shinit <- c(init.vals$xi, rep(0, length(xil)))[1:npsh]
  }
  if (is.null(thetal)) {
    thmat <- as.matrix(rep(1, length(xdat)))
    thetainit <- init.vals$theta
  }else {
    z$trans <- TRUE
    thmat <- cbind(rep(1, length(xdat)), ydat[, thetal])
    thetainit <- c(init.vals$theta, rep(0, length(thetal)))[1:npth]
  }
  if (is.null(etal)) {
    etmat <- as.matrix(rep(1, length(xdat)))
    etainit <- init.vals$eta
  }else {
    z$trans <- TRUE
    etmat <- cbind(rep(1, length(xdat)), ydat[, etal])
    etainit <- c(init.vals$eta, rep(0, length(etal)))[1:npet]
  }
  if (is.null(taul)) {
    tamat <- as.matrix(rep(1, length(xdat)))
    tauinit <- init.vals$tau
  }else {
    z$trans <- TRUE
    tamat <- cbind(rep(1, length(xdat)), ydat[, taul])
    tauinit <- c(init.vals$tau, rep(0, length(taul)))[1:npta]
  }
  if (is.null(eta2l)) {
    e2mat <- as.matrix(rep(1, length(xdat)))
    eta2init <- init.vals$eta2
  }else {
    z$trans <- TRUE
    e2mat <- cbind(rep(1, length(xdat)), ydat[, eta2l])
    eta2init <- c(init.vals$eta2, rep(0, length(eta2l)))[1:npe2]
  }
  
  init <- c(muinit,siginit,shinit)
  if (!theta_zero) {init <- c(init,thetainit)} # add theta init (in case)
  init <- c(init,etainit)                      # add eta init   (always)
  if (!eta2_zero)  {init <- c(init,eta2init)}  # add eta2 init  (in case)
  if (!tau_zero)   {init <- c(init,tauinit)}   # add tau init   (in case)
  
  # functions to calculate neg log-likelihood and gradient:
  gev.lik <- function(a) {
    # computes neg log lik of d-gev model
    mu <- mutlink$linkinv(mumat %*% (a[1:npmu]))
    sigma <- sigma0link$linkinv(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- xilink$linkinv(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    theta <-if(!theta_zero){thetalink$linkinv(thmat %*% (a[seq(npmu + npsc + npsh + 1, length = npth)]))}else{0}
    eta <- etalink$linkinv(etmat %*% (a[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
    eta2  <-if(!eta2_zero){eta2link$linkinv( e2mat %*% (a[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))}else{eta}
    tau   <-if(!tau_zero){taulink$linkinv(  tamat %*% (a[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))}else{0}
    
    ds.t <- ds+theta
    sigma.d <-     sigma/(ds.t^eta2)+tau
    mu.d    <- mu*(sigma/(ds.t^eta)+tau)
    
    y = (xdat - mu.d) / sigma.d
    y <- 1 + xi * y
    
    
    if(!theta_zero) {if(any(theta < 0)) {return(10^6)} } # check definition condition for theta
    if(any(eta <= 0) || any(sigma.d <= 0) || any(y <= 0)) return(10^6)
    if(!tau_zero)   {if(any(tau < 0))    {return(10^6)} } # check definition condition for tau.
    if(!eta2_zero) {if(any(eta2 < 0))    {return(10^6)} } # check definition condition for eta2.
    
    return(sum(log(sigma.d)) + sum(y^(-1/xi)) +  sum(log(y) * (1/xi + 1)))
  }
  
  gev.grad <-  function(a){
    # computes gradient of neg log lik of d-gev model
    mu <- mutlink$linkinv(mumat %*% (a[1:npmu]))
    sigma <- sigma0link$linkinv(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- xilink$linkinv(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    theta <-if(!theta_zero){thetalink$linkinv(thmat %*% (a[seq(npmu + npsc + npsh + 1, length = npth)]))}else{0}
    eta <- etalink$linkinv(etmat %*% (a[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
    eta2  <-if(!eta2_zero){eta2link$linkinv( e2mat %*% (a[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))}else{eta}
    tau   <-if(!tau_zero){taulink$linkinv(  tamat %*% (a[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))}else{0}
    
    sigma.d <-     sigma/((ds + theta)^eta2) + tau
    mu.d    <- mu * (sigma/((ds + theta)^eta) + tau)
    y <-  1 + xi * (xdat - (mu.d))/(sigma.d)
    
    
    dnll.mu.out <- -(xi * (sigma/((ds + theta)^eta) + tau)/(sigma.d)/(y) * (1/xi + 1) + (y)^((-1/xi) - 1) * ((-1/xi) * (xi * (sigma/((ds + theta)^eta) + tau)/(sigma.d))))
    dnll.mu <- apply(mumat,2,function(x)sum(dnll.mu.out*mutlink$mu.eta(mumat %*% (a[1:npmu]))*x))
    dnll.sigma.out <- 1/((ds + theta)^eta2)/(sigma.d) - (y)^((-1/xi) - 1) * 
      ((-1/xi) * (xi * (mu * (1/((ds + theta)^eta)))/(sigma.d) + xi * (xdat - (mu.d)) * (1/((ds + theta)^eta2))/(sigma.d)^2)) - 
      (xi * (mu * (1/((ds + theta)^eta)))/(sigma.d) + xi * (xdat - (mu.d)) * (1/((ds + theta)^eta2))/(sigma.d)^2)/(y) * (1/xi + 1)
    dnll.sigma <- apply(sigmat,2,function(x)sum(dnll.sigma.out*sigma0link$mu.eta(sigmat %*% (a[seq(npmu + 1, length = npsc)]))*x))
    dnll.xi.out <- (y)^((-1/xi) - 1) * ((-1/xi) * ((xdat - (mu.d))/(sigma.d))) + (y)^(-1/xi) * (log((y)) * (1/xi^2)) + 
      ((xdat - (mu.d))/(sigma.d)/(y) * (1/xi + 1) - log(y) * (1/xi^2))
    dnll.xi <- apply(shmat,2,function(x)sum(dnll.xi.out*xilink$mu.eta(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))*x))
    if(!theta_zero){
      dnll.theta.out <- (y)^((-1/xi) - 1) * ((-1/xi) * (xi * (mu * (sigma * ((ds + theta)^(eta - 1) * eta)/((ds + theta)^eta)^2))/
                                                          (sigma.d) + xi * (xdat - (mu.d)) * (sigma * ((ds + theta)^(eta2 - 1) * eta2)/((ds + theta)^eta2)^2)/(sigma.d)^2)) - 
        sigma * ((ds + theta)^(eta2 - 1) * eta2)/((ds + theta)^eta2)^2/(sigma.d) + 
        (xi * (mu * (sigma * ((ds + theta)^(eta - 1) * eta)/((ds + theta)^eta)^2))/(sigma.d) + xi * (xdat - (mu.d)) * 
           (sigma * ((ds + theta)^(eta2 - 1) * eta2)/((ds + theta)^eta2)^2)/(sigma.d)^2)/(y) * (1/xi + 1)
      dnll.theta <-  apply(thmat,2,function(x)sum(dnll.theta.out*thetalink$mu.eta(thmat %*% (a[seq(npmu + npsc + npsh + 1, length = npth)]))*x))
    }
    if(eta2_zero){
      dnll.eta.out <- (y)^((-1/xi) - 1) * ((-1/xi) * (xi * (mu * (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2))/(sigma.d) + xi * (xdat - (mu.d)) * 
                                                        (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2)/(sigma.d)^2)) - sigma * 
        ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2/(sigma.d) + 
        (xi * (mu * (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2))/(sigma.d) + xi * (xdat - (mu.d)) * 
           (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2)/(sigma.d)^2)/(y) * (1/xi + 1)
      dnll.eta <- apply(etmat,2,function(x)sum(dnll.eta.out*etalink$mu.eta(etmat %*% (a[seq(npmu + npsc + npsh + npth + 1, length = npet)]))*x))
    }else{
      dnll.eta.out <- (y)^((-1/xi) - 1) * ((-1/xi) * (xi * (mu * (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2))/(sigma.d))) + xi * 
        (mu * (sigma * ((ds + theta)^eta * log((ds + theta)))/((ds + theta)^eta)^2))/(sigma.d)/(y) * (1/xi + 1)
      dnll.eta <- apply(etmat,2,function(x)sum(dnll.eta.out*etalink$mu.eta(etmat %*% (a[seq(npmu + npsc + npsh + npth + 1, length = npet)]))*x))
      dnll.eta2.out <- (y)^((-1/xi) - 1) * ((-1/xi) * (xi * (xdat - (mu.d)) * (sigma * ((ds + theta)^eta2 * log((ds + theta)))/((ds + theta)^eta2)^2)/(sigma.d)^2)) - 
        sigma * ((ds + theta)^eta2 * log((ds + theta)))/((ds + theta)^eta2)^2/(sigma.d) + xi * (xdat - (mu.d)) * 
        (sigma * ((ds + theta)^eta2 * log((ds + theta)))/((ds + theta)^eta2)^2)/(sigma.d)^2/(y) * (1/xi + 1)
      dnll.eta2 <- apply(e2mat,2,function(x)sum(dnll.eta2.out*eta2link$mu.eta( e2mat %*% (a[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))*x))
    }
    if(!tau_zero){
      dnll.tau.out <- 1/(sigma.d) - (y)^((-1/xi) - 1) * ((-1/xi) * (xi * mu/(sigma.d) + xi * (xdat - (mu.d))/(sigma.d)^2)) - 
        (xi * mu/(sigma.d) + xi * (xdat - (mu.d))/(sigma.d)^2)/(y) * (1/xi + 1)
      dnll.tau <- apply(tamat,2,function(x)sum(dnll.tau.out*taulink$mu.eta(tamat %*% (a[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))*x))
    }
    
    grad.nll <- c(dnll.mu,dnll.sigma,dnll.xi,if(!theta_zero){dnll.theta},dnll.eta,if(!eta2_zero){dnll.eta2},if(!tau_zero){dnll.tau})
    
    if(any(theta < 0)) {return(rep(0,length(grad.nll)))}  # check definition condition for theta
    if(any(eta <= 0) || any(sigma.d <= 0) || any(y <= 0)) return(rep(0,length(grad.nll)))
    if(any(tau < 0))    {return(rep(0,length(grad.nll)))}  # check definition condition for tau.
    if(any(eta2 <= 0))    {return(rep(0,length(grad.nll)))}  # check definition condition for eta2.
    
    return(grad.nll)
  }
  
  # finding minimum of log-likelihood:
  x <- optim(init, gev.lik, hessian = TRUE, method = method, gr = gev.grad,
             control = list(maxit = maxit, ...))
  
  # saving output parameters:
  z$conv <- x$convergence
  mut <- mutlink$linkinv(mumat %*% (x$par[1:npmu]))
  sc0 <- sigma0link$linkinv(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- xilink$linkinv(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  if(!theta_zero){ #When user does NOT set theta parameter to zero (default)
    theta <- thetalink$linkinv(thmat %*% (x$par[seq(npmu + npsc + npsh + 1, length = npth)]))
  }else{ #When user requests theta_parameter to be zero
    theta <- thetalink$linkinv(thmat %*% (0))
  }
  eta <- etalink$linkinv(etmat %*% (x$par[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
  
  if(!eta2_zero){  #When user wants to use eta2 parameter
    eta2 <- eta2link$linkinv(e2mat %*% (x$par[seq(npmu + npsc + npsh + npth + npet + 1,length = npe2)]))
    #transform eta2 to eta~~eta2 oriented
    eta2 <- eta2 - eta
  }else{ #When user requests not to have eta2 parameter (default)
    eta2 <- 0#eta
  }
  
  if(!tau_zero){  #When user does NOT set tau parameter to zero (not default)
    tau <- taulink$linkinv(tamat %*% (x$par[seq(npmu + npsc + npsh + npth + npet + npe2 + 1,length = npta)]))
  }else{ #When user requests tau parameter to be zero
    tau <- taulink$linkinv(tamat %*% (0))
  }
  
  z$nllh <- x$value
  # normalize data to standard Gumbel:
  sc.d  <-      sc0/((ds+theta)^(eta+eta2))+tau  # new
  mut.d <- mut*(sc0/((ds+theta)^eta )+tau) # new
  
  #z$data <-  - log(as.vector((1 + xi * (xdat/sc.d-mut))^(-1/xi)))  # old
  z$data <-  - log(as.vector((1 + xi * ((xdat-mut.d)/sc.d))^(-1/xi))) # new
  z$mle <- x$par
  
  #z$mle$eta2 = z$mle$eta2-z$mle$eta 
  # do not transform eta2 in mle, because the original estimated parameters should be here. 
  # for the transformed parameters, use gev.d.params
  
  test <- try(              # catch error
    z$cov <- solve(x$hessian) # invert hessian to get estimation on var-covar-matrix
    ,silent = TRUE )
  if("try-error" %in% class(test)){
    warning("Hessian could not be inverted. NAs were produced.")
    z$cov <- matrix(NA,length(z$mle),length(z$mle))
  }
  z$se <- sqrt(diag(z$cov)) # sqrt(digonal entries) = standart error of mle's
  z$vals <- cbind(mut, sc0, xi, theta, eta, eta2, tau)
  colnames(z$vals) <- c("mut","sigma0","xi","theta","eta","eta2","tau")
  z$init.vals <- unlist(init.vals)
  # transform eta2 to zero-oriented
  z$init.vals[6] = z$init.vals[6] + z$init.vals[4]
  z$ds <- ds
  z$theta_zero <- theta_zero #Indicates if theta parameter was set to zero by user.
  z$tau_zero <- tau_zero     #Indicates if tau parameter was set to zero by user. (default)
  z$eta2_zero <- eta2_zero   #Indicates if eta2 parameter was set to zero by user. (default)
  if(show) {
    if(z$trans) { # for nonstationary fit
      print(z[c(2, 4)]) # print model, link (3) , conv
      # print names of link functions:
      cat('$link\n')
      #if(!tau_zero){
      print(c(z$link$mutlink$name,z$link$sigma0link$name,z$link$xilink$name,z$link$thetalink$name,z$link$etalink$name,z$link$eta2link$name,z$link$taulink$name))
      #} else{
      #  print(c(z$link$mutlink$name,z$link$sigma0link$name,z$link$xilink$name,z$link$thetalink$name,z$link$etalink$name))
      #}
      cat('\n')
    }else{print(z[4])} # for stationary fit print only conv
    if(!z$conv){ # if fit converged
      print(z[c(5, 7, 9)]) # print nll, mle, se
    }
  }
  class(z) <- "gev.d.fit"
  invisible(z)
}


#### gev.d.init ####

# function to get initial values for gev.d.fit:
# obtain initial values
# by fitting every duration separately

# possible ways to improve:
# take given initial values into account, if there are any
# xi -> mean vs. median ... how do we improve that?
# mu_tilde -> is not very good for small sample sizes yet
# improved initial value for eta, by fitting both mu~d and sigma~d in log-log scale

#' @title get initial values for gev.d.fit
#' @description obtain initial values by fitting every duration separately
#' @param xdat vector of maxima for different durations
#' @param ds vector of durations belonging to maxima in xdat
#' @param link list of 5, link functions for parameters, created with \code{\link{make.link}}
#' @return list of initial values for mu_tilde, sigma_0, xi, theta, eta, eta2, tau
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom ismev gev.fit
#' @keywords internal

gev.d.init <- function(xdat,ds,link){
  durs <- unique(ds)
  mles <- matrix(NA, nrow=length(durs), ncol= 3)
  for(i in 1:length(durs)){
    test <- try(fit <- ismev::gev.fit(xdat[ds==durs[i]],show = FALSE),silent = TRUE)
    if("try-error" %in% class(test) | fit$conv!=0){mles[i,] <- rep(NA,3)}else{mles[i,] <- fit$mle}
  }
  if(all(is.na(mles))){stop('Initial values could not be computed for this dataset.')}
  # get values for sig0 and eta (also mu_0) from linear model in log-log scale
  lmsig <- lm(log(mles[,2])~log(durs))
  lmmu <- lm(log(mles[,1])~log(durs))
  
  # sig0 <- exp Intercept
  siginit <- link$sigma0link$linkfun(exp(lmsig$coefficients[[1]]))
  # eta <- mean of negativ slopes
  etainit <- link$etalink$linkfun(mean(c(-lmsig$coefficients[[2]],-lmmu$coefficients[[2]])))
  eta2init <- 0.0 #etainit + 0.0
  # mean of mu_d/sig_d
  # could try:
  # mu0/sig0 = exp(lmmu$coefficients[[1]])/exp(lmsig$coefficients[[1]])
  muinit <- link$mutlink$linkfun(median(c(mles[,1]/mles[,2]),na.rm = TRUE))
  # mean of shape parameters
  shinit <- link$xilink$linkfun(median(mles[,3],na.rm = TRUE))
  thetainit <- link$thetalink$linkfun(0)
  tauinit <- link$taulink$linkfun(0)
  
  return(list(mu=muinit,sigma=siginit,xi=shinit,theta=thetainit,eta=etainit, eta2=eta2init, tau=tauinit))
}

#### gev.d.lik ####

#' d-GEV Likelihood
#'
#' Computes (log-) likelihood of d-GEV model
#' @param xdat numeric vector containing observations
#' @param ds numeric vector containing corresponding durations (1/60 corresponds to 1 minute, 1 corresponds to 1 hour)
#' @param mut,sigma0,xi,theta,eta,eta2,tau numeric vectors containing corresponding estimates for each of the parameters
#' @param log Logical; if TRUE, the log likelihood is returned.
#'
#' @return single value containing (log) likelihood
#' @export
#'
#' @examples
#' # compute log-likelihood of observation values not included in fit
#' train.set <- example[example$d!=2,]
#' test.set <- example[example$d==2,]
#' fit <- gev.d.fit(train.set$dat,train.set$d,mutl = c(1,2),sigma0l = 1
#'           ,ydat = as.matrix(train.set[c('cov1','cov2')]))
#' params <- gev.d.params(fit,ydat = as.matrix(test.set[c('cov1','cov2')]))
#' gev.d.lik(xdat = test.set$dat,ds = test.set$d,mut = params[,1],sigma0 = params[,2],xi = params[,3]
#'           ,theta = params[,4],eta = params[,5],log=TRUE)
gev.d.lik <- function(xdat,ds,mut,sigma0,xi,theta,eta,log=FALSE,tau=0,eta2=0) {
  eta2 = eta+eta2
  if(any(xi==0)){stop('Function is not defined for shape parameter of zero.')}
  if(any(! c(length(ds),length(mut),length(sigma0),length(xi),length(theta),length(eta),length(eta2),length(tau)) %in%
         c(1,length(xdat)))){
    stop('Input vectors differ in length, but must have the same length.')
  }
  
  ds.t <- ds+theta
  sigma.d <-     sigma0/(ds.t^eta2)+tau
  mu.d    <- mut*(sigma0/(ds.t^eta)+tau)
  y = (xdat - mu.d) / sigma.d # new
  y <- 1 + xi * y
  
  #sigma.d <- sigma0/(ds.t^eta) + tau # old
  #y <- xdat/sigma.d - mut            # old
  #y <- 1 + xi * y                    # old
  
  if(log){
    return(-sum(log(sigma.d) + y^(-1/xi) + log(y) * (1/xi + 1)))
  }else{
    return(-prod(sigma.d * exp(y^(-1/xi)) * y ^ (1/xi + 1)))
  }
  
}

#### gev.d.diag ####

#' Diagnostic Plots for d-gev Models
#'
#' @description  Produces diagnostic plots for d-gev models using
#' the output of the function \code{\link{gev.d.fit}}. Values for different durations can be plotted in
#' different colors of with different symbols.
#' @param fit object returned by \code{\link{gev.d.fit}}
#' @param subset an optional vector specifying a subset of observations to be used in the plot
#' @param cols optional either one value or vector of same length as \code{unique(fit$ds)} to
#' specify the colors of plotting points.
#' The default uses the \code{rainbow} function.
#' @param pch optional either one value or vector of same length as \code{unique(fit$ds)} containing
#' integers or symbols to specify the plotting points.
#' @param which string containing 'both', 'pp' or 'qq' to specify, which plots should be produced.
#' @param mfrow vector specifying layout of plots. If both plots should be produced separately,
#' set to \code{c(1,1)}.
#' @param legend logical indicating if legends should be plotted
#' @param title character vector of length 2, giving the titles for the pp- and the qq-plot
#' @param emp.lab,mod.lab character string containing names for empirical and model axis
#' @param ci logical indicating whether 0.95 confidence intervals should be plotted
#' @param ... additional parameters passed on to the plotting function
#'
#' @export
#' @importFrom graphics plot abline par title
#' @importFrom grDevices rainbow
#' @importFrom evd rgev
#'
#' @examples
#' data('example',package ='IDF')
#'
#' fit <- gev.d.fit(xdat=example$dat,ds = example$d,ydat=as.matrix(example[,c('cov1','cov2')])
#'                  ,mutl=c(1,2),sigma0l=1)
#' # diagnostic plots for complete data
#' gev.d.diag(fit,pch=1,ci = TRUE)
#' # diagnostic plots for subset of data (e.g. one station)
#' gev.d.diag(fit,subset = example$cov1==1,pch=1,ci = TRUE)
gev.d.diag <- function(fit,subset=NULL,cols=NULL,pch=NULL,which='both',mfrow=c(1,2),legend=TRUE,
                       title=c('Residual Probability Plot','Residual Quantile Plot'),
                       emp.lab='Empirical',mod.lab='Model',ci=FALSE,...){
  # check parameter:
  if(!is.element(which,c('both','pp','qq'))) stop("Parameter 'which'= ",which,
                                                  " but only 'both','pp' or 'qq' are allowed.")
  # subset data
  df <- data.frame(data=fit$data,ds=fit$ds)
  if(!is.null(subset)){
    if(dim(df)[1]!=length(subset)){stop("Length of 'subset' does not match length of data
                                        'xdat' used for fitting.")}
    df <- df[subset,]
  }
  # get single durations
  durs <- sort(unique(df$ds))
  # rescale durations to assign colors
  df$cval <- sapply(df$ds,function(d){which(durs==d)})
  
  # sort data
  df <- df[order(df$data),]
  
  # plotting position
  n <- length(df$data)
  px <- (1:n)/(n + 1)
  
  # get 95% confidence intervals
  if(ci){
    samp <- rgev(n * 99, loc = 0, scale = 1, shape = 0)
    samp <- matrix(samp, n, 99)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci.vals <- t(samp[c(3, 97), ])
  }else{ci.vals <- NULL}
  # create plots:
  if(which=='both') par(mfrow=mfrow) # 2 subplots
  # colors and symbols
  if(is.null(cols))cols <- rainbow(length(durs))
  if(length(cols)<length(durs)) cols <- rep(cols, length.out=length(durs)) 
  if(is.null(pch))pch <- 0:(length(durs)-1)#df$cval
  if(length(pch)<length(durs)) pch <- rep(pch, length.out=length(durs))
  
  if(which=='both'|which=='pp'){
    # pp
    plot(px, exp( - exp( - df$data)), xlab = emp.lab, ylab = mod.lab,col=cols[df$cval]
         ,pch=pch[df$cval],ylim=range(exp( - exp(-c(ci.vals,df$data))),na.rm = TRUE),...)
    abline(0, 1, col = 1,lwd=1)
    if(ci){
      lines(px,exp( - exp( - ci.vals[,1])),lty=2)
      lines(px,exp( - exp( - ci.vals[,2])),lty=2)
    }
    title(title[1])
    if(legend){legend('bottomright',legend = round(durs,digits = 2),pch=pch[1:length(durs)],
                      col = cols[1:length(durs)],title = 'Duration [h]',ncol = 2)}
  }
  if(which=='both'|which=='qq'){
    # qq
    plot( - log( - log(px)), df$data, ylab = emp.lab, xlab = mod.lab,col=cols[df$cval]
          ,pch=pch[df$cval],ylim=range(c(ci.vals,df$data),na.rm = TRUE),...)
    abline(0, 1, col = 1,lwd=1)
    if(ci){
      lines(-log(-log(px)),ci.vals[,1],lty=2)
      lines(-log(-log(px)),ci.vals[,2],lty=2)
    }
    title(title[2])
    if(legend){legend('bottomright',legend = round(durs,digits = 2),pch=pch[1:length(durs)],
                      col = cols[1:length(durs)],title = 'Duration [h]',ncol = 2)}
  }
  if(which=='both') par(mfrow=c(1,1)) # reset par
}


#### gev.d.params ####

#' Calculate gev(d) parameters from \code{gev.d.fit} output
#'
#' @description function to calculate mut, sigma0, xi, theta, eta, eta2, tau
#' (modified location, scale offset, shape, duration offset, duration exponent, second duration exponent, intensity offset)
#' from results of \code{\link{gev.d.fit}} with covariates or link functions other than identity.
#' @param fit fit object returned by \code{\link{gev.d.fit}} or \code{\link{gev.fit}}
#' @param ydat A matrix containing the covariates in the same order as used in \code{gev.d.fit}.
#' @seealso \code{\link{IDF-package}}
#' @return data.frame containing mu_tilde, sigma0, xi, theta, eta, eta2, tau (or mu, sigma, xi for gev.fit objects)
#' @export
#'
#' @examples
#' data('example',package = 'IDF')
#' fit <- gev.d.fit(example$dat,example$d,ydat = as.matrix(example[,c("cov1","cov2")])
#'                  ,mutl = c(1,2),sigma0l = 1)
#' gev.d.params(fit = fit,ydat = cbind(c(0.9,1),c(0.5,1)))


gev.d.params <- function(fit,ydat=NULL){
  if(!inherits(fit,c("gev.fit","gev.d.fit")))stop("'fit' must be an object returned by 'gev.d.fit' or 'gev.fit'.")
  if(!is.null(ydat)){
    # check covariates matrix
    if(!is.matrix(ydat))stop("'ydat' must be of class matrix.")
    n.par <- max(sapply(fit$model,function(x){return(ifelse(is.null(x),0,max(x)))}))
    if(n.par>ncol(ydat))stop("Covariates-Matrix 'ydat' has ",ncol(ydat), " columns, but ", n.par," are required.")
  }else{if(!fit$trans){# no model -> no covariates matrix
    ydat <- matrix(1)
  }else{stop("To calculate parameter estimates, covariates matrix 'ydat' must be provided.")}
  }
  
  # number of parameters
  npmu <- length(fit$model[[1]]) + 1
  npsc <- length(fit$model[[2]]) + 1
  npsh <- length(fit$model[[3]]) + 1
  if(inherits(fit,"gev.d.fit")){
    if(!fit$theta_zero){
      npth <- length(fit$model[[4]]) + 1 #Including theta parameter (default)]
    }else{
      npth <- 0
    }#With no theta parameter, asked by user
    npet <- length(fit$model[[5]]) + 1
    if(!fit$eta2_zero){
      npe2 <- length(fit$model[[6]]) + 1 #Including eta2 parameter (not default)]
    }else{
      npe2 <- 0
    }
    if(!fit$tau_zero){
      npta <- length(fit$model[[7]]) + 1 #Including tau parameter (not default)]
    }else{
      npta <- 0
    }#With no tau parameter, asked by user
  }
  
  # inverse link functions
  if(inherits(fit,"gev.d.fit")){
    mulink <- fit$link$mutlink$linkinv
    siglink <- fit$link$sigma0link$linkinv
    shlink <- fit$link$xilink$linkinv
    if(!fit$theta_zero) thetalink <- fit$link$thetalink$linkinv
    etalink <- fit$link$etalink$linkinv
    if(!fit$eta2_zero) eta2link <- fit$link$eta2link$linkinv
    if(!fit$tau_zero) taulink <- fit$link$taulink$linkinv
  }else{
    mulink <- eval(parse(text=fit$link))[[1]]
    siglink <- eval(parse(text=fit$link))[[2]]
    shlink <- eval(parse(text=fit$link))[[3]]
  }
  
  # covariates matrices
  mumat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[1]]],dim(ydat)[1],npmu-1))
  sigmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[2]]],dim(ydat)[1],npsc-1))
  shmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[3]]],dim(ydat)[1],npsh-1))
  if(inherits(fit,"gev.d.fit")){
    if(!fit$theta_zero){thmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[4]]],dim(ydat)[1],npth-1))}
    etmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[5]]],dim(ydat)[1],npet-1))
    if(!fit$eta2_zero) {e2mat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[6]]],dim(ydat)[1],npe2-1))}
    if(!fit$tau_zero)  {tamat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[7]]],dim(ydat)[1],npta-1))}
  }
  
  # calculate parameters
  mut <- mulink(mumat %*% (fit$mle[1:npmu]))
  sc0 <- siglink(sigmat %*% (fit$mle[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (fit$mle[seq(npmu + npsc + 1, length = npsh)]))
  if(inherits(fit,"gev.d.fit")){
    if(!fit$theta_zero){
      theta <- thetalink(thmat %*% (fit$mle[seq(npmu + npsc + npsh + 1, length = npth)]))
    }else{
      theta <- rep(0,dim(ydat)[1])
    }
    eta <- etalink(etmat %*% (fit$mle[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
    if(!fit$eta2_zero){
      eta2 <- eta2link(e2mat %*% (fit$mle[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))
      # transform eta2 from eta2~~eta to eta2~~0
      eta2 <- eta2-eta
    }else{
      eta2 <- rep(0,dim(ydat)[1]) #eta
    }
    if(!fit$tau_zero){
      tau <- taulink(tamat %*% (fit$mle[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))
    }else{
      tau <- rep(0,dim(ydat)[1])
    }
    return(data.frame(mut=mut,sigma0=sc0,xi=xi,theta=theta,eta=eta, eta2=eta2, tau=tau))
  }else{
    return(data.frame(mu=mut,sig=sc0,xi=xi))
  }
}


#### example data ####
#' Sampled data for duration-dependent GEV
#' 
#' @description
#' Randomly sampled data set used for running the example code, containing:
#' \itemize{
#'   \item \code{$xdat}: 'annual' maxima values
#'   \item \code{$ds}: corresponding durations
#'   \item \code{$cov1}, \code{$cov2}: covariates}
#' d-GEV parameters used for sampling:
#' \itemize{
#'   \item \eqn{\tilde{\mu} = 4 + 0.2 cov_1 +0.5 cov_2}
#'   \item \eqn{\sigma_0 = 2+0.5 cov_1}
#'   \item \eqn{\xi = 0.5}
#'   \item \eqn{\theta = 0}
#'   \item \eqn{\eta = 0.5}
#'   \item \eqn{\eta_2 = 0.5}
#'   \item \eqn{\tau = 0}}
#'
#'
#' @docType data
#' @keywords datasets
#' @name example
#' @usage data('example',package ='IDF')
#' @format A data frame with 330 rows and 4 variables
#' 
#' ##############################################################
#' qgev...
#' ####################################################################
#' # This file contains the functions  dgev.d, pgev.d, qgev.d, rgev.d for the duration-dependent-gev.

#### dgev.d() ####

#' d-GEV probability density function
#'
#' @description Probability density function of duration-dependent GEV distribution
#' @param q vector of quantiles
#' @param mut,sigma0,xi numeric value, giving modified location \eqn{\tilde{\mu}}, scale offset \eqn{\sigma_0} and 
#' shape parameter \eqn{\xi}.
#' @param theta numeric value, giving duration offset \eqn{\theta} (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent \eqn{\eta} (defining slope of the IDF curve)
#' @param eta2 numeric value, giving a second duration exponent \eqn{\eta_{2}} (needed for multiscaling). Default: \eqn{\eta_2=0}.
#' @param tau numeric value, giving intensity offset \eqn{\tau} (defining flattening of the IDF curve). Default: \eqn{\tau=0}.
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{dgev}}
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @return list containing vectors of density values for given quantiles.
#' The first element of the list are the density values for the first given duration etc.
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{qgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd dgev 
#'
#' @examples
#' x <- seq(4,20,0.1)
#' # calculate probability density for one duration
#' dgev.d(q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1,d=1)
#' 
#' # calculate probability density for different durations
#' ds <- 1:4
#' dens <- lapply(ds,dgev.d,q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1)
#' 
#' plot(x,dens[[1]],type='l',ylim = c(0,0.21),ylab = 'Probability Density')
#' for(i in 2:4){
#'   lines(x,dens[[i]],lty=i)
#' }
#' legend('topright',title = 'Duration',legend = 1:4,lty=1:4)
dgev.d <- function(q,mut,sigma0,xi,theta,eta,d,eta2=0,tau=0,...) {
  #if(is.null(eta2)){eta2=eta}
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta),length(eta2),length(tau))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta, tau is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative,resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(q)))}else{
    
    #sigma.d <-sigma0/(d+theta)^eta+ tau        # old
    sigma.d <-      sigma0/(d+theta)^(eta+eta2) +tau
    mu.d    <- mut*(sigma0/(d+theta)^eta  +tau)
    
    return(dgev(q,loc=mu.d,scale=sigma.d,shape=xi,...))}
}


#### pgev.d() ####

#' d-GEV cumulative distribution function
#'
#' @description Cumulative probability distribution function of duration-dependent GEV distribution
#' @param q vector of quantiles
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param eta2 numeric value, giving a second duration exponent \eqn{\eta_{2}} (needed for multiscaling). Default: \eqn{\eta_2=0}.
#' @param tau numeric value, giving intensity offset \eqn{\tau} (defining flattening of the IDF curve). Default: \eqn{\tau=0}.
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{pgev}}
#' 
#' @details The duration dependent GEV distribution is defined after 
#' [Koutsoyiannis et al., 1998]:
#' \deqn{G(x)= \exp[-\left( 1+\xi(x/\sigma(d)-\mu_t) \right)^{-1/\xi}] } 
#' with the duration dependent scale \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta} and 
#' modified location parameter \eqn{\mu_t=\mu/\sigma(d)}.
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @return list containing vectors of probability values for given quantiles. 
#' The first element of the list are the probability values for the first given duration etc.
#' 
#' @seealso \code{\link{dgev.d}}, \code{\link{qgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd pgev 
#'
#' @examples
#' x <- seq(4,20,0.1)
#' prob <- pgev.d(q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1,d=1)
pgev.d <- function(q,mut,sigma0,xi,theta,eta,d,tau=0,eta2=0, ...) {
  #if(is.null(eta2)){eta2=eta}
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta),length(eta2),length(tau))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta, tau is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative,resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(q)))}else{
    
    #sigma.d <-sigma0/(d+theta)^eta+tau # old
    sigma.d <-      sigma0/(d+theta)^(eta+eta2) +tau
    mu.d    <- mut*(sigma0/(d+theta)^eta  +tau)
    
    return(pgev(q,loc=mu.d,scale=sigma.d,shape=xi,...))}
}


#### qgev.d() ####

#' d-GEV quantile function 
#'
#' @description Quantile function of duration-dependent GEV distribution (inverse of the cumulative probability distribution function)
#' @param p vector of probabilities
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve for short durations)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param eta2 numeric value, giving a second duration exponent \eqn{\eta_{2}} (needed for multiscaling). Default: \eqn{\eta_2=0}.
#' @param tau numeric value, giving intensity offset \eqn{\tau} (defining flattening of the IDF curve). Default: \eqn{\tau=0}.
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{qgev}}
#' 
#' @details The duration dependent GEV distribution is defined after 
#' [Koutsoyiannis et al., 1998]:
#' \deqn{ G(x)= \exp[-\left( 1+\xi(x/\sigma(d)-\mu_t) \right)^{-1/\xi}] } 
#' with the duration dependent scale \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta} and 
#' modified location parameter \eqn{\mu_t=\mu/\sigma(d)}.
#' 
#' @return list containing vectors of quantile values for given probabilities. 
#' The first element of the list are the q. values for the first given duration etc.
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{dgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd qgev 
#'
#' @examples
#' p <- c(0.5,0.9,0.99)
#' # calulate quantiles for one duration
#' qgev.d(p=p,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3, d=1)
#' 
#' # calculate quantiles for sequence of durations
#' ds <- 2^seq(0,4,0.1)
#' qs <- lapply(ds,qgev.d,p=p,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3)
#' qs <- simplify2array(qs)
#' 
#' plot(ds,qs[1,],ylim=c(3,20),type='l',log = 'xy',ylab='Intensity',xlab = 'Duration')
#' for(i in 2:3){
#'   lines(ds,qs[i,],lty=i)
#' }
#' legend('topright',title = 'p-quantile',
#'        legend = p,lty=1:3,bty = 'n')
qgev.d <- function(p,mut,sigma0,xi,theta,eta,d,tau=0,eta2=0, ...) {
  #if (is.null(eta2)){eta2=eta}
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta), length(eta2), length(tau))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta, eta2, tau is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative, resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(p)))}else{
    #sigma.d <-sigma0/(d+theta)^eta
    #sigma.d <-sigma0/(d+theta)^eta+tau
    
    sigma.d <-      sigma0/(d+theta)^(eta+eta2) +tau
    mu.d    <- mut*(sigma0/(d+theta)^eta  +tau)
    
    return(qgev(p,loc=as.numeric(mu.d)
                ,scale=as.numeric(sigma.d),shape=as.numeric(xi),...))}
  
  #return(qgev(p,loc=as.numeric(mut*sigma.d)                          # old
  #            ,scale=as.numeric(sigma.d),shape=as.numeric(xi),...))} # old
}


#### rgev.d() ####

#' Generation of random variables from d-GEV
#'
#' @description Generation of random variables following duration-dependent GEV.
#' @param n number of random variables per duration
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param eta2 numeric value, giving a second duration exponent \eqn{\eta_{2}} (needed for multiscaling). Default: \eqn{\eta_2=0}.
#' @param tau numeric value, giving intensity offset \eqn{\tau} (defining flattening of the IDF curve). Default: \eqn{\tau=0}.
#' @param d positive numeric value, giving duration
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}
#' 
#' @return list containing vectors of random variables.  
#' The first element of the list are the random values for the first given duration etc.
#' Note that the random variables for different durations are nor ordered (contrary to precipitation maxima of different durations).
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{qgev.d}}, \code{\link{dgev.d}}
#' 
#' @export
#' @importFrom evd rgev 
#'
#' @examples
#' # random sample for one duration
#' rgev.d(n=100,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3,d=1)
#' 
#' # compare randomn samples for different durations
#' ds <- c(1,4)
#' samp <- lapply(ds,rgev.d,n=100,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3)
#' 
#' hist(samp[[1]],breaks = 10,col=rgb(1,0,0,0.5),freq = FALSE
#'      ,ylim=c(0,0.3),xlim=c(3,20),xlab='x',main = 'Random d-GEV samples')
#' hist(samp[[2]],breaks = 10,add=TRUE,col=rgb(0,0,1,0.5),freq = FALSE)
#' legend('topright',fill = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),
#' legend = paste('d=',1:2,'h'),title = 'Duration')
rgev.d <- function(n,mut,sigma0,xi,theta,eta,d,tau=0,eta2=0) {
  #if (is.null(eta2)){eta2=eta}
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta),length(eta2),length(tau))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta, tau is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative, resulting from a negative theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,n))}else{
    
    #sigma.d <-sigma0/(d+theta)^eta+tau # old
    sigma.d <-      sigma0/(d+theta)^(eta+eta2) +tau
    mu.d    <- mut*(sigma0/(d+theta)^eta  +tau)
    
    return(rgev(n,loc=mu.d,scale=sigma.d,shape=xi))}
}
#' 
#' 
