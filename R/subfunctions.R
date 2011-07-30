require(MASS)
require(Matrix)
require(pamr)

softmax <-function(x, gap = FALSE) {
  d <- dim(x)
  maxdist <- x[, 1]
  pclass <- rep(1, d[1])
  for(i in seq(2, d[2])) {
    l <- x[, i] > maxdist
    pclass[l] <- i
    maxdist[l] <- x[l, i]
  }
  dd <- dimnames(x)[[2]]
  if(gap) {
    x <- abs(maxdist - x)
    x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
    gaps <- do.call("pmin", data.frame(x))
  }
  pclass <- if(is.null(dd) || !length(dd))
    pclass
  else
    factor(pclass, levels = seq(d[2]), labels = dd)
  if(gap)
    list(class = pclass, gaps = gaps)
  else
    pclass
}

diag.disc <-function(x, centroids, prior, weight) {
  if(! missing(weight)) {
    posid <- (weight > 0)
    if(any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid,  , drop = FALSE] * weight
      x <- x[posid,  , drop = FALSE] * weight
    }
    else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd <- t(x) %*% centroids
  dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  names(dd0) <- NULL
  scale(dd, dd0, FALSE)
}


soft.shrink <-function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- sign(delta) * dif * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

nsc <- function(x, y = NULL, xtest = NULL, proby = NULL, ytest = NULL, prob.ytest = 
        NULL, threshold = NULL, n.threshold = 30, hetero = NULL, scale.sd = 
        TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, 
        prior = table(y)/length(y), remove.zeros = TRUE, sign.contrast = "both",
           problem.type=c("class", "surv.km","surv.latent"))
{
        this.call <- match.call()

        argy <- ytest
        if(is.null(ytest)) {
                argy <- y
        }
   		if(!is.null(y) & !is.null(proby)){
                stop("Can't specify both y and proby")
        }
        if(!is.null(ytest) & !is.null(prob.ytest)) {
                stop("Can't specify both ytest and prob.ytest")
        }
        if(is.null(y)) {
                y <- apply(proby, 1, which.is.max)
        }
        n.class <- table(y)
        if(min(n.class) == 1) {
                cat("Warning: a class contains only 1 sample")
        }
        if(is.null(xtest)) {
                xtest <- x
                ytest <- y
                prob.ytest <- proby
        }
        norm.cent <- NULL
        if(!is.null(hetero)) {
                norm.cent <- apply(x[, y == hetero], 1, mean)
                x <- abs(t(scale(t(x), center = norm.cent, scale = FALSE)))
                if(!missing(xtest)) {
                        xtest <- abs(t(scale(t(xtest), center = norm.cent, 
                                scale = FALSE)))
                }
        }
        n <- sum(n.class)
        ntest <- ncol(xtest)
        K <- length(prior)
        p <- nrow(x)
        if(is.null(proby)) {
                Y <- model.matrix( ~ factor(y) - 1, data = list(y = y))
        }
        if(!is.null(proby)) {
                Y <- proby
        }
        dimnames(Y) <- list(NULL, names(n.class))

        centroids <- scale(x %*% Y, FALSE, n.class)
        sd <- rep(1, p)
        if(scale.sd) {
                xdif <- x - centroids %*% t(Y)
                sd <- (xdif^2) %*% rep(1/(n - K), n)
                sd <- drop(sqrt(sd))
                offset <- quantile(sd, offset.percent/100)
                sd <- sd + offset
        }
        centroid.overall <- drop(x %*% rep(1/n, n))
        if(is.null(threshold.scale)) {
                threshold.scale <- rep(1, K)
                names(threshold.scale) <- names(n.class)
        }
        if(is.null(se.scale))
                se.scale <- sqrt(1/n.class - 1/n)
        delta <- (centroids - centroid.overall)/sd
        delta <- scale(delta, FALSE, threshold.scale * se.scale)
        if(sign.contrast == "positive") {
                delta <- delta * (delta > 0)
        }
        if(sign.contrast == "negative") {
                delta <- delta * (delta < 0)
        }

        if(!is.null(threshold)) {
                n.threshold <- length(threshold)
        }
        else {
                threshold <- seq(0, max(abs(delta)), length = n.threshold)
        }
        nonzero <- seq(n.threshold)
        errors <- threshold
        yhat <- as.list(seq(n.threshold))
        prob <- array(0, c(ntest, K, n.threshold))
        for(ii in 1:n.threshold) {
                cat(ii)
                delta.shrunk <- soft.shrink(delta, threshold[ii])
                delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * 
                        se.scale))
                nonzero[ii] <- attr(delta.shrunk, "nonzero")
                posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
                dd <- diag.disc((xtest - centroid.overall)/sd, delta.shrunk, 
                        prior, weight = posid)
                yhat[[ii]] <- softmax(dd)
                dd <- safe.exp(dd)
                prob[,  , ii] <- dd/drop(dd %*% rep(1, K))
                if(!is.null(ytest)) {
                        errors[ii] <- sum(yhat[[ii]] != ytest)
                }
                if(!is.null(prob.ytest)) {

                        temp <- c(yhat[[ii]], names(table(y)))
                        Yhat <- model.matrix( ~ factor(temp) - 1, data = list(y
                                 = temp))
                        Yhat <- Yhat[1:length(yhat[[ii]]),  ]
                     
                        errors[ii] <- length(yhat[[ii]]) - sum(Yhat * prob.ytest)
                }
               
        }
        thresh.names <- format(round(threshold, 3))
        names(yhat) <- thresh.names
        attr(yhat, "row.names") <- paste(seq(ntest))
        class(yhat) <- "data.frame"
        if(remove.zeros)
                n.threshold <- match(0, nonzero, n.threshold)
        dimnames(prob) <- list(paste(seq(ntest)), names(n.class), thresh.names)
        object <- list(y = argy, proby = prob.ytest, yhat = yhat[, seq(
                n.threshold)], prob = prob[,  , seq(n.threshold)], centroids = 
                centroids, centroid.overall = centroid.overall, sd = sd, 
                threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(
                n.threshold)], threshold.scale = threshold.scale, se.scale = 
                se.scale, scale.sd=scale.sd, call = this.call, hetero = hetero, norm.cent = 
                norm.cent, prior = prior, offset = offset, sign.contrast = 
                sign.contrast)
        if(!is.null(ytest) | !is.null(prob.ytest))
                object$errors <- errors[seq(n.threshold)]
        class(object) <- "nsc"
        object
}


safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}

permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}


 balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
     
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
cat(i)
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }
 
which.is.max <- function(x)
{
        y <- seq(length(x))[x == max(x)]
        if(length(y) > 1)
                y <- sample(y, 1)
        y
}

my.nsccv <- function(x, y=NULL, proby=NULL, nfold = min(table(y)), folds = NULL, threshold =
        NULL, threshold.scale = NULL, survival.time=NULL, censoring.status=NULL, ngroup.survival=NULL,prior, object, ...)
{
        this.call <- match.call()

        argy <- y
          if(is.null(y)){ y <- as.factor(apply(proby,1,which.is.max))}
        
        n <- length(y)

if(is.null(nfold) & is.null(survival.time)) {nfold <- min(table(y))}
if(is.null(nfold) & !is.null(survival.time)) {nfold <- 10}


 if(is.null(survival.time)){
        if(is.null(folds)) {
          #folds <- balanced.folds(y)  ## yuping deleted
          folds <-balanced.folds(y, nfold)  ## yuping added
        }
       }


        if(!is.null(survival.time)){
        if(is.null(folds)) {
                folds <- split(sample(1:n), rep(1:nfold, length = n))
        }
        }
         
nfold<- length(folds)

        if(missing(prior)) {
                if(missing(object))
                        prior <- table(y)/n
                else prior <- object$prior
        }
    
        if(missing(threshold)) {
                if(missing(object))
                        stop("Must either supply threshold argument, or an nsc object"
                                )
                else {
                        threshold <- object$threshold
                        threshold.scale <- object$threshold.scale
                        se.scale <- object$se.scale
                }
        }
       
        n.threshold <- length(threshold)        ### Set up the data structures
        yhat <- rep(list(y), n.threshold)
        names(yhat) <- paste(seq(n.threshold))
        yhat <- data.frame(yhat)
        n.class <- table(y)
        prob <- array(1, c(n, length(n.class), n.threshold))
        size <- double(n.threshold)
        hetero <-object$hetero
        cv.objects=vector("list",nfold)
        for(ii in 1:nfold) {
                cat("Fold", ii, ":")
                a <- nsc(x[,  - folds[[ii]]], y=argy[ - folds[[ii]]], x[, folds[[ii
                        ]], drop = FALSE], proby=proby[-folds[[ii]],],
                         threshold = threshold, threshold.scale
                         = threshold.scale, se.scale = se.scale, prior = prior,
                          hetero=hetero,
                        ..., remove.zeros = FALSE)
                size <- size + a$nonzero
                prob[folds[[ii]],  ,  ] <- a$prob
                yhat[folds[[ii]],  ] <- a$yhat
                cat("\n")
        cv.objects[[ii]]=a
        }
        if(missing(object))
                size <- round(size/nfold)
        else size <- object$nonzero
        error <- rep(NA, n.threshold)
        loglik <- error
        pvalue.survival <- error
        
        pvalue.survival.func <- function(group, survival.time, censoring.status,ngroup.survival){
            temp <- coxph(Surv(survival.time, censoring.status)~as.factor(group))
            loglik <- 2*(temp$loglik[2]-temp$loglik[1])
            return(1-pchisq(loglik, ngroup.survival-1))
          }
        
        if(!is.null(proby)){proby.temp <-proby}
        else if(!is.null(survival.time)){proby.temp <- pamr.surv.to.class2(survival.time,
                                       censoring.status, n.class=ngroup.survival)$prob
                                       }
        
        for(i in 1:n.threshold) {
      
                if(is.null(survival.time) & is.null(proby)){error[i] <- sum(yhat[, i] != y)/n}
                if(!is.null(survival.time)){
                    
                    temp <- c(yhat[,i],names(table(y)))
                    Yhat <- model.matrix( ~ factor(temp) - 1,
                                       data = list(y = temp))
                     Yhat <- Yhat[1:length(yhat[[ii]]),]
                     error[i] <- (length(yhat[,i])-sum(Yhat*proby.temp))/n
                  }
            
                
                if(is.null(survival.time)){
                  loglik[i] <- sum(log(prob[,  , i][cbind(seq(1, n), unclass(y))]))/                        n}
                
                if(!is.null(survival.time)){
                  pvalue.survival[i]<- pvalue.survival.func(yhat[,i], survival.time,censoring.status, ngroup.survival)
                }
        }

obj<- list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds, cv.objects=cv.objects, pvalue.survival=pvalue.survival,
                call = this.call)
        class(obj) <- "nsccv"
        obj
}



my.pamr.cv <-
function(fit, data, nfold = NULL, folds = NULL ,...)
{
        x <- data$x[fit$gene.subset, fit$sample.subset]

        if( !is.null(data$y) & !is.null(data$proby)){
           stop("Must have exactly one of y and  proby  present in the data object")
         }
        
        y <- NULL
        proby <- NULL
        
        if(!is.null(fit$y)){
           y<-  factor(fit$y[fit$sample.subset])
         }
        
        if(!is.null(fit$proby)){
           proby<-  fit$proby[fit$sample.subset,]
         }
        
        this.call <- match.call()
         junk <- my.nsccv(x, y=y, proby=proby, object = fit, nfold=nfold, folds=folds, survival.time=data$survival.time, censoring.status = data$censoring.status, ngroup.survival=fit$ngroup.survival, problem.type=fit$problem.type, ...)

        junk$call <- this.call
        
        junk$sample.subset <- fit$sample.subset
        class(junk)="pamrcved"
        junk
}



my.pamr.plotcv <- function(fit) {
  par(mar = c(5, 5, 5, 1))
  par(mfrow = c(2, 1))
  n <- nrow(fit$yhat)
  y <- fit$y
  if(!is.null(fit$newy)) {
    y <- fit$newy[fit$sample.subset]
  }
  nc <- length(table(y))
  nfolds <- length(fit$folds)
  err <- matrix(NA, ncol = ncol(fit$yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(fit$yhat), nrow = n)
  ni <- rep(NA, nfolds)
  for(i in 1:nfolds) {
    ii <- fit$folds[[i]]
    ni[i] <- length(fit$folds[[i]])
    err[i,  ] <- apply(temp[ii,  ] != fit$yhat[ii,  ], 2, sum)/ni[i]
  }
  se <- sqrt(apply(err, 2, var)/nfolds)
  plot(fit$threshold, fit$error, ylim = c(-0.1, 0.8), xlab = 
       "Value of threshold  ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = fit$threshold, lab = paste(fit$size), srt = 90, adj = 0)
  mtext("Number of genes", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  lines(fit$threshold, fit$error, col = 2)
  o <- fit$err == min(fit$err)
  points(fit$threshold[o], fit$error[o], pch = "x")
  error.bars(fit$threshold, fit$err - se, fit$err + se)
  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(fit$threshold
                                                 ))
  for(i in 1:(length(fit$threshold) - 1)) {
    s <- pamr.confusion(fit, fit$threshold[i], extra = FALSE)
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }
  plot(fit$threshold, err2[1,  ], ylim = c(-0.1, 1.1), xlab = 
       "Value of threshold ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = fit$threshold, lab = paste(fit$size), srt = 90, adj = 0)     
                                        #       mtext("Number of genes", 3, 4,cex=1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  for(i in 1:nrow(err2)) {
    lines(fit$threshold, err2[i,  ], col = i + 1)
  }
  legend(0, 0.9, dimnames(table(y))[[1]], col = (2:(nc + 1)), lty = 1)
  par(mfrow = c(1, 1))

  return(se);
}

my.pamr.cv.se <- function(fit){
  n <- nrow(fit$yhat)
  y <- fit$y
  if(!is.null(fit$newy)) {
    y <- fit$newy[fit$sample.subset]
  }
  nc <- length(table(y))
  nfolds <- length(fit$folds)
  err <- matrix(NA, ncol = ncol(fit$yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(fit$yhat), nrow = n)
  ni <- rep(NA, nfolds)
  for(i in 1:nfolds) {
    ii <- fit$folds[[i]]
    ni[i] <- length(fit$folds[[i]])
    err[i,  ] <- apply(temp[ii,  ] != fit$yhat[ii,  ], 2, sum)/ni[i]
  }
  se <- sqrt(apply(err, 2, var)/nfolds)
  return(se);
}

error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

multiply.func <- function(xlist, par){
	tmp = xlist[[1]] * par[,1];
	for(i in 2:ncol(par)){
		tmp = tmp + xlist[[i]] * par[,i];
	}
	x.return = matrix(unlist(tmp), nrow=nrow(par));
}

pred.err = function(threshold, train.obj, data.test, type="class"){
	predict_class = pamr.predict(train.obj, data.test$x, threshold=threshold, type=type);
	err = 0
	for(i in 1:length(predict_class)){
		if(predict_class[i]!=data.test$y[i]){err = 1 +err};
	}
	return(errate = err/length(predict_class))
}


lda_project = function(x, y, n.tp=2){
	x = matrix(x, nrow=length(y))
	tt = lda(x, y)$scaling
	return(tt[,1])
}


sqrt_norm = function(x){return(sqrt(sum(x^2)));}
