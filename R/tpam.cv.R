tpam.cv = function(fit, data, nfold = NULL, folds=NULL){
	this.call = match.call()
	y = data$y
	n = length(y)
	if(is.null(nfold)) {nfold <- min(table(y))}
    if(is.null(folds)) {
          folds <-balanced.folds(y, nfold) 
    }
            
	nfold<- length(folds)
	if(is.null(fit$prior)){
		prior = table(y)/n
	}else{
		prior = fit$prior
	}
    threshold = fit$threshold
    threshold.scale = fit$threshold.scale
    se.scale = fit$se.scale
    n.threshold <- length(threshold)        ### Set up the data structures
    yhat <- rep(list(y), n.threshold)
    names(yhat) <- paste(seq(n.threshold))
    yhat <- data.frame(yhat)
    n.class <- table(y)
    prob <- array(1, c(n, length(n.class), n.threshold))
    size <- double(n.threshold)
 	hetero <- fit$hetero
    cv.objects=vector("list",nfold)
    for(ii in 1:nfold) {
        cat("Fold", ii, ":")
        par.train = matrix(nrow=nrow(data$x[[1]]), ncol=length(data$x));
        y.train.tmp = data$y[-folds[[ii]]];
        y.test.tmp = data$y[folds[[ii]]];
        x.train.tmp = matrix(unlist(data$x[[1]][, - folds[[ii]]]), nrow=nrow(data$x[[1]]));
        xlist.train = list();
		xlist.train[[1]] = data$x[[1]][, -folds[[ii]]];
		xlist.test = list();
		xlist.test[[1]] = data$x[[1]][, folds[[ii]]];
		for(kkk in 2:length(data$x)){
			xtmp.kkk = matrix(unlist(data$x[[kkk]][, -folds[[ii]]]), nrow=nrow(data$x[[1]]));
			x.train.tmp = cbind(x.train.tmp, xtmp.kkk);
			xlist.train[[kkk]] = data$x[[kkk]][, -folds[[ii]]];
			xlist.test[[kkk]] = data$x[[kkk]][, folds[[ii]]];	
		}        
        
		par.train = t(apply(x.train.tmp, 1, lda_project, y.train.tmp, length(data$x)));
				
		wx.train.tmp = multiply.func(xlist.train,  par.train);
		
		wx.test.tmp = multiply.func(xlist.test,  par.train);
 		
 		data.temp = list(x = wx.train.tmp, y = y.train.tmp);
 		data.test.temp = list(x = wx.test.tmp, y = y.test.tmp);    
     	a <- nsc(x = wx.train.tmp, y=y.train.tmp, wx.test.tmp, proby=fit$proby[-folds[[ii]],],
                         threshold = threshold, threshold.scale
                         = threshold.scale, se.scale = se.scale, prior = prior,
                          hetero=hetero, remove.zeros = FALSE)
                size <- size + a$nonzero
                prob[folds[[ii]],  ,  ] <- a$prob
                yhat[folds[[ii]],  ] <- a$yhat
                cat("\n")
        cv.objects[[ii]]=a
        }
        
        if(missing(fit))
                size <- round(size/nfold)
        else size <- fit$nonzero
        error <- rep(NA, n.threshold)
        loglik <- error
   for(i in 1:n.threshold) {
  		error[i] <- sum(yhat[, i] != y)/n
  		loglik[i] <- sum(log(prob[,  , i][cbind(seq(1, n), unclass(y))]))/n
   }
  obj<- list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds, cv.objects=cv.objects, call = this.call)
        class(obj) <- "nsccv"
        return(obj)	
}