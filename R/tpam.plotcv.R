tpam.plotcv = function(fit){
	 n <- nrow(fit$yhat)
    y <- fit$y
    if (!is.null(fit$newy)) {
        y <- fit$newy[fit$sample.subset]
    }
    nc <- length(table(y))
    nfolds <- length(fit$folds)
    err <- matrix(NA, ncol = ncol(fit$yhat), nrow = nfolds)
    temp <- matrix(y, ncol = ncol(fit$yhat), nrow = n)
    ni <- rep(NA, nfolds)
    for (i in 1:nfolds) {
        ii <- fit$folds[[i]]
        ni[i] <- length(fit$folds[[i]])
        err[i, ] <- apply(temp[ii, ] != fit$yhat[ii, ], 2, sum)/ni[i]
    }
    se <- sqrt(apply(err, 2, var)/nfolds)
    plot(fit$threshold, fit$error, ylim = c(-0.1, 0.8), xlab = "Value of threshold  ", 
        ylab = "Misclassification Error", type = "n", yaxt = "n")
    axis(3, at = fit$threshold, lab = paste(fit$size), srt = 90, 
        adj = 0)
    mtext("Number of genes", 3, 4, cex = 1.2)
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
    lines(fit$threshold, fit$error, col = 2)
    o <- fit$err == min(fit$err)
    points(fit$threshold[o], fit$error[o], pch = "x")
    error.bars(fit$threshold, fit$err - se, fit$err + se)
}