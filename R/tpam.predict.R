tpam.predict = function(fit, newdata, threshold, type=c("class", "posterior", "centroid", "nonzero"), prior = fit$prior, threshold.scale = fit$threshold.scale){
	
	newx = newdata$x
	obj = pamr.predict(fit, newx, threshold, type= type, prior = prior, threshold.scale = threshold.scale)
	return(obj)	
} 