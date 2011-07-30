tpam.project = function(data, data.test){
	#source("./subfunctions.R")
	X = as.matrix(data$x[[1]])
	for(i in 2:length(data$x)){
		X = cbind(X, data$x[[i]])
	}
	par = t(apply(X, 1, lda_project, data$y, length(data$x)))
	
	wx.train = multiply.func(data$x, par)
	wx.test = multiply.func(data.test$x, par)
	
	wdata.train = list(x = wx.train, y=data$y, genenames=data$genenames, geneid = data$geneid)
	wdata.test = list(x = wx.test, y=data.test$y, genenames =data.test$genenames, geneid = data.test$geneid)
	return(proj.obj = list(data.train = wdata.train, data.test = wdata.test))		
}