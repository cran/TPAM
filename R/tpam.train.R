tpam.train = function(data, data.test){
	proj.obj = tpam.project(data, data.test);
	data.train = proj.obj$data.train;
	fit.obj = pamr.train(data.train)
	return(train.obj = list(proj.obj = proj.obj, fit.obj=fit.obj))
}