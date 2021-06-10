#### source functions to implement CuSum Alarm
#### Author: Meng Wang
#### Email: mengw1@stanford.edu
#### Date: 2021 June

#############################
####    online alarm     ####
#############################
library(xts)
timerange1 = "20160106 0000/20160106 2359"
min.nm = format(timeBasedSeq(timerange1), "%H:%M:%S")
hour.nm = min.nm[(1:length(min.nm)) %% 60 == 1]
	
online.alarming.fn = function(peo.id.1, dir.hr, dir.step, track.par=12, gap.thres = 14, watch.type='fitbit', tune.ind=FALSE, covid.day=NULL, plot.shift.day=0, dat.saved=FALSE) {
	
	#### to create result subfolders
	dir.create('../output/', showWarnings = FALSE, recursive = TRUE)
	dir.fig = "../output/figure/"
    dir.tab = "../output/table/"
    dir.result = "../output/result/"
    dir.clean = "../output/clean/"
    dir.note = "../output/note/"
    dir.create(dir.fig, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir.tab, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir.result, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir.clean, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir.note, showWarnings = FALSE, recursive = TRUE)

	if(dat.saved == FALSE) {
		print('obtaining smoothed RHR data')
		#### to load data
		formatted.dat = dat.format.fn(dir.hr, dir.step, watch.type)
		dat.hr = formatted.dat$dat.hr
		dat.step = formatted.dat$dat.step
				
		save(dat.hr, file=paste(dir.clean, peo.id.1, "_HR.Rdata", sep=""))
		save(dat.step, file=paste(dir.clean, peo.id.1, "_Step.Rdata", sep=""))
		
		#### to obtain RHR
		X.rhr = rhr.fn(dat.hr, dat.step, smth.k.par = 10, rest.min.par = 10, resol.tm.par = 60)		
		save(X.rhr, file=paste(dir.clean, peo.id.1, "_smoothed_RHR.Rdata", sep=""))
	
	} else {
	  print('loading the smoothed RHR data')
	  load(paste(dir.clean, peo.id.1, "_smoothed_RHR.Rdata", sep=""))		
	}
	
	#### to split data to short-term subsets
	if (is.null(X.rhr)) {
		stop.note = 'Total collected RHR data is less than one day. Please wear the smartwatch more often.'
		last.day = dat.hr[nrow(dat.hr),"date"]
		write.table(stop.note, file=paste(dir.note ,peo.id.1, "_not_enough_historical_data_till_", last.day, ".txt", sep=""), row.names = FALSE, col.names = FALSE, quote=FALSE)
	} else {
		hr.dates = colnames(X.rhr)
		day.gap = diff(as.Date(hr.dates))
		gap.ind = (1:length(hr.dates))[day.gap >= gap.thres]
		day.chunk = cbind(c(1, pmin(gap.ind+1, length(hr.dates))), c(pmin(gap.ind, length(hr.dates)), length(hr.dates)) ) 
		
		if (tune.ind == TRUE) {
		   par.seq = seq(0.9, 0.99, by=0.01)
		} else {
		   par.seq = 0.95
		}
		for (par.ind in 1:length(par.seq)) {
			par.0 = par.seq[par.ind]
			#print(par.0)
			
			eval.all = rep(NA, length(hr.dates))
			names(eval.all) = paste(hr.dates, '21:00:00')
			
			online.result.all = c()
			offline.result.all = c()
			
			for (i.chunk in 1:nrow(day.chunk)) {
				hr.dates.keep = hr.dates[day.chunk[i.chunk, 1]:day.chunk[i.chunk, 2]]
				if (length(hr.dates.keep) <= 28){
					nm = paste(hr.dates.keep, '21:00:00')
					eval.all[nm] = 'Baseline'
				} else{
					dir.fig.train = paste(dir.fig, peo.id.1, '_figure_under_par', par.0, "/", sep="")
					dir.tab.train = paste(dir.tab, peo.id.1, '_table_under_par_',par.0, "/", sep="")
					dir.result.train = paste(dir.result, peo.id.1, '_result_under_par_',par.0, "/", sep="")			
					dir.create(dir.fig.train, showWarnings = FALSE, recursive = TRUE)
					dir.create(dir.tab.train, showWarnings = FALSE, recursive = TRUE)
					dir.create(dir.result.train, showWarnings = FALSE, recursive = TRUE)
							           
		            #### to group data in 2*28days			
					tm.ls = unique(c(hr.dates.keep[seq(1, length(hr.dates.keep), by=28)], as.character(as.Date(hr.dates.keep[length(hr.dates.keep)])+1) ))
					
					#### to combine test stats
					print(paste(length(tm.ls) - 2, 'chunck(s) of data to evaluate'))
					for (i in 1:(length(tm.ls) - 2)) {
						print(paste('evaluating data in chunk', i))
						day.start = tm.ls[i]
						day.end = tm.ls[i+2]
						X.rhr.1 = X.rhr[, as.Date(colnames(X.rhr)) >= as.Date(day.start) & as.Date(colnames(X.rhr)) <= (as.Date(day.end) - 1)]
						x = X.rhr.1[, ncol(X.rhr.1)]
						x[(grep("21:00:00", names(x)) + 1):length(x)] = NA
						X.rhr.1[, ncol(X.rhr.1)] = x
						id = paste(peo.id.1, "_from_", day.start, "_to_", day.end, "", sep="")
						
						stats.result = stats.track.fn(id, X.rhr.1, base.num=28, test.r.fixed.1=FALSE, test.r.1=NA, res.quan.1=par.0, pval.thres=0.01)
				
						if (!is.null(stats.result)) {				    						
							res.t = stats.result$res.t
							cusum.t = stats.result$test.t
							cusum.test.r = stats.result$test.r.seq
							cusum.pval = stats.result$test.pval
				        } else {
							res.t = NULL
							cusum.t = NULL
							cusum.test.r = NULL
							cusum.pval = NULL
				        }
											
					    if (i == 1) {
					    	res.t.comb = res.t
					    	cusum.t.comb = cusum.t
					    	cusum.test.r.comb = cusum.test.r
					    	cusum.pval.comb = cusum.pval
					       } else {
					    	day.update = tm.ls[i+1]
					    	res.t.comb = c(res.t.comb, res.t[as.numeric(as.Date(day.update) - as.Date(names(res.t))) <= 0])
					    	cusum.t.comb = c(cusum.t.comb, cusum.t[as.numeric(as.Date(day.update) - as.Date(names(cusum.t))) <= 0])
					    	cusum.test.r.comb = c(cusum.test.r.comb, cusum.test.r[as.numeric(as.Date(day.update) - as.Date(names(cusum.test.r))) <= 0])
					    	cusum.pval.comb = c(cusum.pval.comb, cusum.pval[as.numeric(as.Date(day.update) - as.Date(names(cusum.pval))) <= 0])
					    }
					}
				
					cusum.t.comb.1 = cusum.t.comb[!is.na(cusum.t.comb)]
					cusum.t.ind.comb.1 = base.test.ind.fn(cusum.t.comb.1)
					cusum.pval.comb.1 = cusum.pval.comb[names(cusum.t.comb.1)]
					online.result.comb = cusum.detection.fn(id=peo.id.1, test.t=cusum.t.comb.1, test.t.ind=cusum.t.ind.comb.1, test.pval=cusum.pval.comb.1, pval.thres=0.01, hr.track=track.par)
					
					offline.result.comb = rhr.diff.detection.fn(peo.id.1, res.t.comb, alpha=0.05)
					
					save(peo.id.1, par.0, tm.ls, online.result.comb, offline.result.comb, res.t.comb, cusum.t.comb, cusum.test.r.comb, cusum.pval.comb, file=paste(dir.result.train, peo.id.1, "_result_under_par", par.0, '_in_chunk', i.chunk, ".RData", sep="")) 
					
					online.result.all = rbind(online.result.all, online.result.comb$alarm.sum)
					offline.result.all = rbind(offline.result.all, offline.result.comb)
					
					## to obtain alarming table
					report.tm = 21 
					nm.dum = seq(as.Date(names(cusum.t.comb)[1]), as.Date(names(cusum.t.comb)[length(cusum.t.comb)]), by=1)
				    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
				    alarm.temp = rep(NA, length(nm.dum))
				    names(alarm.temp) = nm.dum
				    alarm.temp[names(cusum.t.comb)[!is.na(names(cusum.t.comb))]] = 0
				    alarm.bump = online.result.comb$alarm.cate
				    alarm.temp[names(alarm.bump)] = alarm.bump
				    last.day = substring(names(alarm.temp[length(alarm.temp)]), 1, 10)
				    alarm.temp = alarm.temp[names(alarm.temp) <= paste(last.day," ", report.tm, ":00:00", sep="")]
	                bk.tm = paste(unique(as.Date(names(alarm.temp))), " ", report.tm, ":00:00", sep="")
	                alarm.bin = cut(1:length(alarm.temp), breaks=c(0, (1:length(alarm.temp))[names(alarm.temp) %in% bk.tm] ))
	                alarm.temp.1 = alarm.temp
	                alarm.temp.1[is.na(alarm.temp)] = -1
	                eval.0 = tapply(alarm.temp.1, alarm.bin, max, na.rm=TRUE)                           
	                names(eval.0) = bk.tm
	                #table(eval.0)
	                eval = eval.0
	                eval[eval == 0] = "Green"
	                eval[eval == 1] = "Yellow"
	                eval[eval == 2] = "Red"
	                eval[eval == -1] = "N/A"
	                eval[1:28] = 'Baseline'					
					eval.all[names(eval)] = eval
								
					#### plot
					pdf(file=paste(dir.fig.train, peo.id.1,"_figure_par", par.0, '_chunk',i.chunk, ".pdf", sep=""), width=17, height=8)
					result.plot.fn(id=paste(peo.id.1, "from data Part", i.chunk, "under threshold =", round(par.0,2)), base.days.int0=NA, res.t=res.t.comb, cusum.t=cusum.t.comb, cusum.test.r=cusum.test.r.comb, offline.result=offline.result.comb, alarm.cate.tm=eval.0, day.info=covid.day, shift.day=0)
					dev.off()	
	                
	                if(plot.shift.day != 0) {
						pdf(file=paste(dir.fig.train, peo.id.1,"_shifted_figure_par", par.0, '_chunk',i.chunk, ".pdf", sep=""), width=17, height=8)
						result.plot.fn(id=paste(peo.id.1, "from data Part", i.chunk,"under threshold =", round(par.0,2)), base.days.int0=NA, res.t=res.t.comb, cusum.t=cusum.t.comb, cusum.test.r=cusum.test.r.comb, offline.result=offline.result.comb, alarm.cate.tm=eval.0, day.info=covid.day, shift.day=plot.shift.day)
						dev.off()
					}									
					
	          }
		   } 			   			
			
			write.csv(online.result.all, file=paste(dir.tab.train, peo.id.1, "_online_result_under_par", par.0, ".csv", sep=""), row.names=FALSE)
	
			write.csv(offline.result.all, file=paste(dir.tab.train, peo.id.1, "_offline_result_under_par", par.0, ".csv", sep=""), row.names=FALSE)

		    tab = cbind(names(eval.all), eval.all)
		    colnames(tab) = c("evaluation time", "alarm for previous 24 hours")
		    head(tab)
		    tab  = tab[order(tab[,1]), ]
		    write.csv(tab, file=paste(dir.tab.train, peo.id.1, "_eval_under", par.0, ".csv", sep=""), row.names=FALSE) 
		    
		    
	    }
	  				
	}
}
		
########################################################
####    to obtain residuals and CuSum statistics    ####
####            based on sliding window             ####
########################################################
rhr.fn = function (dat.hr.1, dat.step.1, smth.k.par = 10, rest.min.par = 10, resol.tm.par = 60) {	
	day.nm = unique(dat.hr.1[,"date"])
	length(day.nm)

	hr.mx = matrix(NA, length(min.nm), length(day.nm))
	rownames(hr.mx) = min.nm
	colnames(hr.mx) = day.nm
	for (i in 1:length(day.nm)) {
		day.1 = day.nm[i]
		day.hr = dat.hr.1[dat.hr.1[,"date"] == day.1, "hr"]
		names(day.hr) = dat.hr.1[dat.hr.1[,"date"] == day.1, "time"]
		hr.mx[names(day.hr), day.1] = day.hr
	}
	head(hr.mx)
		
	hr.mx = apply(hr.mx, 2, as.numeric)
	rownames(hr.mx) = min.nm
	colnames(hr.mx) = day.nm
				
    day.train.0 = colnames(hr.mx)
	dat.train.0 = hr.mx[, day.train.0]
	X = dat.train.0 
	rownames(X) = rownames(dat.train.0)
	x = c(X)
	K = smth.k.par
	x.conv = c(rep(NA, K-1), rollmean(x, k=K, align="right"))
	names(x.conv) = paste(rep(colnames(X), each=nrow(X)), rep(rownames(X),ncol(X)), sep=" ")
	X.conv = matrix(x.conv, nrow=nrow(X))
	rownames(X.conv) = rownames(X)
	colnames(X.conv) = colnames(X)

    X.conv.f = X.conv
	X.conv.f[is.na(X)] = NA
	for (j in 1:length(day.train.0)) {
		day.1 = day.train.0[j]
		day.step.1 = matrix(dat.step.1[dat.step.1[,"date"] == day.1, c("time", "steps")], ncol=2)
		ind = day.step.1[,2]
		names(ind) = day.step.1[,1]
		ind.1 = names(ind)[ind > 0]
		ind.ext = unique(c(unlist(sapply(ind.1, function(tm) rownames(X.conv)[grep(tm, rownames(X.conv))+0:rest.min.par] ))))
		ind.ext = ind.ext[!is.na(ind.ext)]
		X.conv.f[ind.ext ,j] = NA
		X.conv.f[names(ind)[ind > 0], j] = NA
	}
    
    flag = colSums(!is.na(X.conv.f))
    if (sum(flag > 0) <= 1) {
    	#print('collectd RHR data is less than one day')
        X.conv.f = NULL
    } else {
        X.conv.f = X.conv.f[, flag > 0]    	
    }
    return(X.conv.f)
}
	
	
	
cusum.detection.fn = function (id, test.t, test.t.ind, test.pval, pval.thres, hr.track) {

    ####
    risk.score = 1 - test.pval 
    test.sum = cbind(test.t, test.t.ind, risk.score)
    
	ind = rep(0, length(test.t.ind))
	names(ind) = names(test.t.ind)
	ind[test.t.ind > 0] = 1
	ind.mx = with(rle(ind), data.frame(number = values, start = cumsum(lengths) - lengths + 1, end = cumsum(lengths))[order(values),])
    ind.mx = ind.mx[ind.mx$number == 1, ]
    
        
    score.thres = 1 - pval.thres
    risk.score[is.na(risk.score)] = 0
    test.alarming.tm = c()
    run.start = c()
    run.end = c()
    run.max = c()
    run.num = c()
    run.max.val = c()
    run.len = c()
    alarm.bin.cate.comb = c()
    run.alarm.red.comb = c()

    #i = (1:nrow(ind.mx))[rownames(ind.mx) == "2020-09-23 02:00:00"]
    for (i in 1:nrow(ind.mx)) {
    	x.risk = risk.score[ind.mx[i, 2]:ind.mx[i,3]]
    	x.test = test.t[ind.mx[i, 2]:ind.mx[i,3]]

    	if (max(x.risk, na.rm=TRUE) > score.thres) {
    		alarm.ind = min((1:length(x.risk))[x.risk > score.thres])
    		alarm.tm = names(x.risk)[alarm.ind]
    		tm.ext = seq(as.POSIXct(names(x.risk[1]), tz='GMT'), as.POSIXct(names(x.risk[length(x.risk)]), tz='GMT'), by="1 hours")
    		x.test.ext = rep(NA, length(tm.ext))
    		names(x.test.ext) = tm.ext
    		x.test.ext[names(x.test)] = x.test
    		alarm.bin.ind = unique(c(seq((1:length(x.test.ext))[names(x.test.ext) == alarm.tm], length(x.test.ext), by=hr.track), length(x.test.ext)))   
            alarm.bin.tm = names(x.test.ext)[alarm.bin.ind]
            test.bin = cut(1:length(x.test.ext), breaks=c(0,alarm.bin.ind))
            avg.bin = tapply(x.test.ext, test.bin, mean, na.rm=TRUE)
            alarm.bin.tm = alarm.bin.tm[!is.na(avg.bin)]  
            avg.bin = avg.bin[!is.na(avg.bin)]
                      
            alarm.bin.cate = rep(0, length(avg.bin)-1)
            names(alarm.bin.cate) = alarm.bin.tm[-1]
            if(length(avg.bin) >= 2) {
               alarm.bin.cate[diff(avg.bin) > 0] = 1	
               if (length(avg.bin) >= 3) {
               	    for (j in 3:length(avg.bin)) {
	            		avg.bin.1 = avg.bin[(j-2):j]
	            		if (sum(diff(avg.bin.1) >= 0) == 2) {
	            			alarm.bin.cate[j-1] = 2
	            		}
            	    }            	
               }
            }	

            alarm.bin.cate.comb = c(alarm.bin.cate.comb, alarm.bin.cate)
    		
    		if (sum(alarm.bin.cate == 2) > 0) {
    			test.alarming.tm = c(test.alarming.tm, alarm.tm)
	    		run.start = c(run.start, names(x.risk)[1]) 		
	    	    run.end = c(run.end, names(x.risk)[length(x.risk)]) 		
	    		run.max = c(run.max, names(which.max(x.test)))
	    		run.num = c(run.num, which.max(x.test))
	    		run.max.val = c(run.max.val, max(x.test, na.rm=TRUE))
	    		run.len = c(run.len, length(x.risk))	
	    		run.alarm.red = paste(names(alarm.bin.cate)[alarm.bin.cate == 2], collapse="; ")
	    		run.alarm.red.comb = c(run.alarm.red.comb, run.alarm.red)
    		}
	
    	}
    }
    
    if (!is.null(test.alarming.tm)) {
	    alarm.tab = cbind(test.alarming.tm, run.alarm.red.comb, run.start, run.end, run.max, run.num, round(run.max.val,2), run.len)
	    alarm.tab = as.matrix(alarm.tab, ncol=10)
	    colnames(alarm.tab) = c("alarming.time", "red.alarm.time", "anormaly_starting_time", "anormaly_ending_time",  "anormaly_max_time", "duration_hours_to_max", "max_test_stats", "total_hours")	    
	    alarm.sum = cbind(rep(id, nrow(alarm.tab)), alarm.tab[,  "anormaly_starting_time"], alarm.tab[, "anormaly_ending_time"], alarm.tab[,  "alarming.time"], alarm.tab[,"red.alarm.time"], test.sum[alarm.tab[, "alarming.time"], 2], test.sum[alarm.tab[, "anormaly_max_time"], 2], alarm.tab[,"total_hours"], alarm.tab[,"max_test_stats"],  round(as.numeric(test.sum[alarm.tab[, "alarming.time"], 1]), 2), round(1 - as.numeric(test.sum[alarm.tab[, "alarming.time"], 3]), 2) )
	    colnames(alarm.sum) = c("ParticipantID", "anormaly starting time", "anormaly ending time", "initial alarming time", "red alarming time",  "alarming hour", "duration hours to max", "total hours", "max CuSum statistics", "CuSum statistics at alarm", "p-value")
	    rownames(alarm.sum) = NULL
	    
    } else {
    	alarm.sum = matrix(c(id, rep("N/A", 10)), ncol=11)
        colnames(alarm.sum) = c("ParticipantID", "anormaly starting time", "anormaly ending time", "initial alarming time", "red alarming time", "alarming hour", "duration hours to max", "total hours", "max CuSum statistics", "CuSum statistics at alarm", "p-value")
    }
    
    return(list(alarm.sum=alarm.sum,  alarm.cate=alarm.bin.cate.comb))
}



#### preformatting the dataset
dat.format.fn = function(dir.hr, dir.step, watch.type) {
	dat.hr = read.csv(dir.hr)
	dat.hr = as.matrix(dat.hr)
	#head(dat.hr)
	
	# to summarize hear rate in one minute resolution
	hr.med = tapply(as.numeric(dat.hr[, "heartrate"]), substring(dat.hr[,"datetime"], 1, 16), median, na.rm=TRUE)
	dat.hr.1 = cbind(paste(names(hr.med), ":00", sep=""), hr.med )
	dat.hr.1 = cbind(substring(dat.hr.1[,1], 1, 10), substring(dat.hr.1[,1], 12), dat.hr.1[,2])
	colnames(dat.hr.1) = c("date", "time", 'hr')
	dat.hr.1[,"hr"] = as.numeric(dat.hr.1[,"hr"])
	#head(dat.hr.1)
	
	dat.step.0 = matrix(0, nrow(dat.hr.1), ncol=4)
	colnames(dat.step.0) = c("date", "time", 'datetime', 'steps')
	dat.step.0[,'date'] = dat.hr.1[,'date']
	dat.step.0[,'time'] = dat.hr.1[,'time']
	dat.step.0[,'datetime'] = paste(dat.hr.1[,'date'], dat.hr.1[,'time'])
    head(dat.step.0)
	dat.step = read.csv(dir.step)
	dat.step = as.matrix(dat.step)
	head(dat.step)
	
	if (watch.type == 'apple') {
		for (i in 1:nrow(dat.step)) {
		    min.intv = seq(as.POSIXct(dat.step[i, 'start_datetime'], tz='GMT', tryFormats ="%Y-%m-%d %H:%M"), as.POSIXct(dat.step[i, 'end_datetime'], tz='GMT',  tryFormats ="%Y-%m-%d %H:%M"), by="1 min")
		    min.intv = gsub(' GMT','', min.intv)
		    dat.step.0[dat.step.0[,'datetime'] %in% min.intv, 'steps'] = as.numeric(dat.step[i, 'steps'])		
		}			
	}


    if (watch.type == 'fitbit') {
		dat.step.1 = dat.step[, c("datetime", "steps")]
		dat.step.1 = cbind(substring(dat.step.1[,1], 1, 10), substring(dat.step.1[,1], 12), dat.step.1[,1:2])
		colnames(dat.step.1) = c("date", "time", 'datetime', 'steps')
		tm.intersect = intersect(dat.step.1[, 'datetime'], dat.step.0[, 'datetime'])
		dat.step.0[match(tm.intersect, dat.step.0[, 'datetime']), 'steps'] = dat.step.1[match(tm.intersect, dat.step.1[, 'datetime']), 'steps']
		#head(dat.step.1)    	
    }


    return(list(dat.hr=dat.hr.1, dat.step=dat.step.0))
}

#### to obtain baseline residuals 
base.stand.res.fn = function(X.conv.f.base, resol.tm=60, test.r.fixed=FALSE, test.r=NULL, res.quan=0.9){
	cnt.tm = nrow(X.conv.f.base)/resol.tm	
	X.bin = sapply(1:cnt.tm, function(k) {
		X.sub = X.conv.f.base[(1+(k-1)*resol.tm):(k*resol.tm), ]
		c(X.sub)
	} )
	dim(X.bin)
	
	x.cen = apply(X.bin, 2, mean, na.rm=TRUE)
	x.sd = apply(X.bin, 2, sd, na.rm=TRUE)
	
	z.mx = matrix(NA, cnt.tm, ncol(X.conv.f.base))
	for (j in 1:ncol(X.conv.f.base)) {
		x = X.conv.f.base[,j]
		x.bin = matrix(x, resol.tm, cnt.tm)
		x.avg = colMeans(x.bin, na.rm=TRUE)
	    z.mx[,j] = (x.avg - x.cen)/x.sd		
	}
	
	z.abs.thres=4
	z.v = c(z.mx)
	z.v = pmin(z.v, z.abs.thres)
	z.v = pmax(z.v, - z.abs.thres)
	names(z.v) = paste(rep(colnames(X.conv.f.base), each=cnt.tm), rep(rownames(X.conv.f.base)[seq(1, nrow(X.conv.f.base), by=resol.tm)], ncol(X.conv.f.base)), sep=" ")
	
	if (test.r.fixed ==FALSE) {
		test.r = quantile(z.v, res.quan, na.rm=TRUE)/2
	}	
	cusum.v = cusum.fn(z.v, test.r)
	return(list(z.v=z.v, test.r=test.r, cusum.v=cusum.v))
}


#### to update residuals along with time
update.stand.res.fn = function(X.conv.f.base, x.conv.f, resol.tm=60, z.abs.thres=4){
	cnt.tm = nrow(X.conv.f.base)/resol.tm
	X.bin = sapply(1:cnt.tm, function(k) {
		X.sub = X.conv.f.base[(1+(k-1)*resol.tm):(k*resol.tm), ]
		c(X.sub)
	} )
	dim(X.bin)
	
	x.cen = apply(X.bin, 2, mean, na.rm=TRUE)
	x.sd = apply(X.bin, 2, sd, na.rm=TRUE)
	
	x.bin = matrix(x.conv.f, resol.tm, cnt.tm)
	x.avg = colMeans(x.bin, na.rm=TRUE)
	z.v = (x.avg - x.cen)/x.sd
	z.v = pmin(z.v, z.abs.thres)
	z.v = pmax(z.v, - z.abs.thres)		
	names(z.v) = paste(rep(substring(names(x.conv.f)[1], 1, 10), cnt.tm), rownames(X.conv.f.base)[seq(1, nrow(X.conv.f.base), by=resol.tm)], sep=" ")
	
	return(z.v)
}

### to obtain initial CuSum stats
cusum.fn = function(z.v, test.r=1) {
	z.t = z.v[!is.na(z.v)]
	day.diff = as.numeric(diff(as.Date(names(z.t))))
	day.sel = c(0, which(day.diff > 7), length(z.t))
	
	w.t = c()
	for(k in 1:(length(day.sel)-1)) {
		w.t.1 = as.numeric(day.sel[k+1] - day.sel[k])
		w.t.1[1] = 0
		z.t.1 = z.t[(day.sel[k]+1):day.sel[k+1]]
		for (i in 1:length(z.t.1)) {
			w.t.1[i+1] = max(0, w.t.1[i]+z.t.1[i]-test.r) 
		}
		w.t.1 = w.t.1[-1]
		names(w.t.1) = names(z.t.1)
		w.t = c(w.t, w.t.1)		
	}

	return(cusum.stats = w.t)
}



#### to update CuSum stats along with time
udpate.cusum.fn = function (x.res, cusum.t, test.r=1) {
	z.t = x.res[!is.na(x.res)]
	w.t = as.numeric(length(z.t))
			
	if (as.numeric(as.Date(names(z.t)[1]) - as.Date(names(cusum.t)[length(cusum.t)])  ) > 7 ) {
		w.t[1] = 0	
	} else {
		day.sel = sort(sort(unique(as.Date(names(cusum.t))), decreasing=TRUE)[1:3], decreasing=FALSE)
		avg = sapply(day.sel, function(ind) mean(cusum.t[as.Date(names(cusum.t)) == ind], na.rm=TRUE) )			
		if ( sum(diff(avg) < 0) == 2) {
			w.t[1] = 0
		} else {
		    w.t[1] = cusum.t[length(cusum.t)]			
		}
	}
	
	for (i in 1:length(z.t)) {
		w.t[i+1] = max(0, w.t[i]+z.t[i]-test.r) 
	}
	w.t = w.t[-1]
	names(w.t) = names(z.t)		

	return(cusum.stats = w.t)
}


#### to obtain emprical null distribution for test stats	
null.ecdf.fn = function(base.test) {
	 ind.0 = (1:(length(base.test)-1))[base.test[1:(length(base.test)-1)] == 0]
	 ind.0 = ind.0[base.test[ind.0 + 1] > 0]
	 max.step = max(diff(ind.0))
	 null.ecdf = list(1:max.step)			
	 for (step.k in 1:max.step) {
		 max.test = sapply(1:length(ind.0), function(j) max(base.test[ind.0[j]:(ind.0[j] + step.k)]))
		 max.test = max.test[!is.na(max.test)]
	     null.ecdf[[step.k]] = ecdf(max.test)
	 }
     return(null.ecdf)
 }
 
 
#### run index for CuSum stats
base.test.ind.fn = function(base.test) {
	 base.test.ind = rep(0, length(base.test))
	 base.test.ind[1] = ifelse(base.test[1] == 0, 0, 1)
	 for (i in 2:length(base.test.ind)) {
	 	 if (base.test[i] > 0) {
	 	 	 base.test.ind[i] = base.test.ind[i-1] + 1
	 	 }
	 }
	 names(base.test.ind) = names(base.test)
	 return(base.test.ind)			 	
}			

#### to obtain residuals and CuSum stats compared to sliding window baseline

stats.track.fn = function (id, X.rhr.1, base.num=28, test.r.fixed.1=FALSE, test.r.1=NA, res.quan.1=0.9, pval.thres=0.01, resol.tm.1=60, red.alarm.thres=24) {
	
	# to obtain initial baseline residuals, CuSum stats, null distribution and p-value
	X.conv.f.base = X.rhr.1[, 1:base.num]
	base.info = base.stand.res.fn(X.conv.f.base, resol.tm.1, test.r.fixed.1, test.r=test.r.1, res.quan.1)
    test.r.base = base.info$test.r
	base.res = base.info$z.v		
    test.r.seq = rep(test.r.base, length(base.res))
    names(test.r.seq) = names(base.res)   	
	base.test = base.info$cusum.v
	null.ecdf.test = null.ecdf.fn(base.test)			
    base.test.ind = base.test.ind.fn(base.test)
   
    base.days = colnames(X.conv.f.base)
    base.days.0 = base.days
	day.test = colnames(X.rhr.1)[(base.num+1):ncol(X.rhr.1)]
	res.t = base.res
	
	score.thres = 1 - pval.thres
	test.t = base.test
	test.t.ind = base.test.ind
	risk.score = rep(NA, length(base.test))
	names(risk.score) = base.test
	for (j in 1:length(base.test.ind)) {
		if (base.test.ind[j] >= 2 & base.test.ind[j] <= length(null.ecdf.test)) {
			risk.score[j] =  null.ecdf.test[[base.test.ind[j]]](base.test[j] - 0.001)
		}
		if (base.test.ind[j] > length(null.ecdf.test)) {
			risk.score[j] = 1
		}				
	}
	names(risk.score) = names(base.test)	
	
	day.score.seq = rep(0, length(base.days))
	names(day.score.seq) = base.days
	track.tab = c()	
	alarm.day = c()
	# to obtain residuals, CuSum stats, and p-value along with time
	for (k in 1:length(day.test)) { 
		day.1 = day.test[k]	
		x.conv.f =  X.rhr.1[, day.1] 
		names(x.conv.f) = paste(rep(day.1, length(x.conv.f)), names(x.conv.f), sep=" ")
			
		# to obtain standardized residuals based on baseline
	    x.res = update.stand.res.fn(X.conv.f.base, x.conv.f, resol.tm=resol.tm.1, z.abs.thres=4)

        if (sum(!is.na(x.res)) >= 2 ) { 
            res.t = c(res.t, x.res)
	        x.test.r = rep(test.r.base, length(x.res))
	        names(x.test.r) = names(x.res)
            test.r.seq = c(test.r.seq, x.test.r)
            
			# to obtain cusum stats
 			x.test = udpate.cusum.fn(x.res, cusum.t=test.t, test.r=test.r.base)			
			test.t = c(test.t, x.test)
			
			x.test.ind = rep(0, length(x.test))
			names(x.test.ind) = names(x.test)
			x.test.ind[1] = ifelse(x.test[1] == 0, 0, 1)
			if (test.t.ind[length(test.t.ind)] > 0 & x.test.ind[1] == 1) {
				x.test.ind[1] = 1 + test.t.ind[length(test.t.ind)]
			}
			for (i in 2:length(x.test)) {
				if (x.test[i] > 0) {
				  x.test.ind[i]	= x.test.ind[i-1] + 1					
				}				
			}
			test.t.ind = c(test.t.ind, x.test.ind)

			x.risk.score = rep(NA, length(x.test))
			names(x.risk.score) = names(x.test)
			for (j in 1:length(x.test.ind)) {
				if (x.test.ind[j] >= 2 & x.test.ind[j] <= length(null.ecdf.test)) {						x.risk.score[j] = null.ecdf.test[[x.test.ind[j]]]( test.t[(1:length(test.t))[names(test.t.ind) == names(x.test[j])] ]  - 0.001)
				}
				if (x.test.ind[j] > length(null.ecdf.test)) {
					x.risk.score[j] = 1
				}				
			}
			risk.score = c(risk.score, x.risk.score)	
									
						
			if (sum(x.risk.score > score.thres, na.rm=TRUE) == 0) {
                
                day.score = 0
                names(day.score) = day.1
                day.score.seq = c(day.score.seq, day.score)
                		    			    			       
			} else {
				day.score = 1
                names(day.score) = day.1
				x.latest.run.ind = test.t.ind[(length(test.t.ind) -x.test.ind[length(x.test.ind)]+1):length(test.t.ind)]
				x.test.current.avg = mean(x.test[intersect(names(x.latest.run.ind), names(x.test))], na.rm=TRUE)	            
				day.dum = unique(substring(names(test.t.ind), 1, 10))
			    past.day = day.dum[length(day.dum)-1]
			    past.tm = names(x.latest.run.ind[as.Date(names(x.latest.run.ind)) == past.day])
			    if (length(past.tm) > 0) {
				    x.test.past.avg = mean(test.t[past.tm], na.rm=TRUE)
				    if (x.test.ind[length(x.test.ind)] >= red.alarm.thres) {
				    	if (x.test.current.avg > x.test.past.avg) {
				    		day.score = 2
				    		names(day.score) = day.1
				    		
				    		if (day.score.seq[length(day.score.seq)] == 2) {
				    			track.tab.dum = cbind(rep(id, length(x.latest.run.ind)), names(x.latest.run.ind), x.latest.run.ind, round(test.t[names(x.latest.run.ind)],3), round(1-risk.score[names(x.latest.run.ind)],3) )
				    			colnames(track.tab.dum) = c("ParticipantID", "track alarming time", "alarming hour",  "CuSum statistics at alarm", "p-value")
				    			rownames(track.tab.dum) = NULL
				    			track.tab = rbind(track.tab, track.tab.dum)
				    		    track.tab = track.tab[!duplicated(track.tab), ]
				    			alarm.day = c(alarm.day, day.1)
				    		}
				    	}
				    }
			    }
			    day.score.seq = c(day.score.seq, day.score) 
			   			     			    

			} 	
			
				if (day.score < 1) {
				    base.days = c(base.days, day.1)
				    X.conv.f.base = cbind(X.conv.f.base, x.conv.f)
				    colnames(X.conv.f.base)[ncol(X.conv.f.base)] = day.1
				    X.conv.f.base = X.conv.f.base[, max(1, ncol(X.conv.f.base) - base.num+1):ncol(X.conv.f.base)]
					base.info = base.stand.res.fn(X.conv.f.base, resol.tm.1, test.r.fixed.1, test.r=test.r.1, res.quan.1)
	
					base.res = base.info$z.v
					base.test = base.info$cusum.v
					test.r.base = base.info$test.r				
					null.ecdf.test = null.ecdf.fn(base.test)						
				    base.test.ind = base.test.ind.fn(base.test)
			    }
					
		} 
	
	}
		
	return(list(base.days.0=base.days.0, res.t = res.t, test.t=test.t, test.t.ind=test.t.ind, test.pval=1-risk.score, test.r.seq= test.r.seq, day.score.seq=day.score.seq, alarm.day=alarm.day, track.tab=track.tab))
}



#########################################
####    RHR-Diff offline detection   ####
#########################################
       
# reference: https://github.com/mwgrassgreen/RankScan
# paper: Arias-Castro, Ery, Rui M Castro, Ervin Tánczos, and Meng Wang. 2018. "Distribution-free detection of structured anomalies: Permutation and rank-based scans." Journal of the American Statistical Association 113 (522): 789-801.


#### rank scan detection
rhr.diff.detection.fn = function (id, res.t, alpha=0.05, B=1000) {
	
	res.t[is.na(res.t)] = 0
	nm.dum = seq(as.Date(names(res.t)[1]), as.Date(names(res.t)[length(res.t)]), by=1)
    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
    res.t.1 = rep(0, length(nm.dum))
    names(res.t.1) = nm.dum
    res.t.1[names(res.t)] = res.t
    
	z.v = res.t.1
	n = length(z.v)
	Q = floor(log2(n))
	null.result = crt.val.rank(n, Q, alpha, B) 
	d = null.result$crt.val
	null.B = null.result$test.null
	#d = 43120.513817 # n=2^15
	a = rank(z.v)
	names(a) = NULL
	q = floor(log2(length(z.v)))
	result = rank_scan_est(a, q, d)
	pval = sapply(abs(result[,3]), function(x) sum(null.B >= x)/(1+B))	
	detect = cbind(names(z.v)[result[,1]], names(z.v)[result[,2]], round(result[,3],3), round(pval, 3))
	if (ncol(detect) < 3) {
		detect.pos = matrix(c(id, rep("N/A", 4)), ncol=5)
	} else {
		detect = detect[!duplicated(detect[,3]), ]
		detect = matrix(detect, ncol=4)	
        hr.diff = difftime(detect[,2], detect[,1], unit="hours")
		detect = detect[as.numeric(hr.diff) >= 24, ] # to remove the detected interval < 24 hours
		detect = matrix(detect, ncol=4)
		if (sum(detect[,3] > 0) == 0) {
			detect.pos = matrix(c(id, rep("N/A", 4)), ncol=5)
		} else {
			detect.pos = matrix(detect[detect[,3] > 0, ],ncol=4)
			detect.pos = matrix(cbind(rep(id, nrow(detect.pos)), detect.pos ), ncol=5)			
		}
	}
	colnames(detect.pos) = c("ParticipantID", "staring time of detected interval", "ending time of detected interval", "rank scan test statistics", "p-value")
	return(detect.pos)
}



# to get critival value for rank scan test
# only scan intervals in length 2^(1:(log2(N)-1))

GetMaxAtk <- function(X, mu.X, k){
            # k--the length of signal interval candidate	
	        N <- length(X)
	        Intv <- c(rep(1, k), rep(0, N-k))
	        return(max(convolve(X - mu.X, Intv)[1:(N-k+1)])/sqrt(k))
}

GetRankScan <- function(Rank.X, Q){
	        N <- length(Rank.X)
	        # scan intervals in dyalic lengths from 2 to [N/2]
	        #Q <- floor(log2(N))
	        max.at.k <- Q
	        for (q in 1:(Q-1)){
	        	max.at.k[q] <- GetMaxAtk(Rank.X, (N+1)/2, 2^q)
	        }
	        stats <- max(max.at.k)
            return(stats)
}


crt.val.rank = function (N, Q, alpha, B) {
	rank.scan.null <- rep(NA, B)
	for (b in 1:B){
		Rank.X.null <- sample(1:N)			
		rank.scan.null[b] <- GetRankScan(Rank.X.null, Q)
	}
	
	crt.val <- quantile(rank.scan.null, 1 - alpha)
    return(list(crt.val=crt.val, test.null=rank.scan.null))
}



#### to get identified anomalous intervals from rank scans

rank_scan_est = function(a, q, d){
	
	# a is the rank sequence of the orignial data X
	# q is the log_2 of maximum length of the candidate intervals
	# d is the threshold for normalized sum of ranks in one interval i.e., 1/sqrt(|S|) sum_{v in S} (R_v - (N+1)/2) in the notation of (Arias-Castro et al 2018)
	# reference for reporting the identified intervals: Jeng, X. J., Cai, T. T., and Li, H. (2010), “Optimal Sparse Segment Identification With Application in Copy Number Variation Analysis” 

	
	b = length(a)
	x = rep(0, b*(2^q+1)) # store all the intervals in length 2^(1:q)
	for (i in 1:(b-1)) {
		#print(i)
		# i is the starting index of a candidate interval
		for (j in pmin(b, i+2^(1:q)-1)) {
		    # j is its ending index
		    #print(j)
			x[(i-1)*2^q + j] = sum(a[i:j] - (b+1)/2)/sqrt(j-i+1)
			#print(x[(i-1)*2^q + j])
		}
	}
	
	
	k= which(abs(x)>d);
	i = ceiling(k/(2^q+1)); 
	j = k - (i-1)*2^q;
	list = cbind(i,j, x[k]);
	
	start = rep(0,1);
	end = rep(0,1);
	Rank_scan = rep(0,1);
	t=1;
	
	while (length(list)> 3) {
	ind = which(abs(list[,3]) == max(abs(list[,3])));
	len.ind = length(ind)
	start[t:(t+len.ind-1)] = list[ind,1];
	end[t:(t+len.ind-1)] = list[ind,2];
	Rank_scan[t:(t+len.ind-1)] = list[ind,3];
	II = c()
	for (l in 1:len.ind) {
		s = t+l-1
		II = c(II, which(list[,1]<=end[s] & list[,2]>=start[s]))
	}
	#II = which(list[,1]<=end[t] & list[,2]>=start[t]);
	list = list[-II,];
	t = t+len.ind; 
	} 
	
	if(length(list)==3) {
	start[t] = list[1];
	end[t] = list[2];
	Rank_scan[t] = list[3];
	}
	
	peaks = cbind(start, end, Rank_scan); 
	return(peaks)
}



#########################################
####    plot of detection result     ####
#########################################

result.plot.fn = function (id, base.days.int0, res.t, cusum.t, cusum.test.r, offline.result, alarm.cate.tm, day.info, shift.day) {
	nm = names(day.info)
	if(!is.null(day.info)) {
      day.info = as.character(as.Date(as.matrix(day.info)) - shift.day)
    }
           
    par(mfrow=c(2,1))
    par(mai=c(1.5, 1, 1, 0.1))		
    
    nm.dum = seq(as.Date(names(res.t)[1]), as.Date(names(res.t)[length(res.t)]), by=1)
    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
    res.t.1 = rep(0, length(nm.dum))
    names(res.t.1) = nm.dum
    res.t.1[names(res.t)] = res.t
	plot(res.t.1, type="l", xaxt="n", xlab="", ylab="residuals", cex.lab=2, cex.main=2, main=paste("offline detection for", id), ylim=c(-4, 4)) 
	pos = (1:length(res.t.1))[!duplicated(substring(names(res.t.1), 1, 10))]
	names(pos) = unique(substring(names(res.t.1), 1, 10))
    test.r.thres =  cusum.test.r
    test.r.thres.1 =  test.r.thres[names(res.t.1)]
    plot(stepfun(1:(length(res.t.1)), c(test.r.thres.1,test.r.thres.1[length(test.r.thres.1)] )), pch=".", col='darkgreen', add=TRUE,lwd=2)
    abline(h=c(0),lty=c(2), lwd=2, col='darkgreen')  
	abline(v=(1:length(res.t.1))[!duplicated(substring(names(res.t.1), 1, 10))], lty=2, col="grey")	
	line.pos=-2
	z.v = res.t.1
	detect = offline.result[, c("staring time of detected interval", "ending time of detected interval")]
    detect = matrix(detect, ncol=2)
	if (!is.null(detect) & detect[1,1] != "N/A"){
		for (i in 1:nrow(detect)) {
				arrows((1:length(z.v))[names(z.v) == detect[i,1] ], line.pos, (1:length(z.v))[names(z.v) == detect[i,2] ], line.pos, col="red", lwd=2, length=0.1)
		        arrows((1:length(z.v))[names(z.v) == detect[i,2] ], line.pos, (1:length(z.v))[names(z.v) == detect[i,1] ], line.pos, col="red", lwd=2, length=0.1)		

	    }	
	}
	
	names(pos) = as.character(as.Date(names(pos)) - shift.day)
	axis(1, pos,  names(pos), las=2, cex.axis=1)
    if (!is.null(day.info)) {
       if (sum(names(pos) %in% day.info[1]) > 0) {
         axis(1, pos[day.info[1]],  day.info[1], col.axis="red", las=2, cex.axis=1)       	
       }
       if (sum(names(pos) %in% day.info[2]) > 0) {
         axis(1, pos[day.info[2]],  day.info[2], col.axis="purple", las=2, cex.axis=1)		       	
       }	
	}


    ##### 
    nm.dum = seq(as.Date(names(cusum.t)[1]), as.Date(names(cusum.t)[length(cusum.t)]), by=1)
    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
    cusum.t.1 = rep(NA, length(nm.dum))
    names(cusum.t.1) = nm.dum
    cusum.t.1[names(cusum.t)] = cusum.t
    
	test.t = cusum.t.1
	test.alarming.tm = alarm.cate.tm
	plot(test.t, type="l", xaxt="n", xlab="", ylab="CuSum statistics", cex.lab=2, cex.main=2, main=paste("online detection for", id))
	points(test.t, pch=".", cex=3)
	abline(v=(1:length(test.t))[!duplicated(substring(names(test.t), 1, 10))], lty=2, col="grey")
	abline(v=(1:length(test.t))[names(test.t) %in% names(test.alarming.tm[test.alarming.tm == 1])], lty=2, col="yellow")
	abline(v=(1:length(test.t))[names(test.t) %in% names(test.alarming.tm[test.alarming.tm == 2])], lty=1, col="red")
	pos = (1:length(test.t))[!duplicated(substring(names(test.t), 1, 10))]
	names(pos) = unique(substring(names(test.t), 1, 10))
    names(pos) = as.character(as.Date(names(pos)) - shift.day)
	axis(1, pos,  names(pos), las=2, cex.axis=1)
    if (!is.null(day.info)) {
       if (sum(names(pos) %in% day.info[1]) > 0) {
         axis(1, pos[day.info[1]],  day.info[1], col.axis="red", las=2, cex.axis=1)       	
       }
       if (sum(names(pos) %in% day.info[2]) > 0) {
         axis(1, pos[day.info[2]],  day.info[2], col.axis="purple", las=2, cex.axis=1)		       	
       }	
	}


}

