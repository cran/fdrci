fdr_od <-
function(obsp,permp,pnm,ntests,thres,cl=.95,c1=NA){
      z_ = qnorm(1 - (1 - cl)/2) # two-tailed test
	pcount = rep(NA,length(permp))
	for(p_ in 1:length(permp)){
	   permp[[p_]][,pnm] = ifelse(permp[[p_]][,pnm]<thres,1,0)
	   pcount[p_] = sum(permp[[p_]][,pnm],na.rm=TRUE)
	}
	# over-dispersion parameter is observed variance of p in permuted data / expected
	p = mean(pcount,na.rm=TRUE)/ntests # estimate p
	e_vr = ntests*p*(1 - p)
	o_vr = var(pcount,na.rm=TRUE)
	if( is.na( c1 ) ) {
	   c1 = o_vr / e_vr
	   if( !is.na( c1 ) ) if( c1 < 1) c1 = 1
	} 

	nperm = length(permp)
	mo = ntests
	ro = sum(obsp < thres)
	vp = sum(pcount)
	vp1 = vp
  	rslt = rep(NA,4)
  	if(ro > 0){
		if(vp == 0) vp = 1
		mean.vp = vp / nperm
        	fdr0 = mean.vp / ro
	num = mo - ro
	denom = mo - (vp/nperm)
	if( num >= denom ){
		pi0 = 1
	} else {
        		pi0 = (mo - ro)/(mo - (vp/nperm))
	}
        	if( pi0 < 0.5 ) pi0 = 0.5    # updated 10/14/16 to limit influence of pi0
        	fdr = fdr0 * pi0    # updated calculation of fdr to be robust to ro = mtests
		
		 # variance of FDR
        	mp = nperm * mo
        	t1 = 1 / vp
        	denom = mp - vp
        	t2 = 1 / denom
        	t3 = 1 / ro
        	denom = ntests - ro
        	if( denom < 1 ) denom = 1  # updated 10/14/16 to avoid inf t4
        	t4 = 1 /  denom
        	s2fdr = (t1 + t2 + t3 + t4) * c1
        	ul = exp(log(fdr) + z_ * sqrt(s2fdr))
        	ll = exp(log(fdr) - z_ * sqrt(s2fdr))
	
		rslt = c(fdr,ll,ul,pi0)
		rslt = ifelse(rslt > 1, 1, rslt) # FDR > 1 does not make sense, thus set to 1 in this case
		rslt = c(rslt,c1,ro,vp1)
	} 
	return(rslt)
}

