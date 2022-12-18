
# simulate some data
mu0=100
sigma0=15
delta=2*sigma0
theta=2
d=simuleer1(mu0,sigma0,delta,theta,n0=50,n1=100)   # 2 steekproeven (cases en controls)
rm(mu0,sigma0,delta,theta)

# visualize the data
#par(mfrow=c(2,1))
dev.new()
dotplot(d)
dev.new()
histogramplot(d, z=5)

# NONPARAMETRIC ANALYSIS
nonparres=nonparametricanalysis(d, bootstrap=T, nbootstraps=10000, cdplot=F, Rplot=F, Kplot=F, Bplots=F) 
#nonparres$AUC
#rm(nonparres)

# BAYESIAN ANALYSIS
bayesres=bayesiananalysis(d, thinning=10, updates=10000, iterations=10000, nchains=4)
#bayesres$samenvatting        # may give a large table
#plot(bayesres$mcmcresults)   # check convergence, may give a lot of graphs !!!
bayesres1=bayesianresults(bayesres, cumdist=TRUE, kolmogorov=FALSE)
#bayesres1$AUC
qq2=bayesres$mcmcresult[[1]][,which(colnames(bayesres$mcmcresult[[1]])=="auc")]
for (j in 2:length(bayesres$mcmcresult)) {qq2=c(qq2,bayesres$mcmcresult[[j]][,which(colnames(bayesres$mcmcresult[[1]])=="auc")])}
dev.new()
hist(qq2, xlab="AUC", main="")
rm(qq2,j)
#rm(bayesres1,qq2,j)

# PARAMETRIC ANALYSIS, compare models
analyse2(d)
analyse3(d, bootstrap=F, popsampling=T, nbootstraps=1000, model="weibull")
stats=steekproevenverdeling3(d, niter=10000, model="gaussian")
dev.new()
hist(stats[,1], xlab="estimated optimal cutoff point", main="sampling distribution", sub="means & SDs sampled from distribution with estimated means,SDs & SEs", cex.sub=0.7)
rm(stats)


#rm(d, ff4, bayesres)


#########################################################################################################################################################################################


logit=function(x){log(x/(1-x))}                         # the logit-transform of x (scalar/vector/matrix)
expit=function(x){1/(1+exp(-x))}                        # the anti-logit-transform of x (scalar/vector/matrix)
simuleer1 = function(mu0,sigma0,delta,theta,n0,n1) {    # simulate two samples of sizes n0, n1 from two normal distributions with means (mu0,mu0+delta) and SD (sigma0,sigma0*theta)
   mu1=(mu0+delta)
   sigma1=sigma0*theta
   y0=sort(rnorm(n0,mu0,sigma0))
   y1=sort(rnorm(n1,mu1,sigma1))
   return(list(y0=y0, y1=y1))                           # return a list with the two samples of values
}
empfunxie=function(x,y) {                               # calculate the empirical cumulative distribution of a sample of values y at a vector of values x
	p =sapply(x,function(x,y){sum(y<=x)/length(y)},y)
	nn=sapply(x,function(x,y){sum(y==x)},y)
    logitp=log(p/(1-p))
    selogitp = sqrt(1/(length(y)*p*(1-p)))
    return(cbind(p, logitp, selogitp, nn))              # return a matrix with nrows=length(x), and four columns: (i) the cumulative probabilities, (ii) the associated logit-transforms, 
}                                                       # (iii) the standard errors of the logit-transformed cumulative probabilities and (iv) the number of values in y at the various x-values
analyse1=function(d) {     # nonparametric analysis of two samples of the same biomarker, presented as a list with two vectors y0 and y1
   ksi=min(res$y0,res$y1)/10
   unieke_y=c(min(d$y0,d$y1)-ksi,sort(unique(c(d$y0,d$y1))),max(d$y0,d$y1)+ksi)   # determine the unique biomarker-values in the 2 samples
   G=empfunxie(unieke_y, d$y1)                                                    # determine both empirical cumulative distribution functions
   H=empfunxie(unieke_y, d$y0)
   for (i in 2:nrow(G)) {                                                         # if the cumulative probabilities equal zero or 1, impute in their standard errors the values 
      if(H[i,3]==Inf) {H[i,3]=H[(i-1),3]}                                         #   of the closet probabilities unequal to zero or 1
      if(G[i,3]==Inf) {G[i,3]=G[(i-1),3]}
   }
   for (i in (nrow(G)-1):1) {
      if(H[i,3]==Inf) {H[i,3]=H[(i+1),3]}
      if(G[i,3]==Inf) {G[i,3]=G[(i+1),3]}
   }
   dif=H[,1]-G[,1]                                                                # calculate the differences between the two cumulative distribution functions at all unique biomarker-values
   sedif=sqrt(H[,1]*(1-H[,1])/length(d$y0) + G[,1]*(1-G[,1])/length(d$y1))        #  and the associated standard errors
   ff1=1-H[,1]         # 1-specificity to calculate AUC
   ff2=1-G[,1]         # sensitivity
   rorde=order(ff1,ff2)   # sort with respect to 1-specificity and sensitivity
   ff1=ff1[rorde]
   ff2=ff2[rorde]
   auc=sum((ff1[2:length(ff1)]-ff1[1:(length(ff1)-1)]) * ff2[1:(length(ff2)-1)] + (ff1[2:length(ff1)]-ff1[1:(length(ff1)-1)])*(ff2[2:length(ff2)]-ff2[1:(length(ff1)-1)])/2)
   seauc=sqrt((auc*(1-auc) + (length(d$y1)-1)*((auc/(2-auc))-auc^2) + (length(d$y0)-1)*((2*auc^2/(1+auc))-auc^2)) / (length(d$y1)*length(d$y0)))
   xx5=c(auc,seauc)
   names(xx5)=c("AUC", "SE")
   return(list(unieke_y=unieke_y, G=G, H=H, dif=dif, sedif=sedif, n0=length(d$y0), n1=length(d$y1), y0=d$y0, y1=d$y1, auc=xx5)) # return an extended list with the original data, the unique values, 
}                                                                                                                               #  the empirical cum.dists, differences, AUC and SEs
histogramplot = function(res, z=10) {           # a histogram of the two samples with mirrored y-axis for the two samples
   ksi=min(res$y0,res$y1)/10
   #z=5
   unieke_y=c(min(res$y0,res$y1)-ksi,sort(unique(c(res$y0,d$y1))),max(res$y0,res$y1)+ksi)
   ff1=min(unieke_y)+(0:z)*(max(unieke_y) - min(unieke_y))/z
   ff2=sapply(ff1[-1],function(x,yy){sum(yy<=x)},yy=res$y0)
   ff3=sapply(ff1[-1],function(x,yy){sum(yy<=x)},yy=res$y1)
   ff2[2:z]=ff2[2:z]-ff2[1:(z-1)]
   ff3[2:z]=ff3[2:z]-ff3[1:(z-1)]
   plot(ff2, type="h", lwd=10, col=3, ylim=c(-max(c(ff2,ff3)),max(c(ff2,ff3))), xaxt="n", xlab="biomarker y", yaxt="n", ylab="frequency", lend=1)
   points(-ff3, type="h", lwd=10, col=2, lend=1)
   abline(h=0, col=1, lty=1)
   axis(1, at=1:z, labels=paste(round(ff1[1:z],0),"-",round(ff1[2:(z+1)],0), sep=""))
   axis(2, at=seq(-z*(1+trunc(max(c(ff2,ff3))/z)),z*(1+trunc(max(c(ff2,ff3))/z)),z), labels=abs(seq(-z*(1+trunc(max(c(ff2,ff3))/z)),z*(1+trunc(max(c(ff2,ff3))/z)),z)))
   legend(x=8, y=max(c(ff2,ff3)), legend=c("controls", "cases"), pch=16, lwd=10, col=c(3,2), bty="n")
}
dotplot=function(d) {   # a simple dot-plot of the two samples
   plot(rep(0,length(d$y0))+rnorm(length(d$y0),0,0.05), d$y0, xaxt="n", xlab="", ylab="biomarker y", xlim=c(-0.5,1.5), ylim=c(min(d$y0,d$y1),max(d$y0,d$y1)), pch=16, col=3)
   points(rep(1,length(d$y1))+rnorm(length(d$y1),0,0.05), d$y1, pch=16, col=2)
   axis(1, at=c(0,1), labels=c("controls","cases"))
}
cumdistplot = function(res, plot=TRUE) {   # only called from nonparametricanalysis: calculate Youden-index cutoff-point with sens/spec and max(difference)
   if (plot==TRUE) {                       # plot the cumulative distributions, and their 95% CIs (calculated with SEs from analysis1)
      dev.new()
      plot (res$unieke_y, res$H[,1], type="s", lty=1, lwd=2, col=3, xlab="biomarker y", ylab="cumulative distribution")
      lines(res$unieke_y, 1/(1+exp(-(res$H[,2]-1.96*res$H[,3]))), type="s", lty=3, lwd=1, col=3)
      lines(res$unieke_y, 1/(1+exp(-(res$H[,2]+1.96*res$H[,3]))), type="s", lty=3, lwd=1, col=3)
      lines(res$unieke_y, res$G[,1], type="s", lty=1, lwd=2, col=2)
      lines(res$unieke_y, 1/(1+exp(-(res$G[,2]-1.96*res$G[,3]))), type="s", lty=2, lwd=1, col=2)
      lines(res$unieke_y, 1/(1+exp(-(res$G[,2]+1.96*res$G[,3]))), type="s", lty=2, lwd=1, col=2)
      abline(v=mean(res$unieke_y[which(res$dif==max(res$dif))]), lty=2, col=1)
      legend(x=min(res$unieke_y)-1, y=1, legend=c("controls","cases"), lwd=2, col=c(3,2), bty="n")
      text(x=3+median(res$unieke_y[which(res$dif==max(res$dif))]), y=0.05,paste("largest difference = ",round(max(res$dif),4)), adj=0)
      text(x=3+median(res$unieke_y[which(res$dif==max(res$dif))]), y=0,   paste("optimal cutoff y_c = ",round(median(res$unieke_y[which(res$dif==max(res$dif))]),0)), adj=0)
   }
   # best cut-off value y_c of the biomarker Y with associated specificity and sensitivity
   return(list(optimal_cutoff_point=mean(res$unieke_y[which(res$dif==max(res$dif))]), maxdif=max(res$dif),
               specificity=res$H[which(res$dif==max(res$dif)),1], sensitivity=1-res$G[which(res$dif==max(res$dif)),1]))
}
rocplot=function(res, plot=FALSE) {   # only called from nonparametricanalysis: calculate Qpoint, difference at Qpoint and sens/spec
   ff1=1-res$H[,1]
   ff2=1-res$G[,1]
   rorde=order(ff1,ff2)
   ff1=ff1[rorde]
   ff2=ff2[rorde]
   ddd=sqrt((ff1-0)^2 + (ff2-1)^2)
   if (plot==TRUE) {     # plot the roc-curve with 95% CIs
      dev.new()
      plot(ff1, ff2, type="s", xlab="1 - specificity", ylab="sensitivity", xlim=c(0,1), ylim=c(0,1), lwd=2, lty=1, col=1)
      points(ff1[which(ddd==min(ddd))],ff2[which(ddd==min(ddd))],col=2, pch=16, lwd=6)
      abline(a=0, b=1, lty=2, lwd=1, col=1)
      text(x=ff1[which(ddd==min(ddd))], y=ff2[which(ddd==min(ddd))]+0.05, "Q-point", srt=45)
      text(x=0.5, y=0.20, paste("AUC = ", round(res$auc[1],4),"(SE = ", round(res$auc[2],4),")",sep=""))
      text(x=0.5, y=0.15, paste("biomarker cut-off at Q-point =",round(res$unieke_y[rorde][which(ddd==min(ddd))],0),"with"))
      text(x=0.5, y=0.10, paste("sensitivity =", round(ff2[which(ddd==min(ddd))],2), "and specificity =", 1-ff1[which(ddd==min(ddd))]))
      lines(ff1, ff2+1.96*sqrt(ff2*(1-ff2)/length(res$y1)), type="s", lty=2, lwd=1,col=1)
      lines(ff1, ff2-1.96*sqrt(ff2*(1-ff2)/length(res$y1)), type="s", lty=2, lwd=1,col=1)
      lines(ff1+1.96*sqrt(ff1*(1-ff1)/length(res$y0)), ff2, type="s", lty=3, lwd=1,col=5)
      lines(ff1-1.96*sqrt(ff1*(1-ff1)/length(res$y0)), ff2, type="s", lty=3, lwd=1,col=5)
      legend(x=0.5, y=0.075, legend=c("95% CI given specificity", "95% CI given sensitivity"), col=c(1,5), lty=c(2,3), bty="n", cex=0.7)
   }
   # q-point and associated specificity and sensitivity
   return(list(optimal_cutoff_point=res$unieke_y[rorde][which(ddd==min(ddd))], cdistdif=res$dif[which(ddd==min(ddd))],    # mediaan van de min(...)'s??
               specificity=1-ff1[which(ddd==min(ddd))], sensitivity=ff2[which(ddd==min(ddd))]))
}
kolmogorovplot = function(res, plot=FALSE) {
   difference= res$dif - max(res$dif)
   v=which(res$dif==max(res$dif))[1]
   Var_v = res$H[v,1]*(1-res$H[v,1])/res$n0 + res$G[v,1]*(1-res$G[v,1])/res$n1
   vvv=c()
   for (jj in 1:length(difference)) {
      vvv[jj]=res$H[jj,1] * (1-res$H[jj,1])/res$n0+res$G[jj,1] * (1-res$G[jj,1])/res$n1 + 
                 2*(abs(res$H[jj,1]-res$H[v,1])*min(res$H[jj,1],res$H[v,1])/res$n0+abs(res$G[jj,1]-res$G[v,1])*min(res$G[jj,1],res$G[v,1])/res$n1)
   }
   tt1=which(difference+1.96*sqrt(vvv)<0)
   if (plot==TRUE) {
      dev.new()
      plot(res$unieke_y, res$dif, xlab="biomarker y", ylab="Kolmogorov distance", type="b", lty=1, col=1, lwd=1, pch=1, cex=0.5,
            ylim=c(min(max(res$dif)+difference-1.96*sqrt(vvv)), max(max(res$dif)+1.96*sqrt(vvv))))
      lines(res$unieke_y, res$dif, col=1, lty=1, lwd=1)
      abline(h=max(res$dif) - 1.96*res$sedif[res$dif==max(res$dif)], lty=3)
      abline(h=0, lty=1)
      lines(res$unieke_y, max(res$dif)+difference, lty=2, lwd=1,col=1)
      lines(res$unieke_y, max(res$dif)+difference-1.96*sqrt(vvv), lty=2, lwd=1,col=1)
      lines(res$unieke_y, max(res$dif)+difference+1.96*sqrt(vvv), lty=2, lwd=1,col=1)
      abline(v=res$unieke_y[max(which(which(difference+1.96*sqrt(vvv)<0) < v))+1], col=1, lty=3)
      abline(v=res$unieke_y[(min(tt1[which(tt1>v)])-1)], col=1, lty=3)
   }
   return(list(optimal_cutoff_point=res$unieke_y[v],
          lower_limit_optimal_cutoff_point=res$unieke_y[max(which(which(difference+1.96*sqrt(vvv)<0) < v))+1], 
          upper_limit_optimal_cutoff_point=res$unieke_y[(min(tt1[which(tt1>v)])-1)] )  )
}
eenbootstrapfunctie = function(d, nbootstraps)  {
   ycc1=maxdif=sens1=spec1=ycc2=mdif2=sens2=spec2=aucc=c()
   for (ii in 1:nbootstraps) {
      nrs0=nrs1=yy0=yy1=unieke_yy=GG=HH=diff=fff1=fff2=rordee=NA
      nrs0=sample(1:length(d$y0), length(d$y0), replace=TRUE)
      nrs1=sample(1:length(d$y1), length(d$y1), replace=TRUE)
      yy0=d$y0[nrs0]
      yy1=d$y1[nrs1]
      unieke_yy=c(min(yy0,yy1)-1,sort(unique(yy0,yy1)),max(yy0,yy1)+1)
      GG=empfunxie(unieke_yy, yy1)
      HH=empfunxie(unieke_yy, yy0)
      diff=abs(HH[,1]-GG[,1])
      ddd=sqrt((1-HH[,1]-0)^2 + (GG[,1])^2)
      ycc1[ii]=mean(unieke_yy[which(diff==max(diff))])        # mediaan?
      maxdif[ii]=(HH[,1] - GG[,1])[which(diff==max(diff))]    # mediaan over de max?
      sens1[ii]=1-GG[which(diff==max(diff)),1]
      spec1[ii]=HH[which(diff==max(diff)),1]
      ycc2[ii]=mean(unieke_yy[which(ddd==min(ddd))])
      mdif2[ii]=diff[which(ddd==min(ddd))]
      sens2[ii]=1-GG[which(ddd==min(ddd)),1]
      spec2[ii]=HH[which(ddd==min(ddd)),1]
      fff1=1-HH[,1]
      fff2=1-GG[,1]
      rordee=order(fff1,fff2)
      fff1=fff1[rordee]
      fff2=fff2[rordee]
      aucc[ii]=sum((fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)]) * fff2[1:(length(fff2)-1)] + (fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)])*(fff2[2:length(fff2)]-fff2[1:(length(fff1)-1)])/2)
   }
   rm(nrs0,nrs1,yy0,yy1,unieke_yy,GG,HH,diff,fff1,fff2,rordee,ii)
   return(cbind(ycc1,maxdif,sens1,spec1,ycc2,mdif2,sens2,spec2,aucc))
}
nonparametricanalysis = function(d, bootstrap=T, nbootstraps=10000, cdplot=FALSE, Rplot=FALSE, Kplot=FALSE, Bplots=FALSE) {
   dd=analyse1(d)
   xxx1=cumdistplot(dd, plot=cdplot)   #cutoff, sens, spec, maxdif
   xxx2=rocplot(dd, plot=Rplot)       #cutoff, sens, spec
   xxx3=kolmogorovplot(dd,plot=Kplot)   # cutoff, lower & upper limit, same as with 
   xxx4=NA
   if (bootstrap==FALSE) {
      xxx4=rbind(
         paste("optimal cutoff-point according to the Youden index = ",round(xxx1[[1]],0),"(95% CI:",round(xxx3[[2]],0)," -", round(xxx3[[3]],0),")"),
         paste("largest diference between the cumulative distributions =", round(xxx1[[2]],2)),
         paste("with sensitivity =",round(xxx1[[4]],3), " and specificity =",round(xxx1[[3]],3)),
         paste("",""),
         paste("optimal cutoff-point according to the Q-point =",round(xxx2[[1]],0)),
         paste("with diference between the cumulative distributions =", round(xxx2[[2]],2)),
         paste("with sensitivity =",xxx2[[4]], " and specificity = ",xxx2[[3]]),
         paste("",""),
         paste("AUC =", round(dd$auc[1], 4), " SE=", round(dd$auc[2],4))
      )
      print(xxx4,quote=F)
   }
   ff4=NA
   if (bootstrap==TRUE | Bplots==TRUE) {
      ff4 = eenbootstrapfunctie(dd, nbootstraps) 
      ff5=rbind(
          c(mean(ff4[,1]), quantile(ff4[,1], probs=c(0.025,0.975))),
          c(mean(ff4[,2]), quantile(ff4[,2], probs=c(0.025,0.975))),
          c(mean(ff4[,3]), quantile(ff4[,3], probs=c(0.025,0.975))),
          c(mean(ff4[,4]), quantile(ff4[,4], probs=c(0.025,0.975)))
       )
       ff6=rbind(
          c(mean(ff4[,5]), quantile(ff4[,5], probs=c(0.025,0.975))),
          c(mean(ff4[,6]), quantile(ff4[,6], probs=c(0.025,0.975))),
          c(mean(ff4[,7]), quantile(ff4[,7], probs=c(0.025,0.975))),
          c(mean(ff4[,8]), quantile(ff4[,8], probs=c(0.025,0.975)))
       )
      ff7=c(mean(ff4[,9]), quantile(ff4[,9], probs=c(0.025,0.975)))
      colnames(ff5)=colnames(ff6)=names(ff7)=c("mean","95% CI lower limit","95% CI upper limit")
      row.names(ff5)=row.names(ff6)=c("cutoff-point","difference between cum. distributions","sensitivity","specificity")
      if (Bplots==TRUE) {
         dev.new()
         layout(matrix(c(1,3,2,4),ncol=2,nrow=2))
         hist(ff4[,1], xlab="best cutoff point according to the Youden Index", main="")
         hist(ff4[,2], xlab="max. difference between the cum. distributions", main="")
         hist(ff4[,3], xlab="sensitivity", main="")
         hist(ff4[,4], xlab="specificity", main="")
         dev.new()
         layout(matrix(c(1,3,2,4),ncol=2,nrow=2))
         hist(ff4[,5], xlab="best cutoff point according to the Q-point", main="")
         hist(ff4[,6], xlab="difference between the cum. distributions", main="")
         hist(ff4[,7], xlab="sensitivity", main="")
         hist(ff4[,8], xlab="specificity", main="")
      }
      legerij=rep("",3)
      xxx5=rbind(ff5,legerij,ff6,legerij,ff7)
      row.names(xxx5)=c("Optimal Cutoff according to the Youden-Index","Max. Difference of the cum. distributions","Sensitivity","Specificity","",
                        "Cutoff according to the Q-point","Difference of the cum. distributions","Sensitivity","Specificity","","AUC")
      print(xxx5,quote=FALSE)
      return(list(Youdenresults=ff5, Qpointresults=ff6, AUC=ff7, bootstrapresults=ff4))
   }
}
modelstring12 = "
   model {
      for (j in 1:T) {
         for (i in 1:N0) {
            event0[i,j] ~ dpois(incidence0[i,j])
            incidence0[i,j] <- Y0[i,j] * deltacumHaz0[j]
         }
         for (i in 1:N1) {
            event1[i,j] ~ dpois(incidence1[i,j])
            incidence1[i,j] <- Y1[i,j] * deltacumHaz1[j]
         }
         deltacumHaz0[j] ~ dgamma(mu0[j], c)
         deltacumHaz1[j] ~ dgamma(mu1[j], c)
         mu0[j] <- dL0.star[j] * c
         mu1[j] <- dL0.star[j] * c
         S0[j] <- exp(-sum(deltacumHaz0[1:j]))
         S1[j] <- exp(-sum(deltacumHaz1[1:j]))
      }
      c <- 0.001
      r <- 0.1
      for (j in 1:T) {
         dL0.star[j] <- r * (t[j+1] - t[j])
         KolmDist[j] <- abs(S0[j]-S1[j])
      }
      for (j in 1:(T-1)) {
         surf[j] <- ((1-S0[(j+1)]) - (1-S0[j])) * (S1[j]) + ((1-S0[(j+1)]) - (1-S0[j])) * (S1[(j+1)]-S1[j])/2
      }
      maxKD <- max(KolmDist)
      # voor welke j is KolmDist = maxKD  ???  doe dit maar buiten JAGS
      auc  <- sum(surf)  + (1-S0[(1)]) * S1[1]/2 + S0[T] * (1-S1[T])/2
   }
"
bayesiananalysis=function(d, thinning=10, updates=10000, iterations=10000, nchains=4) {
   library(rjags)
   T = length(unique(c(d$y0,d$y1)))      # aantal unieke failure times
   ksi=0.1
   t=c(sort(unique(c(d$y0,d$y1))), max(c(d$y0,d$y1))+ksi)
   eps = 0.000001 
   event0=Y0=matrix(NA, nrow=length(d$y0), ncol=T)
   for (i in 1:length(d$y0)) {
      for (j in 1:T) {
         Y0[i,j] <- 1*((d$y0[i] - t[j] + eps)>0)
         event0[i,j] <- Y0[i,j] * 1*((t[j+1] - d$y0[i] -eps)>0)
      }
   }
   event1=Y1=matrix(NA, nrow=length(d$y1), ncol=T)
   for (i in 1:length(d$y1)) {
      for (j in 1:T) {
         Y1[i,j] <- 1*((d$y1[i] - t[j] + eps)>0)
         event1[i,j] <- Y1[i,j] * 1*((t[j+1] - d$y1[i] -eps)>0)
      }
   }
   data_jags12 = list(N0=length(d$y0), N1=length(d$y1), T=T, t=t, Y0=Y0,event0=event0, Y1=Y1,event1=event1)
   model12=jags.model(textConnection(modelstring12), data=data_jags12, n.chains = nchains, n.adapt=0)
   update(model12, n.iter=updates)
   output12=coda.samples(model12,variable.names=c("S0", "S1", "KolmDist", "maxKD", "auc"), n.iter= iterations, thin=thinning)
   sumoutput12=summary(output12)
   return(list(t=t, samenvatting=sumoutput12, mcmcresult=output12))
}
bayesianresults=function(bayesres, cumdist=FALSE, kolmogorov=FALSE) {
   sumoutput12=bayesres$samenvatting
   output12=bayesres$mcmcresult
   t=bayesres$t
   qq1=sapply(row.names(sumoutput12[[1]]), function(x){substr(x,1,2)})
   ksi=0.1
   if (cumdist==TRUE) {
      St0=sumoutput12[[1]][which(qq1=="S0"),]
      St1=sumoutput12[[1]][which(qq1=="S1"),]
      Stlimits0=sumoutput12[[2]][which(qq1=="S0"),]
      Stlimits1=sumoutput12[[2]][which(qq1=="S1"),]
      dev.new()
      plot(1,1, type="n", xlab="biomarker y", ylab="cumulative distribution", xlim=c(min(t), max(t)), ylim=c(0,1))
      lines(c(min(t)-ksi,t), c(0,empfunxie(t,d$y0)[,1]), type="s", col=1, lty=1, lwd=2)
      lines(c(min(t)-ksi,t), c(0,1-St0[,1],1), col=3, lwd=2, lty=1)
      lines(c(min(t)-ksi,t), c(0,1-Stlimits0[,1],1), col=3, lwd=1, lty=2)
      lines(c(min(t)-ksi,t), c(0,1-Stlimits0[,5],1), col=3, lwd=1, lty=2)
      lines(c(min(t)-ksi,t), c(0,empfunxie(t,d$y1)[,1]), type="s", col=1, lty=1, lwd=2)
      lines(c(min(t)-ksi,t), c(0,1-St1[,1],1), col=2, lwd=2, lty=1)
      lines(c(min(t)-ksi,t), c(0,1-Stlimits1[,1],1), col=2, lwd=1, lty=2)
      lines(c(min(t)-ksi,t), c(0,1-Stlimits1[,5],1), col=2, lwd=1, lty=2)
      abline(v=t[which(abs(St0-St1) == max(abs(St0-St1)))], lty=2, col=1, lwd=1)
      text(x=3+t[which(abs(St0-St1)==max(abs(St0-St1)))], y=0.05, paste("largest difference =", round(max(abs(St0-St1)),4)), adj=0)
      text(x=3+t[which(abs(St0-St1)==max(abs(St0-St1)))], y=0,    paste("optimal cutoff y_c = ",round(t[which(abs(St0-St1)==max(abs(St0-St1)))],0)), adj=0)
   }
   maxKD=optimumy=sens1=spec1=sens2=spec2=mdif2=optimum2=matrix(NA, nrow=length(output12), ncol=nrow(output12[[1]]))
   for (j in 1:length(output12))   {           # 4 chains
      qq2=output12[[j]]
      for (i in 1:nrow(qq2))  {     # 100000 samples
         dddd=NA
         maxKD[j,i]=qq2[i,which(qq1=="ma")]   # maxKD
         optimumy[j,i]=median(t[which(qq2[i,which(qq1=="Ko")]==maxKD[j,i])])   # KolmDist
         sens1[j,i]=qq2[i,which(qq1=="S1")][which(qq2[i,which(qq1=="Ko")]==maxKD[j,i])]
         spec1[j,i]=1-qq2[i,which(qq1=="S0")][which(qq2[i,which(qq1=="Ko")]==maxKD[j,i])]
         dddd=(qq2[i,which(qq1=="S1")]-1)^2 + (qq2[i,which(qq1=="S0")]-0)^2
         optimum2[j,i]=median(t[which(dddd == min(dddd))])
         mdif2[j,i]= abs((1-qq2[i,which(qq1=="S0")][which(dddd == min(dddd))]) - (1-qq2[i,which(qq1=="S1")][which(dddd == min(dddd))]))
         sens2[j,i]=qq2[i,which(qq1=="S1")][which(dddd == min(dddd))]
         spec2[j,i]=1-qq2[i,which(qq1=="S0")][which(dddd == min(dddd))]
      }
   }
   if (kolmogorov==TRUE) {
      ff1=sumoutput12[[1]][which(qq1=="Ko"),] 
      ff2=sumoutput12[[2]][which(qq1=="Ko"),]
      dev.new()
      plot(1,1, type="n", col=1, lwd=2, lty=2, xlab="biomarker y", ylab="Kolmogorov Distance", 
         xlim=c(min(t),max(t)), ylim=c(0,max(ff2[,5])))
      lines(t[1:(length(t)-1)], ff1[,1], type="s", col=1, lwd=2, lty=1)
      lines(t[1:(length(t)-1)], ff2[,1], type="s", col=1, lwd=1, lty=2)
      lines(t[1:(length(t)-1)], ff2[,5], type="s", col=1, lwd=1, lty=2)
      abline(v=mean(as.numeric(optimumy)), lty=3, col=1, lwd=1)
      abline(v=quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[1], lty=3, col=1, lwd=1)
      abline(v=quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[3], lty=3, col=1, lwd=1)
      text(x=round(mean(as.numeric(optimumy)),1),y=0.05,round(mean(as.numeric(optimumy)),0))
      text(x=round(quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[1],1),y=0.01,round(quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[1],0))
      text(x=round(quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[3],1),y=0.01,round(quantile(as.numeric(optimumy), probs=c(0.025,0.5,0.975))[3],0))
   }
   sig2=sig3=sig4=sig5=sig6=sig7=sig8=c()
   i=1
   for (i in 1:length(output12)) {
      qq2=qq3=NA
      qq2=optimumy[i,]
      qq3=ar(qq2,aic=TRUE)
      sig2[i]=qq3$var.pred
      if (qq3$order>0) {sig2[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=sens1[i,]
      qq3=ar(qq2,aic=TRUE)
      sig3[i]=qq3$var.pred
      if (qq3$order>0) {sig3[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=spec1[i,]
      qq3=ar(qq2,aic=TRUE)
      sig4[i]=qq3$var.pred
      if (qq3$order>0) {sig4[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=optimum2[i,]
      qq3=ar(qq2,aic=TRUE)
      sig5[i]=qq3$var.pred
      if (qq3$order>0) {sig5[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=mdif2[i,]
      qq3=ar(qq2,aic=TRUE)
      sig6[i]=qq3$var.pred
      if (qq3$order>0) {sig6[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=sens2[i,]
      qq3=ar(qq2,aic=TRUE)
      sig7[i]=qq3$var.pred
      if (qq3$order>0) {sig7[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
      qq2=qq3=NA
      qq2=spec2[i,]
      qq3=ar(qq2,aic=TRUE)
      sig8[i]=qq3$var.pred
      if (qq3$order>0) {sig8[i]=qq3$var.pred/(1-sum(qq3$ar)^2)}
   }
   se.tv.opti1 = sqrt(mean(sig2)/(length(output12)* nrow(output12[[1]])))
   se.tv.sens1 = sqrt(mean(sig3)/(length(output12)* nrow(output12[[1]])))
   se.tv.spec1 = sqrt(mean(sig4)/(length(output12)* nrow(output12[[1]])))
   se.tv.opti2 = sqrt(mean(sig5)/(length(output12)* nrow(output12[[1]]))) 
   se.tv.mdif2 = sqrt(mean(sig6)/(length(output12)* nrow(output12[[1]])))
   se.tv.sens2 = sqrt(mean(sig7)/(length(output12)* nrow(output12[[1]])))
   se.tv.spec2 = sqrt(mean(sig8)/(length(output12)* nrow(output12[[1]])))
   MaxKolDistance = c(ff1=sumoutput12[[1]][which(qq1=="ma"),],ff2=sumoutput12[[2]][which(qq1=="ma"),])
   AUC = c(ff1=sumoutput12[[1]][which(qq1=="au"),],ff2=sumoutput12[[2]][which(qq1=="au"),])
   Sens1=c(mean(as.numeric(sens1)),sd(as.numeric(sens1)), 
         sd(as.numeric(sens1))/sqrt(length(output12)*nrow(output12[[1]])),
         se.tv.sens1, quantile(as.numeric(sens1), probs=c(0.025,0.25,0.5,0.75,0.975)))
   Spec1=c(mean(as.numeric(spec1)),sd(as.numeric(spec1)), 
         sd(as.numeric(spec1))/sqrt(length(output12)*nrow(output12[[1]])),
         se.tv.spec1, quantile(as.numeric(spec1), probs=c(0.025,0.25,0.5,0.75,0.975)))
   OptimalCutoff=c(mean(as.numeric(optimumy)),sd(as.numeric(optimumy)), 
      sd(as.numeric(optimumy))/sqrt(length(output12)*nrow(output12[[1]])),
      se.tv.opti1, quantile(as.numeric(optimumy), probs=c(0.025,0.25,0.5,0.75,0.975)))
   Qpointcutoff=c(mean(as.numeric(optimum2)), sd(as.numeric(optimum2)), sd(as.numeric(optimum2))/sqrt(length(output12)*nrow(output12[[1]])),
      se.tv.opti2, quantile(as.numeric(optimum2), probs=c(0.025,0.25,0.5,0.75,0.975)))
   Qpointdist=c(mean(as.numeric(mdif2)), sd(as.numeric(mdif2)), sd(as.numeric(mdif2))/sqrt(length(output12)*nrow(output12[[1]])),
      se.tv.mdif2, quantile(as.numeric(mdif2), probs=c(0.025,0.25,0.5,0.75,0.975)))
   Sens2=c(mean(as.numeric(sens2)), sd(as.numeric(sens2)), sd(as.numeric(sens2))/sqrt(length(output12)*nrow(output12[[1]])),
      se.tv.sens2, quantile(as.numeric(sens2), probs=c(0.025,0.25,0.5,0.75,0.975)))
   Spec2=c(mean(as.numeric(spec2)), sd(as.numeric(spec2)), sd(as.numeric(spec2))/sqrt(length(output12)*nrow(output12[[1]])),
      se.tv.spec2, quantile(as.numeric(spec2), probs=c(0.025,0.25,0.5,0.75,0.975)))
   names(AUC)=names(MaxKolDistance)=names(OptimalCutoff)=names(Sens1)=names(Spec1)=c("Mean","SD", "Naive.SE", "Time-series.SE","Pct2.5%","Pct25%","Pct50%","Pct75%","Pct97.5%")
   names(Qpointcutoff)=names(Qpointdist)=names(Sens2)=names(Spec2)=c("Mean","SD", "Naive.SE", "Time-series.SE","Pct2.5%","Pct25%","Pct50%","Pct75%","Pct97.5%")
   legerij=rep("",9)
   xxxx5=rbind(OptimalCutoff, MaxKolDistance, Sens1, Spec1, legerij, Qpointcutoff, Qpointdist, Sens2, Spec2, legerij, AUC)
   row.names(xxxx5)=c("Optimal Cutoff according to the Youden Index","Max. Difference of the cum. distributions","Sensitivity", "Specificity","",
                      "Cutoff according to the Q-point", "Difference of the cum. distributions","Sensitivity", "Specificity", "", "AUC")
   print(xxxx5, quote=FALSE)
   return(list(OptimalCutoff=OptimalCutoff, MaxKolDistance=MaxKolDistance, Sens1=Sens1, Spec1=Spec1, 
               QpointCutoff=Qpointcutoff,   QpointDistance=Qpointdist,     Sens2=Sens2, Spec2=Spec2, AUC=AUC,
               optimumy=optimumy, maxKD=maxKD, sens1=sens1, spec1=spec1, optimum2=optimum2, sens2=sens2, spec2=spec2, mdif2=mdif2))
}
tweebootstrap=function(d, nbootstraps, model="gaussian") {
   YI=QP=maxdif=QPdif=Sens1=Spec1=Sens2=Spec2=AUC=c()
   for (iter in 1:nbootstraps) {
      nrs0=nrs1=yy0=yy1=dif=ddd=NA
      nrs0=sample(1:length(d$y0), length(d$y0), replace=TRUE)
      nrs1=sample(1:length(d$y1), length(d$y1), replace=TRUE)
      yy0=d$y0[nrs0]
      yy1=d$y1[nrs1]
      vv0=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)}, lower=.Machine$double.eps^0.25, upper=100, x=yy0)
      shape_k0=vv0$root
      scale_lambda0=(sum(yy0^shape_k0)/length(yy0))^(1/shape_k0)
      vv1=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)}, lower=.Machine$double.eps^0.25, upper=100, x=yy1)
      shape_k1=vv1$root
      scale_lambda1=(sum(yy1^shape_k1)/length(yy1))^(1/shape_k1)
      vv20=optim(c(1,median(yy0)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=yy0, control=list(fnscale=-1))
      cauchyscale0=vv20$par[1]
      cauchyloc0=vv20$par[2]
      vv21=optim(c(1,median(yy1)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=yy1, control=list(fnscale=-1))
      cauchyscale1=vv21$par[1]
      cauchyloc1=vv21$par[2]
      uniekey = sort(unique(c(yy0,yy1)))
      if (model=="gaussian") {pp0=pnorm((uniekey-mean(yy0))/sd(yy0))
                              pp1=pnorm((uniekey-mean(yy1))/sd(yy1))}
      if (model=="logistic") {pp0=plogis(uniekey, location=mean(yy0), scale=sd(yy0)*sqrt(3)/pi)
                              pp1=plogis(uniekey, location=mean(yy1), scale=sd(yy1)*sqrt(3)/pi)}
      if (model=="lognormal") {pp0=plnorm(uniekey, mean=mean(log(yy0)), sd=sd(log(yy0)))
                               pp1=plnorm(uniekey, mean=mean(log(yy1)), sd=sd(log(yy1)))}
      if (model=="t") {pp0=pt((uniekey-mean(yy0))/sd(yy0), df=length(yy0)-1)
                       pp1=pt((uniekey-mean(yy1))/sd(yy1), df=length(yy1)-1)}
      if (model=="gamma") {pp0=pgamma(uniekey, shape=mean(yy0)^2/var(yy0), scale=var(yy0)/mean(yy0))
                           pp1=pgamma(uniekey, shape=mean(yy1)^2/var(yy1), scale=var(yy1)/mean(yy1))}
      if (model=="weibull") {pp0=pweibull(uniekey, shape=shape_k0, scale=scale_lambda0)
                             pp1=pweibull(uniekey, shape=shape_k1, scale=scale_lambda1)}
      if (model=="cauchy") {pp0=pcauchy(uniekey, location=cauchyloc0, scale=cauchyscale0)
                            pp1=pcauchy(uniekey, location=cauchyloc1, scale=cauchyscale1)}
      dif=abs(pp0-pp1)
      ddd=sqrt((1-pp0)^2 + (pp1)^2)
      YI[iter]=uniekey[which(dif==max(dif))]
      QP[iter]=uniekey[which(ddd==min(ddd))] 
      maxdif[iter]=max(dif)
      QPdif[iter]=dif[which(ddd==min(ddd))]
      Sens1[iter]=pp0[which(dif==max(dif))]
      Spec1[iter]=1-pp1[which(dif==max(dif))]
      Sens2[iter]=pp0[which(ddd==min(ddd))]
      Spec2[iter]=1-pp1[which(ddd==min(ddd))]
      fff1=1-pp0
      fff2=1-pp1
      rangorde=order(fff1,fff2)
      fff1=fff1[rangorde]
      fff2=fff2[rangorde]
      AUC[iter]=sum((fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)]) * fff2[1:(length(fff2)-1)] + (fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)])*(fff2[2:length(fff2)]-fff2[1:(length(fff1)-1)])/2)
   }
   xxx=cbind(YI,maxdif,Sens1,Spec1,QP,QPdif,Sens2,Spec2,AUC)
   colnames(xxx)=c("Youden index", "max. difference", "Sensitivity", "Specificity", "Q-point","difference","Sensitivity","Specicity","AUC")
   return(xxx)
}
analyse2 = function(d) {
   # nog toevoegen: cauchy-verdeling ?  Niet: chikw-verdeling (lijkt teveel op gamma-verdeling), log logistic (lijkt erg op lognormaal en is niet makkelijk(er))
   m0=mean(d$y0)
   s0=sd(d$y0)
   vv0=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)},
      lower=.Machine$double.eps^0.25, upper=100, x=d$y0)
   shape_k0=vv0$root
   scale_lambda0=(sum(d$y0^shape_k0)/length(d$y0))^(1/shape_k0)
   m1=mean(d$y1)
   s1=sd(d$y1)
   vv1=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)},
      lower=.Machine$double.eps^0.25, upper=100, x=d$y1)
   shape_k1=vv1$root
   scale_lambda1=(sum(d$y1^shape_k1)/length(d$y1))^(1/shape_k1)
   vv20=optim(c(1,median(d$y0)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=d$y0, control=list(fnscale=-1))
   cauchyscale0=vv20$par[1]
   cauchyloc0=vv20$par[2]
   vv21=optim(c(1,median(d$y1)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=d$y1, control=list(fnscale=-1))
   cauchyscale1=vv21$par[1]
   cauchyloc1=vv21$par[2]

   unieke_y=sort(unique(c(d$y0,d$y1)))

   exp0=pexp(unieke_y, rate=1/m0)
   wei0=pweibull(unieke_y, shape=shape_k0, scale=scale_lambda0)
   nor0=pnorm(unieke_y, mean=m0, sd=s0)
   log0=plogis(unieke_y, location=m0, scale=s0*sqrt(3)/pi)
   lno0=plnorm((unieke_y), mean=mean(log(d$y0)), sd=sd(log(d$y0)))
   tve0=pt((unieke_y-mean(d$y0))/sd(d$y0), df=length(d$y0)-1)
   gam0=pgamma(unieke_y, shape=mean(d$y0)^2/var(d$y0), scale=var(d$y0)/mean(d$y0))
   cau0=pcauchy(unieke_y, location=cauchyloc0, scale=cauchyscale0)

   exp1=pexp(unieke_y, rate=1/m1)
   wei1=pweibull(unieke_y, shape=shape_k1, scale=scale_lambda1)
   nor1=pnorm(unieke_y, mean=m1, sd=s1)
   log1=plogis(unieke_y, location=m1, scale=s1*sqrt(3)/pi)
   lno1=plnorm((unieke_y), mean=mean(log(d$y1)), sd=sd(log(d$y1)))
   tve1=pt((unieke_y-mean(d$y1))/sd(d$y1), df=length(d$y1)-1)
   gam1=pgamma(unieke_y, shape=mean(d$y1)^2/var(d$y1), scale=var(d$y1)/mean(d$y1))
   cau1=pcauchy(unieke_y, location=cauchyloc1, scale=cauchyscale1)

   dev.new()
   plot(1,1, type="n", xlim=c(min(unieke_y),max(unieke_y)), ylim=c(0,1), xlab="biomarker y", ylab="cumulative distribution")
   lines(unieke_y, empfunxie(unieke_y, d$y0)[,1], lty=1)
   lines(unieke_y, wei0, col=2, lty=1)
   lines(unieke_y, nor0, col=3, lty=1)
   lines(unieke_y, log0, col=4, lty=1)
   lines(unieke_y, lno0, col=5, lty=1)
   lines(unieke_y, tve0, col=6, lty=1)
   lines(unieke_y, gam0, col=7, lty=1)
   lines(unieke_y, cau0, col=8, lty=1)
   lines(unieke_y, empfunxie(unieke_y, d$y1)[,1], lty=1)
   lines(unieke_y, wei1, col=2, lty=2)
   lines(unieke_y, nor1, col=3, lty=2)
   lines(unieke_y, log1, col=4, lty=2)
   lines(unieke_y, lno1, col=5, lty=2)
   lines(unieke_y, tve1, col=6, lty=2)
   lines(unieke_y, gam1, col=7, lty=2)
   lines(unieke_y, cau1, col=8, lty=2)
   legend(x=min(unieke_y), y=1, legend=c("Weibull","Normal","Logistic","log Normal","t","gamma","Cauchy"), col=2:8, lty=1, bty="n", adj=0)

   dif2=abs(wei0-wei1)
   dif3=abs(nor0-nor1)
   dif4=abs(log0-log1)
   dif5=abs(lno0-lno1)
   dif6=abs(tve0-tve1)
   dif7=abs(gam0-gam1)
   dif8=abs(cau0-cau1)

   res2_0=empfunxie(unieke_y, d$y0)[,1]-wei0 
   res3_0=empfunxie(unieke_y, d$y0)[,1]-nor0 
   res4_0=empfunxie(unieke_y, d$y0)[,1]-log0 
   res5_0=empfunxie(unieke_y, d$y0)[,1]-lno0 
   res6_0=empfunxie(unieke_y, d$y0)[,1]-tve0 
   res7_0=empfunxie(unieke_y, d$y0)[,1]-gam0 
   res8_0=empfunxie(unieke_y, d$y0)[,1]-cau0 
   res2_1=empfunxie(unieke_y, d$y1)[,1]-wei1 
   res3_1=empfunxie(unieke_y, d$y1)[,1]-nor1 
   res4_1=empfunxie(unieke_y, d$y1)[,1]-log1 
   res5_1=empfunxie(unieke_y, d$y1)[,1]-lno1 
   res6_1=empfunxie(unieke_y, d$y1)[,1]-tve1
   res7_1=empfunxie(unieke_y, d$y1)[,1]-gam1 
   res8_1=empfunxie(unieke_y, d$y1)[,1]-cau1 
   reskw2=sum(res2_0^2+res2_1^2)^0.5
   reskw3=sum(res3_0^2+res3_1^2)^0.5
   reskw4=sum(res4_0^2+res4_1^2)^0.5
   reskw5=sum(res5_0^2+res5_1^2)^0.5
   reskw6=sum(res6_0^2+res6_1^2)^0.5
   reskw7=sum(res7_0^2+res7_1^2)^0.5
   reskw8=sum(res8_0^2+res8_1^2)^0.5

   llike1=sum(dexp(d$y0, rate=1/m0, log=T))+sum(dexp(d$y1, rate=1/m1, log=T))
   llike2=sum(dweibull(d$y0, shape=shape_k0, scale=scale_lambda0, log=TRUE))+sum(dweibull(d$y1, shape=shape_k1, scale=scale_lambda1, log=TRUE))
   llike3=sum(dnorm(d$y0, mean=m0,sd=s0,log=TRUE))+sum(dnorm(d$y1,mean=m1,sd=s1,log=TRUE))
   llike4=sum(dlogis(d$y0, location=m0, scale=s0*sqrt(3)/pi,log=TRUE))+sum(dlogis(d$y1, location=m1, scale=s1*sqrt(3)/pi,log=TRUE))
   llike5=sum(dlnorm((d$y0), mean=mean(log(d$y0)),sd=sd(log(d$y0)),log=TRUE))+sum(dlnorm((d$y1),mean=mean(log(d$y1)),sd=sd(log(d$y1)),log=TRUE))
   llike6=sum(log(dt((d$y0-mean(d$y0))/sd(d$y0), df=length(d$y0)-1,log=FALSE)/sd(d$y0)))+sum(log(dt((d$y1-mean(d$y1))/sd(d$y1), df=length(d$y1)-1,log=FALSE)/sd(d$y1)))
   llike7=sum(dgamma(d$y0, shape=mean(d$y0)^2/var(d$y0), scale=var(d$y0)/mean(d$y0),log=TRUE))+sum(dgamma(d$y1,shape=mean(d$y1)^2/var(d$y1), scale=var(d$y1)/mean(d$y1),log=TRUE))
   llike8=sum(dcauchy(d$y0, location=cauchyloc0, scale=cauchyscale0, log=TRUE))+sum(dcauchy(d$y1, location=cauchyloc1, scale=cauchyscale1, log=TRUE))

   xxx=cbind(
      c(unieke_y[which(abs(dif2)==max(abs(dif2)))],unieke_y[which(abs(dif3)==max(abs(dif3)))],unieke_y[which(abs(dif4)==max(abs(dif4)))],
        unieke_y[which(abs(dif5)==max(abs(dif5)))], unieke_y[which(abs(dif6)==max(abs(dif6)))], unieke_y[which(abs(dif7)==max(abs(dif7)))],
        unieke_y[which(abs(dif8)==max(abs(dif8)))]),
      c(max(abs(dif2)), max(abs(dif3)), max(abs(dif4)), max(abs(dif5)), max(abs(dif6)), max(abs(dif7)), max(abs(dif8))),
      c(reskw2,reskw3,reskw4,reskw5, reskw6, reskw7, reskw8),
      c(llike2,llike3,llike4,llike5, llike6, llike7, llike8)
   )
   colnames(xxx)=c("optimal cutoff-point","maximal distance","root sum of squared residuals","log likelihood")
   row.names(xxx)=c("weibull", "gaussian", "logistic", "lognormal", "t", "gamma", "Cauchy")
   return(xxx)
}
steekproevenverdeling3 = function(d, niter=10000, model="gaussian") {
   library(mvtnorm)
   par0=par1=rep(NA,2)
   Sigma0=Sigma1=matrix(0,ncol=2,nrow=2)
   if (model!="weibull" & model!="cauchy" & model!="lognormal") {
      par0=c(mean(d$y0), log(sd(d$y0)))
      Sigma0[1,1]=var(d$y0)/length(d$y0)
      Sigma0[2,2]=var(d$y0)/(2*length(d$y0)*var(d$y1))
      par1=c(mean(d$y1), log(sd(d$y1)))
      Sigma1[1,1]=var(d$y1)/length(d$y1)
      Sigma1[2,2]=var(d$y1)/(2*length(d$y1)*var(d$y1))
   }
   if (model=="weibull") {  # par[1]=log(a) en a=shape & a>0; par[2]=log(sigma) en sigma=scale & sigma>0
      vv10=optim(c(0.5,1), function(par,x){length(x)*(par[1]-par[2])+sum((exp(par[1])-1)*(log(x)-par[2]))-sum((x/exp(par[2]))^exp(par[1]))}, x=d$y0, control=list(fnscale=-1), hessian=TRUE)
      vv11=optim(c(0.5,1), function(par,x){length(x)*(par[1]-par[2])+sum((exp(par[1])-1)*(log(x)-par[2]))-sum((x/exp(par[2]))^exp(par[1]))}, x=d$y1, control=list(fnscale=-1), hessian=TRUE)
      par0=vv10$par
      Sigma0=solve(-vv10$hessian)
      par1=vv11$par
      Sigma1=solve(-vv11$hessian)
   }
   if (model=="lognormal") {   # par[1]=mean van de log(y)  en par[2]=log(sd van de log(y))
      par0=c(mean(log(d$y0)), log(sd(log(d$y0))))
      Sigma0[1,1]=var(log(d$y0))/length(d$y0)
      Sigma0[2,2]=var(log(d$y0))/(2*length(d$y0)*var(log(d$y0)))
      par1=c(mean(log(d$y1)), log(sd(log(d$y1))))
      Sigma1[1,1]=var(log(d$y1))/length(d$y1)
      Sigma1[2,2]=var(log(d$y1))/(2*length(d$y1)*var(log(d$y1)))
   }
   if (model=="cauchy") {     # par[1]=log(s) en s=scale & s>0; par[2]=l is location & -inf < l < inf
      vv20=optim(c(1,median(d$y0)), function(par,x) {-length(x)*(par[1]+log(pi))-sum(log(1+((x-par[2])/exp(par[1]))^2))}, x=d$y0, control=list(fnscale=-1), hessian=TRUE)
      vv21=optim(c(1,median(d$y1)), function(par,x) {-length(x)*(par[1]+log(pi))-sum(log(1+((x-par[2])/exp(par[1]))^2))}, x=d$y1, control=list(fnscale=-1), hessian=TRUE)
      par0=vv20$par
      Sigma0=solve(-vv20$hessian)
      par1=vv21$par
      Sigma1=solve(-vv21$hessian)
   }
   stats=npstats=matrix(NA, ncol=9, nrow=niter)
   bewaardata=list()
   for (iter in 1:niter) {
      parms0=rmvnorm(1, par0, Sigma0)
      parms1=rmvnorm(1, par1, Sigma1)
      if (model=="weibull") {
         yy0=rweibull(length(d$y0), shape=exp(parms0[1]), scale=exp(parms0[2]))
         yy1=rweibull(length(d$y1), shape=exp(parms1[1]), scale=exp(parms1[2]))
      }
      if (model=="gaussian") {
         yy0=rnorm(length(d$y0), mean=parms0[1], sd=exp(parms0[2]))
         yy1=rnorm(length(d$y1), mean=parms1[1], sd=exp(parms1[2]))
      }
      if (model=="logistic") {
         yy0=rlogis(length(d$y0), location=parms0[1], scale=exp(parms0[2])*sqrt(3)/pi)
         yy1=rlogis(length(d$y1), location=parms1[1], scale=exp(parms1[2])*sqrt(3)/pi)
      }
      if (model=="t") {
         yy0=rt(length(d$y0), df=length(d$y0)-1)*exp(parms0[2])+parms0[1]
         yy1=rt(length(d$y1), df=length(d$y1)-1)*exp(parms1[2])+parms1[1]
      }
      if (model=="lognormal") {
         yy0=rlnorm(length(d$y0), meanlog=parms0[1], sdlog=exp(parms0[2]))
         yy1=rlnorm(length(d$y1), meanlog=parms1[1], sdlog=exp(parms1[2]))
      }
      if (model=="gamma") {
         yy0=rgamma(length(d$y0), shape=parms0[1]^2/exp(2*parms0[2]), scale=exp(2*parms0[2])/parms0[1])
         yy1=rgamma(length(d$y1), shape=parms1[1]^2/exp(2*parms1[2]), scale=exp(2*parms1[2])/parms1[1])
      }
      if (model=="cauchy") {
         yy0=rcauchy(length(d$y0), location=parms0[2], scale=exp(parms0[1]))
         yy1=rcauchy(length(d$y1), location=parms1[2], scale=exp(parms1[1]))
      }
      uniekey=c(min(yy0,yy1)-1,sort(unique(yy0,yy1)),max(yy0,yy1)+1)
      #GG=empfunxie(uniekey, yy1)
      #HH=empfunxie(uniekey, yy0)
      #diff=abs(HH[,1]-GG[,1])
      #dddd=sqrt(((1-HH[,1])-0)^2 + ((1-GG[,1])-1)^2)
      #fff1=1-HH[,1]
      #fff2=1-GG[,1]
      #rangorde=order(fff1,fff2)
      #fff1=fff1[rangorde]
      #fff2=fff2[rangorde]
      #AUC=sum((fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)]) * fff2[1:(length(fff2)-1)] + (fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)])*(fff2[2:length(fff2)]-fff2[1:(length(fff1)-1)])/2)
      #npstats[iter,]=c(uniekey[diff==max(diff)],max(diff),HH[(diff==max(diff)),1],1-GG[(diff==max(diff)),1],
      #                 uniekey[dddd==min(dddd)],diff[dddd==min(dddd)],HH[(dddd==min(dddd)),1],1-GG[(dddd==min(dddd)),1], AUC)
      bewaardata[[iter]]=list(y0=yy0, y1=yy1)
      if (model=="weibull") {pp0=pweibull(uniekey, shape=shape_k0, scale=scale_lambda0)
                            pp1=pweibull(uniekey, shape=shape_k1, scale=scale_lambda1)}
      if (model=="gaussian") {pp0=pnorm((uniekey-mean(yy0))/sd(yy0))
                              pp1=pnorm((uniekey-mean(yy1))/sd(yy1))}
      if (model=="logistic") {pp0=plogis(uniekey, location=mean(yy0), scale=sd(yy0)*sqrt(3)/pi)
                              pp1=plogis(uniekey, location=mean(yy1), scale=sd(yy1)*sqrt(3)/pi)}
      if (model=="t") {pp0=pt((uniekey-mean(yy0))/sd(yy0), df=length(yy0)-1)
                       pp1=pt((uniekey-mean(yy1))/sd(yy1), df=length(yy1)-1)}
      if (model=="lognormal") {pp0=plnorm(uniekey, mean=mean(log(yy0)), sd=sd(log(yy0)))
                               pp1=plnorm(uniekey, mean=mean(log(yy1)), sd=sd(log(yy1)))}
      if (model=="gamma") {pp0=pgamma(uniekey, shape=mean(yy0)^2/var(yy0), scale=var(yy0)/mean(yy0))
                           pp1=pgamma(uniekey, shape=mean(yy1)^2/var(yy1), scale=var(yy1)/mean(yy1))}
      if (model=="cauchy") {pp0=pcauchy(uniekey, location=cauchyloc0, scale=cauchyscale0)
                            pp1=pcauchy(uniekey, location=cauchyloc1, scale=cauchyscale1)}
      dif=abs(pp0-pp1)
      ddd=sqrt(((1-pp0)-0)^2 + ((1-pp1)-1)^2)
      fff1=1-pp0
      fff2=1-pp1
      rangorde=order(fff1,fff2)
      fff1=fff1[rangorde]
      fff2=fff2[rangorde]
      AUC=sum((fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)]) * fff2[1:(length(fff2)-1)] + (fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)])*(fff2[2:length(fff2)]-fff2[1:(length(fff1)-1)])/2)
      stats[iter,]=c(uniekey[dif==max(dif)], max(dif), 1-pp1[dif==max(dif)], pp0[dif==max(dif)],
                     uniekey[ddd==min(ddd)], dif[ddd==min(ddd)], 1-pp1[ddd==min(ddd)], pp0[ddd==min(ddd)],AUC)
   }
   return(stats)#bewaardata)
}
analyse3 = function(d, bootstrap=T, popsampling=F, nbootstraps=10000, model="gaussian") {
   YI=QP=maxdif=QPdif=Sens1=Spec1=Sens2=Spec2=AUC=c()
   yy0=d$y0
   yy1=d$y1
   if (model=="weibull") {
      vv0=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)}, lower=.Machine$double.eps^0.25, upper=100, x=yy0)
      shape_k0=vv0$root
      scale_lambda0=(sum(yy0^shape_k0)/length(yy0))^(1/shape_k0)
      vv1=uniroot(function(k,x) {((sum(x^k * log(x)))/sum(x^k)) - 1/k - sum(log(x))/length(x)}, lower=.Machine$double.eps^0.25, upper=100, x=yy1)
      shape_k1=vv1$root
      scale_lambda1=(sum(yy1^shape_k1)/length(yy1))^(1/shape_k1)
   }
   if (model=="cauchy") {
      vv20=optim(c(1,median(yy0)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=yy0, control=list(fnscale=-1))
      cauchyscale0=vv20$par[1]
      cauchyloc0=vv20$par[2]
      vv21=optim(c(1,median(yy1)), function(par,x) {-length(x)*log(par[1]*pi)-sum(log(1+((x-par[2])/par[1])^2))} ,  x=yy1, control=list(fnscale=-1))
      cauchyscale1=vv21$par[1]
      cauchyloc1=vv21$par[2]
   }
   uniekey = sort(unique(c(yy0,yy1)))
   if (model=="gaussian") {pp0=pnorm((uniekey-mean(yy0))/sd(yy0))
                           pp1=pnorm((uniekey-mean(yy1))/sd(yy1))}
   if (model=="logistic") {pp0=plogis(uniekey, location=mean(yy0), scale=sd(yy0)*sqrt(3)/pi)
                           pp1=plogis(uniekey, location=mean(yy1), scale=sd(yy1)*sqrt(3)/pi)}
   if (model=="lognormal") {pp0=plnorm(uniekey, mean=mean(log(yy0)), sd=sd(log(yy0)))
                            pp1=plnorm(uniekey, mean=mean(log(yy1)), sd=sd(log(yy1)))}
   if (model=="t") {pp0=pt((uniekey-mean(yy0))/sd(yy0), df=length(yy0)-1)
                    pp1=pt((uniekey-mean(yy1))/sd(yy1), df=length(yy1)-1)}
   if (model=="gamma") {pp0=pgamma(uniekey, shape=mean(yy0)^2/var(yy0), scale=var(yy0)/mean(yy0))
                        pp1=pgamma(uniekey, shape=mean(yy1)^2/var(yy1), scale=var(yy1)/mean(yy1))}
   if (model=="weibull") {pp0=pweibull(uniekey, shape=shape_k0, scale=scale_lambda0)
                          pp1=pweibull(uniekey, shape=shape_k1, scale=scale_lambda1)}
   if (model=="cauchy") {pp0=pcauchy(uniekey, location=cauchyloc0, scale=cauchyscale0)
                         pp1=pcauchy(uniekey, location=cauchyloc1, scale=cauchyscale1)}
   dif=abs(pp0-pp1)
   ddd=sqrt((1-pp0)^2 + (pp1)^2)
   YI=uniekey[which(dif==max(dif))]
   QP=uniekey[which(ddd==min(ddd))] 
   maxdif=max(dif)
   QPdif=dif[which(ddd==min(ddd))]
   Sens1=pp0[which(dif==max(dif))]
   Spec1=1-pp1[which(dif==max(dif))]
   Sens2=pp0[which(ddd==min(ddd))]
   Spec2=1-pp1[which(ddd==min(ddd))]
   fff1=1-pp0
   fff2=1-pp1
   rangorde=order(fff1,fff2)
   fff1=fff1[rangorde]
   fff2=fff2[rangorde]
   AUC=sum((fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)]) * fff2[1:(length(fff2)-1)] + (fff1[2:length(fff1)]-fff1[1:(length(fff1)-1)])*(fff2[2:length(fff2)]-fff2[1:(length(fff1)-1)])/2)
   if (bootstrap==TRUE & popsampling==FALSE) {xxx=tweebootstrap(d, nbootstraps, model="gaussian")}
   if (popsampling==TRUE) {stats=steekproevenverdeling3(d, nbootstraps, model="gaussian")}
   dev.new()
   plot(uniekey,  empfunxie(uniekey,yy0)[,1], type="l", xlim=c(min(uniekey),max(uniekey)), ylim=c(0,1), 
      sub=paste(model,"distributions"), xlab="biomarker y", ylab="cumulative distribution function", lwd=2, lty=1, col=1, cex.sub=0.7)
   lines(uniekey, empfunxie(uniekey,yy1)[,1], lty=1, col=1, lwd=2)
   lines(uniekey, pp0, lty=1, col=3, lwd=2)
   pp0x=pp0
   pp0x[pp0>0.99]=0.99
   pp0x[pp0<0.01]=0.01
   pp1x=pp1
   pp1x[pp1>0.99]=0.99
   pp1x[pp1<0.01]=0.01
   lines(uniekey,expit(logit(pp0)-1.96*sqrt(1/(pp0x*(1-pp0x)*length(yy0)))), lty=2, lwd=1, col=3)
   lines(uniekey,expit(logit(pp0)+1.96*sqrt(1/(pp0x*(1-pp0x)*length(yy0)))), lty=2, lwd=1, col=3)
   lines(uniekey, pp1, type="l", lwd=2, lty=1, col=2)
   lines(uniekey,expit(logit(pp1)-1.96*sqrt(1/(pp1x*(1-pp1x)*length(yy1)))), lty=2, lwd=1, col=2)
   lines(uniekey,expit(logit(pp1)+1.96*sqrt(1/(pp1x*(1-pp1x)*length(yy1)))), lty=2, lwd=1, col=2)
   text(x=3+median(uniekey[which(dif==max(dif))]), y=0.05,paste("largest difference = ",round(max(dif),4)), adj=0)
   text(x=3+median(uniekey[which(dif==max(dif))]), y=0,   paste("optimal cutoff y_c = ",round(median(uniekey[which(dif==max(dif))]),0)), adj=0)
   abline(v=YI, lty=2, col=1, lwd=1)
   xxxx=as.matrix(c(YI,maxdif,Sens1,Spec1, QP, QPdif, Sens2,Spec2,AUC),ncol=1,nrow=9)
   row.names(xxxx)=c("Youden Index","max. difference between cum. distributions","Sensitivity","Specificity",
                     "Q-point","Difference between cum. distributions","Sensitivity","Specificity","AUC")
   if (bootstrap==TRUE) {xxxx=cbind(c(YI,maxdif,Sens1,Spec1, QP, QPdif, Sens2,Spec2,AUC),
                                  t(apply(xxx,2,function(x){quantile(x,probs=c(0.50, 0.025,0.975))})))
      colnames(xxxx)=c("estimate", "median", "95% CI lower limit", "95% CI higher limit")
      row.names(xxxx)=c("Youden Index","max. difference between cum. distributions","Sensitivity","Specificity",
                        "Q-point","Difference between cum. distributions","Sensitivity","Specificity","AUC")
   }
   if (popsampling ==TRUE) {
      yyyy=cbind(xxxx,cbind( apply(stats,2,median), t(apply(stats,2,function(x){quantile(x,probs=c(0.025,0.975))}))))
      colnames(yyyy)=c("estimate", "median", "95% CI lower limit", "95% CI higher limit")
      row.names(yyyy)=c("Youden Index","max. difference between cum. distributions","Sensitivity","Specificity",
                        "Q-point","Difference between cum. distributions","Sensitivity","Specificity","AUC")
      xxxx=yyyy
   }
   return(xxxx)
}


#########################################################################################################################################################################################






