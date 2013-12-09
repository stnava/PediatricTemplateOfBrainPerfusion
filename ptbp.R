library(ANTsR)
library(visreg)
data('aal',package='ANTsR')
library(pheatmap)
library(e1071)
thickinds<-c(1:36,39:40,43:70,79:90)
demog<-read.csv("demo_thick_cbf.csv")
demog<-demog[,3:228]
jhu<-read.csv("JHU_whitematter.csv")
demog<-cbind(demog,jhu[,4:ncol(jhu)])
popul<-demog[,1:22]
thickness<-demog[,grep("Thick",names(demog))]
thickness <- thickness[,thickinds]
colnames( thickness )<-c( paste("Thick",aal$label_name[thickinds]) )
CBF<-demog[,grep("CBF",names(demog))]
CBF <- CBF[,thickinds]
gcbf<-apply(as.matrix(CBF),FUN=mean,MARGIN=1)
colnames( CBF )<-c( paste("CBF",aal$label_name[thickinds]) )
FA<-demog[,grep("FA",names(demog))]
MD<-demog[,grep("MD_",names(demog))]
boldin<-read.csv("AAL_BoldConn.csv")
boldin<-boldin[,c(4,95:ncol(boldin))]
gbold<-apply(boldin,FUN=mean,MARGIN=1)
globmess<-read.csv("global_measures.csv")
globmess<-globmess[,2:ncol(globmess)]
usesubs<-rep(TRUE,nrow(demog))
usesubs[ c(155,142) ] <- FALSE      # for thickness
usesubs[ c(155,142,which( is.na(gbold) ),which( is.na(gcbf) ),which( gcbf < 10  ), which( is.na(jhu$RD_JHU.43) ) )] <- FALSE  # for cbf
tempdf<-data.frame( cbind( popul, usesubs ) )
write.csv(tempdf,"good_v_bad_cbf_subjects.csv",row.names=F)
uids <- unique( demog$SubID )
for ( i in 1:length(uids) )
  {
  ww<-which( demog$SubID == uids[i] )
  if ( length( ww ) > 1 ) usesubs[ ww[2:length(ww)] ] <- FALSE
  }
popul<-popul[usesubs,]
CBF<-CBF[usesubs,]
FA<-FA[usesubs,]
MD<-MD[usesubs,]
bold<-boldin[usesubs,]
thickness<-thickness[usesubs,]
globmess<-globmess[usesubs,]
inclevs<-levels( popul$Income )
inclevs[1]<-NA
inclevs[2:6]<-110.0
inclevs[7:9]<-130.0
inclevs[10:11]<-150.0
inclevs[12]<-170.0
inclevs[13]<-190.0
inclevs[c(14:17,19,20,32)]<-210.0
inclevs[18]<-35.0
inclevs[21:22]<-50.0
inclevs[23]<-60.0
inclevs[24:25]<-70.0
inclevs[26]<-85.0
inclevs[27:28]<-90.0
inclevs[c(29,33)]<-NA
inclevs[34]<-15.0
inclevs[30:31]<-20.0
inclevs<-as.numeric( inclevs )
levels( popul$Income )<-( inclevs )
popul$Income<-as.numeric( as.character( popul$Income ) )
myincome<-c( impute( cbind(popul$blank,popul$Income ) ) )
myincome2<-myincome
myincome2[ myincome >= 120 ]<-6
myincome2[ myincome <= 120 ]<-5
myincome2[ myincome <= 90 ]<- 4
myincome2[ myincome <= 70 ]<- 3
myincome2[ myincome <= 50 ]<- 2
myincome2[ myincome <= 20 ]<- 1
popul$AgeAtScan <- as.numeric(as.character(popul$AgeAtScan))
myiq<-( c( impute( cbind( popul$blank, as.numeric(as.character(popul$FullScaleIQ)) ) ) ) )
myiq2<-( c( impute( cbind( popul$blank, as.numeric(as.character(popul$Verbal.IQ)) ) ) ) )
myiq3<-( c( impute( cbind( popul$blank, as.numeric(as.character(popul$Performance.IQ)) ) ) ) )
mylad<-c( impute( cbind( popul$blank, as.numeric(as.character(popul$Teen.Ladder.SES.score ) )   ) ) )
myladc<-c( impute( cbind( popul$blank, as.numeric(as.character(popul$Teen.Ladder.Community.Score ) )   ) ) )
####################################################################################################
####################################################################################################
# for ( i in 1:ncol(CBF) ) {
#   CBF[,i] <- CBF[,i] / thickness[,i] 
#   }
FIQ      <- myiq
VIQ      <- myiq2
PIQ      <- myiq3
Ladder  <- ( mylad )
RIncome <- rank(myincome)
BV<-globmess$BrainVolume
myglobal<-globmess[,2:ncol(globmess)]
brainpreds<-as.matrix(cbind(thickness,CBF,FA,MD,myglobal))/BV
brainpreds<-as.matrix(cbind(thickness,CBF,FA,myglobal,bold))/BV
# + I(popul$AgeAtScan^2)
mylm<-bigLMStats( lm( (brainpreds) ~ VIQ + PIQ + RIncome + Ladder + popul$Sex * I(popul$AgeAtScan^1)   ) )
for ( i in 1:nrow(mylm$beta.pval) ) {
  mylm$beta.pval[i,]<-p.adjust( mylm$beta.pval[i,] , method="BH" )
  print( paste(row.names(mylm$beta.pval)[i] , min(mylm$beta.pval[i,] )  ) )
}
colnames( mylm$beta.pval )<-colnames(brainpreds)
mygroups<-c( 1:nrow(mylm$beta.pval) )
braingroups<-c( rep( max(mygroups)+1, ncol(thickness) )
               , rep( max(mygroups)+2, ncol(CBF) )
               , rep( max(mygroups)+3, ncol(FA) ) 
               , rep( max(mygroups)+4, ncol(myglobal) ) 
               , rep( max(mygroups)+5, ncol(bold) ) ) #, rep( max(mygroups)+5, ncol(myglobal) ) )
mygroups<-c( mygroups , braingroups )
mybpcor<-cor( brainpreds )
mybpcor[  mybpcor  < 0.8 ] <- 0
myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Force", outfile="./PediatricTemplateOfBrainPerfusion/results.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0.85, zoom = T , mygroup=mygroups )
                                        # myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp2.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0. , mygroup=mygroups )
#
mypopulp<-cbind(PIQ,VIQ, RIncome, I(popul$AgeAtScan^1) ) # , I(popul$AgeAtScan^2))
maxtest<-6   # kfolds
testind<-2   # which predictor
maxnruns<-10 # number of sub-samples in model building 
myspars<-c( 0.02 )  # select about 1 to 2% of predictors 
nv2<-ncol(mypopulp) # cca sparsevectors
spar2<-( -1.0/nv2 * 0.9  )
whichpreds<-rep(0,ncol(brainpreds))
ntests <- c(1:maxtest)
cvpredictions<-matrix( rep(NA,nrow(brainpreds)*ncol(mypopulp)) , ncol=ncol(mypopulp) )
selector<-(rep(ntests, nrow(popul) ))
selector<-sample(selector, nrow(popul) )
selector<-selector[1:nrow(popul)]
for ( myrun in 1:maxnruns ) {
  for ( whichtotest in 1:(maxtest/2) ) {
    select1<-( selector != whichtotest & selector < (maxtest/2) )
    select1<-sample(select1) 
    select2<-!select1
    colnames(mypopulp)[ncol(mypopulp)]<-"Age"
    mydemogp1<-subset(mypopulp,     select1)
    mydemogp2<-subset(mypopulp,     select2)
    brainpreds1<-subset(brainpreds, select1)
    brainpreds2<-subset(brainpreds, select2 )
    myresults<-c()
    for ( myspar in myspars  ) {
      myz<-( -1 )
      sccan<-sparseDecom2( inmatrix= list(brainpreds1,mydemogp1) , nvecs=nv2, its=55, mycoption=1, sparseness=c(myspar, spar2 ), z=myz )
      mytrainbrain2<- brainpreds1[,] %*% as.matrix( sccan$eig1 )
      mytestbrain2 <- brainpreds2[,] %*% as.matrix( sccan$eig1 )
      mytrainbrain<-cbind(mytrainbrain2)
      mytestbrain<-cbind(mytestbrain2) 
      library(randomForest)
      for ( ind in 1:ncol(mypopulp) ) {
        mydf1        <- data.frame( mycognition=mydemogp1[,ind], imgs = scale( mytrainbrain ), sex=subset(popul$Sex,select1)  ) 
        mydf2        <- data.frame(  imgs = scale( mytestbrain ) , sex=subset(popul$Sex,select2) ) 
        my.rf        <- randomForest( mycognition ~ . , data=mydf1, ntree = 5000 )
        mypred <- predict( my.rf , newdata = mydf2 )
        plot(mypred,mydemogp2[,ind])
        mycor<-cor.test(mypred,mydemogp2[,ind])
        myres<-paste(myspar,mycor$est ,mycor$p.value, mean(abs(mypred-mydemogp2[,ind])) )
        myresults<-c(myresults,myres)
        if ( ind == testind ) {
          ww<-which( abs(sccan$eig1[,which.max(my.rf$importance[1:ncol(mypopulp)] )]) > 0 )
          whichpreds[ww]<-whichpreds[ww]+1
        }
        cvpredictions[select2,ind]<-mypred
      }
    }
                                        # print( myresults )
  }
  par(mfrow=c(1,ncol(mypopulp)))
  for ( ind in 1:ncol(mypopulp) ) {
    plot(cvpredictions[,ind],mypopulp[,ind])
    mycor<-cor.test(cvpredictions[,ind],mypopulp[,ind])
    myres<-paste(myspar,mycor$est ,mycor$p.value, mean(abs(cvpredictions[,ind]-mypopulp[,ind])) )
    print(myres)
  }
  print(paste("Run",myrun))
} # myrun

dev.new()
par(mfrow=c(1,1))
ww<-which( whichpreds > max(whichpreds)*0.1 )
brainpredsSub<-brainpreds[,ww]
select1 <- selector <  (maxtest/2) 
select2 <- selector >= (maxtest/2) 
mydemogp1<-subset(mypopulp,     select1 )
mydemogp2<-subset(mypopulp,     select2 )
globinds<-c(2:ncol(globmess))
globinds<-c(3,5)
brainpreds1<-cbind( subset(brainpredsSub, select1 ) )
brainpreds2<-cbind( subset(brainpredsSub, select2 ) )
mydf2        <- data.frame(  imgs = scale( brainpreds2 )  )
mydf1        <- data.frame( mycognition=mydemogp1[,testind], imgs = scale( brainpreds1 ) )
my.rf        <- randomForest( mycognition ~ . , data=mydf1, ntree = 5000 )
mypred <- predict( my.rf , newdata = mydf2 )
plot(mypred,mydemogp2[,testind])
mycor<-cor.test(mypred,mydemogp2[,testind])
myres<-paste(myspar,mycor$est ,mycor$p.value, mean(abs(mypred-mydemogp2[,testind])) )
print( my.rf$importance/sum(my.rf$importance) )
print(myres)

stop("ASSS")
#################################################
############### now the test data ###############
#################################################
testdatabrain<-scale(brainpreds2) %*% scale(as.matrix(sccan$eig1))
testdatademog<-scale(mydemogp2) %*% scale(t(sccan$eig2))
par(mfrow=c(1,3))
for ( i in 1:3 ) {
  proj1<-sccan$projections[,i]
  proj2<-sccan$projections2[,i] 
  proj1<-testdatabrain[,i]
  proj2<-testdatademog[,i]
  ccc<-round( cor.test( proj1 , proj2 )$est * 100 ) / 100
  mytit<-paste("Eig",as.character(i),"pval",sccan$ccasummary[1,i+1],"cor",ccc )
  print(mytit)
  ww<-which( abs(sccan$eig2[,i]) > 0 )
  print( colnames(mydemogp)[ww] ) 
  ww<-which( abs(sccan$eig1[,i]) > 0 )
  print( colnames(brainpreds)[ww] )
  plot( proj1 ,proj2 , main=mytit ) 
}



stop("AAL above")
# not sure how to achieve this yet
# fused<-joinEigenanatomy( myothk, NA,  eanat$eigenanatomyimages, graphdensity )
#
mycthk<-cbind( residuals( lm(as.matrix(  (demog[,thickinds+thickoff]) )~cthk ) ),  cthk ) 
mybthk<-cbind( residuals( lm(as.matrix(  (demog[,thickinds+ cbfoff ]) )~gcbf ) ), gcbf ) 
mycthk<-cbind( as.matrix(  (demog[,thickinds+thickoff]) )/cthk ,  cthk ) 
mybthk<-cbind( as.matrix(  (demog[,thickinds+ cbfoff ]) )/gcbf , gcbf ) 
nv<-40; spar<-(0.02)
ceanat<-sparseDecom( inmatrix= as.matrix(mycthk) , nvecs=nv, its=3, mycoption=0, sparseness=spar, z=myz, inmask="NA" )
beanat<-sparseDecom( inmatrix= as.matrix(mybthk) , nvecs=nv, its=3, mycoption=0, sparseness=spar, z=myz, inmask="NA" )
cfused<-joinEigenanatomy( mycthk, NA,  ceanat$eigenanatomyimages, 0.1 )
bfused<-joinEigenanatomy( mybthk, NA,  beanat$eigenanatomyimages, 0.1 )
bmyimages <- t( bfused$fusedlist )
cmyimages <- t( cfused$fusedlist )
myocthk<-cbind( as.matrix(  (demog[,thickinds+thickoff ]) ) , cthk ) 
myobthk<-cbind( as.matrix(  (demog[,thickinds+ cbfoff ]) ) , gcbf )
mybthk<-residuals( lm( myobthk ~ myocthk ) ) 
mycthk<-residuals( lm( myocthk ~ myobthk ) ) 
brainpreds<-cbind( mycthk %*% as.matrix( cmyimages ),  mybthk %*% as.matrix( bmyimages ) )
colnames(brainpreds)<-c(paste("TEanat",1:ncol(cmyimages),sep=''),(paste("CEanat",c(1:ncol(bmyimages)),sep='')))
mylm<-bigLMStats( lm( (brainpreds) ~ IQ  + RIncome + Ladder + demog$Sex * I(demog$AgeAtScan^1)  ) )
for ( i in 1:nrow(mylm$beta.pval) ) {
  mylm$beta.pval[i,]<-p.adjust( mylm$beta.pval[i,] , method="BH" )
  print( paste(row.names(mylm$beta.pval)[i] , min(mylm$beta.pval[i,] )  ) )
}
colnames( mylm$beta.pval )<-colnames(brainpreds)
mybpcor<-cor( brainpreds )
mygroups<-c( 1:nrow(mylm$beta.pval) )
braingroups<-c( rep( max(mygroups)+1, ncol(cmyimages) ),
                rep( max(mygroups)+2, ncol(bmyimages) ) )
mygroups<-c( mygroups , braingroups )
pheatmap( mybpcor,cluster_rows=F,cluster_cols=F)
myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp2.html", logvals=T, correlateMyOutcomes = mybpcor, corthresh = 0.5 ) 
myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp3.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0.5 ) 
myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Force", outfile="/Users/stnava/Downloads/temp.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0.5 , mygroup=mygroups, zoom=T) 



stop("Ass")



mythk<-as.matrix( cbind( scale(demog[,thickinds+thickoff]), scale(demog[,thickinds+cbfoff]) ) )
myrthk<-residuals( lm( mythk ~   demog$Sex * demog$AgeAtScan + cthk + gcbf  ) ) 
myrthk<-residuals( lm( mythk ~  cthk + gcbf  ) ) 
nv<-100; spar<-(0.01)
eanat<-sparseDecom( inmatrix= as.matrix(myothk) , nvecs=nv, its=3, mycoption=0, sparseness=spar, z=myz, inmask="NA" )
gds<-c(1:100)/200
gdvecs <-c()
for ( graphdensity in gds ) {
  myimages <- eanat$eigenanatomyimages
  fused<-joinEigenanatomy( myrthk, NA,  eanat$eigenanatomyimages, graphdensity )
  gdvecs<-c(gdvecs,dim(fused$fusedlist)[1])
}
plot(gds,gdvecs,type='l')
fused<-joinEigenanatomy( myrthk, NA,  eanat$eigenanatomyimages, 0.05 )
myimages <- t( fused$fusedlist )
# myimages<-eanat$eigenanatomyimages
tempinc<-cbind(  as.numeric( myiq )  ,  as.numeric( mylad ) , rank(myincome)  ) #,
myglobalpreds<-cbind( cthk, gcbf )
mylm<-bigLMStats( lm( myglobalpreds ~ tempinc   + demog$Sex * demog$AgeAtScan  + I(demog$AgeAtScan^2)  ) )
print( mylm$beta.pval )
brainpreds<-myothk %*% as.matrix( myimages )
mylm<-bigLMStats( lm( brainpreds ~ tempinc   + demog$Sex * demog$AgeAtScan  + I(demog$AgeAtScan^2) + gcbf + cthk ) )
print( paste( min(p.adjust( mylm$beta.pval[1,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[1,] , method='BH') )  ) )
print( paste( min(p.adjust( mylm$beta.pval[2,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[2,] , method='BH') )  ) )
print( paste( min(p.adjust( mylm$beta.pval[3,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[3,] , method='BH') )  ) )
ind<-3 ; myimages[which( abs(myimages[,ind]) > 0 ),ind]
print(which( abs(myimages[,ind]) > 0 ))
aal$label_name[thickinds][ which(abs( myimages[,ind]) > 0 ) %% length(thickinds) ]
brainpreds<-myothk %*% as.matrix( myimages )
mylm<-bigLMStats( lm( brainpreds ~ tempinc   + demog$Sex * demog$AgeAtScan  + I(demog$AgeAtScan^2)  ) )
colnames(brainpreds)<-paste("Vox",c(1:ncol(brainpreds)),sep='')
colnames( mylm$beta.pval )<-colnames(brainpreds)
mybpcor<-cor( brainpreds )
mybpcor[  mybpcor  < 0.9 ] <- 0
myd3<-regressionNetworkViz( mylm , sigthresh=0.01, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp2.html", logvals=T, correlateMyOutcomes = mybpcor, corthresh = 0.9 ) 
myd3<-regressionNetworkViz( mylm , sigthresh=0.01, whichviz="Force", outfile="/Users/stnava/Downloads/temp.html", logvals=T, correlateMyOutcomes = mybpcor, corthresh = 0.9 ) 


# now residualize & decompose again 
brainpreds<-myrthk %*% as.matrix( myimages )
mylm<-bigLMStats( lm( brainpreds ~ tempinc    + I(demog$AgeAtScan^2)  + demog$Sex * demog$AgeAtScan ) )
print( paste( min(p.adjust( mylm$beta.pval[1,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[1,] , method='BH') )  ) )
print( paste( min(p.adjust( mylm$beta.pval[2,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[2,] , method='BH') )  ) )
print( paste( min(p.adjust( mylm$beta.pval[3,] , method='BH') ), which.min( p.adjust( mylm$beta.pval[3,] , method='BH') )  ) )
print(dim(myimages))
ind<-5 ; myimages[which( abs(myimages[,ind]) > 0 ),ind]
print(which( abs(myimages[,ind]) > 0 ))
aal$label_name[thickinds][ which(abs( myimages[,ind]) > 0 ) %% length(thickinds) ]
pheatmap( cor(brainpreds) ,cluster_rows=F,cluster_cols=F)
################ now some cool visualization with d3 ##########################################
colnames(brainpreds)<-paste("Vox",c(1:ncol(brainpreds)),sep='')
colnames( mylm$beta.pval )<-colnames(brainpreds)
regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp2.html") 
regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Force", outfile="/Users/stnava/Downloads/temp.html") 
################################################################################################################
stop("UseEanatAbove")
nv<-3
if ( ! exists("np") ) np<-25
sparval<-( -0.05 ) * length( thickinds ) / length( mystudyinds )
sccan<-sparseDecom2( inmatrix=list( mythk, as.matrix(tempinc) ) , nvecs=nv, robust=0,
                    its=29, mycoption=1 ,  perms=np, sparseness=c( sparval , 0.3 ), z=-1, ell1=11 ) # cbf 
ww<-which( abs(sccan$eig1[,1]) > 0 )
print( ww )
print( sccan$eig1[ww,1] ) 
print( aal$label_name[thickinds][ ww %% length(thickinds) ] )
print("&")
ww<-which( abs(sccan$eig1[,2]) > 0 )
print( ww )
print( sccan$eig1[ww,2] )
print( aal$label_name[thickinds][ ww %% length(thickinds) ] )
print("&")
ww<-which( abs(sccan$eig1[,3]) > 0 )
print( ww )
 print(sccan$eig1[ww,3] )
print( aal$label_name[thickinds][ ww %% length(thickinds) ] )
mypro1<- sccan$projections  # mythk %*% as.matrix(sccan$eig1)
mypro2<- sccan$projections2 # tempinc %*% as.matrix(sccan$eig2)
print(sccan$eig2)
if ( length( dev.list() ) == 0 ) dev.new() # dev.off()
for ( i in 1:nv ) {
  ccc<-round( cor.test(mypro1[,i] ,mypro2[,i])$est * 100 ) / 100
  plot( mypro1[,i], mypro2[,i], main=paste("Eig",as.character(i),"pval",sccan$ccasummary[1,i+1],"cor",ccc ) ) 
  Sys.sleep(3)
}
stop("Unfinished work below")
mythk<-as.matrix( cbind( demog[,mystudyinds] ) )
# mythk<-residuals( lm( mythk ~  demog$Sex * demog$AgeAtScan + I(demog$AgeAtScan^2) + cthk + gcbf ) ) 
mysvd<-svd(mythk)
mycor<-cor( mythk, (mysvd$u) )
pheatmap( mycor ,cluster_rows = F, cluster_cols = F )
mymeans<-apply(abs(mycor),FUN=mean,M=2)
mycbfmean1<- mean(abs(mycor[(length(thickinds)+1):(2*length(thickinds)),1]))
mythkmean1<- mean(abs(mycor[     1               :(1*length(thickinds)),1]))
mycbfmeans<-apply(abs(mycor[(length(thickinds)+1):(2*length(thickinds)),]),FUN=mean,M=2)
mythkmeans<-apply(abs(mycor[1               :(1*length(thickinds)),]),FUN=mean,M=2)
plot(mycbfmeans,type='l')
points(mythkmeans,type='l',col='red')
stop("Unfinished work below")
par(mfrow=c(1,2))
# demog$Full.4.IQ[ demog$Full.4.IQ == 0 ]<-mean( demog$Full.4.IQ[ demog$Full.4.IQ > 0 ] ,na.rm=T)
 demog$AgeAtScan<-as.numeric( as.character( demog$AgeAtScan ) )
 outc<-as.numeric( as.character( demog$Teen.Ladder.SES.score ) )
# outc<-as.numeric( as.character( demog$Teen.Ladder.Community.Score ) )
# outc<-demog$Income
# outc<-as.numeric( demog$Full.Scale.IQ )
outc[ is.na( outc ) ]<-median( outc ,na.rm=T)
pv<-rep(NA,lgind )
bvol[ is.na(bvol) ]<-mean(bvol,na.rm=T)
if ( mean(  demog$AgeAtScan ) > 100 ) demog$AgeAtScan<-demog$AgeAtScan/360
nstudies<-4
peakageresults1<-matrix( rep( 0, nstudies*lgind ) , nrow = lgind )
peakageresults2<-matrix( rep( 0, nstudies*lgind  ) , nrow = lgind )
peakageresults<-matrix( rep( 0, nstudies*lgind  ) , nrow = lgind )
ct<-1
# for ( whichstudy in c( "Thickness" , "rCBF" , "CBF", "tCBF" ) ) {
for ( whichstudy in c( "Thickness"  ) ) {
  mysel<-( bvol > 0.5 & demog$AgeAtScan < 25 ) # get rid of outlier data 
  data('aal',package='ANTsR')
  aalimg<-antsImageRead('template/Labels/aal.nii.gz',3)
  aalvec<-aalimg > 0
  aalvals<-aalimg[ aalvec ]
  aalimg[ aalimg > 90 ]<-0 
  for ( ind in c(1:lgind) ) {
  if ( whichstudy == "rCBF" ) {
  # relative CBF
  rform<-paste("(CBF_AAL.",gind[ind],"/gcbf) ~   1 + Thickness_AAL.AAL",gind[ind],"",sep='')
  rform<-paste("(CBF_AAL.",gind[ind],"/gcbf) ~   1 ",sep='')
  rcbf<-winsor( residuals( lm( rform , data = demog , subset=mysel) ) , tr=0.05 )
  form1<-paste(" rcbf ~   1 + Sex + bvol ",sep='')
  form2<-paste(" rcbf ~   1 + Sex + bvol + AgeAtScan + I(AgeAtScan^2) ",sep='')
  form2b<-paste("rcbf ~   1 + Sex + bvol + I(AgeAtScan^2)  ",sep='')
  }
  if ( whichstudy == "CBF" ) {
  # raw CBF 
  rform<-paste("(CBF_AAL.",gind[ind],") ~   1 ",sep='')
  rcbf<-winsor( residuals( lm( rform , data = demog , subset=mysel) ) , tr=0.05 )
  form1<-paste(" rcbf ~   1 + Sex + bvol ",sep='')
  form2<-paste(" rcbf ~   1 + Sex + bvol + AgeAtScan + I(AgeAtScan^2)   ",sep='')
  form2b<-paste("rcbf ~   1 + Sex + bvol + I(AgeAtScan^2)   ",sep='')
  }
  if ( whichstudy == "tCBF" ) {
  # thickness adjusted CBF 
  rform<-paste("(CBF_AAL.",gind[ind],") ~   1  + Thickness_AAL.AAL",gind[ind],"",sep='')
  rcbf<-winsor( residuals( lm( rform , data = demog , subset=mysel) ) , tr=0.05 )
  form1<-paste(" rcbf ~   1 + Sex + bvol ",sep='')
  form2<-paste(" rcbf ~   1 + Sex + bvol + AgeAtScan + I(AgeAtScan^2)   ",sep='')
  form2b<-paste("rcbf ~   1 + Sex + bvol + I(AgeAtScan^2)   ",sep='')
  }
  if ( whichstudy == "Thickness" ) {
  # cortical thickness
  form1<-paste("Thickness_AAL.AAL",gind[ind]," ~   1 + Sex ",sep='')
  form2<-paste("Thickness_AAL.AAL",gind[ind]," ~   1 + Sex + AgeAtScan + I(AgeAtScan^2) ",sep='')
  form2b<-paste("Thickness_AAL.AAL",gind[ind]," ~  1 + Sex + I(AgeAtScan^2) ",sep='')
  form1<-paste("Thickness_AAL.AAL",gind[ind]," ~   1 + Sex+ AgeAtScan + I(AgeAtScan^2) ",sep='')
  form2<-paste("Thickness_AAL.AAL",gind[ind]," ~   1 + Sex+ AgeAtScan + I(AgeAtScan^2) + Income ",sep='')
  form2b<-form2
  }
  # select the model
  tempdf<-data.frame( demog, bvol=bvol )
  preddf<-subset( tempdf,  mysel )
  mdl1<-lm( as.formula(form1) , data = preddf  ) 
  mdl2<-lm( as.formula(form2) , data = preddf  )
  mdl2b<-lm( as.formula(form2b) , data = preddf  )
  print(paste("AAL:", ind , aal$label_name[ind] ))
  print( summary( mdl2 ) )
}
  adaffdafadfdaf
  {
  pv[ ind ] <-  anova( mdl1, mdl2 )$P[2]
  isquad<-TRUE
  if ( summary(mdl2)$fstatistic[1] > summary(mdl2b)$fstatistic[1] ) { isquad<-FALSE } #; mdl2<-mdl2b }
  newages <-  c(1:(20*5))/5 
  myagepred<-data.frame( AgeAtScan = newages , bvol = rep( mean(bvol) , length(newages) ), Sex = rep( preddf$Sex[1], Income = rep( mean(demog$Income) , length(newages) ) , length(newages) )  )
  mycbfpred<-predict( mdl2 , newdata=myagepred )
  peakage<-myagepred$AgeAtScan[ which.max( mycbfpred ) ]
  if ( peakage < 1 | peakage > 19.9 ) peakage<-myagepred$AgeAtScan[ which.min( mycbfpred ) ]
  peaks<-paste( "Peak-Age-max",  myagepred$AgeAtScan[ which.max( mycbfpred ) ],"Peak-Age-min",  myagepred$AgeAtScan[ which.min( mycbfpred ) ],"peakage",  peakage , isquad , qv[ind] )
  peakageresults1[ ind , ct ]<-myagepred$AgeAtScan[ which.max( mycbfpred ) ]
  peakageresults2[ ind , ct ]<-myagepred$AgeAtScan[ which.min( mycbfpred ) ]
  peakageresults[ ind , ct ] <- peakage
  aalimg[ aalimg == as.numeric(ind)  ] <- peakage # 100 * ( 1 - pv[ ind ] ) #
  myoutstring <- paste( aal$label_name[gind[ind]], pv[ind] , peaks )
  print( summary( mdl2 ) )
  tagger<-""
  if ( pv[ ind ] < 0.01 ) tagger<-"*"
  plot(  myagepred$AgeAtScan , mycbfpred , type='l' ,main= peaks)
  mytit<-paste(whichstudy,aal$label_name[gind[ind]],tagger ) 
  visreg( mdl2 , xvar="AgeAtScan", main=mytit)
  pdf(paste( paste("figs/",whichstudy,"_",sep=''),aal$label_name[gind[ind]],"_ct_%03d",".pdf",sep='') , onefile = FALSE )
  visreg( mdl2 , xvar="AgeAtScan", main=mytit )
#  plot(  myagepred$AgeAtScan , mycbfpred , type='l' ,main= peaks)
  dev.off()
}
#
qv<-p.adjust( pv , method='BH' )
print("")
print( paste(whichstudy, "Final Corrected q-values" , sum( qv < 0.0500001 ) ) )
print("")
for ( ind in c(1:lgind) ) {
  if ( qv[ ind ] < 0.05 ) print( paste( aal$label_name[gind[ind]], ind, qv[ind], "peak", peakageresults[ ind , ct ] ) )
}
# pheatmap( peakageresults1 , cluster_rows = F, cluster_cols = F)
antsImageWrite( aalimg, paste('aal_peak_age',whichstudy,'.nii.gz',sep='') )
ct<-ct+1  
}




######################################
ccahello <- function( sccanx , namer="Study" ) {
  eps<-1.e-8
  for ( ind in 1:nv ) {
    print(paste("Sccanvec",ind,"pvalue",sccanx$ccasummary[1,ind+1],"corr",sccanx$ccasummary[2,ind+1]))
    print( paste( colnames(iqin)[ abs(sccanx$eig2[,ind]) > eps ] ) )
    print( ( sccanx$eig2[,ind] )[ abs(sccanx$eig2[,ind]) > eps] / sum(  abs(sccanx$eig2[,ind])  ) )
    print("Imaging Predictors")
    myanat<-aal$label_name[1:90][ abs(sccanx$eig1[,ind]) > eps ]
    wt2<-( sccanx$eig1[,ind] )[ abs(sccanx$eig1[,ind]) > eps ] 
    for ( j in 1:length(myanat) )  {
      selector<-myanat[j]
      print(  paste( selector , wt2[j]/sum(abs(wt2))  , sum(abs(wt2[1:j]))/sum(abs(wt2)) , which( aal$label_name == myanat[j] ) ) )
    }
  }
  print( paste( "Adjusted p-values:",namer)  )
  print( p.adjust( sccanx$ccasummary[1,2:(ncol(sccanx$ccasummary)-1)] , method="BH" ) )
######################################
}

####### STANDARDIZE THE DATA #######
thk<-as.matrix( preddf[,9:98] )
cbf<-thk # as.matrix( preddf[,99:188] )
thk<-residuals( lm( thk ~ bvol                  , data  = preddf ) ) 
cbf<-residuals( lm( cbf ~ bvol , data = preddf  ) )
iqin<-as.matrix( preddf[,c(6,478:483)] )
iq<-matrix( as.numeric(iqin),nrow=nrow(cbf))
iq<-impute( iq )
colnames( iq ) <- colnames( iqin )

##########################################
#### FORM TESTING AND TRAINING SAMPLES ###
##########################################
mysel<-rnorm( nrow( preddf ) ) < 0.5
## TRAIN set
traindf<-subset( preddf, mysel )
cbf1<-subset( cbf, mysel )
thk1<-subset( thk, mysel )
iq1<-subset( iq , mysel )
## test set
testdf<-subset( preddf, !mysel )
cbf2<-subset( cbf, !mysel )
thk2<-subset( thk, !mysel )
iq2<-subset( iq ,  !mysel )
####
nperm<-100
myrob<-0
nv<-ncol(iq)/2
myspar<-c( 0.1 , -1/ncol(iq) )
myspar<-c( 0.1 , -0.05/nv )
#######
sccanc <- sparseDecom2( inmatrix=list( cbf1 , iq1 ), inmask = c( NA , NA ) , mycoption = 0, sparseness=myspar, nvecs=nv, its=20, smooth=0, cthresh = c(0, 0), robust=myrob , perms=nperm )
m1 <- as.matrix( sccanc$eig1 ) # ;# m1[ m1 < 0.0001 ]<-0
m2 <- as.matrix( sccanc$eig2 ) # ; m2[ m2 < 0.0001 ]<-0
proj1 <-  cbf2  %*% m1
proj2 <-  iq2  %*% m2
ccahello( sccanc, 'CBF' )
for ( i in 1:ncol(proj1) ) print( cor.test( proj1[,i] , proj2[,i] )$p.value ) 
####################################################################################
####################################################################################
sccant <- sparseDecom2( inmatrix=list( thk1 , iq1 ), inmask = c( NA , NA ) , mycoption = 0, sparseness=myspar, nvecs=nv, its=20, smooth=0, cthresh = c(0, 0), robust=myrob , perms=nperm )
m1 <- as.matrix( sccant$eig1 )
m2 <- as.matrix( sccant$eig2 )
# m1[ m1 < 0.01 ]<-0
# m2[ m2 < 0.01 ]<-0
proj1 <-  thk2  %*% m1
proj2 <-  iq2  %*% m2
sccant$eig1<-m1
sccant$eig2<-m2
ccahello( sccant )
for ( i in 1:ncol(proj1) ) print( cor.test( proj1[,i] , proj2[,i] )$p.value ) 


