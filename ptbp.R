library(ANTsR)
data('aal',package='ANTsR')
library(pheatmap)
library(e1071)
thickinds<-c(1:36,39:40,43:70,79:90)
thickinds<-c(1:70,79:90)
demog<-read.csv("AAL_demo_thickness.csv")
cbf<-read.csv("AAL_meancbf.csv")
demog<-cbind( demog, cbf[,4:ncol(cbf)] )
thickoff<-23
cthk<-apply( demog[  , thickinds+thickoff ], MARGIN=1, FUN=mean )
cbfoff<-23+102
gcbf<-apply( demog[  , thickinds+cbfoff ], MARGIN=1, FUN=mean )
demog$Sex[ demog$Sex == "F " ]<-"F"
usesubs<-rep(TRUE,nrow(demog))
usesubs[ c(155,142) ] <- FALSE      # for thickness
usesubs[ c(155,142,which( is.na(gcbf) ),which( gcbf < 10  )) ] <- FALSE  # for cbf
tempdf<-data.frame( cbind( demog, usesubs ) )
write.csv(tempdf,"good_v_bad_cbf_subjects.csv",quote=F,row.names=F)
uids <- unique( demog$SubID )
for ( i in 1:length(uids) )
  {
  ww<-which( demog$SubID == uids[i] )
  if ( length( ww ) > 1 ) usesubs[ ww[2:length(ww)] ] <- FALSE
  }
demog<-demog[usesubs,]
cthk<-apply( demog[  , thickinds+thickoff ], MARGIN=1, FUN=mean )
gcbf<-apply( demog[  , thickinds+cbfoff ], MARGIN=1, FUN=mean )
inclevs<-levels( demog$Income )
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
levels( demog$Income )<-( inclevs )
demog$Income<-as.numeric( as.character( demog$Income ) )
demog<-data.frame( demog , myincome=(myincome) , gcbf=gcbf, tgcbf=(gcbf/cthk) )
residagainstthickness<-FALSE 
if ( residagainstthickness ) demog[  , thickinds+cbfoff ]<-as.matrix( residuals( lm( as.matrix(demog[  , thickinds+cbfoff ]) ~ cthk ) ) )
mdl<-lm( gcbf ~ AgeAtScan * Sex + myincome  , data=mydf )
print( summary( mdl  ) )
pdf("~/Downloads/cbf_v_income.pdf")
visreg( mdl , "myincome", main="CBF vs Income" )
dev.off()
####################################################################################################
demog$AgeAtScan <- as.numeric(as.character(demog$AgeAtScan))
myincome<-c( impute( cbind(demog$blank,demog$Income ) ) )
myiq<-( c( impute( cbind( demog$blank, as.numeric(as.character(demog$FullScaleIQ)) ) ) ) )
myiq2<-( c( impute( cbind( demog$blank, as.numeric(as.character(demog$Verbal.IQ)) ) ) ) )
mylad<-c( impute( cbind( demog$blank, as.numeric(as.character(demog$Teen.Ladder.SES.score ) )   ) ) )
myladc<-c( impute( cbind( demog$blank, as.numeric(as.character(demog$Teen.Ladder.Community.Score ) )   ) ) )
####################################################################################################
myoffset<-cbfoff
myoffset<-thickoff
if ( FALSE ) {
# look at k = 5 !
for ( k in 1:length(thickinds) ) 
  {
  mdl<-lm( demog[,thickinds[k]+myoffset] ~ 1 + Sex*AgeAtScan +    I(AgeAtScan^2) + myiq  , data = demog  )
  dd<-stepAIC( mdl , direction = c("both") , trace =  0 )
  print(summary( lm( formula(dd) , data=demog ) ) )
  print( paste( k, aal$label_name[thickinds][k] ) )
  Sys.sleep(3)
  }
}
####################################################################################################
myincome2<-myincome
myincome2[ myincome >= 120 ]<-6
myincome2[ myincome <= 120 ]<-5
myincome2[ myincome <= 90 ]<- 4
myincome2[ myincome <= 70 ]<- 3
myincome2[ myincome <= 50 ]<- 2
myincome2[ myincome <= 20 ]<- 1
bvol<-apply( demog[  , thickinds+myoffset ], MARGIN=1, FUN=mean )
mythk<-as.matrix( cbind( demog[,thickinds+myoffset] ) )
mythk<-residuals( lm( mythk ~  demog$Sex * demog$AgeAtScan + cthk ) ) 
tempinc<-cbind(  as.numeric( myiq ) ,as.numeric( myiq2 ) , as.numeric( myincome ) , as.numeric(myincome2), as.numeric( mylad ) , myladc  )
nv<-3
if ( ! exists("np") ) np<-2500
sparval<-( 0.1 )
# if ( myoffset == thickoff ) sparval<-0.05
sccan<-sparseDecom2( inmatrix=list( mythk, as.matrix(tempinc) ) , nvecs=nv, robust=0,
                    its=40, mycoption=1 ,  perms=np, sparseness=c( sparval , 0.3 ), z=1, ell1=11 ) # cbf
print( aal$label_name[thickinds][ abs(sccan$eig1[,1]) > 0 ] )
print("&")
print( aal$label_name[thickinds][ abs(sccan$eig1[,2]) > 0 ] )
print("&")
if ( nv > 2 ) print( aal$label_name[thickinds][ abs(sccan$eig1[,3]) > 0 ] )
mypro1<- mythk %*% as.matrix(sccan$eig1)
mypro2<- tempinc %*% as.matrix(sccan$eig2)
print(sccan$eig2)
if ( length( dev.list() ) == 0 ) dev.new() # dev.off()
for ( i in 1:nv ) {
  ccc<-round( cor.test(mypro1[,i] ,mypro2[,i])$est * 100 ) / 100
  plot( mypro1[,i], mypro2[,i], main=paste("Eig",as.character(i),"pval",sccan$ccasummary[1,i+1],"cor",ccc ) ) 
  Sys.sleep(3)
}
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


