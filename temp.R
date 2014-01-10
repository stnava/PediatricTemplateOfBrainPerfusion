library(ggplot2)
library(randomForest)
library(ANTsR)
library(visreg)
library(psych)
library(pheatmap)
library(e1071)
votedbrainregions1<-rep(0,293) # know this is n-variables , name later
votedbrainregions2<-rep(0,293)
collectresults<-matrix( rep(NA,2000), ncol=4 )
for ( kk in 1:nrow(collectresults) ) {
    mylist<-ls()
    keeplist<-c( "collectresults","kk","votedbrainregions1","votedbrainregions2")
    for ( myk in keeplist ) mylist<-mylist[ mylist != myk ]
    rm(list=mylist)
    gc()
doboot<-T
data('aal',package='ANTsR')
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
colnames( CBF )<-c( paste("CBF",aal$label_name[thickinds]) )
FA<-winsor( demog[,grep("FA",names(demog))] , trim = 0.02 )
MD<-winsor( demog[,grep("MD_",names(demog))], trim = 0.02 )
RD<-winsor( demog[,grep("RD_",names(demog))], trim = 0.02 )
globmess<-read.csv("global_measures.csv")
globmess<-globmess[,2:ncol(globmess)]
usebold <- TRUE
if ( usebold ) {
  boldin<-read.csv("AAL_BoldConn.csv")
#  boldin<-boldin[,c(95:ncol(boldin))]
  boldin<-boldin[,c(5:94)] # degree 
  boldin<-boldin[,thickinds]
  colnames( boldin )<-c( paste("Bold",aal$label_name[thickinds]) )
}
usesubs<-rep(TRUE,nrow(demog))
usesubs[ c(155,142) ] <- FALSE      # for thickness
usesubs[ c(155,142,which( is.na(globmess$GrayCBF) ),which( globmess$GrayCBF < 10  ), which( is.na(jhu$RD_JHU.43) ) )] <- FALSE  # for cbf
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
RD<-RD[usesubs,]
if ( usebold ) {
  bold<-boldin[usesubs,]
  gbold<-apply(as.matrix(bold),FUN=mean,MARGIN=1)
  bold<-winsor( impute(bold), trim=0.05 )
  gbold<-apply( bold ,FUN=mean,MARGIN=1)
}
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
FIQ      <- myiq
VIQ      <- myiq2
PIQ      <- myiq3
Ladder  <- ( mylad )
RIncome <- rank(myincome)
gbold<-apply(as.matrix(bold),FUN=mean,MARGIN=1)
gcbf<-apply(as.matrix(CBF),FUN=mean,MARGIN=1)
BV<-globmess$BrainVolume
globmess<-cbind( globmess, gbold )
myglobal<-globmess
mygroups<-c( 1:6 )  # nrow(mylm$beta.pval) )
if ( ! usebold ) {
brainpreds<-as.matrix(cbind(thickness,CBF,FA,MD))/BV
braingroups<-c( rep( max(mygroups)+1, ncol(thickness) )
               , rep( max(mygroups)+2, ncol(CBF) )
               , rep( max(mygroups)+3, ncol(FA) ) 
               , rep( max(mygroups)+4, ncol(MD) ) ) #, rep( max(mygroups)+5, ncol(myglobal) ) )
} else {
brainpreds<-as.matrix(cbind(thickness/myglobal$GrayThickness,CBF/myglobal$GrayCBF,MD/myglobal$WhiteMD,bold/gbold,myglobal))
brainpreds<-as.matrix(cbind(thickness/myglobal$GrayThickness,CBF/BV,MD/BV,bold/BV,myglobal))
brainpreds<-as.matrix(cbind(thickness/BV,                    CBF/BV,MD/BV,bold/BV,myglobal))
braingroups<-c( rep( max(mygroups)+1, ncol(thickness) )  
               , rep( max(mygroups)+2, ncol(CBF)  )  
               , rep( max(mygroups)+3, ncol(FA) ) 
               , rep( max(mygroups)+4, ncol(bold)) , rep( max(mygroups)+5, ncol(myglobal) ) ) 
# brainpreds<-bold # as.matrix(residuals( lm( bold ~ gbold) ) )
# braingroups<-c( rep( max(mygroups)+1, ncol(brainpreds) )  )
}
whichmodality<-157:204
whichmodality<-1:ncol(brainpreds)
mylm<-bigLMStats( lm( (brainpreds) ~  PIQ + VIQ   + Ladder  + popul$Sex * I(popul$AgeAtScan^1)  ) ) # + I(popul$AgeAtScan^2)  ) )
for ( i in 1:nrow(mylm$beta.pval) ) {
  locpv<-p.adjust( mylm$beta.pval[i,whichmodality] , method="BH" )
  print( paste(row.names(mylm$beta.pval)[i] , min( locpv )  ) )
  siglocpv<-locpv[ locpv <= 0.05 ]
  signifregions<-names( siglocpv )
  signifregions<-signifregions[ order( siglocpv ) ]
  siglocpv<-siglocpv[ order( siglocpv ) ]
  if ( min( siglocpv ) < 0.05 ) {
      print( paste( signifregions , collapse=" / ") )
      signifregions<-factor(signifregions,levels=(signifregions),ordered=TRUE)
      dfn<-data.frame( OutcomeNames=signifregions , LogQVals=log(siglocpv))
      tit<-paste("Log Q-values for Significant Regions:",row.names(mylm$beta.pval)[i])
      pdf(paste(row.names(mylm$beta.pval)[i],'_qvals.pdf',sep=''),width=7,height=4)
      print( ggplot(dfn, aes(x=OutcomeNames,y=LogQVals) ) + ggtitle(tit)
            + theme(plot.title = element_text(size = rel(2)))
            + geom_boxplot() +  theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=1)) )
      dev.off()
  }
}

# make some heatmaps of the data 
pdf('heatmap_thick.pdf',width=7)
pheatmap( cor( thickness ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
pdf('heatmap_cbf.pdf',width=7)
pheatmap( cor( CBF ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
pdf('heatmap_MD.pdf',width=7)
pheatmap( cor( MD ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
pdf('heatmap_bold.pdf',width=7)
pheatmap( cor( bold ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
brainpredssub<-as.matrix(cbind(thickness/BV,CBF/BV,MD/BV,bold/BV ))
pdf('heatmap_all1.pdf',width=7)
pheatmap( cor( brainpredssub ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
brainpredssub<-as.matrix(cbind(thickness/myglobal$GrayThickness,CBF/globmess$GrayCBF,MD/globmess$WhiteMD,bold/gbold,myglobal))
pdf('heatmap_all2.pdf',width=7)
pheatmap( cor( brainpredssub ) , cluster_rows = F, cluster_cols = F,show_colnames = T,show_rownames = T, fontsize=6 )
dev.off()
    
colnames(votedbrainregions1)<-names(brainpreds)
colnames(votedbrainregions2)<-names(brainpreds)
colnames( mylm$beta.pval )<-colnames(brainpreds)
mybpcor<-cor( brainpreds )
mybpcor[  mybpcor  < 0.8 ] <- 0
mygroups<-c( 1:nrow(mylm$beta.pval) )
mygroups<-c( mygroups , braingroups )
myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Force", outfile="./PediatricTemplateOfBrainPerfusion/results.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0.85, zoom = T , mygroup=mygroups )
lddf<-data.frame( thick=brainpreds[,15],  PIQ=PIQ , VIQ=VIQ   , Ladder=Ladder ,Sex=popul$Sex , age=popul$AgeAtScan )
lddf
mdl<-lm( thick ~ .,data=lddf )
# print( summary( mdl ) )
# visreg(mdl)
stop("Regression")
######
#  myd3<-regressionNetworkViz( mylm , sigthresh=0.05, whichviz="Sankey", outfile="/Users/stnava/Downloads/temp2.html", logvals=T, correlateMyOutcomes = NA, corthresh = 0. , mygroup=mygroups )
#######
dVIQPIQ<-(PIQ-VIQ)-min(PIQ-VIQ)
mypopulp<-data.frame( PIQ=PIQ, VIQ=VIQ , Gen=as.factor(popul$Sex) ,Ladder=Ladder   ,  I(popul$AgeAtScan)) #,, I(popul$AgeAtScan^2) )
nv2<-ncol(mypopulp) # cca
myz<-( -1 )
myspars<-c( 0.01 )
spar2<-( -0.5 )
selector<-(rep(c(1:2), nrow(popul) ))
selector<-selector[1:nrow(popul)]
testind<-1
if ( doboot ) {
  selector<-sample(selector)
  select1<-selector != 1 
  select2<-selector == 1
  mydemogp1<-subset(mypopulp,     select1)
  mydemogp2<-subset(mypopulp,     select2)
  mycompv <- t.test( mydemogp1[,testind], mydemogp2[,testind] )$p.value 
  mycompv2 <- t.test( mydemogp1[,ncol(mydemogp2)], mydemogp2[,ncol(mydemogp2)] )$p.value 
  while ( (  mycompv2 < 0.5 ) ) {
    selector<-sample(selector)
    select1<-selector != 1 
    select2<-selector == 1
    mydemogp1<-subset(mypopulp,     select1)
    mydemogp2<-subset(mypopulp,     select2)
    #mycompv <- t.test( mydemogp1[,testind], mydemogp2[,testind] )$p.value 
    mycompv2 <- t.test( mydemogp1[,ncol(mydemogp2)], mydemogp2[,ncol(mydemogp2)] )$p.value 
  }
}
select1<-selector != 1 
select2<-selector == 1
mydemogp1<-subset(mypopulp,     select1)
mydemogp2<-subset(mypopulp,     select2)
brainpreds1<-subset(brainpreds, select1)
bold1<-subset(bold, select1)
brainpreds2<-subset(brainpreds, select2 )
myresults<-c()
# par(mfrow=c(1,ncol(mypopulp)))
for ( myspar in myspars  ) {
  mymat<-data.matrix(mydemogp1)
  mymats<-list((brainpreds1),(mymat))
  sccan<-sparseDecom2boot( inmatrix= mymats , perms=0,nvecs=nv2, its=10, mycoption=1, sparseness=c(myspar, spar2 ), z=myz , robust=0 , nboot = 60 , nsamp = (2.0/3.0) )
  eigpreds<-sccan$eig1
  sccan2<-sparseDecom2(    inmatrix= mymats, perms=0,nvecs=nv2, its=10, mycoption=1, sparseness=c(myspar, spar2 ), z=myz , robust=0 )
  eigpreds2<-sccan2$eig1
  mytrainbrain2<- brainpreds1[,] %*% as.matrix( eigpreds2[,] )
  mytestbrain2 <- brainpreds2[,] %*% as.matrix( eigpreds2[,] )
  mytrainbrain<-  brainpreds1[,] %*% as.matrix( eigpreds[,] )
  mytestbrain<-   brainpreds2[,] %*% as.matrix( eigpreds[,] )
  for ( ind in testind ) {
    mydf1        <- data.frame( imgs = scale( mytrainbrain ), mycognition=(mydemogp1[,ind]) )
    mydf2        <- data.frame( imgs = scale( mytestbrain  ) ) 
    my.rf        <- randomForest( mycognition ~ . , data=mydf1 ) #, ntree = 5000 )
    mypred <- predict( my.rf , newdata = mydf2 )
#    plot(mypred,mydemogp2[,ind])
    if (typeof(mydemogp2[,ind]) != "integer" ) {
      mycor<-cor.test(mypred,mydemogp2[,ind])
      myp1<-mycor$p.value
      myerr1<-mean(abs(mypred-as.numeric(mydemogp2[,ind])))
    } else {
      mycor<-list(est=NA)
      myerr1<-sum( mydemogp2[,ind] == mypred )/length(mypred)*100
      myerr<-paste("M/F",myerr1,"%")
    }
    myres<-paste(names(mypopulp)[ind],myspar,mycor$est ,mycor$p.value, myerr1 )
    myresults<-c(myresults,myres)
    if ( ind == testind ) {
      my.rf$importance<-( my.rf$importance/sum(my.rf$importance) )
      print( my.rf$importance )
      wmax<-which.max(my.rf$importance[1:ncol(mypopulp)] )
      for(wmax in 1:nv2 ) {
        print(paste("EigRegions",wmax))
        ww<-which( abs(eigpreds[,wmax]) > 0 )
        print( colnames( brainpreds1 )[ww] )
        votedbrainregions1[ww]<-votedbrainregions1[ww]+1
      }
    }
  }

    for ( ind in testind ) {
    mydf1        <- data.frame( imgs = scale( mytrainbrain2 ), mycognition=(mydemogp1[,ind]) )
    mydf2        <- data.frame( imgs = scale( mytestbrain2  ) ) 
    my.rf        <- randomForest( mycognition ~ . , data=mydf1) #, ntree = 5000 )
    mypred <- predict( my.rf , newdata = mydf2 )
#    plot(mypred,mydemogp2[,ind])
    if (typeof(mydemogp2[,ind]) != "integer" ) {
      mycor<-cor.test(mypred,mydemogp2[,ind])
      myp2<-mycor$p.value
      myerr2<-mean(abs(mypred-as.numeric(mydemogp2[,ind])))
    } else {
      mycor<-list(est=NA)
      myerr2<-sum( mydemogp2[,ind] == mypred )/length(mypred)*100
      myerr<-paste("M/F",myerr2,"%")
    }
    myres<-paste(names(mypopulp)[ind],myspar,mycor$est ,mycor$p.value, myerr2 )
    myresults<-c(myresults,myres)
    if ( ind == testind ) {
      my.rf$importance<-( my.rf$importance/sum(my.rf$importance) )
      print( my.rf$importance )
      wmax<-which.max(my.rf$importance[1:ncol(mypopulp)] )
      for(wmax in 1:nv2 ) {
        print(paste("EigRegions",wmax))
        ww<-which( abs(eigpreds2[,wmax]) > 0 )
        print( colnames( brainpreds1 )[ww] )
        votedbrainregions2[ww]<-votedbrainregions2[ww]+1
    }
    }
  }

}
print( myresults )
collectresults[kk,]<-c(myerr1,myerr2,myp1,myp2)
if ( kk > 5 ) {
    wh1<-collectresults[,3] > 0.05 & collectresults[,4] > 0.05
    wh2<-collectresults[,3] < 0.05 & collectresults[,4] > 0.05
    wh3<-collectresults[,3] > 0.05 & collectresults[,4] < 0.05
    wh4<-collectresults[,3] < 0.05 & collectresults[,4] < 0.05
    plot(collectresults[wh1,1],collectresults[wh1,2],xlim=c(8,12), ylim=c(8,12))
    points(collectresults[wh2,1],collectresults[wh2,2],col='blue')
    points(collectresults[wh3,1],collectresults[wh3,2],col='green')
    points(collectresults[wh4,1],collectresults[wh4,2],col='red')
    print( t.test(collectresults[,1],collectresults[,2],paired=TRUE) )
    print(paste("ConsistencyMeasure",sum(votedbrainregions1==0)/sum(votedbrainregions1>1),sum(votedbrainregions1==0),"ConsistencyMeasure",sum(votedbrainregions2==0)/sum(votedbrainregions2>1),sum(votedbrainregions2==0)))
}
}
stop("TODO: write out figures and demog table from this doc")

