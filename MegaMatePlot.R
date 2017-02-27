#################### FUNCTIONS for plotting megamate allele frequency data ################

cumPosition<-function(myData,genomeIdx){
  cumPos=with(myData, data.frame("Chr"=Chr,"Pos"=Pos,"Position"=rep(0,dim(myData)[1])))
  levels(cumPos$Chr)<-c(levels(cumPos$Chr),levels(genomeIdx[,1])[5])
  for (i in 1:dim(genomeIdx)[1]){
    chr<-genomeIdx[i,1]
    cumPos[cumPos$Chr==chr,3]<-cumPos$Pos[cumPos$Chr==chr]+genomeIdx[genomeIdx$V1==chr,3]
  }
  myData<-cbind(myData[,1:3],cumPos[3],myData[,4:dim(myData)[2]])
  return(myData)
}



cleanData<-function(myData,dataName="",MADs=5) {
  #write to log file some general data
  myMessage<-paste0("All loci: ",dim(myData)[1],"\n")
  readDepth<-grep("_readDepth",names(myData))
  med<-median(myData[,readDepth])
  myMessage<-paste0(myMessage,"Median read depth: ",med," \n")
  mad<-mad(myData[,readDepth])
  myMessage<-paste0(myMessage,"Median average deviation of read depth: ",mad," \n")
  #remove sequences with readDepth<2
  remove<-myData[,readDepth]<2
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because readDepth<2\n")
  #remove sequences that are extreme outliers in the distribution
  remove<-myData[,readDepth]>med+MADs*mad
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because more than ", med+MADs*mad," reads (",MADs,"xMADs)\n")
  remove<-myData[,readDepth]<med-MADs*mad
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because less than ", med-MADs*mad," reads (",MADs,"xMADs)\n")
  #write some final data to logFile
  myMessage<-paste0(myMessage,"Final number of loci: ",dim(myData)[1],"\n")
  logFile<-file(description=paste0("../finalData/log_",dataName),open="a")
  writeLines(myMessage,logFile)
  close(logFile)
  return(myData)
}



removeNonvariable<-function(myData) {
  #remove sequences that are not polymorphic
  CBfreq<-grep("_CBfreq",names(myData))
  remove<-(rowMeans(myData[,CBfreq]==1)==1 | rowMeans(myData[,CBfreq]==0)==1)
  myData<-myData[!remove,]
  return(myData)
}



plotFreq <- function (myData, extraData, useLoci=c(2)) {
  CBfreq<-grep("_CBfreq",names(myData))
  plot(myData$Position,myData[,CBfreq[1]], type='p',xlab="Position(bp)",ylab="Hawaii allele frequency",
       main=extraData$mainTitle,ylim=c(0,1),xlim=c(0,max(myData$Position)),pch=16,cex=0.4,col="#00000077")
  points(myData$Position,myData[,CBfreq[2]],type='p',pch=16,cex=0.4,col="#FF000077")
  text(extraData$midpoint,0.98,extraData$chrNum,cex=1, col="black")
  abline(v=genomeIdx[,3],col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="dark red",lwd=0.6)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=1, col="red",adj=extraData$labelAdj[useLoci])
  legend("bottomright",legend=c("1F3","mel-26"),col=c(1,2),pch=16,cex=0.9)
}



plotSmoothedDiff <- function (myData, extraData, sm=1001, useLoci=c(2)) {
  CBfreq<-grep("_CBfreq",names(myData))
  CBdiff<-myData[,CBfreq[2]]-myData[,CBfreq[1]]
  smdiff<-smootheByChr(CBdiff,myData$Chr,sm)
  plot(myData$Position,smdiff, type="n",xlab="Position(bp)", ylab="Hawaii allele freq difference", main=extraData$mainTitle, 
       ylim=c(-1,1), xlim=c(0,max(myData$Position)))
  plotByChr(myData, smdiff, chrList=extraData$chrs,lwd=4)
  text(extraData$midpoint,0.98,extraData$chrNum,cex=1, col="black")
  abline(v=genomeIdx[,3],col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=1, col="red",adj=extraData$labelAdj[useLoci])
  abline(h=0,lty=5,col="dark gray")
}



smootheByChr <- function( data2smoothe, chrData, smWin) {
  smoothedData<-c(sgolayfilt(data2smoothe[chrData=="CHROMOSOME_I"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_II"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_III"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_IV"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_V"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_X"],p=3,n=smWin))
  return(smoothedData)
}



plotByChr<-function(myData, yData, chrList, chrColumn="Chr", ...) {
  colours<-brewer.pal(6,"Dark2")
  for (c in chrList) {
    i=which(myData[,chrColumn]==c)
    lines(myData[i,"Position"], yData[i], col=colours[match(c,chrList)],...)
  }
}



plotSmPvals<-function(myData, log10pVals, extraData, sm=1001, useLoci=c(2)) {
  smPvals<-smootheByChr(log10pVals, myData$Chr, sm)
  plot(myData$Position ,smPvals, type="n",xlab="Position(bp)", ylab="-log10(p value)", main=extraData$mainTitle, 
       ylim=c(0,max(smPvals)), xlim=c(0,max(myData$Position)))
  plotByChr(myData, smPvals, chrList=extraData$chrs, lwd=4)
  text(extraData$midpoint, 0.98*max(smPvals), extraData$chrNum, cex=1, col="black")
  abline(v=genomeIdx[,3],col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=1, col="red",adj=extraData$labelAdj[useLoci])
  abline(h=-log10(0.05),lty=5,col="dark red")
  text(90000000,-log10(0.05),"FDR=0.05", col="dark red",pos=3)
}
