clean <- function(file,advanced=FALSE,slider=0,dummy=NULL,keeprows=NULL,force=NULL,shiny=FALSE,Wmin=NA,Wmax=NA,experiod=0,dates=FALSE){
    if (is.null(file)){
        return(NULL)
    }
    if(shiny==TRUE){
        file=data.frame(file, stringsAsFactors = FALSE)
        list2env(file,envir=environment())
        qvdata=read.table(datapath,skip=3,sep="|",dec=",")
    }else{
        qvdata=read.table(file,skip=2,sep="|",dec=",")
    }
        qvdata=qvdata[,c(2,3,5,7,4)]
        names(qvdata)=c("Date","Time","Quality","W","Q")
        qvdata$Time=as.character(qvdata$Time)
        qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
        qvdata$Quality=gsub('\\s+', '',qvdata$Quality)
        qvdata$W=0.01*qvdata$W
        qvdata=qvdata[with(qvdata,order(W)),]
        if(length(keeprows)!=0){
            qvdata=qvdata[keeprows,]
        }
        if(advanced==TRUE){
            years=as.numeric(format(qvdata$Date, "%Y"))
            qvdata=qvdata[which(years<=slider[2] & years >= slider[1]),]

        }
        if(dates==TRUE){
            datesA=qvdata$Date
            qvdata=qvdata[which(datesA<=experiod[2] & datesA >= experiod[1]),]


        }
        if(sum(unlist(lapply(dummy,length)))!=0){
            dummydata=as.data.frame(dummy)
            dummydata=round(dummydata,3)
            dummydata$Date=Sys.Date()
            dummydata$Time=format(Sys.time(),"%H:%M:%S")
            dummydata$Quality="dummy"
            dummydata=dummydata[,c("Date","Time","Quality","W","Q")]
            qvdata=rbind(qvdata,dummydata)

        }
        if(sum(unlist(lapply(force,length)))!=0){
            forcedata=as.data.frame(force)
            forcedata=round(forcedata,3)
            forcedata$Date=Sys.Date()
            forcedata$Time=format(Sys.time(),"%H:%M:%S")
            forcedata$Quality="forcepoint"
            forcedata=forcedata[,c("Date","Time","Quality","W","Q")]
            qvdata=rbind(qvdata,forcedata)
        }
        #order again with new data from dummy or force
        qvdata=qvdata[with(qvdata,order(W)),]
        if(is.na(Wmin)) Wmin=min(qvdata$W)
        if(is.na(Wmax)) Wmax=max(qvdata$W)
        qvdata=subset(qvdata,W >= Wmin & W <=Wmax )
        wq=as.matrix(qvdata[,c("W","Q")])

    return(list("wq"=wq,"qvdata"=qvdata))
}
