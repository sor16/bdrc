library(roxygen2)
    #' Clean data
    #'
    #' This function takes in the data and cleans it, if the advanced checkbox has bin checked it conditions and fixes the data depending on the input of the user ec.
clean <- function(file,advanced=FALSE,slider=0,dummy=NULL,keeprows=NULL,force=NULL,shiny=FALSE,Wmin=NA,Wmax=NA, experiod=FALSE,dates=0){
    if (is.null(file)){
        return(NULL)
    }
    if(shiny==TRUE){
        file=data.frame(file, stringsAsFactors = FALSE)
        list2env(file,envir=environment())
        qvdata=read.table(datapath,skip=2,sep="|",dec=",")
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

        if(experiod==TRUE){
            qvdata=qvdata[which(qvdata$Date<=dates[1] | qvdata$Date >= dates[2]),]
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
