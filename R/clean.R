library(roxygen2)
#'Cleans the input file
#'
#'This function takes in the data and cleans it, if advanced is TRUE and experiod is TRUE it conditions and fixes
#'the data depending on the input of the user.
#'@param file Is a string that contains the name of a txt file that contains stage and flow data from a certain river
#'@param advanced TRUE or FALSE depending if you want to use the advanced settings
#'@param slider Input should be a vector with two integers. The integers span the date constrainment for the data.
#'So every datapoint that is not inside that date range will be discarded.
#'@param dummy
#'@param keeprows
#'@param force
#'@param shiny
#'@param Wmin
#'@param Wmax
#'@param experiod TRUE or FALSE depending on whether you want to exclude a date range from the data
#'@param dates Input is a vector with two Date values. The dates span the range that you want to exclude from the data
#'@return If all the parameters are either default or used as described the output will be a matrix with the clean stage
#'and flow values as well as a data frame with the clean data values. That is cointaining date, time, quality, stage and flow.

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
        qvdata$Q=gsub('\\s+', '',qvdata$Q)
        qvdata=qvdata[qvdata$W!=0,]
        qvdata$Q=as.numeric(as.character(gsub(",",".",qvdata$Q)))
        qvdata$W=0.01*qvdata$W
        qvdata=qvdata[with(qvdata,order(W)),]
        qvdata_before=qvdata

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

    return(list("wq"=wq,"qvdata"=qvdata,"qvdata_before"=qvdata_before))
}
