library(roxygen2)
#'Cleans the input file
#'
#'This function takes in the data and cleans it, if advanced is TRUE and exclude is TRUE it conditions and fixes
#'the data depending on the input of the user.
#'@param file Is a string that contains the name of a txt file that contains stage and flow data from a certain river
#'@param advanced Logical,depending if you want to use the advanced settings
#'@param includedates A vector with two integers. The integers represent the year range the user wants to extract from the data. Parameter advanced needs to be true.
#'So every datapoint that is not inside that date range will be discarded.
#'@param dummy A list with information on the dummy point, with elements W and Q, stage and discharge respectively. Parameter advanced needs to be true.
#'@param keeprows A logical vector which indicates whether to keep a datapoint or not. Parameter advanced needs to be true.
#'@param force A list with information on the force point, with elements W and Q, stage and discharge respectively. Parameter advanced needs to be true.
#'@param shiny Logical, whether the function should read the data as done in shiny or not.
#'@param Wmin Numeric, minimum stage of rating curve.
#'@param Wmax Numeric, maximum stage of rating curve.
#'@param exclude Logical depending on whether the user wants to exclude a date range from the data or not.
#'@param excludedates Vector with two Date values of the form %Y-%m-%d. The dates span the range which to exclude from the data. Parameters advanced and exclude need to be true.
#'@return If all the parameters are either default or used as described the output will be a list containing a matrix with the cleaned stage
#'and flow values as well as a data frame with the cleaned data cointaining date, time, quality, stage and flow.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}


clean <- function(file,advanced=TRUE,includedates=c(1950,as.numeric(format(Sys.Date(), "%Y"))),dummy=NULL,keeprows=NULL,force=NULL,shiny=FALSE,Wmin=NA,Wmax=NA, exclude=TRUE,excludedates=c(Sys.Date()-1,Sys.Date()-1)){
    require(xlsx)

    if (is.null(file)){
        return(NULL)
    }
    name=file
    if(shiny==TRUE){
        list2env(file,envir=environment())
        #file=data.frame(file, stringsAsFactors = FALSE)

        if(type=='text/plain'){
            qvdata=read.table(datapath,skip=2,sep="|",dec=",")
            qvdata=qvdata[,c(2,3,5,7,4)]

        }else if(type=="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"){
            qvdata=read.xlsx(datapath,sheetIndex=1)

        }else{return(NULL)}

    }else{
        if(gsub(".*\\.", "",file)=='txt'){
        qvdata=read.table(file,skip=2,sep="|",dec=",")
        qvdata=qvdata[,c(2,3,5,7,4)]
        }else if(gsub(".*\\.", "",file)=='xlsx'){
            qvdata=read.xlsx(file,sheetIndex=1)
        }else{return(NULL)}
    }
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
            qvdata=qvdata[which(years<=includedates[2] & years >= includedates[1]),]

            if(exclude==TRUE){
                qvdata=qvdata[which(qvdata$Date<=excludedates[1] | qvdata$Date >= excludedates[2]),]
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
        }
        #order again with new data from dummy or force
        qvdata=qvdata[with(qvdata,order(W)),]
        if(is.na(Wmin)) Wmin=min(qvdata$W)
        if(is.na(Wmax)) Wmax=max(qvdata$W)
        qvdata=subset(qvdata,W >= Wmin & W <=Wmax )
        wq=as.matrix(qvdata[,c("W","Q")])

    return(list("wq"=wq,"qvdata"=qvdata,"qvdata_before"=qvdata_before))
}
