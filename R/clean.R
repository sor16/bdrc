#'Data cleaning
#'
#'Takes in stage-discharge data and cleans it, and subsets it according to user inputs.
#'@param file A string that contains the name of a txt file that contains stage and flow data from a certain river
#'@param advanced Logical,depending if you want to use the advanced settings
#'@param includedates A vector with two integers. The integers represent the year range the user wants to extract from the data.
#'Every datapoint that is not inside that date range will be discarded. Parameter advanced needs to be TRUE.
#'@param dummy A list with information on the dummy point, with elements W and Q, stage and discharge respectively. Parameter advanced needs to be TRUE.
#'@param keeprows A logical vector of the same length as the data which indicates whether to keep a datapoint or not. Parameter advanced needs to be TRUE.
#'@param force A list with information on the force point, with elements W and Q, stage and discharge respectively. Parameter advanced needs to be TRUE.
#'@param shiny Logical, whether the function should read the data as done in shiny or not. Only TRUE when used in a shiny interface.
#'@param Wmin Numeric, minimum stage of rating curve.
#'@param Wmax Numeric, maximum stage of rating curve.
#'@param exclude Logical depending on whether the user wants to exclude a date range from the data or not.
#'@param excludedates Vector with two Date values of the form %Y-%m-%d. The dates span the range which to exclude from the data. Parameters advanced and exclude need to be true.
#'@return Returns a list of objects that ar the input to model1BH and model1BH, especially the cleaned input data. Other outputs are a matrix of stage and discharge values
#' and the cleaned data before it changes due to the interactive elements such as dummypoint,forcepoint,keeprows etc.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}


clean <- function(file,advanced=TRUE,includedates=c(1950,as.numeric(format(Sys.Date(), "%Y"))),dummy=NULL,keeprows=NULL,force=NULL,shiny=FALSE,Wmin=NA,Wmax=NA, exclude=TRUE,excludedates=c(Sys.Date()-1,Sys.Date()-1)){
    suppressPackageStartupMessages(require(xlsx))

    if (is.null(file)){
        return(NULL)
    }
    name=file
    if(shiny==TRUE){
        list2env(file,envir=environment())

        if(type=='text/plain'){
            observedData=read.table(datapath,skip=2,sep="|",dec=",")
            observedData=observedData[,c(2,3,5,7,4)]
            names(observedData)=c("Date","Time","Quality","W","Q")
            observedData$Date=as.Date(gsub("\\.","-",observedData$Date),"%d-%m-%Y")
        }else if(type=="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"){
            observedData=read.xlsx(datapath,sheetIndex=1)
            names(observedData)=c("Date","Time","Quality","W","Q")

        }else{return(NULL)}

    }else{
        if(gsub(".*\\.", "",file)=='txt'){
        observedData=read.table(file,skip=2,sep="|",dec=",")
        observedData=observedData[,c(2,3,5,7,4)]
        #added
        names(observedData)=c("Date","Time","Quality","W","Q")
        observedData$Date=as.Date(gsub("\\.","-",observedData$Date),"%d-%m-%Y")
        }else if(gsub(".*\\.", "",file)=='xlsx'){
            observedData=read.xlsx(file,sheetIndex=1)
            #added
            names(observedData)=c("Date","Time","Quality","W","Q")
            observedData$Date=as.Date(gsub("\\.","-",observedData$Date),"%d-%m-%Y")
        }else{return(NULL)}
    }
        observedData$Time=as.character(observedData$Time)
        observedData$Q=gsub('\\s+', '',observedData$Q)
        observedData=observedData[observedData$W!=0,]
        observedData$Q=as.numeric(as.character(gsub(",",".",observedData$Q)))
        observedData$W=0.01*observedData$W
        observedData=observedData[with(observedData,order(W)),]
        observedData_before=observedData

        if(advanced==TRUE){
            if(length(keeprows)!=0){
                observedData=observedData[keeprows,]
            }
            years=as.numeric(format(observedData$Date, "%Y"))
            observedData=observedData[which(years<=includedates[2] & years >= includedates[1]),]

            if(exclude==TRUE){
                observedData=observedData[which(observedData$Date<=excludedates[1] | observedData$Date >= excludedates[2]),]
            }


            if(sum(unlist(lapply(dummy,length)))!=0){
                dummydata=as.data.frame(dummy)
                dummydata=round(dummydata,3)
                dummydata$Date=Sys.Date()
                dummydata$Time=format(Sys.time(),"%H:%M:%S")
                dummydata$Quality="dummy"
                dummydata=dummydata[,c("Date","Time","Quality","W","Q")]
                observedData=rbind(observedData,dummydata)

            }
            if(sum(unlist(lapply(force,length)))!=0){
                forcedata=as.data.frame(force)
                forcedata=round(forcedata,3)
                forcedata$Date=Sys.Date()
                forcedata$Time=format(Sys.time(),"%H:%M:%S")
                forcedata$Quality="forcepoint"
                forcedata=forcedata[,c("Date","Time","Quality","W","Q")]
                observedData=rbind(observedData,forcedata)
            }
        }
        #order again with new data from dummy or force
        observedData=observedData[with(observedData,order(W)),]
        if(is.na(Wmin)) Wmin=min(observedData$W)
        if(is.na(Wmax)) Wmax=max(observedData$W)
        observedData=subset(observedData,W >= Wmin & W <=Wmax )
        wq=as.matrix(observedData[,c("W","Q")])

    return(list("wq"=wq,"observedData"=observedData,"observedData_before"=observedData_before))
}
