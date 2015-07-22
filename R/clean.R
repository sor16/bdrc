clean <- function(file,advanced=FALSE,slider=0,dummy=NULL,keeprows=NULL){
    if (is.null(file)){
        return(NULL)
    }
    file=data.frame(file, stringsAsFactors = FALSE)
    list2env(file,envir=environment())
    if (type =="text/plain"){
        qvdata=read.table(datapath,skip=3,sep="|",dec=",")
        qvdata=qvdata[,c(2,3,5,7,4)]
        names(qvdata)=c("Date","Time","Quality","W","Q")
        qvdata$Time=as.character(qvdata$Time)
        qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
        qvdata$Quality=gsub('\\s+', '',qvdata$Quality)
        qvdata$W=0.01*qvdata$W
        qvdata=qvdata[with(qvdata,order(W)),]
        if(!is.null(length(keeprows))){
            qvdata=qvdata[keeprows,]
        }
        if(advanced==TRUE){
            years=as.numeric(format(qvdata$Date, "%Y"))
            qvdata=qvdata[which(years<=slider[2] & years >= slider[1]),]

        }
        if(!is.null(length(dummy))){
            dummydata=as.data.frame(dummy)
            dummydata=round(dummydata,3)
            dummydata$Date=Sys.Date()
            dummydata$Time=format(Sys.time(),"%H:%M:%S")
            dummydata$Quality="dummy"
            dummydata=dummydata[,c("Date","Time","Quality","W","Q")]
            qvdata=rbind(qvdata,dummydata)

        }
        #order again with new data from dummy or force
        qvdata=qvdata[with(qvdata,order(W)),]
        wq=as.matrix(qvdata[,c("W","Q")])
    }
    return(list("wq"=wq,"qvdata"=qvdata))
}
