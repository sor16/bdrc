clean <- function(file,advanced=FALSE,slider=0,dummy=NULL,toggle=NULL){
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
        if(advanced==TRUE){
            years=as.numeric(format(qvdata$Date, "%Y"))
            qvdata=qvdata[which(years<=slider[2] & years >= slider[1]),]

        }
        if(!is.null(dummy)){
            dummydata=as.data.frame(dummy)
            dummydata$Date=Sys.Date()
            dummydata$Time=format(Sys.time(),"%H:%M:%S")
            dummydata$Quality="dummy"
            dummydata=dummydata[,c(3,4,5,1,2)]
            qvdata=rbind(qvdata,dummydata)

        }
        if(!is.null(toggle)){
            toggle=as.data.frame(toggle)
            qvdata=qvdata[-which(rowSums(toggle)==rowSums(qvdata[,c("W","Q")])),]
        }
        qvdata=qvdata[with(qvdata,order(W)),]
        wq=as.matrix(qvdata[,c("W","Q")])
    }
    return(list("wq"=wq,"qvdata"=qvdata))
}
