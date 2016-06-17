library(rzmq);      stopifnot(packageVersion("rzmq") >= "0.7.7")
library(jsonlite);  stopifnot(packageVersion("jsonlite") >= "0.9.22")

context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://localhost:5556")
i <- 0

sequences <- list(sequences=list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG"))
sequences.json <- charToRaw(toJSON(sequences, auto_unbox=TRUE))

for(i in 1:4){
   send.socket(socket, sequences.json, serialize=FALSE)
   tbl.json <- fromJSON(rawToChar(receive.socket(socket, unserialize=FALSE)))
   tbl.tmp <- as.data.frame(sapply(tbl.json, rbind))
   tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
   print(tbl)
   }
