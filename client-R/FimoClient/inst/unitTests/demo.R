library(pbdZMQ)
library(jsonlite);  stopifnot(packageVersion("jsonlite") >= "0.9.22")

context = init.context()
socket = init.socket(context,"ZMQ_REQ")
host <- "localhost"
port <- 5000
server.url <- sprintf("tcp://%s:%d", host, port)
connect.socket(socket, server.url)
i <- 0

sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG")
pvalThreshold <- 0.0001
msg <- list(sequences=sequences, pvalThreshold=pvalThreshold)

for(i in 1:2){
   printf("---- sending request to %s", server.url)
   print(msg)
   sequences.json <- charToRaw(toJSON(msg, auto_unbox=TRUE))
   send.socket(socket, sequences.json, serialize=FALSE)
   msg.raw <- receive.socket(socket, unserialize=FALSE)
   s <- rawToChar(msg.raw)
   tbl.as.list <- fromJSON(s)
   tbl.tmp <- as.data.frame(sapply(tbl.as.list, rbind))
   tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
   colnames(tbl)[1] <- "motif"
   print(tbl)
   }
