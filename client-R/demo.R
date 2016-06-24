library(rzmq);      stopifnot(packageVersion("rzmq") >= "0.7.7")
library(jsonlite);  stopifnot(packageVersion("jsonlite") >= "0.9.22")

context = init.context()
socket = init.socket(context,"ZMQ_REQ")
connect.socket(socket,"tcp://localhost:5558")
i <- 0

sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG")
sequences.json <- charToRaw(toJSON(sequences, auto_unbox=TRUE))
#sequences.json <- as.character(toJSON(sequences, auto_unbox=TRUE))
#sequences.json <- toJSON(sequences, auto_unbox=TRUE)

for(i in 1:1000){
   send.socket(socket, sequences.json, serialize=FALSE)
   s <- receive.string(socket)
   print(s)
   #tbl.as.list <- fromJSON(s)
   #tbl.tmp <- as.data.frame(sapply(tbl.as.list, rbind))
      # each cell in the dataframe is a list; each a scalar.  squash those 
      # cells to scalars all the way down
   #tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
   #printf("response %d", i)
   #print(tbl)
   }
