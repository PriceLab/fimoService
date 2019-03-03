printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.FimoClient <- setClass("FimoClientClass",
                        slots = c(host="character",
                                  port="integer",
                                  state="environment",
                                  quiet="logical")
                            )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("requestMatch", signature="obj", function(obj, sequences, pvalThreshold) standardGeneric("requestMatch"))
setGeneric("requestMatchForRegions", signature="obj", function(obj, tbl.regions, genomeName, pvalThreshold)
              standardGeneric("requestMatchForRegions"))
#------------------------------------------------------------------------------------------------------------------------
FimoClient <- function(host, port, quiet=FALSE)
{
    port <- as.integer(port)

    if(!quiet)
       printf("initializing zmq context")

    socketContext <- init.context()
    socket <- init.socket(socketContext, "ZMQ_REQ")
    uri <- sprintf("tcp://%s:%d", host, port)
    connect.socket(socket, uri)

    if(!quiet)
       printf("attempting socket connection to %s", uri)


    if(!quiet)
       printf("socket connected")

    state <- new.env(parent=emptyenv())
    state$socketContext <- socketContext
    state$socket <- socket

   .FimoClient(host=host, port=port, state=state, quiet=quiet)

} # FimoClient, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("requestMatch", "FimoClientClass",

    function(obj, sequences, pvalThreshold){
       stopifnot(is.list(sequences))
       xyz <- "Fimoclient::requestMatch"
       if(is.null(names(sequences))){
           artificial.names <- sprintf("seq%04d", 1:length(sequences))
           names(sequences) <- artificial.names
           }
       msg <- list(sequences=sequences, pvalThreshold=pvalThreshold)

       sequences.json.raw <- charToRaw(toJSON(msg, auto_unbox=TRUE))
       #sequences.json.raw <- charToRaw(toJSON(sequences, auto_unbox=TRUE))

       if(!obj@quiet)
           printf("sending %d sequences to server", length(sequences))

       send.socket(obj@state$socket, sequences.json.raw, serialize=FALSE)
       msg.raw <- receive.socket(obj@state$socket, unserialize=FALSE)
       s <- rawToChar(msg.raw)

       if(!obj@quiet){
          printf("response recieved from server")
          printf("nchar(s): %d", nchar(s))
          }

       tbl.as.list <- fromJSON(s)
       tbl.tmp <- as.data.frame(sapply(tbl.as.list, rbind))
       if(nrow(tbl.tmp) == 0 || ncol(tbl.tmp) == 0){
          printf("no motif matches for any sequence")
          return(data.frame())
          }
         # each cell in the dataframe is a list; each a scalar.  squash those
         # cells to scalars all the way down
       tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
       colnames(tbl)[1] <- "motif"
       tbl
       }) # request

#------------------------------------------------------------------------------------------------------------------------
setMethod("requestMatchForRegions", "FimoClientClass",

    function(obj, tbl.regions, genomeName, pvalThreshold){
       stopifnot(is.data.frame(tbl.regions))
       stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.regions)))
       stopifnot(genomeName %in% c("hg38"))
       require(BSgenome.Hsapiens.UCSC.hg38)
       sequences <- as.list(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                tbl.regions$chrom, tbl.regions$start, tbl.regions$end)))
       stopifnot(is.list(sequences))
       xyz <- "Fimoclient::requestMatch"
       if(is.null(names(sequences))){
           artificial.names <- sprintf("seq%04d", 1:length(sequences))
           names(sequences) <- artificial.names
           }
       msg <- list(sequences=sequences, pvalThreshold=pvalThreshold)
       sequences.json.raw <- charToRaw(toJSON(msg, auto_unbox=TRUE))
       #sequences.json.raw <- charToRaw(toJSON(sequences, auto_unbox=TRUE))

       if(!obj@quiet)
           printf("sending %d sequences to server", length(sequences))

       send.socket(obj@state$socket, sequences.json.raw, serialize=FALSE)
       msg.raw <- receive.socket(obj@state$socket, unserialize=FALSE)
       s <- rawToChar(msg.raw)

       if(!obj@quiet){
          printf("response recieved from server")
          printf("nchar(s): %d", nchar(s))
          }

       tbl.as.list <- fromJSON(s)
       tbl.tmp <- as.data.frame(sapply(tbl.as.list, rbind))
       if(nrow(tbl.tmp) == 0 || ncol(tbl.tmp) == 0){
          printf("no motif matches for any sequence")
          return(data.frame())
          }
         # each cell in the dataframe is a list; each a scalar.  squash those
         # cells to scalars all the way down
       tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
       colnames(tbl)[1] <- "motif"
       tbl
       }) # request

#------------------------------------------------------------------------------------------------------------------------
