printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.FimoClient <- setClass("FimoClientClass",
                        slots = c(host="character",
                                  port="integer",
                                  state="environment",
                                  #socketContext="R6",
                                  #socket="R6",
                                  #socketContext="externalptr",
                                  #socket="externalptr",
                                  quiet="logical")
                            )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("requestMatch", signature="obj", function(obj, sequences) standardGeneric("requestMatch"))
#------------------------------------------------------------------------------------------------------------------------
FimoClient <- function(host, port, quiet=FALSE)
{
    port <- as.integer(port)

    if(!quiet)
       printf("initializing zmq context")

    socketContext = zmq$Context()
    socket = socketContext$socket("ZMQ_REQ")
    uri <- sprintf("tcp://%s:%d", host, port)

    if(!quiet)
       printf("attempting socket connection to %s", uri)

    socket$connect(uri)

    if(!quiet)
       printf("socket connected")

    #context = init.context()
    #socket = init.socket(context,"ZMQ_REQ")
    #uri <- sprintf("tcp://%s:%d", host, port)
    #socketContext = init.context()
    #result <- connect.socket(socket, uri)

    state <- new.env(parent=emptyenv())
    state$socketContext <- socketContext
    state$socket <- socket
    
   .FimoClient(host=host, port=port, state=state, quiet=quiet)

} # FimoClient, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("requestMatch", "FimoClientClass",

    function(obj, sequences){
       stopifnot(is.list(sequences))
       if(is.null(names(sequences))){
           artificial.names <- sprintf("seq%04d", 1:length(sequences))
           names(sequences) <- artificial.names
           }
       sequences.json.raw <- charToRaw(toJSON(sequences, auto_unbox=TRUE))
       if(!obj@quiet)
           printf("sending %d sequences to server", length(sequences))

       obj@state$socket$send(sequences.json.raw, serialize=FALSE)
       msg.raw <- obj@state$socket$receive(unserialize=FALSE)
       #browser()
       s <- rawToChar(msg.raw)
       #send.socket(obj@state$socket, sequences.json, serialize=FALSE)
       #response <- receive.socket(obj@state$socket, unserialize=FALSE)
       #s <- receive.string(obj@state$socket)
       #printf("--- socket begore receive")
       #print(obj@state$socket)
       #s <- rawToChar(obj@state$socket$receive(unserialize=FALSE))

       if(!obj@quiet)
           printf("response recieved from server")
       #printf("nchar(s): %d, %s", nchar(s), s)
       tbl.as.list <- fromJSON(s)
       tbl.tmp <- as.data.frame(sapply(tbl.as.list, rbind))
         # each cell in the dataframe is a list; each a scalar.  squash those 
         # cells to scalars all the way down
       tbl <- as.data.frame(lapply(tbl.tmp, unlist), stringsAsFactors=FALSE)
       tbl
       }) # request

#------------------------------------------------------------------------------------------------------------------------
