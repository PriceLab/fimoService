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
# added for Fimo 5.0.3, in which the (apparently) the tsv table is in a different format so that, even after accomodating
# the comments at the end of the tsv file (in python, via pandas.read_csv(f, comment="#")) the old and never-quite-understood
# toJSON/unlist/do.call sequence no longer works.   here is more transparent approach.
.jsonToDataFrame <- function(jsonString)
{
  # browser()
  # xzy <- "entering .jsonToDataFrame, with jsonString"

  x2 <- fromJSON(jsonString)
  x3 <-lapply(x2, function(element) unlist(element, use.names=FALSE))

  tbl <- as.data.frame(do.call(cbind, x3), stringsAsFactors=FALSE)

  colnames(tbl)[grep("p-value", colnames(tbl))] <- "pValue"
  colnames(tbl)[grep("q-value", colnames(tbl))] <- "qValue"
  colnames(tbl)[grep("motif_id", colnames(tbl))] <- "motif"

  tbl$motif <- as.character(tbl$motif)
  tbl$strand <- as.character(tbl$strand)
  tbl$sequence_name <- as.character(tbl$sequence_name)
  tbl$matched_sequence <- as.character(tbl$matched_sequence)
  tbl$start <- as.numeric(tbl$start)
  tbl$stop <- as.numeric(tbl$stop)
  tbl$score <- as.numeric(tbl$score)
  tbl$pValue <- as.numeric(tbl$pValue)
  tbl$qValue <- as.numeric(tbl$qValue)

  return(tbl)

} # .jsonToDataFrame
#------------------------------------------------------------------------------------------------------------------------
setMethod("requestMatch", "FimoClientClass",

    function(obj, sequences, pvalThreshold){
       stopifnot(is.list(sequences))
       if(is.null(names(sequences))){
           artificial.names <- sprintf("seq%04d", 1:length(sequences))
           names(sequences) <- artificial.names
           }
       msg <- list(sequences=sequences, pvalThreshold=pvalThreshold)

       sequences.json.raw <- charToRaw(toJSON(msg, auto_unbox=TRUE))

       if(!obj@quiet)
           printf("sending %d sequences to server", length(sequences))

       send.socket(obj@state$socket, sequences.json.raw, serialize=FALSE)
       msg.raw <- receive.socket(obj@state$socket, unserialize=FALSE)
       character.response <- rawToChar(msg.raw)

       if(!obj@quiet){
          printf("response recieved from server")
          printf("nchar(character.response): %d", nchar(character.response))
          }

       if(nchar(character.response) == 0)
           return(data.frame())

       tbl.fimo <- .jsonToDataFrame(character.response)
       return(tbl.fimo)
       }) # request

#------------------------------------------------------------------------------------------------------------------------
setMethod("requestMatchForRegions", "FimoClientClass",

    function(obj, tbl.regions, genomeName, pvalThreshold){
       stopifnot(is.data.frame(tbl.regions))
       stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.regions)))

       stopifnot(genomeName %in% c("hg38", "tair10"))

       if(genomeName == "hg38"){
          require(BSgenome.Hsapiens.UCSC.hg38)
          sequences <- as.list(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                   tbl.regions$chrom, tbl.regions$start, tbl.regions$end)))
          names(sequences) <- sprintf("%s:%d-%d", tbl.regions$chrom, tbl.regions$start, tbl.regions$end)
          } # hg38
       if(genomeName == "tair10"){
          require(BSgenome.Athaliana.TAIR.TAIR9)
          chrom.name <- tbl.regions$chrom
             # tair9 (aka tair10?) need the chromsome ("sequence") names: Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC
          if(all(nchar(chrom.name) < 2))
             tbl.regions$chrom <- sprintf("chr%s", chrom.name)
          tbl.regions$chrom <- sub("chr", "Chr", tbl.regions$chrom)
          sequences <- as.list(as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9,
                                                   tbl.regions$chrom, tbl.regions$start, tbl.regions$end)))
          names(sequences) <- sprintf("%s:%d-%d", tbl.regions$chrom, tbl.regions$start, tbl.regions$end)
          } # tair10

       if(!obj@quiet)
           printf("requesting matches for %d sequences", length(sequences))
       tbl.fimo <- requestMatch(obj, sequences, pvalThreshold)
       tbl.locs <- .expandChromLocStrings(tbl.fimo$sequence_name)

       tbl.fimo$start <- tbl.fimo$start + tbl.locs$start
       tbl.fimo$end <- tbl.fimo$stop + tbl.locs$start
       tbl.fimo$chrom <- tbl.locs$chrom
       coi <- c("chrom", "start", "end", "motif", "strand", "score", "pValue", "qValue", "matched_sequence")
       tbl.out <- tbl.fimo[, coi]
       if(!obj@quiet)
           printf("   got %d hits", nrow(tbl.out))
       tbl.out
       }) # request

#------------------------------------------------------------------------------------------------------------------------
.parseChromLocString <- function(chromLocString)
{
    chromLocString <- gsub(",", "", chromLocString);
    tokens.0 <- strsplit(chromLocString, ":", fixed=TRUE)[[1]]
    stopifnot(length(tokens.0) == 2)
    chrom <- tokens.0[1]
    if(!grepl("chr", chrom))
        chrom <- sprintf("chr%s", chrom)

    tokens.1 <- strsplit(tokens.0[2], "-")[[1]]
    stopifnot(length(tokens.1) == 2)
    start <- as.integer(tokens.1[1])
    end <- as.integer(tokens.1[2])

    return(list(chrom=chrom, start=start, end=end))

} # .parseChromLocString
#------------------------------------------------------------------------------------------------------------------------
.expandChromLocStrings <- function(chromLocStrings)
{
   locs <- lapply(chromLocStrings, .parseChromLocString)
   tbl.locs <- as.data.frame(do.call(rbind, locs))

     #  all the new columns are lists().  throw away that structure:

   tbl.locs$chrom <- as.character(tbl.locs$chrom)
   tbl.locs$start <- as.numeric(tbl.locs$start)
   tbl.locs$end <- as.numeric(tbl.locs$end)

   tbl.locs

} # .expandChromLocStrings
#------------------------------------------------------------------------------------------------------------------------
