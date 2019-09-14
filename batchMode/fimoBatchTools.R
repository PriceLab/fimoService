library(RUnit)
library(GenomicRanges)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_createFastaFileForFimo()
  test_runFimo()
  test_fixMotifNamesTruncatedAt100characters()
  test_expandFimoTable()
  test_fimoBatch()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
createFastaFileForFimo <- function(tbl.regions, fastaFileName)
{
   sequences <- with(tbl.regions, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end))

   if(is(sequences, "DNAString"))
       sequences <- DNAStringSet(sequences)
   names(sequences) <- with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end))
   message(sprintf("creating fasta file '%s' for %d sequences, for FIMO", fastaFileName, length(sequences)))

   writeXStringSet(sequences, fastaFileName)

} # createFastaFileForFimo
#------------------------------------------------------------------------------------------------------------------------
test_createFastaFileForFimo <- function()
{
   message(noquote(sprintf("--- test_createFastaFileForFimo")))

     # a < 1kb region in the promoter of GATA2, where TBX15 hits may be found
   tbl.regions <- data.frame(chrom="chr3", start=128497569, end=128498329, stringsAsFactors=FALSE)
   createFastaFileForFimo(tbl.regions, "smallTest.fa")
   checkTrue(file.exists("smallTest.fa"))
   x <- readDNAStringSet("smallTest.fa")
   checkEquals(length(x), 1)

   tbl.regions.2 <- data.frame(chrom=c("chr3", "chr4"),
                               start=c(128497569, 128498569),
                               end=  c(128498329, 128498589),
                               stringsAsFactors=FALSE)
   createFastaFileForFimo(tbl.regions.2, "smallTest2.fa")
   x <- readDNAStringSet("smallTest2.fa")
   checkEquals(length(x), 2)

} # test_createFastaFileForFimo
#------------------------------------------------------------------------------------------------------------------------
runFimo <- function(fastaFileName, resultsDirectory, threshold=5e-4)
{
   printf("--- running FIMO")
   FIMO <- file.path(Sys.getenv("HOME"), "meme", "bin", "fimo")  # true on both hagfish & khaleesi
   MOTIFS <- "~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme"
   #cmd <- sprintf("%s --oc %s --thresh -%f --text --verbosity 1 %s %s",
   #               FIMO, resultsDirectory, threshold, MOTIFS, fastaFileName)

   base.name <- basename(fastaFileName)
   tsv.name <- sub(".fa", ".tsv", base.name, fixed=TRUE)
   tsv.path <- file.path(resultsDirectory, tsv.name)

   cmd <- sprintf("%s --thresh %f --verbosity 1 --text %s %s > %s",
                  FIMO, threshold, MOTIFS, fastaFileName, tsv.path)

   print(cmd)
   system(cmd)
   return(tsv.path)
   # /users/pshannon/meme/bin/fimo --oc . --thresh 1e-4 ~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme chr11-small.fa

} # runFimo
#------------------------------------------------------------------------------------------------------------------------
test_runFimo <- function()
{
   message(noquote(sprintf("--- test_runFimo")))

   fastaFileName <- "smallTest.fa"  # created in test_createFastaFileforFimo
   checkTrue(file.exists(fastaFileName))
   resultsDirectory <- tempdir()
   resultsFile <- runFimo(fastaFileName, resultsDirectory, threshold=1e-2)
   checkTrue(file.exists(resultsFile))

} # test_runFimo
#------------------------------------------------------------------------------------------------------------------------
runFimoGATA2.big <- function()
{
      # GATA2 and the full span of all enhancers
   tbl.regions <- data.frame(chrom="chr3", start=128013674, end=128712841, stringsAsFactors=FALSE)
      # only 263kb around GATA2's enhancers
   tbl.regions <- data.frame(chrom="chr3", start=128383794, end=128647775, stringsAsFactors=FALSE)
     # 1kb region for faster testing
   start.loc <- 128383794
   end.loc <- start.loc + 1000
   tbl.regions <- data.frame(chrom="chr3", start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   with(tbl.regions, printf("span: %d", end-start))

   resultsDirectory <- "tmp-fimo-out-3"

   if(!dir.exists(resultsDirectory))
     dir.create (resultsDirectory)

   fastaFilename <- file.path(resultsDirectory, "test.fa")
   createFastaFileForFimo(tbl.regions, fastaFilename)
   checkTrue(file.exists(fastaFilename))
   checkTrue(file.size(fastaFilename) > with(tbl.regions, end-start))  # bigger than the base count


      #  4 minutes for 1e-4, 263k    201090 lines
      #  5 minutes for 1e-3, 263k   1510766 lnes
      # 18 minutes for 1e-2, 263k  12033436 lines

   system.time(runFimo(fastaFilename, resultsDirectory, threshold=1e-5))
   fimo.results.file <- file.path(resultsDirectory, "test.tsv")
   checkTrue(file.exists(fimo.results.file))

   tbl <- read.table(fimo.results.file, sep="\t", as.is=TRUE, nrow=-1, header=TRUE)  # two chopped names
   tbl.fixed <- fixMotifNamesTruncatedAt100characters(tbl)
   fimo.results.file.fixed <- file.path(resultsDirectory, "fimo-fixed.tsv")
   write.table(tbl.fixed, fimo.results.file.fixed, sep="\t", quote=FALSE, row.names=FALSE)

} # runFimoGATA2.big
#------------------------------------------------------------------------------------------------------------------------
fixMotifNamesTruncatedAt100characters <- function(tbl)
{
   char.100.lines <- which(nchar(tbl$motif_id) == 100)
   printf("found %d 100 character lines", length(char.100.lines))
   if(length(char.100.lines) == 0)
      invisible(tbl)

   motifDb.indices <- lapply(tbl$motif_id[char.100.lines], function(id) grep(id, names(MotifDb)))
   motifDb.names <- names(MotifDb)[as.integer(motifDb.indices)]
   tbl$motif_id[char.100.lines] <- motifDb.names

   invisible(tbl)

} # fixMotifNamesTruncatedAt100characters
#------------------------------------------------------------------------------------------------------------------------
test_fixMotifNamesTruncatedAt100characters <- function()
{
   printf("--- test_fixMotifNamesTruncatedAt100characters")
   #  a couple of instances of one long motif name:
   #  128 Mmusculus;Rnorvegicus;Xlaevis;Stropicalis;Ggallus;Hsapiens;Btaurus;Ocuniculus-jaspar2018-NFYA-MA0060           NA chr3:128383794-128647775 104174 104189      + 16.1461 1.29e-06      NA GCCAGCCAATCAACGC
   #  129 Mmusculus;Rnorvegicus;Xlaevis;Stropicalis;Ggallus;Hsapiens;Btaurus;Ocuniculus-jaspar2018-NFYA-MA0060           NA chr3:128383794-128647775 104210 104225      - 16.9888 3.79e-07      NA CGCAGCCAATGGGAGG

   start.loc <- 128383794 + 104170
   end.loc <- start.loc + 30
   tbl.regions <- data.frame(chrom="chr3", start=start.loc, end=end.loc, stringsAsFactors=FALSE)

   resultsDirectory <- "tmp-fimo-out-3"

   if(!dir.exists(resultsDirectory))
     dir.create (resultsDirectory)

   fastaFilename <- file.path(resultsDirectory, "test.fa")
   createFastaFileForFimo(tbl.regions, fastaFilename)
   checkTrue(file.exists(fastaFilename))
   checkTrue(file.size(fastaFilename) > with(tbl.regions, end-start))  # bigger than the base count


      #  4 minutes for 1e-4, 263k    201090 lines
      #  5 minutes for 1e-3, 263k   1510766 lnes
      # 18 minutes for 1e-2, 263k  12033436 lines

   system.time(runFimo(fastaFilename, resultsDirectory, threshold=1e-5))
   fimo.results.file <- file.path(resultsDirectory, "test.tsv")
   checkTrue(file.exists(fimo.results.file))

   tbl <- read.table(fimo.results.file, sep="\t", as.is=TRUE, nrow=-1, header=TRUE)  # two chopped names
   checkEquals(tbl$motif_id[1],
               "Mmusculus;Rnorvegicus;Xlaevis;Stropicalis;Ggallus;Hsapiens;Btaurus;Ocuniculus-jaspar2018-NFYA-MA0060")
   tbl.fixed <- fixMotifNamesTruncatedAt100characters(tbl)
   checkEquals(tbl.fixed$motif_id[1],
               "Mmusculus;Rnorvegicus;Xlaevis;Stropicalis;Ggallus;Hsapiens;Btaurus;Ocuniculus-jaspar2018-NFYA-MA0060.1")

} # test_fixMotifNamesTruncatedAt100characters
#------------------------------------------------------------------------------------------------------------------------
parseLocStrings <- function(locStrings)
{
   match <- regexpr("(?<chromosome>chr.*):(?<startPos>\\d+)-(?<endPos>\\d+)", locStrings, perl=TRUE)
   match.count <- length(match)
   if(match.count != length(locStrings)){
      stop("fimoBatchTools.parseLocStrings encountered unexpected sequence_name: %s", locStrings[1])
      }

   columns <-     attr(match, "capture.names")
   starts <-  as.data.frame(attr(match, "capture.start"))
   lengths <- as.data.frame(attr(match, "capture.length"))
   chroms <- substring(locStrings, starts$chromosome, starts$chromosome + lengths$chromosome -1)
   start.locs <- as.integer(substring(locStrings, starts$startPos, starts$startPos + lengths$startPos -1))
   end.locs <- as.integer(substring(locStrings, starts$endPos, starts$endPos + lengths$endPos -1))
   data.frame(chrom=chroms, start=start.locs, end=end.locs, stringsAsFactors=FALSE)

} # parseLocStrings
#------------------------------------------------------------------------------------------------------------------------
# add chrom coloumn.  adjust starts and stops.  rearrange columns
expandFimoTable <- function(tbl)
{
     # our convention is that fimo is called with a sequence_name which is the chromloc of the target sequence
     # here we depend upon that in order to add concrete chrom, start, end  locations to each match

   stopifnot("sequence_name" %in% colnames(tbl))  # this is

   tbl.locs <- parseLocStrings(tbl$sequence_name)
   number.of.distinct.chromosomes <- length(unique(tbl.locs$chrom))
   if(number.of.distinct.chromosomes > 1){
      stop("fimoBatchTools.expandFimoTable found %d chromosomes in table, only one permited",
           number.of.distinct.chromosomes)
      }
   tbl$chrom <- tbl.locs$chrom
   tbl$start <- tbl$start + tbl.locs$start
   tbl$end   <- tbl$stop + tbl.locs$start

   tfs <- unlist(lapply(tbl$motif_id, function(motifName) mcols(MotifDb[motifName])$geneSymbol))

   tbl$tf <- tfs
   coi <- c("chrom", "start", "end", "tf", "strand", "score", "p.value", "matched_sequence", "motif_id")
   tbl.final <- tbl[, coi]
   new.order <- order(tbl.final$start, decreasing=FALSE)
   tbl.final <- tbl.final[new.order,]

} # expandFimoTable
#------------------------------------------------------------------------------------------------------------------------
test_expandFimoTable <- function(tbl)
{
   printf("--- test_expandFimoTable")
   tbl.directFromFimo <- get(load("tbl.fimoBeforeColumnsFixes.RData"))
   tbl <- expandFimoTable(tbl.directFromFimo)

   checkEquals(nrow(tbl.directFromFimo), nrow(tbl))
   checkEquals(nrow(unique(tbl.directFromFimo)), nrow(tbl))
   checkTrue(all(tbl$chrom == "chr3"))
   checkEquals(sort(unique(tbl$tf)), c("CEBPZ", "FOXI1", "NFYA", "NFYB", "NFYC", "PBX3", "ZBTB14"))
   checkEquals(colnames(tbl), c("chrom", "start", "end", "tf", "strand", "score", "p.value", "matched_sequence", "motif_id"))

} # test_expandFimoTable
#------------------------------------------------------------------------------------------------------------------------
fimoBatch <- function(tbl.regions, matchThreshold)
{
   resultsDirectory <- tempdir()

   if(!dir.exists(resultsDirectory))
     dir.create (resultsDirectory)

   fastaFilename <- file.path(resultsDirectory, "forFimo.fa")
   createFastaFileForFimo(tbl.regions, fastaFilename)

   system.time(runFimo(fastaFilename, resultsDirectory, threshold=matchThreshold))
   fimo.results.file <- file.path(resultsDirectory, "forFimo.tsv")
   checkTrue(file.exists(fimo.results.file))

   tbl <- read.table(fimo.results.file, sep="\t", as.is=TRUE, nrow=-1, header=TRUE)  # two chopped names
   tbl.fixed <- fixMotifNamesTruncatedAt100characters(tbl)
   tbl.expanded <- expandFimoTable(tbl.fixed)
   invisible(tbl.expanded)

} # fimoBatch
#------------------------------------------------------------------------------------------------------------------------
test_fimoBatch <- function()
{
   printf("--- test_fimoBatch")
   chomosome <- "chr3"
   start.loc <- 128383794 + 104170
   start.loc.2 <- 28000000
   end.loc <- start.loc + 30
   end.loc.2 <- start.loc.2 + 1000
   tbl.regions <- data.frame(chrom="chr3", start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   tbl.match <- fimoBatch(tbl.regions, 1e-5)
   checkEquals(dim(tbl.match), c(13, 9))

     # now send two regions
   tbl.regions.2 <- data.frame(chrom=c("chr3", "chr3"),
                               start=c(start.loc, start.loc.2),
                               end=c(end.loc,     end.loc.2),
                               stringsAsFactors=FALSE)
   tbl.match.2 <- fimoBatch(tbl.regions.2, 1e-6)
   checkEquals(dim(tbl.match.2), c(7, 9))
     # make sure that both regions were used
   checkEquals(length(which(tbl.match.2$start < start.loc)), 4)
   checkEquals(length(which(tbl.match.2$start > start.loc)), 3)

} # test_fimoBatch
#------------------------------------------------------------------------------------------------------------------------
# if(!interactive()) test_runFimoBig()
