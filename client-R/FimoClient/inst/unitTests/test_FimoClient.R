library(FimoClient)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 600159
if(!exists("fimoServerStarted")){
   fimoServerStarted <- TRUE
   meme.file <- system.file(package="FimoClient", "extdata", "human.jaspar2018.meme")
   stopifnot(file.exists(meme.file))
   cmd <- sprintf("make -f ~/github/fimoService/server/makefile PORT=%d MOTIFS=%s", FIMO_PORT, meme.file)
   print(cmd)
   system(cmd)
   printf("--- sleeping 5, making sure fimo server is awake")
   Sys.sleep(5)
   }

if(!exists("fc"))
    fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)


#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_.jsonToDataFrame()
   test_.expandChromLocStrings()
   test_matchTert()
   test_matchByRegion()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   checkEquals(is(fc), "FimoClientClass")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_matchTert <- function()
{
   printf("--- test_matchTert")

   sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG")

   pval <- 0.001

   tbl <- requestMatch(fc, sequences, pvalThreshold=pval)
   dim(tbl)
   checkTrue("data.frame" %in% is(tbl))
   checkEquals(ncol(tbl), 9)
   checkTrue(nrow(tbl) > 10)

   pval <- 0.00000001
   tbl <- requestMatch(fc, sequences, pvalThreshold=pval)
   checkEquals(nrow(tbl), 0)

} # test_matchTert
#------------------------------------------------------------------------------------------------------------------------
test_hg38.matchByRegion <- function()
{
   printf("--- test_hg38.matchByRegion")
        # part of proximal promoter of gata2
   tbl.regions <- data.frame(chrom="chr3",
                             start=c(128488353-1000, 128488353+1000),
                             end=c(128493411-1000, 128493411+1000),
                             stringsAsFactors=FALSE)
   tbl.regions <- data.frame(chrom="chr3", start=128486353, end=128493411, stringsAsFactors=FALSE)
   tbl <- requestMatchForRegions(fc, tbl.regions, "hg38", 1e-6)
   checkEquals(ncol(tbl), 9)
   checkEquals(colnames(tbl), c("chrom", "start", "end", "motif", "strand", "score", "pValue", "qValue", "matched_sequence"))
   checkTrue(nrow(tbl) > 50)
   checkTrue(all(tbl$start >= 128486353))
   checkTrue(all(tbl$end <= 128493411))
   checkTrue(all(tbl$chrom == "chr3"))

} # test_hg38.matchByRegion
#------------------------------------------------------------------------------------------------------------------------
test_tair10.matchByRegion <- function()
{
   printf("--- test_tair10.matchByRegion")
        # part of proximal promiter of WBC19, AT3G55130
   fc <- FimoClient("localhost", 60001)
   tbl.regions <- data.frame(chrom="3", start=20434092, end=20438604, stringsAsFactors=FALSE)
   tbl <- requestMatchForRegions(fc, tbl.regions, "tair10", 1e-6)
   checkEquals(ncol(tbl), 9)
   checkEquals(colnames(tbl), c("chrom", "start", "end", "motif", "strand", "score", "pValue", "qValue", "matched_sequence"))
   checkTrue(nrow(tbl) > 50)
   checkTrue(all(tbl$start >= 128486353))
   checkTrue(all(tbl$end <= 128493411))
   checkTrue(all(tbl$chrom == "chr3"))

} # test_tair10.matchByRegion
#------------------------------------------------------------------------------------------------------------------------
# new explicit data.frame creation from json after upgrading - painfully! - to fimo 5.0.4
test_.jsonToDataFrame <- function()
{
   printf("--- test_.jsonToDataFrame")

   tbl <- requestMatch(fc, list(klf1="ATCGATCGAAGGGTGAGGCATCGATCG"), 10e-4)
   dim(tbl)
   checkEquals(ncol(tbl), 9)
   checkTrue(nrow(tbl) > 30)
   checkTrue("Hsapiens-jaspar2018-KLF5-MA0599.1" %in% tbl$motif)
   checkTrue(all(tbl$sequence_name == "klf1"))
   checkTrue("GCCTCACCCT" %in% tbl$matched_sequence)
   checkTrue(tbl$score[1] > 10)
   checkTrue(tbl$qValue[1] < 0.0025)

   tbl <- requestMatch(fc, list(klf1 = "ATCGATCGAAGGGTGAGGCATCGATCG", fli1 = "ACTACAGGAAGTGGCAGCAGCAGCAGCAG"), pvalThreshold=1e-4)
   checkTrue(nrow(tbl) > 20)

   tbl <- requestMatch(fc, list(bogus=LETTERS), pvalThreshold=10e-4)
   checkEquals(dim(tbl), c(0, 9))

   #   produces this, after rawToChar(msg.raw)
   #  "{\"motif_id\":{\"0\":\"Hsapiens-SwissRegulon-KLF1.SwissRegulon\"},\"motif_alt_id\":{\"0\":null},\"sequence_name\":{\"0\":\"klf1\"},\"start\":{\"0\":9},\"stop\":{\"0\":19},\"strand\":{\"0\":\"+\"},\"score\":{\"0\":13.7528},\"p-value\":{\"0\":0.0000137},\"q-value\":{\"0\":0.000466},\"matched_sequence\":{\"0\":\"AAGGGTGAGGC\"}}"
   # requestMatch(fc, list(klf1 = "ATCGATCGAAGGGTGAGGCATCGATCG", fli1 = "ACTACAGGAAGTGGCAGCAGCAGCAGCAG"), produces this
   #  "{\"motif_id\":{\"0\":\"Hsapiens-jaspar2018-FLI1-MA0475.1\",\"1\":\"Hsapiens-SwissRegulon-KLF1.SwissRegulon\"},\"motif_alt_id\":{\"0\":null,\"1\":null},\"sequence_name\":{\"0\":\"fli1\",\"1\":\"klf1\"},\"start\":{\"0\":4,\"1\":9},\"stop\":{\"0\":14,\"1\":19},\"strand\":{\"0\":\"+\",\"1\":\"+\"},\"score\":{\"0\":17.6067,\"1\":13.7528},\"p-value\":{\"0\":0.000000205,\"1\":0.0000137},\"q-value\":{\"0\":0.0000147,\"1\":0.000987},\"matched_sequence\":{\"0\":\"ACAGGAAGTGG\",\"1\":\"AAGGGTGAGGC\"}}"

} # study_JSON.table.parsing.problem
#------------------------------------------------------------------------------------------------------------------------
test_.expandChromLocStrings <- function()
{
   printf("--- test_.expandChromLocStrings")
   chromLocs <- c("chr3:128489353-128494411", "chr3:128487353-128492411")
   tbl.locs <- FimoClient:::.expandChromLocStrings(chromLocs)
   checkEquals(dim(tbl.locs), c(2, 3))
   checkTrue(is.character(tbl.locs$chrom))
   checkTrue(all(tbl.locs$chrom == "chr3"))
   checkEquals(tbl.locs$start, c(128489353, 128487353))
   checkEquals(tbl.locs$end, c(128494411, 128492411))

} # test_.expandChromLocStrings
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
