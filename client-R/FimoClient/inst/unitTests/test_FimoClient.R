library(FimoClient)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fc"))
    fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_.jsonToDataFrame()
     # depends on server restart:
     # (cd ~/github/fimoService/server; make -f makefile.pshannon unitTests)
     # to load the meme file
   #test_rreb1()

      # these next two tests do not tell me what meme file should
      # be loaded into the FimoServer.  thus disabled until I make
      # time to figure that out

   #test_request.small.100x()
   #test_request.large()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   checkEquals(is(fc), "FimoClientClass")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_request.small.100x <- function()
{
   printf("--- test_request.small.100x")
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)
   sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG")

   pval <- 0.1

   tbl <- requestMatch(fc, sequences, pvalThreshold=pval)
   dim(tbl)
   checkTrue("data.frame" %in% is(tbl))
   checkEquals(ncol(tbl), 9)
   checkTrue(nrow(tbl) > 30)

   pval <- 0.000001
   tbl <- requestMatch(fc, sequences, pvalThreshold=pval)
   checkEquals(dim(tbl), c(0,0))

} # test_request.small.100x
#------------------------------------------------------------------------------------------------------------------------
test_request.large <- function()
{
   printf("--- test_request.large")
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)

   count <- 1000
   sequences <- as.list(rep("CCCCTTCCGGG", count))
   names(sequences) <-  sprintf("tert_mut.%03d", 1:count)
   pvalThreshold <- 0.001

   max <- 3
   for(i in 1:max){
      printf("--- request %d", i)
      tbl <- requestMatch(fc, sequences, pvalThreshold)
      checkTrue("data.frame" %in% is(tbl))
      expected.row.count <- 4 * count
      checkEquals(dim(tbl), c(expected.row.count, 9))
      } # for i

} # test_request_large
#------------------------------------------------------------------------------------------------------------------------
test_rreb1 <- function()
{
   printf("--- test_MA0073.1")

   sequence <- list(test="CTTGGCCCCAGCACCCCCCGCCCCGAGGCCCGG")
     # FimoServer started with "../pfms/rreb1.human.meme",
     # created in ~/github/fimoService/pfms/prepateMatrixFile.R
     #   pfms <- query(MotifDb, c("sapiens", "RREB1"))
     #  export(pfms, "rreb1.human.meme", 'meme')

   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   tbl.fimo <- requestMatch(fc, sequence, pvalThreshold=0.0001)
   checkTrue(all(tbl.fimo$motif %in% c("Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D",
                                       "Hsapiens-SwissRegulon-RREB1.SwissRegulon",
                                       "Hsapiens-SwissRegulon-RREB1.SwissRegulon",
                                       "Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D")))
   checkTrue(all(tbl.fimo$strand == "-"))
   checkEquals(which(tbl.fimo$score > 10), c(1, 2))
   checkEquals(which(tbl.fimo$score < 10), c(3, 4))
     #                                      motif sequence.name start stop strand    score  p.value  q.value       matched.sequence
     # 1 Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D          test     5   26      - 13.44440 6.78e-06 0.000163 TCGGGGCGGGGGGTGCTGGGGC
     # 2 Hsapiens-SwissRegulon-RREB1.SwissRegulon          test     5   26      - 12.77780 8.34e-06 0.000200 TCGGGGCGGGGGGTGCTGGGGC
     # 3 Hsapiens-SwissRegulon-RREB1.SwissRegulon          test     3   24      -  7.34444 8.09e-05 0.000971 GGGGCGGGGGGTGCTGGGGCCA
     # 4 Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D          test     3   24      -  7.66667 8.50e-05 0.001020 GGGGCGGGGGGTGCTGGGGCCA

   # compare to trenas MotifMatcher
   library(trena)
   pfms=query(MotifDb, c("sapiens", "RREB1"), c("jaspar2018", "hocomoco", "swissregulon"))
   mm <- MotifMatcher("hg38", as.list(pfms), quiet=TRUE)
   tbl.regions <- data.frame(chrom="chr14", start=92513603, end=92513635, stringsAsFactors=FALSE)
   tbl.mm <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=70)

    #                                  motifName chrom motifStart motifEnd strand motifScore motifRelativeScore                  match chromStart chromEnd                  seq status          shortMotif
    # 7 Hsapiens-SwissRegulon-RREB1.SwissRegulon chr14   92513605 92513626      -   10.91554          0.7570944 TGGCCCCAGCACCCCCCGCCCC   92513603 92513635 CTTGGCCCCAGCACCCC...     wt  RREB1.SwissRegulon
    # 5 Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D chr14   92513605 92513626      -   10.84615          0.7601078 TGGCCCCAGCACCCCCCGCCCC   92513603 92513635 CTTGGCCCCAGCACCCC...     wt RREB1_HUMAN.H10MO.D
    # 1       Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92513606 92513625      +   10.72727          0.7564103   GGCCCCAGCACCCCCCGCCC   92513603 92513635 CTTGGCCCCAGCACCCC...     wt            MA0073.1
    # 3       Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92513608 92513627      +   10.72727          0.7564103   CCCCAGCACCCCCCGCCCCG   92513603 92513635 CTTGGCCCCAGCACCCC...     wt            MA0073.1
    # 8 Hsapiens-SwissRegulon-RREB1.SwissRegulon chr14   92513607 92513628      -   10.67795          0.7406150 GCCCCAGCACCCCCCGCCCCGA   92513603 92513635 CTTGGCCCCAGCACCCC...     wt  RREB1.SwissRegulon
    # 6 Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D chr14   92513607 92513628      -   10.61538          0.7439353 GCCCCAGCACCCCCCGCCCCGA   92513603 92513635 CTTGGCCCCAGCACCCC...     wt RREB1_HUMAN.H10MO.D
    # 4       Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92513609 92513628      +   10.09091          0.7115385   CCCAGCACCCCCCGCCCCGA   92513603 92513635 CTTGGCCCCAGCACCCC...     wt            MA0073.1
    # 2       Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92513607 92513626      +   10.00000          0.7051282   GCCCCAGCACCCCCCGCCCC   92513603 92513635 CTTGGCCCCAGCACCCC...     wt            MA0073.1

   checkTrue(all(tbl.fimo$motif %in% tbl.mm$motifName))

} # test_rreb1
#------------------------------------------------------------------------------------------------------------------------
test_rreb1_regions <- function()
{
   printf("--- test_rreb1_regions")
   tbl.regions <- data.frame(chrom="chr22", start=36253061, end=36253097, stringsAsFactors=FALSE)
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   tbl.fimo <- requestMatchForRegions(fc, tbl.regions, "hg38", pvalThreshold=0.000001)
   checkTrue(length(grep("IRF1", tbl.fimo$motif)) > 4)

} # test_rreb1_regions
#------------------------------------------------------------------------------------------------------------------------
# new explicit data.frame creation from json after upgrading - painfully! - to fimo 5.0.4
test_.jsonToDataFrame <- function()
{
   printf("--- test_.jsonToDataFrame")
    
   tbl <- requestMatch(fc, list(klf1="ATCGATCGAAGGGTGAGGCATCGATCG"), 10e-4)
   checkEquals(dim(tbl), c(1, 9))
   checkEquals(tbl$motif, "Hsapiens-SwissRegulon-KLF1.SwissRegulon")
   checkEquals(tbl$sequence_name, "klf1")
   checkEquals(tbl$matched_sequence, "AAGGGTGAGGC")
   checkTrue(tbl$score > 13)
   checkTrue(tbl$qValue < 0.0005)
               
   tbl <- requestMatch(fc, list(klf1 = "ATCGATCGAAGGGTGAGGCATCGATCG", fli1 = "ACTACAGGAAGTGGCAGCAGCAGCAGCAG"), pvalThreshold=10e-4)
   checkEquals(dim(tbl), c(2, 9))

   tbl <- requestMatch(fc, list(bogus=LETTERS), pvalThreshold=10e-4)
   checkEquals(dim(tbl), c(0, 9))
   
   #   produces this, after rawToChar(msg.raw)
   #  "{\"motif_id\":{\"0\":\"Hsapiens-SwissRegulon-KLF1.SwissRegulon\"},\"motif_alt_id\":{\"0\":null},\"sequence_name\":{\"0\":\"klf1\"},\"start\":{\"0\":9},\"stop\":{\"0\":19},\"strand\":{\"0\":\"+\"},\"score\":{\"0\":13.7528},\"p-value\":{\"0\":0.0000137},\"q-value\":{\"0\":0.000466},\"matched_sequence\":{\"0\":\"AAGGGTGAGGC\"}}"
   # requestMatch(fc, list(klf1 = "ATCGATCGAAGGGTGAGGCATCGATCG", fli1 = "ACTACAGGAAGTGGCAGCAGCAGCAGCAG"), produces this
   #  "{\"motif_id\":{\"0\":\"Hsapiens-jaspar2018-FLI1-MA0475.1\",\"1\":\"Hsapiens-SwissRegulon-KLF1.SwissRegulon\"},\"motif_alt_id\":{\"0\":null,\"1\":null},\"sequence_name\":{\"0\":\"fli1\",\"1\":\"klf1\"},\"start\":{\"0\":4,\"1\":9},\"stop\":{\"0\":14,\"1\":19},\"strand\":{\"0\":\"+\",\"1\":\"+\"},\"score\":{\"0\":17.6067,\"1\":13.7528},\"p-value\":{\"0\":0.000000205,\"1\":0.0000137},\"q-value\":{\"0\":0.0000147,\"1\":0.000987},\"matched_sequence\":{\"0\":\"ACAGGAAGTGG\",\"1\":\"AAGGGTGAGGC\"}}"

} # study_JSON.table.parsing.problem
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
