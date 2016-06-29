library(FimoClient)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "whovian"
FIMO_PORT <- 5558
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   #test_constructor()
   test_request.small.100x()
   test_request.large()

} # runTests    
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)    

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_request.small.100x <- function()
{
   printf("--- test_request.small.100x")
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)
   sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut= "CCCCTTCCGGG")

   max <- 3
   for(i in 1:max){
      printf("--- request %d", i)
      tbl <- requestMatch(fc, sequences)
      checkTrue("data.frame" %in% is(tbl))
      checkEquals(dim(tbl), c(4, 9))
      } # for i

} # test_request.small.100x
#------------------------------------------------------------------------------------------------------------------------
test_request.large <- function()
{
   printf("--- test_request.large")
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)

   count <- 1000
   sequences <- as.list(rep("CCCCTTCCGGG", count))
   names(sequences) <-  sprintf("tert_mut.%03d", 1:count)

   max <- 3
   for(i in 1:max){
      printf("--- request %d", i)
      tbl <- requestMatch(fc, sequences)
      checkTrue("data.frame" %in% is(tbl))
      expected.row.count <- 4 * count
      checkEquals(dim(tbl), c(expected.row.count, 9))
      } # for i

} # test_request_large
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
