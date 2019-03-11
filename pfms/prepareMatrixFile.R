library(MotifDb)
pfms <- query(MotifDb, c("sapiens"), c("FLI1", "KLF1"))
length(pfms)

curated.names <- c("Hsapiens-jaspar2018-FLI1-MA0475.1", "Hsapiens-HOCOMOCOv10-KLF1_HUMAN.H10MO.C")
pfms <- MotifDb[curated.names]
export(pfms, "human-klf1-fli1.meme", 'meme')
