library(MotifDb)
pfms <- query(MotifDb, c("sapiens", "jaspar2018", "MA0481.1"), c("foxp1"))
pfms <- query(MotifDb, c("athaliana", "jaspar2018"))
length(pfms)

# curated.names <- c("Hsapiens-jaspar2018-FLI1-MA0475.1", "Hsapiens-HOCOMOCOv10-KLF1_HUMAN.H10MO.C")
# pfms <- MotifDb[curated.names]
export(pfms, "athaliana--jaspar2018.meme", 'meme')
