library(MotifDb)

pfms <- query(MotifDb, c("sapiens", "RREB1"))
length(pfms)
export(pfms, "rreb1.human.meme", 'meme')
