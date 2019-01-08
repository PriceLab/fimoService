library(MotifDb)
pfms <- query(MotifDb, "sapiens", c("jaspar2018"))
length(pfms)
export(pfms, "human-jaspar2018.meme", 'meme')
