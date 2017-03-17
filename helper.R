dir = "Rat_data"
pData = "sample_info.txt"
annot = ragene21sttranscriptcluster.db::ragene21sttranscriptcluster.db
hom <- read.delim("homologene.data", header = F,
                  colClasses = c("character", "character", "character", "NULL", "NULL", "NULL"),
                  col.names = c("id", "tax", "entrez", "NULL", "NULL", "NULL"))
taxid = "10116"