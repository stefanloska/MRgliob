library(oligo)
library(Biobase)
library(AnnotationDbi)

# Get expression ####

# Read data ####

# get sample meta data
tab <- read.delim("Rat_data/sample_info.txt", stringsAsFactors = F)
rownames(tab) <- tab$filename
tab

# get microarray signal
library(oligo)
ab <- read.celfiles(filenames = list.celfiles("Rat_data", full.names = T), phenoData = AnnotatedDataFrame(tab))


# Preprocessing ####

# rma - two levels to compare
R <- lapply(c("probeset", "core"), function(target){
  oligo::rma(ab, target = target)
})

# look into id whne the same when differ
p <- as.numeric(rownames(R[[1]]))
c <- as.numeric(rownames(R[[2]]))

k <- 100
p_ <- p[p %in% 17610200:17610400]
c_ <- c[c %in% 17610200:17610400]
plot(c(min(c(p_, c_)), max(c(p_, c_))), c(0, 2), type = "n")
points(p_, rep(1, length(p_)), pch=15, col = 2)
points(c_, rep(1.05, length(c_)), pch=15, col = 3)

# choose core
R <- oligo::rma(ab, target = "core")

# preview
pData(R)
fData(R)
exprs(R)[1:5,]
annotation(R)
experimentData(R)


# Get annotations ####
# 2 different ways tested

# get annotation first way
featureData(R) <- getNetAffx(R, "transcript") # the source is annotation(R) i.e. "pd.ragene.2.1.st"
fData(R)[1:5,] # looks like nothing here but see later

# get annotation second way
library(ragene21sttranscriptcluster.db)
columns(ragene21sttranscriptcluster.db)

xx <- as.list(ragene21sttranscriptclusterENTREZID)
xx[1:5]
length(xx)
dim(R)

x <- select(ragene21sttranscriptcluster.db, keys = rownames(R), keytype = "PROBEID", columns = "ENTREZID")

# compare
which(rownames(R) %in% as.character(17610290:17610310)) # this is where probeset ids start differ from probe ids, compare plot earlier
fData(R)[5083:5087,1:8]
x[5080:5100,]
strsplit(fData(R)["17610327","geneassignment"], " /// ")

# genes missing egid in ragene21sttranscriptcluster.db and pd.ragene.2.1.st compared
xna <- unique(x[is.na(x[,2]),1])
length(xna)
fna <- rownames(R)[is.na(fData(R)$geneassignment)]
length(fna)
miss <- xna[!xna %in% fna]
fna[!fna %in% xna]

fData(R)[miss, "geneassignment"]
fData(R)[miss, "geneassignment"][615] # probe with egid in fData but not in x
# indeed no such egid in ragene21sttranscriptcluster.db:
select(ragene21sttranscriptcluster.db, keys = "682635", keytype = "ENTREZID", columns = "PROBEID")

# extract egids from fData, for those probes that have egis in fData but not in x
miss <- fData(R)[miss, "geneassignment"]
miss <- strsplit(miss, " /// ")
miss <- lapply(miss, function(x){
  ans <- strsplit(x, " // ")
  unique(sapply(ans, function(x) tail(x, 1)))
})
miss <- unique(unlist(miss))
miss <- miss[order(as.numeric(miss))] # looks like they are pseudogenes - check examples in ncbi

# optimal slution - use ragene21sttranscriptcluster.db - easier


# Remove ambigous probesets ####

# choose columns for feature data
library(ragene21sttranscriptcluster.db)
fd <- select(ragene21sttranscriptcluster.db, keys = rownames(R), keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL", "GENENAME"))
length(unique(fd$PROBEID))
dim(fd)
head(fd)
# remove NAs
fd <- fd[!is.na(fd$ENTREZID),]
dim(fd)
head(fd)
# count frequency of probesets
freq <- table(fd$PROBEID)
length(freq)
sum(freq > 1)
# select unique probesets
dambg <- names(freq)[freq == 1]
# subset feature data to unique probesets
fd <- fd[fd$PROBEID %in% dambg, ]
dim(fd)
head(fd)
rownames(fd) <- fd$PROBEID

# subser eSet to unique probes with annotations
dim(R)
R <- R[rownames(fd)]
dim(R)
# replace fData
fData(R) <- fd
# check dims
dim(fd)
dim(R)
dim(fData(R))

# keep lines with at least one observation per group and not constant observations, i.e. sd in each group != 0
(x <- exprs(R)[1,])
sel <- apply(exprs(R), 1, function(x){
  sds <- tapply(x, pData(R)$class, sd, na.rm = T)
  all(sds != 0)
})

R <- R[sel,]

# keep only one probeset per gene; this will be the one with highest variance
sel <- tapply(seq_along(fData(R)$ENTREZID), fData(R)$ENTREZID, c)
sel["100233206"]

table(sapply(sel, length))

(x <- which(fData(R)$ENTREZID == "100233206"))
fData(R)[x,]

sel <- tapply(seq_along(fData(R)$ENTREZID), fData(R)$ENTREZID, function(x){
  if (length(x) == 1) x else{
    vars <- apply(exprs(R)[x,], 1, var, na.rm = T)
    x[which.max(vars)]
  }
})

sel <- as.vector(sel)
sel["100233206"]
sel[1:5]

R <- R[sel,]


# Convert to human ####

# get the data
download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data", "homologene.data")

#read the file
hom <- read.delim("homologene.data", header = F,
                  colClasses = c("character", "character", "character", "NULL", "NULL", "NULL"),
                  col.names = c("id", "tax", "entrez", "NULL", "NULL", "NULL"))

# restric to human and rat
hom <- hom[hom$tax %in% c("9606", "10116"),]

# remove disambiguities (and no pair), i.e. for each id need exactly one human and one rat gene
table(table(hom$id))
sel <- tapply(hom$tax, hom$id, function(x){
  length(x) == 2 & length(unique(x)) == 2
})

i <- names(sel)[!sel][1:2]
hom[hom$id == i[1],]
hom[hom$id == i[2],]

sel <- names(sel)[sel]

hom <- hom[hom$id %in% sel,]

table(table(hom$id))

# convert fData
ids <- hom$id[match(fData(R)$ENTREZID, hom$entrez)]
ens <- hom$entrez[match(ids, hom$id)]

fData(R)$ENTREZID[1:5]
ids[1:5]
ens[1:5]
hom[hom$id %in% ids[1:5],]
sum(is.na(ens))

fData(R)$HUMENTREZID <- ens
dim(R)
R <- R[!is.na(fData(R)$HUMENTREZID),]
dim(R)

rownames(R) <- fData(R)$HUMENTREZID

# To functions ####

read_cell <- function(dir, pData){
  # get sample meta data
  tab <- read.delim(paste(dir, pData, sep = "/"), stringsAsFactors = F)
  rownames(tab) <- tab$filename

  # get microarray signal
  ab <- read.celfiles(filenames = list.celfiles(dir, full.names = T), phenoData = AnnotatedDataFrame(tab))
  ab
}


rma_core <- function(ab){
  R <- oligo::rma(ab, target = "core")
  R
}


annotate <- function(R, annot){
  # choose columns for feature data
  fd <- select(annot, keys = rownames(R), keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL", "GENENAME"))

  fd <- fd[!is.na(fd$ENTREZID),]

  freq <- table(fd$PROBEID)
  dambg <- names(freq)[freq == 1]
  fd <- fd[fd$PROBEID %in% dambg, ]
  rownames(fd) <- fd$PROBEID

  R <- R[rownames(fd)]
  fData(R) <- fd

  R
}


dedegen <- function(R){
  # keep only one probeset per gene; this will be the one with highest variance
  sel <- tapply(seq_along(fData(R)$ENTREZID), fData(R)$ENTREZID, c)

  sel <- tapply(seq_along(fData(R)$ENTREZID), fData(R)$ENTREZID, function(x){
    if (length(x) == 1) x else{
      vars <- apply(exprs(R)[x,], 1, var, na.rm = T)
      x[which.max(vars)]
    }
  })

  sel <- as.vector(sel)

  R <- R[sel,]
  R
}


to_human <- function(R, hom, taxid){
  # restric to human and taxid
  hom <- hom[hom$tax %in% c("9606", taxid),]

  # remove disambiguities (and no pair), i.e. for each id need exactly one human and one taxid gene
  sel <- tapply(hom$tax, hom$id, function(x){
    length(x) == 2 & length(unique(x)) == 2
  })
  sel <- names(sel)[sel]
  hom <- hom[hom$id %in% sel,]

  # convert fData
  ids <- hom$id[match(fData(R)$ENTREZID, hom$entrez)]
  ens <- hom$entrez[match(ids, hom$id)]

  fData(R)$HUMENTREZID <- ens
  R <- R[!is.na(fData(R)$HUMENTREZID),]
  rownames(R) <- fData(R)$HUMENTREZID

  R
}


get_exprs <- function(dir, pData, annot, hom, taxid){
  ab <- read_cell(dir, pData)
  R <- rma_core(ab)
  R <- annotate(R, annot)
  R <- dedegen(R)
  R <- to_human(R, hom, taxid)
  R
}
