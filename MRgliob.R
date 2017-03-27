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

read_cel <- function(data_dir, s_info, skip = NULL){
  # get sample meta data
  tab <- read.delim(paste(data_dir, s_info, sep = "/"), colClasses = c("character", "factor", "character"))
  rownames(tab) <- tab$filename
  # check sample_info vs. cel files
  if (! setequal(tab$filename, list.celfiles(data_dir))) stop("sample info file doesn't fit available cel files")

  # skip what is meant to skip
  tab <- tab[!rownames(tab) %in% skip,]

  # get microarray signal
  ab <- read.celfiles(filenames = paste(data_dir, tab$filename, sep = "/"), phenoData = AnnotatedDataFrame(tab))
  colnames(ab) <- tab$id
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
  # human to human just adds extra column in fData and changes rownames
  if (taxid == "9606") {
    fData(R)$HUMENTREZID <- fData(R)$ENTREZID
    rownames(R) <- fData(R)$HUMENTREZID
    return(R)
  }

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


get_exprs <- function(data_dir, s_info, skip = NULL, annot, hom, taxid){
  ab <- read_cel(data_dir, s_info, skip)
  R <- rma_core(ab)
  R <- annotate(R, annot)
  R <- dedegen(R)
  R <- to_human(R, hom, taxid)
  R
}


hom <- read.delim("homologene.data", header = F,
                  colClasses = c("character", "character", "character", "NULL", "NULL", "NULL"),
                  col.names = c("id", "tax", "entrez", "NULL", "NULL", "NULL"))

Rat <- get_exprs(data_dir = "Rat_data",
                 s_info = "sample_info.txt",
                 annot = ragene21sttranscriptcluster.db::ragene21sttranscriptcluster.db,
                 hom = hom,
                 taxid = "10116")


# PCA ####

e <- exprs(Rat)-rowMeans(exprs(Rat))
s <- svd(e)
dim(s$v)

pars <- par(no.readonly = T)
par(pty = "s")
labs <- s$d^2
labs <- round(labs[1:2]/sum(labs)*100, 0)
labs <- paste(c("PC1: ", "PC2: "), labs, "%", sep = "")
suppressWarnings(plot(s$v[, 1:2], xlab = labs[1], ylab = labs[2], pch = 16, col = as.numeric(pData(Rat)$class) + 1, labels = F, tick = F))
abline(h = 0, v = 0, lty = 3)
text(s$v[, 1], s$v[, 2], pData(Rat)$id, pos = 1)
par(pars)


pca <- function(R, lbs = pData(R)$id){
  # svd
  e <- sweep(exprs(R), 1, rowMeans(exprs(R)))
  s <- svd(e)

  # plot
  pars <- par(no.readonly = T)
  on.exit(par(pars))
  par(pty = "s")
  # calcuklate variance percentage
  labs <- s$d^2
  labs <- round(labs[1:2]/sum(labs)*100, 0)
  labs <- paste(c("PC1: ", "PC2: "), labs, "%", sep = "")
  # draw
  suppressWarnings(plot(s$v[, 1:2], xlab = labs[1], ylab = labs[2], labels = F, tick = F,
                        pch = 16, col = as.numeric(pData(R)$class) + 1))
  abline(h = 0, v = 0, lty = 3)
  text(s$v[, 1], s$v[, 2], lbs, pos = 1)

  invisible(s)
}

pca(Rat)


# GEO data ####

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12657/matrix/GSE12657_series_matrix.txt.gz",
              "GSE12657_series_matrix.txt.gz")

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE12657&format=file",
              "GSE12657_RAW.tar")


# ad hoc way ####
pd <- read.delim("GSE12657_series_matrix.txt.gz", row.names = 1, header = F, skip = 35, nrows = 1)
pd <- rbind(pd, read.delim("GSE12657_series_matrix.txt.gz",row.names = 1, header = F, skip = 56, nrows = 1))
pd <- data.frame(t(pd), stringsAsFactors = F)
colnames(pd) <- c("class", "id")
rownames(pd) <- pd$id
pd$class <- factor(pd$class)
pd <- AnnotatedDataFrame(pd)

ge <- read.delim("GSE12657_series_matrix.txt.gz", row.names = 1, skip = 56, nrow = 12625)
ge <- data.matrix(ge)


# more like parsing way ####
tab <- readLines("GSE12657_series_matrix.txt.gz")
meta <- grepl("^!|^$", tab)

# long way to get pheno data - advantage: extract all fields, disadv: id taken from !Sample_geo_accession instead "ID_REF"
pd <- tab[meta]
pd <- pd[-(1 : grep("^$", pd))]
pd <- pd[seq_len(length(pd) - 2)]
pd <- gsub("\"", "", pd)
pd <- sub("^!", "", pd)
pd <- strsplit(pd , "\t")
pd <- do.call(cbind, pd)
colnames(pd) <- pd[1, ]
pd <- pd[-1, ]
pd <- data.frame(pd, stringsAsFactors = F)
pd <- data.frame(class = factor(pd$Sample_characteristics_ch1), id = pd$Sample_geo_accession, stringsAsFactors = F)
rownames(pd) <- pd$id
pd <- AnnotatedDataFrame(pd)


# short way to get pheno data - just get what we need
pd <- grep("^!Sample_characteristics|^\"ID_REF\"", tab, value = T)
pd <- gsub("\"", "", pd)
pd <- strsplit(pd , "\t")
pd <- do.call(cbind, pd)
pd <- pd[-1, ]
pd <- data.frame(class = factor(pd[, 1]), id = pd[, 2], stringsAsFactors = F)
rownames(pd) <- pd$id
pd <- AnnotatedDataFrame(pd)

ge <- tab[!meta]
ge <- gsub("\"", "", ge)
ge <- strsplit(ge , "\t")
ge <- do.call(rbind, ge)
colnames(ge) <- ge[1, ]
rownames(ge) <- ge[, 1]
ge <- ge[-1, -1]
mode(ge) <- "numeric"

# test the annotation library
library(hgu95av2.db)
fd <- select(hgu95av2.db, keys = c("930_at", "931_at"), keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL", "GENENAME"))

# build eSet object
G <- ExpressionSet(ge, pd)

# add annotations
G <- annotate(G, hgu95av2.db::hgu95av2.db)
G <- dedegen(G)
dim(G)
hom <- read.delim("homologene.data", header = F,
                  colClasses = c("character", "character", "character", "NULL", "NULL", "NULL"),
                  col.names = c("id", "tax", "entrez", "NULL", "NULL", "NULL"))
G <- to_human(G, hom, "9606")
dim(G)
fData(G)[1:5,]
exprs(G)[1:5,1:5]

# To functions ####
read_geo <- function(geo_file){
  tab <- readLines(geo_file)
  meta <- grepl("^!|^$", tab)

  pd <- grep("^!Sample_characteristics|^\"ID_REF\"", tab, value = T)
  pd <- gsub("\"", "", pd)
  pd <- strsplit(pd , "\t")
  pd <- do.call(cbind, pd)
  pd <- pd[-1, ]
  pd <- data.frame(class = factor(pd[, 1]), id = pd[, 2], stringsAsFactors = F)
  rownames(pd) <- pd$id
  pd <- AnnotatedDataFrame(pd)

  ge <- tab[!meta]
  ge <- gsub("\"", "", ge)
  ge <- strsplit(ge , "\t")
  ge <- do.call(rbind, ge)
  colnames(ge) <- ge[1, ]
  rownames(ge) <- ge[, 1]
  ge <- ge[-1, -1]
  mode(ge) <- "numeric"

  G <- ExpressionSet(ge, pd)

  G
}

build_GEO <- function(geo_file, annot, hom, taxid){
  G <- read_geo(geo_file)
  G <- annotate(G, annot)
  G <- dedegen(G)
  G <- to_human(G, hom, taxid)
  G
}

G <- build_GEO("GSE12657_series_matrix.txt.gz", hgu95av2.db::hgu95av2.db, hom, "9606")
pca(G)

