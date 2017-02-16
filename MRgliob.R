# get sample meta data
tab <- read.delim("Rat_data/sample_info.txt", stringsAsFactors = F)
rownames(tab) <- tab$filename
tab

# get microarray signal
library(oligo)
ab <- read.celfiles(filenames = list.celfiles("Rat_data", full.names = T), phenoData = AnnotatedDataFrame(tab))

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
R <- R[[2]]

# preview
pData(R)
fData(R)
exprs(R)[1:5,]
annotation(R)
experimentData(R)

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



