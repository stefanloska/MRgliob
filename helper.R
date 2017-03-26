dir = "Rat_data"
pData = "sample_info.txt"
annot = ragene21sttranscriptcluster.db::ragene21sttranscriptcluster.db
hom <- read.delim("homologene.data", header = F,
                  colClasses = c("character", "character", "character", "NULL", "NULL", "NULL"),
                  col.names = c("id", "tax", "entrez", "NULL", "NULL", "NULL"))
taxid = "10116"



# Review sdv ####

plot.matrix <- function(M) {
  pars <- par(no.readonly = T)
  on.exit(par(pars))
  par(pty = "s")
  plot.default(M, pch = 16, tck = 1)

}



M <- matrix(c(1, 1,
              5, 5,
              2, 4,
              4, 2), 4, 2, byrow = T)


M

plot(M)
S
S <- svd(M)

diag(S$d)


S$u %*% diag(S$d) %*% t(S$v)

M %*% S$v

S$u %*% diag(S$d)


# centered

M <- sweep(M, 2, colMeans(M))

M

plot(M)

S <- svd(M)
S
diag(S$d)


S$u %*% diag(S$d) %*% t(S$v)

M %*% S$v
S$u %*% diag(S$d)

M_ <- S$u %*% diag(S$d)
t(M_) %*% M_
t(S$u %*% diag(S$d)) %*% S$u %*% diag(S$d)
diag(S$d) %*% t(S$u) %*% S$u %*% diag(S$d)
t(S$u) %*% S$u
S$u %*% t(S$u)


S$v %*% t(S$v)
t(S$v) %*% S$v

S <- svd(M, nu = nrow(M), nv = ncol(M))
S

