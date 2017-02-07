# AnnotationDbi::inpIDMapper code extracted for analysis
# found a bug:
# function fails with uninformative error when there are no homology hits
# wanted to know what happens

ids = c("11363", "11364")
srcSpecies = "MUSMU"
destSpecies = "HOMSA"
srcIDType = "EG"
destIDType = "EG"
keepMultGeneMatches = FALSE
keepMultProtMatches = FALSE
keepMultDestIDMatches = TRUE

ids = c("11363", "11364") # fine
AnnotationDbi::inpIDMapper(ids, srcSpecies, destSpecies,
                           srcIDType = srcIDType, destIDType = destIDType,
                           keepMultGeneMatches = keepMultGeneMatches,
                           keepMultProtMatches = keepMultProtMatches,
                           keepMultDestIDMatches = keepMultDestIDMatches)

ids = c("11364") # error
AnnotationDbi::inpIDMapper(ids, srcSpecies, destSpecies,
                           srcIDType = srcIDType, destIDType = destIDType,
                           keepMultGeneMatches = keepMultGeneMatches,
                           keepMultProtMatches = keepMultProtMatches,
                           keepMultDestIDMatches = keepMultDestIDMatches)


if (srcSpecies == destSpecies) {
  stop("The srcSpecies and destSpecies should be different. 'No matter where you go, there you are.' - Bucakaroo Banzai")
}
if (srcIDType == "SYMBOL") {
  stop("I refuse to attempt this on the grounds that SYMBOLS are nearly worthless as IDs, and it would be irresponsible to proceed.  If you must use symbols, you will have to process them one step at a time and double check your work, because they cannot be relied upon to be uniquely mapped onto a single gene.")
}
setupVals = AnnotationDbi:::.getMappingData(srcSpecies)
srcSpcAbrv = setupVals[1]
srcDBAbrv = setupVals[2]
if (!is.na(setupVals[3])) {
  protMap = get(setupVals[3])
}
centralID1 = setupVals[4]
require(paste0("hom.", srcSpcAbrv, ".inp.db"), character.only = TRUE)
homMap = get(paste0("hom.", srcSpcAbrv, ".inp", destSpecies))
if (srcIDType != centralID1) {
  toSrcEGMap = get(paste0("org.", srcSpcAbrv, ".", srcDBAbrv,
                          srcIDType))
}
mapBackVals = AnnotationDbi:::.getMappingData(destSpecies)
destSpcAbrv = mapBackVals[1]
destDBAbrv = mapBackVals[2]
centralID2 = mapBackVals[4]
if (!is.na(mapBackVals[3])) {
  geneMap = get(mapBackVals[3])
}
if (srcIDType == centralID1) {
  genes = ids
  names(genes) = ids
} else {
  ids = AnnotationDbi:::.cleanup(ids)
  genes = mget(as.character(ids), revmap(toSrcEGMap), ifnotfound = NA)
  genes = AnnotationDbi:::.cleanup(genes)
}
genes = AnnotationDbi:::.handleMultipleMatches(genes, keepMultGeneMatches)
genes = AnnotationDbi:::.cleanup(genes)
if (exists("protMap", inherits = FALSE)) { # TRUE
  inpIDs = mget(as.character(genes), protMap, ifnotfound = NA)
  inpIDs = AnnotationDbi:::.reLabel(genes, inpIDs, "inpIDs")
  inpIDs = AnnotationDbi:::.cleanup(inpIDs) # inpIDs = named list()
} else {
  inpIDs = genes
}
destList = lapply(inpIDs, function(x) {
  mget(as.character(x), homMap, ifnotfound = NA)
}) # destList = named list()
destIDs = lapply(destList, function(x) {
  x = x[!is.na(x)]
}) # destIDs = named list()
destIDs = lapply(destIDs, unlist) # destIDs = named list()
dindex = vector()
for (i in seq_len(length(destIDs))) {
  dindex = c(dindex, !is.null(destIDs[[i]]))
} # dindex = logical(0)
dnames = names(destList)[dindex] # dnames = character(0)
destIDs = destIDs[dindex] # destIDs = named list()
destIDs = AnnotationDbi:::.handleMultipleMatches(destIDs, keepMultProtMatches) # destIDs = named list()
destIDs = unlist(destIDs) # destIDs = NULL
if (length(destIDs) == length(dnames)) {
  names(destIDs) = dnames
  ########################
  # Error in names(destIDs) = dnames : attempt to set an attribute on NULL
  ########################
  # destIDs <- structure(destIDs, .Names = dnames)
} else {
  stop("Names are not congruent")
}
finIDs = AnnotationDbi:::.handleMultipleMatches(destIDs, keepMultProtMatches)
uniqIDs = AnnotationDbi:::.cleanup(finIDs)
if (exists("geneMap", inherits = FALSE)) {
  EGIDs = mget(as.character(uniqIDs), revmap(geneMap),
               ifnotfound = NA)
  EGIDs = AnnotationDbi:::.reLabel(uniqIDs, EGIDs, "EGIDs")
  EGIDs = AnnotationDbi:::.cleanup(EGIDs)
} else {
  EGIDs = uniqIDs
}
if (destIDType != centralID2 && toupper(destDBAbrv) != destIDType) {
  resultMap = get(paste0("org.", destSpcAbrv, ".", destDBAbrv,
                         destIDType))
  resultIDs = mget(as.character(EGIDs), resultMap, ifnotfound = NA)
} else {
  resultIDs = EGIDs
}
resultIDs = AnnotationDbi:::.reLabel(EGIDs, resultIDs, "resultIDs")
if (keepMultDestIDMatches == FALSE) {
  resultIDs = AnnotationDbi:::.handleMultipleMatches(resultIDs, keepMultiples = TRUE)
}
resultIDs = AnnotationDbi:::.cleanup(resultIDs)
resultIDs
