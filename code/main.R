source("code/preprocessing.R")
source("code/levenshtein.funs.R")
# options -------
path <- 'labels'
# b1=20 # brightness step 1
# b2=15 # brightness step 2
# c1=60 # contrast step 1
# c2=25 # constrast step 2
# nr=4L # noise radius for reduction
# di=4L # de-speckle iterations
# ----------------


# ----------------
flist <- list.files(path=path, pattern='*.jpeg$', full.names = TRUE)
names <- gsub(x=flist, pattern=paste(path, '/', sep=''), replacement='', fixed=TRUE)

cropped.images <- lapply(flist, crop.images)
pp.images <- lapply(cropped.images, apply.filters)


