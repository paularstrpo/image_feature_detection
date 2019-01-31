source("code/preprocessing.funs.R")
require(tesseract)
require(splitstackshape)
require(readr)

# ----------------------------
# 0) Options 
# ----------------------------
path <- 'labels'
# b1=20 # brightness step 1
# b2=15 # brightness step 2
# c1=60 # contrast step 1
# c2=25 # constrast step 2
# nr=4L # noise radius for reduction
# di=4L # de-speckle iterations
# ----------------------------


# ----------------------------
# 1) Preprocessing
# ----------------------------

flist <- list.files(path=path, pattern='*.jpeg$', full.names = TRUE)
names <- gsub(x=flist, pattern=paste(path, '/', sep=''), replacement='', fixed=TRUE)

cropped.images <- lapply(flist, crop.images)
pp.images <- lapply(cropped.images, apply.filters)
# ----------------------------

# ----------------------------
# 2) Do OCR & Tidy the Result
# ----------------------------

pp.labels <- sapply(pp.images, image_ocr)
label.df <- data.frame(as.matrix(pp.labels))
colnames(label.df) <-  c('raw_ocr_dump')
label.df$label_file <- factor(gsub('.*labels\\/', '', rownames(histo_label_df)  ))
label.df$Image.ID <- factor(gsub('\\..*', '', gsub('.*-', '', rownames(histo_label_df)  )))

#TODO: swap this for tidyr
label.df <- cSplit(label.df, "raw_ocr_dump", sep="\n", direction = "wide")
# ----------------------------

# ----------------------------
# 3) Load & Merge Metadata
# ----------------------------
meta_label_data <- read.csv('metadata/Image_IDs_MSSM_slides.csv')
histo_label_data.table_merge <- merge(as.data.frame(label.df), meta_label_data, by = 'Image.ID')
histo_label_data.table_merge <- as.data.frame(histo_label_data.table_merge)
histo_label_data.table_merge$raw_ocr_dump_02 <- as.character(histo_label_data.table_merge$raw_ocr_dump_02)

# ----------------------------
# 4) Grab & Prep Samples
# ----------------------------
sample_list <- read_excel('Freezerworks_Physical_Sample_List_06162017.csv')
sample_list_biopsy <- sample_list[sample_list$IBD_SpecimenType %in% 'Biopsies',] #restrict to biopsy GUAIDs

# ----------------------------
# 5) Load Funs for Fuzzy Match
# ----------------------------
source("code/levenshtein.funs.R")