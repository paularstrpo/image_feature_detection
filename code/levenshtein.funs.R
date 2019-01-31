# FUNCTIONS
# maximum levenshtein similarity (and minimum distance)

require(tesseract)
require(splitstackshape)
require(RecordLinkage)

Lsim <- function(string, stringVector){
#!
#!similarity: 1 - d(str1,str2) / max(A,B)
  if(!grepl(glob2rx('MS|Ms'), string)){									#all MS (standard-of-care strings) are easily lifted
  similarity = levenshteinSim(string, stringVector);
  out <- stringVector[similarity == max(similarity)]
  						  }else {out <- 'SOC'}
  out <- paste(out, collapse = '__')
  return(out)
#!
}
#!
Lsim2 <- function(string) {
		Lsim(string, as.character(factor(sample_list_biopsy$GUAID)))
						  }
#!
Ldist <- function(string, stringVector){
#!
#!similarity: 1 - d(str1,str2) / max(A,B)
  if(!grepl(glob2rx('MS|Ms'), string)){									#all MS (standard-of-care strings) are easily lifted
  Ldist = levenshteinDist(string, stringVector);
  out <- min(Ldist)
  						  }else {out <- NA}
  return(as.numeric(out))
#!
}
#!
Ldist2 <- function(string) {
		Ldist(string, as.character(factor(sample_list_biopsy$GUAID)))
						  }


    function(label.path, label.format, label.metatable){

        flist <- list.files(path=label.path, pattern = paste('*', label.format, '$', sep=''), full.names = TRUE)
        extracted.labels <- sapply(flist, ocr)
        extracted.labels <- data.frame(as.matrix(extracted.labels))

        
        colnames(extracted.labels) <-  c('raw_ocr_dump')
        extracted.labels$label_file <- factor(gsub('.*labels\\/', '', rownames(extracted.labels)))
        extracted.labels$Image.ID <- factor(gsub('\\..*', '', gsub('.*-', '', rownames(extracted.labels))))
        
        #!split out ocr dump via tidy tools
        extracted.labels.tidy <- cSplit(extracted.labels, "raw_ocr_dump", sep="\n", direction = "wide")

        histo_label_data.table_merge <- merge(as.data.frame(extracted.labels.tidy), label.metatable, by = 'Image.ID')
        histo_label_data.table_merge <- as.data.frame(histo_label_data.table_merge)
        histo_label_data.table_merge$raw_ocr_dump_02 <- as.character(histo_label_data.table_merge$raw_ocr_dump_02)
    }