#!
#!decode histology image labels
# library(jpeg)
# setwd('~/analysis/projects/IBD/Prospective_data/histology_images/')
# img <- readJPEG("labels/label-110468.jpeg", native = TRUE)
# img_array <- readJPEG("labels/label-110468.jpeg", native = FALSE)

# #!verify it was read in correctly
# if(exists("rasterImage")){
#       plot(1:2, type='n')
#       rasterImage(img,1,1,2,2)
# }
#!

#STRATEGY
#!image preprocess for contrast / brightness / crop -> OCR -> fuzzy matching on ground-trouth set of Freezerworks GUAIDs for all
#																non standard-of-care samples
#!result: tractable ambiguity space for labels
#!
########ocr tesseract deconvolution for simple label IDs
library(tesseract) #default to english alphanumeric set
library(splitstackshape)
library(tidyverse)
library(RecordLinkage)
#!
setwd('/Users/user/analysis/projects/Prospective_IBD/image_deconvolution/image_feature_detection')
#setwd('~/analysis/projects/IBD/Prospective_data/histology_images/')
#!
#ocr(flist2[which(grepl('112237', flist2))])
#!enhance contrast -brightness-contrast 10x5   -monochrome
#convert -crop 740x580+50+20 test.jpeg test2.jpeg
#!
#!################################################################################################################
#!################################################################################################################
#!croppoing to visual portion 
for j in `find *.jpeg`; do convert -crop 740x580+50+20 $j enhanced2/enhanced_$j; done
#!now enhance contrast somehwat 
#convert -level 50x100% -type grayscale -depth 8
for j in `find *.jpeg`; do convert -level 50x100% -type grayscale -depth 8 $j step2/pp_$j; done
#!convert to 300 dpi
#convert enhanced_label-112675.jpeg -units "PixelsPerInch" -density 300 -resample "300x" test.jpeg
#!
#for j in `find *.jpeg`; do convert -brightness-contrast 10x5 $j enhanced/enhanced_$j; done
#!################################################################################################################
#!################################################################################################################
#! lift images 
flist <- list.files('labels/enhanced2/step2', pattern = '*jpeg', full.names = T)
#!
#test <- ocr(flist[1])
#test <- sapply(flist, ocr )
#!default english
#apply OCR on histology labels
#
extracted_labels_pp <- sapply(flist, ocr )
#
#!
#!################################################################################################################
#!################################################################################################################
#!manage output 

histo_label_df <- data.frame(as.matrix(extracted_labels_pp))
#histo_label_df <- data.frame(as.matrix(test))
colnames(histo_label_df) <-  c('raw_ocr_dump')
histo_label_df$label_file <- factor(gsub('.*labels\\/', '', rownames(histo_label_df)  ))
histo_label_df$Image.ID <- factor(gsub('\\..*', '', gsub('.*-', '', rownames(histo_label_df)  )))
#!
#!split out ocr dump via tidy tools
histo_label_data.table <- cSplit(histo_label_df, "raw_ocr_dump", sep="\n", direction = "wide")
#!read in meta data
#meta_label_data <- read.csv('metadata/Image_IDs_MSSM_slides.csv')
#!merge with master data 
histo_label_data.table_merge <- merge(as.data.frame(histo_label_data.table), meta_label_data, by = 'Image.ID')
histo_label_data.table_merge <- as.data.frame(histo_label_data.table_merge)
histo_label_data.table_merge$raw_ocr_dump_02 <- as.character(histo_label_data.table_merge$raw_ocr_dump_02)
#!
#!use partial matching in known GUAID  table to severely restrict the matching problem to a tractable load 
#!list of all possible GUAIDs
#!grab freezerworks data
#sample_list <- read.csv('Freezerworks_Physical_Sample_List_06162017.csv')
sample_list_biopsy <- sample_list[sample_list$IBD_SpecimenType %in% 'Biopsies',] #restrict to biopsy GUAIDs
#sample_list_biopsy[agrepl('292*45', sample_list_biopsy$GUAID),]
#sample_list_biopsy[agrepl('EFZ 310961 A', sample_list_biopsy$GUAID, max.distance = list(insertions = .2, cost = .17)),]
#!
#!################################################################################################################
#!################################################################################################################
#!maximum levenshtein similarity (and minimum distance)
#
Lsim = function(string, stringVector){
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
Ldist = function(string, stringVector){
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
#!################################################################################################################
#!################################################################################################################

#!print most similar GUAIDs
#!dominant case of image barcoding being placed into dump_02 slot, but not necessarily true
#!
apply(histo_label_data.table_merge, 2, function(x) sum(grepl('Z' ,x)))
#!
histo_label_data.table_merge$GUAID_imputed <- do.call(rbind, lapply(histo_label_data.table_merge[, 4], function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD <- do.call(rbind, lapply(histo_label_data.table_merge[, 4], function(x) Ldist2(x)))
#!take from dump_01
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_01) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 3]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 3]), function(x) Ldist2(x)))
#!take from dump_03
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_03) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 5]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 5]), function(x) Ldist2(x)))
#!dump_04
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_04) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 6]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 6]), function(x) Ldist2(x)))
#!dump_05
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_05) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 7]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 7]), function(x) Ldist2(x)))
#!dump_07
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_07) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 9]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 9]), function(x) Ldist2(x)))
#!dump_09
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_09) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 11]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 11]), function(x) Ldist2(x)))
#!dump_10
iidx <- grepl('Z', histo_label_data.table_merge$raw_ocr_dump_10) 
histo_label_data.table_merge$GUAID_imputed[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 12]), function(x) Lsim2(x)))
histo_label_data.table_merge$minimum_LD[iidx] <- 
				do.call(rbind, lapply(as.matrix(histo_label_data.table_merge[iidx, 12]), function(x) Ldist2(x)))
#!no more shifts

histo_label_data.table_merge$SOC_or_EZF <- factor(ifelse(grepl('SOC', histo_label_data.table_merge$GUAID_imputed), 
														'standard_of_care', 'Freezerworks'))
#!remove faulty matching attempts
iidx <- histo_label_data.table_merge$Block.ID != ""
histo_label_data.table_merge$GUAID_imputed[iidx] <- 'SOC'
#!
histo_label_data.table_merge$GUAID_imputed_length <- as.numeric(nchar(histo_label_data.table_merge$GUAID_imputed))
#histo_label_data.table_merge <- histo_label_data.table_merge[order(histo_label_data.table_merge$GUAID_imputed_length),]
#!
#!################################################################################################################
#!################################################################################################################
write.csv(histo_label_data.table_merge, file = 'results/Image_IDs_MSSM_slides__matched_with_GUAIDs.csv')
Image_IDs_MSSM_slides__matched_with_GUAIDs <- histo_label_data.table_merge
save(Image_IDs_MSSM_slides__matched_with_GUAIDs, file = 'RData/Image_IDs_MSSM_slides__matched_with_GUAIDs.RData')

#!################################################################################################################
#!################################################################################################################



#!################################################################################################################
#!################################################################################################################
#diagnostic plots 
ggplot(histo_label_data.table_merge, aes(log10(GUAID_imputed_length), minimum_LD)) + geom_point()

reconstruction_matching_amguity <- 
	ggplot(histo_label_data.table_merge, aes(log10(GUAID_imputed_length), color = SOC_or_EZF)) + geom_rug() + geom_density() + 
				geom_text_repel(data = 
						histo_label_data.table_merge[log10(histo_label_data.table_merge$GUAID_imputed_length) > 3,], 
						aes(y = .1,  label = Image.ID), size = 2) + 
					ggtitle( 'optical-character-histology-label-text-imputed GUAID ambiguity')

pdf('figures/reconstruction_matching_amguity.pdf', height = 8, width = 10)
reconstruction_matching_amguity
dev.off()

ggplot(histo_label_data.table_merge, aes(minimum_LD)) + geom_rug() + geom_density()


#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!################################################################################################################
#!
histo_label_data.table_merge <- histo_label_data.table_merge[order(-histo_label_data.table_merge$GUAID_imputed_length),]
#!
sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > 1)
sum(histo_label_data.table_merge$GUAID_imputed_length > 3 & histo_label_data.table_merge$GUAID_imputed_length < 9)

sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > log10(22))
sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > log10(34))
#!
#!
sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > 1 & grepl('EF',histo_label_data.table_merge$raw_ocr_dump_02 ))

sum(grepl('PS',histo_label_data.table_merge$raw_ocr_dump_02 ) | grepl('MS',histo_label_data.table_merge$raw_ocr_dump_01 ))
sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > 1 & grepl('PS',histo_label_data.table_merge$raw_ocr_dump_02 ))
sum(log10(histo_label_data.table_merge$GUAID_imputed_length) > 1 & grepl('EF',histo_label_data.table_merge$raw_ocr_dump_02 ))
#!
test <- histo_label_data.table_merge[log10(histo_label_data.table_merge$GUAID_imputed_length) > 1 & 
									!grepl('EF',histo_label_data.table_merge$raw_ocr_dump_02 )
									#!grepl('EF',histo_label_data.table_merge$raw_ocr_dump_02 ) & 
									#!grepl('Ms',histo_label_data.table_merge$raw_ocr_dump_02 ) & 
									#!grepl('MS',histo_label_data.table_merge$raw_ocr_dump_02 ) 
										,]
#!
Lsim('EFZ 299991 A fF', as.character(factor(sample_list_biopsy$GUAID)))
Lsim2('EFZ 299991 A fF')
#!
#!
# Lsim2(histo_label_data.table_merge$raw_ocr_dump_2[1500])
# Lsim2(histo_label_data.table_merge$raw_ocr_dump_2[1])
#!
#!
#!
#!#################################################
#load('Biopsy_PD_perGUAID_RK_11_19_18.RDa')
GRID_GUAID_table <- Biopsy_PD_perGUAID[, c('GRID', 'GUAID')]
#!
master_image_table <- read.csv('metadata/Image_IDs_MSSM_slides.csv')
#!
#!
master_GUAID_biopsy_sheet <- read.csv('metadata/BIOPSIES_FOR_JANSSEN_master_sheet_master.csv')
#!
master_GUAID_biopsy_sheet_list <- list()
master_GUAID_biopsy_sheet_list[[1]] <- master_GUAID_biopsy_sheet[, c()]
#!
#!ideas for histology deconvolution

#########
#install.packages("devtools") If devtools is not available
#devtools::install_github("bnosac/image", subdir = "image.darknet", build_vignettes = TRUE)
# library(image.darknet)

# #define darknet model

# yolo_tiny_voc <- image_darknet_model(type = 'detect', model = 'tiny-yolo-voc.cfg', 
# 		weights = '/Library/Frameworks/R.framework/Versions/3.5/Resources/library/image.darknet/models/tiny-yolo-voc.weights',  
# 						labels = '/Library/Frameworks/R.framework/Versions/3.5/Resources/library/image.darknet/include/darknet/data/voc.names')
# #!
# x <- image_darknet_detect(file = "labels/label-110468.jpeg", 
#                           object = yolo_tiny_voc,
#                           threshold = 0.19)
           
##########




