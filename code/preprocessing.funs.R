# FUNCTIONS
# image pre-processing

require(magick)
require(magrittr)

crop.images <- function(image, crop.param="740x580+50+20"){

    i <- image_read(path=image) %>% 
         image_crop(crop.param) %>%
         image_convert(colorspace='gray')
    
    return(i)

}

image.contrast <- function(image, brightness=10, contrast=5){
    contrast <- 100 + contrast
    brightness <- 100 + brightness
    i <- image_modulate(image, brightness=brightness) %>% # 1) brightness
         image_contrast(sharpen = contrast)   # 2) contrast
    
    return(i)
}


image.denoise <- function(image, noise.radius=1L, despeckle.iter=1L){
    i <- image_reducenoise(image, radius=noise.radius) %>%  # 3) noise reduction
         image_despeckle(times = despeckle.iter)         # 4) despeckle
    
    return(i)
}


apply.filters <- function(image, b1=20, b2=15, c1=60, c2=25, nr=4L, di=4L){
i <- image.contrast(image, brightness=b1, contrast=c1)
i <- image.denoise(i, noise.radius=nr, despeckle.iter=di)
i <- image.contrast(i, brightness=b2, contrast=c2)
return(i)
}
