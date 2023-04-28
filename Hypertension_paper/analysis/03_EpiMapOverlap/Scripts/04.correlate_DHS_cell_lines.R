##### pearson correlation between DHS HEK and cardiovascular #####
### data comes from running jaccard between DHS profiles

### DHS data ###
sample.meta <- data.frame("CM" = c(1, 0.141514, 0.175783),
                          "VSMC" = c(0.141514, 1, 0.191493),
                          "HEK293" = c( 0.175783, 0.191493, 1))

rownames(sample.meta) <- c('CM','VSMC','HEK293')

### plot ####
library(corrplot)
corrplot(as.matrix(sample.meta))

corrplot(as.matrix(sample.meta), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, col = COL1("Blues", 20), col.lim = c(0, 1))

         