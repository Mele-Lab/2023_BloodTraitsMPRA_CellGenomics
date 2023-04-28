### make plot enrichmnent prioritiztion variants #####
### read data ####

top <- read.delim('~/Downloads/Enrichment_GO/_FINAL_GO.csv', sep = ',')
bottom <- read.delim('~/Downloads/metascape_bottom/Enrichment_GO/_FINAL_GO.csv', sep = ',')

### merge results #####
top$prioritized <- 'top'
bottom$prioritized <- 'bottom'

all_res <- rbind(top,bottom)


### top pathways ####
top_pathways <- (top[order(top$LogP),'Description'][1:10])

labels = c("regulation of MAPK cascade"  ="regulation of \nMAPK cascade",
           "cGMP-PKG signaling pathway" ="cGMP-PKG \nsignaling pathway",
           "positive regulation of MAPK cascade" = "positive regulation \nof MAPK cascade",
           "Aldosterone synthesis and secretion" = "Aldosterone synthesis \nand secretion",
           "regulation of anatomical structure size" = "regulation of anatomical \nstructure size",
           "enzyme-linked receptor protein signaling pathway" = "enzyme-linked receptor \nprotein signaling pathway",
           "vascular process in circulatory system" = "vascular process in\n circulatory system",
           "NABA ECM GLYCOPROTEINS" = "NABA ECM \nGLYCOPROTEINS",
           "cyclic nucleotide metabolic process" = "cyclic nucleotide \nmetabolic process",
           "cell population proliferation" = "cell population \nproliferation")

ggplot(all_res[all_res$Description %in% top_pathways,], aes(y=Enrichment, x=prioritized, fill=-(LogP))) + 
  geom_bar(position="dodge", stat="identity") + 
  facet_grid(~ Description, labeller = labeller(Description = labels)) + 
  theme_classic() + xlab('')+
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    strip.text.x = element_text(size=6)) + scale_fill_gradientn(colours = c("yellow", "orange", "red", "darkred"))
