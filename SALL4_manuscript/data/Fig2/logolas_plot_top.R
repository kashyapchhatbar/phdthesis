library(Logolas)
library(universalmotif)
library(svglite)
files <- list.files(path="top_motif/", pattern="*.tsv", 
                    full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, sep="\t", check.names=FALSE, row.names=1,
                  header=TRUE) # load file
  tt <- motif_rc(read_matrix(x, sep="\t", rownames=TRUE))
  
  bg = c(0.32078681018555705, 0.23702305406056712, 0.2645935416499317, 0.1775965941039441)
  names(bg) = c('T', 'C', 'A', 'G')
  
  
  # svglite(gsub("tsv", "svg", x), pointsize=14)
  png(gsub("tsv", "png", x), pointsize=14)
  logomaker(tt["motif"], type="EDLogo", bg=bg, colors = c( "#AC9D93", "#D40000", "#FF8080", "#AC9D93"))
  dev.off()
  
})
