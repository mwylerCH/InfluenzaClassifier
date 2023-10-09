# script to identify closest reference from distance matrix (gotree)

library(data.table)
library(tidyverse)

FOLDER <- '//wsl.localhost/Ubuntu-22.04/home/mwyler/uscita'

Files <- list.files(path= FOLDER, pattern='.distMatrix', all.files=FALSE, 
           full.names=T)

# make table for output
output <- as.data.frame(matrix(nrow = 1, ncol = 11))
colnames(output) <- c('ID', 'Subtipe', 'Genotype', 
                      'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS')

# run through files
for (MatriceName in Files){
  #MatriceName <- "//wsl.localhost/Ubuntu-22.04/home/mwyler/uscita/combined1_comb.distMatrix"
  #MatriceName <- "//wsl.localhost/Ubuntu-22.04/home/mwyler/uscita/combined4_comb.distMatrix"
  Matrice <- fread(MatriceName, data.table = F, header = F)
  
  # keep only column of interest
  colnames(Matrice) <- c('Seq', Matrice$V1)
  Matrice <- Matrice[, !grepl('NrRef', colnames(Matrice))]
  
  # get new sequence name
  dirtyName <- colnames(Matrice)[!grepl('Seq', colnames(Matrice))]
  dirtyName <- unlist(strsplit(dirtyName, '\\|'))
  output$ID[1] <- dirtyName[grepl('witzerland', dirtyName)]
  
  # get Segment
  segmenti <- c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'MP')
  Segment <- dirtyName[grepl(paste0(segmenti, collapse = '|'), dirtyName)]
  
  # get closest distance
  colnames(Matrice)[2] <- 'distance'
  closestRef <- Matrice[Matrice$distance > 0, ] %>% 
    slice_min(distance, n = 1) %>% 
    pull(Seq)
  
  # clean name Reference
  output[1, Segment] <- gsub('.*\\|(NrRef[[:digit:]]+)\\|.*' ,'\\1', closestRef)
}

output[,4:ncol(output)] <- apply(output[,4:ncol(output)], 1, function(x){gsub('NrRef', '', x)})

# Reference genotype -------------------------

# load reference table
folderScript <- dirname(rstudioapi::getSourceEditorContext()$path)
referenza <- fread(paste0(folderScript, '/GenoTypeTable.csv'), 
                   data.table = F, header = T)
colnames(referenza)[8] <- "NA"

# check which combination is the sample
refFilter <- referenza[referenza$PB2 == output$PB2 &
            referenza$PB1 == output$PB1 &
            referenza$PA == output$PA &
            referenza$HA == output$HA &
            referenza$NP == output$NP &
            referenza$`NA` == output$`NA` &
            referenza$MP == output$MP &
            referenza$NS == output$NS, 
          ]
# print out
if (nrow(refFilter) > 0){
  output[, 2:3] <- refFilter[, 1:2]
  fwrite(output[,], sep = "\t")
  
} else {
  fwrite(output[,], sep = "\t")
  stop("ERROR: Could not Identify adeguate reference.", call. = F)
}
