# subscript to identify closest reference from distance matrix (gotree)
# Michele Wyler, IVI Mittelhäusern

suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(tidyverse)))

args <- commandArgs(TRUE)
# transposed matrix (marker as row)
FOLDER <- args[2]
folderScript <- args[1]

# FOLDER <- '//wsl.localhost/Ubuntu-22.04/home/mwyler/uscita'
#folderScript <- '//wsl.localhost/Ubuntu-22.04/home/mwyler/GenotapeClassifier/'

Files <- list.files(path= FOLDER, pattern='.distMatrix', all.files=FALSE, 
                    full.names=T)

# make table for output
output <- as.data.frame(matrix(nrow = 1, ncol = 11))
colnames(output) <- c('ID', 'Subtipe', 'Genotype', 
                      'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS')

# run through distance matrix files
for (MatriceName in Files){
  # MatriceName <- "//wsl.localhost/Ubuntu-22.04/home/mwyler/uscita/combined1_comb.distMatrix"
  Matrice <- fread(MatriceName, data.table = F, header = F)
  
  # keep only column of interest (sample)
  colnames(Matrice) <- c('Seq', Matrice$V1)
  Matrice <- Matrice[, !grepl('ID[[:digit:]]+SEQ', colnames(Matrice))]
  Matrice <- Matrice[, !grepl('REF[[:digit:]]+Nr', colnames(Matrice))]
  
  # keep only Row of interest (references)
  Matrice <- Matrice[grepl('REF[[:digit:]]+Nr', Matrice$Seq), ]
  
  
  # get new sequence name (common name without segment information)
  dirtyName <- colnames(Matrice)[!grepl('Seq', colnames(Matrice))]
  dirtyName <- unlist(strsplit(dirtyName, '\\|'))
  CleanName <- dirtyName[grepl('[⁠[:alnum:]]{4,}', dirtyName)]
  if (length(CleanName) != 1){
    output$ID[1] <- 'InputSample'
  } else {
    output$ID[1] <- CleanName
  }
  
  # get Segment of the input
  segmenti <- c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'MP')
  Segment <- dirtyName[grepl(paste0('\\b', segmenti, '\\b', collapse = '|'), dirtyName)]
  
  
  # get closest reference distance
  colnames(Matrice)[2] <- 'distance'
  closestRef <- Matrice[Matrice$distance > 0, ] %>% 
    slice_min(distance, n = 1) %>% 
    pull(Seq)
  
  # clean name Reference
  output[1, Segment] <- gsub('REF([[:digit:]]+)Nr' ,'\\1', closestRef)
}

output[1, 'HA'] <- 20
output[,4:ncol(output)] <- apply(output[,4:ncol(output)], 1, function(x){gsub('NrRef', '', x)})

# Reference genotype/subtype -------------------------

# load reference table
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

