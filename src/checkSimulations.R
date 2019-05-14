###############################################################################
# Check simulations contain all pairs
###############################################################################

## Count files with output
find ./ -name outputchr*.txt | wc

for i in {1..100}
do
 echo sim$i `(ls sim$i/outputchr*.txt | wc -l)` | grep -v 22$
done 


## Make files with lines per simulation
ls -d */ -I lista>lista

for j in {1..22}
do
  for i in $(cat lista)
  do
  wc -l "$i"outputchr$j.txt>>lineas_$j.txt
  done
done

# Sexual chromosomes
j=X
for i in $(cat lista)
do
wc -l "$i"outputchr$j.txt>>lineas_$j.txt
done

j=Y
for i in $(cat lista)
do
wc -l "$i"outputchr$j.txt>>lineas_$j.txt
done

# Autosomes
### Check all files have the same pairs
badSims <- lapply(1:22, function(chr) {
  chrlist <- read.table(paste0("lineas_", chr, ".txt"), as.is = TRUE)
})
sapply(badSims, function(x) length(unique(x[,1])))



# Autosomes
## Load pairs per chromosome
makePairsDF <- function(file){
  pairs <- readLines(file)
  pairs <- pairs[grep("\"chr", pairs)]
  pairs <-  gsub('[1] \"', "", pairs, fixed = TRUE)
  pairs <-  gsub('\"', "", pairs, fixed = TRUE)
  
  a <- strsplit(pairs, ": ")
  pairs <- data.frame(chr = sapply(a, `[`, 1), pairs = sapply(a, `[`, 2), 
                      stringsAsFactors = FALSE)
  rownames(pairs) <- pairs$chr
  
  pairs
}

pairs <- makePairsDF("divide.out")

## Check files
badSims <- sapply(1:22, function(chr) {
  chrlist <- read.table(paste0("lineas_", chr, ".txt"), as.is = TRUE)
  
  num <- pairs[paste0("chr", chr), 2]
  
  chrlist[chrlist$V1 != num, 2]
})
badSims <- unlist(badSims)

badSims <- lapply(1:100, function(sim){
  badSims[grep(paste0("sim", sim, "/"), badSims)]
})
badSims <- badSims[lengths(badSims) > 0]

# Female
pairs <- makePairsDF("divideFemale.out")

badSims <- sapply(c(1:22, "X"), function(chr) {
  chrlist <- read.table(paste0("female/lineas_", chr, ".txt"), as.is = TRUE)
  
  num <- pairs[paste0("chr", chr), 2]
  
  chrlist[chrlist$V1 != num, 2]
})
badSims <- unlist(badSims)

badSims <- lapply(1:100, function(sim){
  badSims[grep(paste0("sim", sim, "/"), badSims)]
})
badSims <- badSims[lengths(badSims) > 0]


# Male
pairs <- makePairsDF("divideMale.out")

badSims <- sapply(c(1:22, "X", "Y"), function(chr) {
  chrlist <- read.table(paste0("male/lineas_", chr, ".txt"), as.is = TRUE)
  
  num <- pairs[paste0("chr", chr), 2]
  
  chrlist[chrlist$V1 != num, 2]
})
badSims <- unlist(badSims)

badSims <- lapply(1:100, function(sim){
  badSims[grep(paste0("sim", sim, "/"), badSims)]
})
badSims <- badSims[lengths(badSims) > 0]