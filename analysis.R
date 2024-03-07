library(dplyr)
library(RMySQL)
library(ShortRead)
library(lubridate)
library(parallel)
library(GenomicRanges)
source('lib.R')

seqDataPath <- '454_run_data'
refDBpath   <- '~/data_current/referenceGenomes/blat/hg38.2bit' # Available via UCSC.
dbConn      <- dbConnect(MySQL(), group='core')                 # DB dump included.

d <- bind_rows(lapply(c(395, 309, 300, 117, 104, 105), function(run){
       message(paste0('run: ', run))
       d <- dbGetQuery(dbConn, paste0('select * from runs_samples where run_accession = "', run, '"'))
       o <- readFastq(paste0(seqDataPath, paste0('/run_', run, '.fastq')))
       b <- narrow(o, 1, 8)
       bind_rows(lapply(split(d, 1:nrow(d)), function(x){
         x$sample_name <- toupper(x$sample_name)
         message(paste0('   sample: ', x$sample_name))
         adaptSeq <- substr(x$linker_sequence, 1, 12)
         k  <- o[b@sread == x$barcode_sequence]
         if(length(k) > 0){
           k2 <- subseq(k@sread, 9, width(k))
           
           #browser()
           #i <- vcountPattern(x$primer_sequence, narrow(k2, 1, nchar(x$primer_sequence)+1), max.mismatch = 3, with.indels = TRUE) == 1
           
           names(k2) <- as.character(k@id)
           writeFasta(k2, 'tmpDNAseq.fasta')
           system(paste0('cutadapt -g X', x$primer_sequence, ' -a ', adaptSeq, ' -e 0.15 -O 5 -o tmpDNAseq2.fasta tmpDNAseq.fasta'))
           k3 <- readFasta('tmpDNAseq2.fasta')
           invisible(file.remove(c('tmpDNAseq.fasta', 'tmpDNAseq2.fasta')))
           
           enzyme <- ifelse(grepl('MU', x$sample_name), 'MU', unlist(strsplit(x$sample_name, '\\-'))[6])
           return(tibble(setName = x$sample_name, 
                         readID = paste0(as.character(k3@id), '-', enzyme), 
                         readSeq = as.character(k3@sread),
                         fullReadSeq = as.character(k[k@id == k3@id]@sread),
                         fullReadQual = as.character(k[k@id == k3@id]@quality@quality))) 
         } else {
           return(tibble())
         }
         }))
     }))

# Update set names to match the published set names.
d$setName <- sub('Brady', 'Brady-pseudo', d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p410-d060606-1326ng-Mu-PCR1only-Aug2010', 'Brady-pseudo-CD4zeta-Tcell-p410-d060606-1326ng-MuPCR1only-Aug2010', d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p410-d060606-663ng-Mu-PCR1and2-Aug2010',  'Brady-pseudo-CD4zeta-Tcell-p410-d060606-663ng-MuPCR1and2-Aug2010',  d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p410-d060606-663ng-Mu-PCR1only-Aug2010',  'Brady-pseudo-CD4zeta-Tcell-p410-d060606-663ng-MuPCR1only-Aug2010',  d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p707-d062299-2343ng-Mu-PCR1and2-Aug2010', 'Brady-pseudo-CD4zeta-Tcell-p707-d062299-2343ng-MuPCR1and2-Aug2010', d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p707-d062299-2343ng-Mu-PCR1only-Aug2010', 'Brady-pseudo-CD4zeta-Tcell-p707-d062299-2343ng-MuPCR1only-Aug2010', d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p707-d062299-Mu-PCR1only-Aug2010',        'Brady-pseudo-CD4zeta-Tcell-p707-d062299-Mu-PCR1onlyAug2010',  d$setName, ignore.case = TRUE)
d$setName <- sub('Brady-pseudo-CD4zeta-Tcell-p801-d111898-Mse-PCR1only-Sep2010',       'Brady-pseudo-CD4zeta-Tcell-p801-d111898-Mse-PCR1onlySep2010', d$setName, ignore.case = TRUE)
d$setName <- toupper(d$setName)


# Limit samples to those found in the previous paper, "Decade-Long Safety and Function of Retroviral-Modified Chimeric Antigen Receptor T-cells."
published_sets <- toupper(unique(unlist(lapply(strsplit(readLines('table_S5.txt'), '\\s'), '[', 1))))
d <- subset(d, setName %in% published_sets)


# Convert (and correct) dates to a more usable format.
o <- toupper(sub('^D', '', unlist(lapply(strsplit(unique(d$setName), '\\-'), '[', 6))))
o <- sub('[A-Z]$', '', o)
o <- sub('063904', '06302004', o)
o <- mdy(o)
k <- tibble(setName = unique(d$setName), date = o)
d <- left_join(d, k, by = 'setName')


# Extract subject ids from sequencing set names.
d$subject <- unlist(lapply(strsplit(d$setName, '\\-'), '[', 5))


# Adapter trimmed reads must be at least 25nt long to progress.
d$readSeqLength <- nchar(d$readSeq)
d <- subset(d, readSeqLength >= 25)


# Align the first 25 NTs of reads to the vector and remove reads that align
# well to the vector because these are more than likely internal vector reads.

o <- DNAStringSet(substr(d$readSeq, 1, 25))
names(o) <- d$readID
writeFasta(o, 'vectorAlnTest.fasta')

system(paste0('./blat CD4zetaVector.2bit vectorAlnTest.fasta vectorAlnTest.psl ', 
              ' -tileSize=11 -stepSize=9 -repMatch=2048 -out=psl -t=dna ',
              ' -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))

b <- parseBLAToutput('vectorAlnTest.psl') %>%
     dplyr::filter(alignmentPercentID >= 97, tNumInsert <= 1, qNumInsert <= 1, 
                   tBaseInsert <= 2, qBaseInsert <= 2, qStart <= 3, matches >= 20)

d <- subset(d, ! readID %in% b$qName)
invisible(file.remove(c('vectorAlnTest.fasta', 'vectorAlnTest.psl')))


# Create a column to server as a splitting vector and split the reads across 24 
# cores then align each read block to the reference geneome.

d$n <- ntile(1:nrow(d), 50)
cluster <- makeCluster(24)
clusterExport(cluster, 'refDBpath')

invisible(parLapply(cluster, split(d, d$n), function(x){
  library(Biostrings)
  source('lib.R')
  
  t <- tmpFile()
  
  o <- DNAStringSet(x$readSeq)
  names(o) <- x$readID
  writeXStringSet(o, t)
          
  system(paste0('./blat ', refDBpath, ' ', t, ' ', paste0(t, '.psl') , 
                ' -tileSize=11 -stepSize=9 -repMatch=2048 -out=psl -t=dna ',
                ' -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))
  
  invisible(file.remove(t))
}))

              
# Alignments must have a %seq id >= 97 and start no latter than 3 nt from the start of reads.
# Only minor indels will be accepted.

b <- bind_rows(lapply(list.files(pattern = '*.psl$'), function(x){
       b <- parseBLAToutput(x)
       if(nrow(b) == 0) return(tibble())
  
       dplyr::filter(b, alignmentPercentID >= 97, tNumInsert <= 1, qNumInsert <= 1, 
                     tBaseInsert <= 2, qBaseInsert <= 2, qStart <= 3) %>%
       dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tStart, tEnd, alignmentPercentID)
}))

saveRDS(b, 'blatResult.rds')

invisible(file.remove(list.files(pattern = '*.psl$')))


# Alignments must be at least 25nt long.
b$alnLength <- b$qEnd - b$qStart + 1
b <- subset(b, b$alnLength >= 25)

# Find instances where a read has more than one alignment.
# Alignments are already filtered such that they start near the beginning of reads.
# 3' adapter trimming is weak and alignment to the end of trimmed reads can not be enforced. 
# When more than one read alignment is found, accept the longest provided that the second 
# longest is 10nt shorter otherwise remove the read from the analysis. 

b <- bind_rows(lapply(split(b, b$qName), function(x){
       if(nrow(x) > 1){
         x <- arrange(x, desc(alnLength))
         if((x[1,]$alnLength - x[2,]$alnLength) >= 10){
           return(x[1,])
         } else {
           return(tibble())
         }
       } else {
         return(x)
       }
     }))


# Create position ids and add back metadata columns.
b$posid <- paste0(b$tName, b$strand, ifelse(b$strand == '+', b$tStart, b$tEnd))
b <- left_join(b, select(d, readID, subject), by = c('qName' = 'readID'))

frags <- bind_rows(lapply(split(b, paste(b$subject, b$posid)), function(x){
           tibble(subject = x$subject[1], 
                  readIDs = list(x$qName),
                  seqnames = x$tName[1], 
                  strand = x$strand[1], 
                  start = as.integer(mean(x$tStart)), 
                  end = as.integer(mean(x$tEnd)), 
                  reads = n_distinct(x$qName))
         }))


# Standardize fragments.
frags$n <- 1:nrow(frags)

frags2 <- bind_rows(lapply(split(frags, frags$subject), function(x){
            g <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
            g2 <- standardize_sites(g, counts.col = 'reads', sata.gap = 5)
            as.data.frame(g2)
       })) %>% arrange(n) %>% mutate(posid = paste0(seqnames, strand, ifelse(strand == '+', start, end)))
  

# Unnest reads to recreate read-level data.
frags3 <- tidyr::unnest(frags2, readIDs)
frags3 <- left_join(frags3, select(d, readID, setName), by = c('readIDs' = 'readID'))


# Compare numbers of recovered sites to previous published result.
r1 <- group_by(frags3, setName) %>% summarise(nReads = n_distinct(readIDs), nSites = n_distinct(posid)) %>% ungroup()
r0 <- read.table('table_S5.txt', header = FALSE)
r0 <- r0[, 1:3]
r0$V1 <- toupper(r0$V1)
r0 <- left_join(r0, r1, by = c('V1' = 'setName'))
names(r0) <- c('setName', 'alignableReads1', 'nSites1', 'alignableReads2', 'nSites2')

readr::write_tsv(r0, 'table_S5_comparison.tsv')


# Write out demultiplexed fastq files.
invisible(lapply(split(d, paste(d$subject, d$date)), function(x){
  fileName <- paste0('CD4zeta_', gsub('-', '', paste0(x$subject[1], '_', as.character(x$date[1]))), '_R1.fastq.gz')
  
  R1 <- ShortRead::ShortReadQ(sread = DNAStringSet(x$fullReadSeq), 
                              id = BStringSet(x$readID), 
                              quality = BStringSet(x$fullReadQual))
  
  writeFastq(R1, file = file.path('fastq', fileName), compress = TRUE)
}))

