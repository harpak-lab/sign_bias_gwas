## Clean UKB 
GWAS 
#Remove scientific notation
options(scipen = 999)  
options(warn=1)

#Argument, libraries, and settings
suppressMessages(library("optparse")) #load opt parser library
suppressMessages(library("dplyr"))
suppressMessages(library("dtplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("data.table"))
suppressMessages(library("Biostrings"))

option_list = list(
  
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input GWAS summary file", metavar="character"),
  make_option(c("-b", "--blocks"), type="character", default=NULL, 
              help="LD block file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input GWAS summary file).n", 
       call.=FALSE)
}
if (is.null(opt$blocks)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (LD block file).n", 
       call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", 
       call.=FALSE)
}

# Load GWAS data
print("------------------------")
print("Loading GWAS data...")

sdat <- fread(opt$file, sep="\t")

print("Data loaded.")

# Load UKB variants
variants<-fread(file="Neale.variants.tsv", 
  sep="\t")

variants <- variants %>% select(-AC,-minor_AF,-minor_allele)
sdat <- left_join(sdat, variants, join_by(variant)) 
rm(variants)

# Remove split multi-allelic entries 
sdat <- sdat[ , if (.N == 1) .SD, by = rsid]
sdat <- sdat[ , if (.N == 1) .SD, by = .(chr,pos)]

#subset matched variants
vars_matched <- fread(file="aou_matched_rsids.tsv",
                      sep="\t")
sdat <- sdat[rsid %in% vars_matched$rsid]

# Load LD block file
ld <- read.table(opt$blocks, sep="", 
                 header = TRUE, 
                 check.names=FALSE, 
                 stringsAsFactors = FALSE)
ld$chr <- as.numeric(substr(ld$chr, 4, nchar(ld$chr)))
ld$block_id <- 1:nrow(ld)

# separate variant columns (keep only autosomes)
print("------------------------")
print("Processing data...")

sdat <- sdat %>% filter(!chr %in% c("X", "Y")) %>%
  mutate(chr = as.numeric(chr),
         pos = as.numeric(pos)) %>%
  filter(as.numeric(chr) %in% 1:22) %>%
  filter(!is.na(se)) %>% #remove NA SEs
    
  # polarize effects w.r.t to minor alleles (alt is effect allele in UKB)
  mutate(alt_freq = ifelse(alt == minor_allele, 
                           minor_AF, 1 - minor_AF),
         alt_is_minor = alt_freq <= 0.5,
         m_Beta = ifelse(alt_is_minor, beta, -beta),
         MAF = minor_AF)


# assign ancestral states
print("------------------------")
print("Assigning ancestral states...")

adir<-"/ancestral_seqs/homo_sapiens_ancestor_"
anc <- c()
for (i in 1:22) {
  fasta_file <- paste0(adir, 
                       i, ".fa")
  aseq <- readDNAStringSet(fasta_file)
  positions <- sdat %>%
    filter(as.numeric(chr) == i) %>%
    pull(pos)
  positions<-as.numeric(positions)
  print(paste0("Loaded chr ", i, " anc sequence.."))
  views <- Views(aseq[[1]], start = positions, width = 1)
  anc <- c(anc, toupper(as.character(views)))
  print(paste0("Retrieved chr ", i, " anc sequence positions.."))
}
sdat <- sdat %>%
  mutate(anc_allele = unlist(anc))

rm(aseq)
rm(views)

# Polarize wrt to derived alleles
sdat<- sdat%>%
  mutate(
    alt_is_anc = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'), alt == anc_allele, NA),
    d_beta = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'),
                    ifelse(alt_is_anc, -beta, beta), NA),
    DAF = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'),
                 ifelse(alt_is_anc, 1 - alt_freq, alt_freq), NA)
  )

# Assign SNPs to LD blocks
print("------------------------")
print("Assigning LD blocks...")

setDT(ld)

sdat <- sdat[, block_id :=
       ld[.SD,
          on = .(chr,
                 start <= pos,
                 stop  >= pos),
          mult   = "first",      # pick the first match if overlapping
          block_id,             
          nomatch = NA]          # NA for SNPs outside any block
]

sdat <- sdat[ !is.na(block_id) ] #remove SNPs not assigned blocks

#remove unnecessary columns and QC for info and HWE
print("------------------------")
print("Variant QC....")

sdat <- subset(sdat, info >= 0.8 & p_hwe >= 1e-10)
valid <- c("A","C","G","T")
sdat <- sdat[
  nchar(ref) == 1 & nchar(alt) == 1 &
    ref %in% valid  & alt %in% valid
]

#Get minor allele counts
sdat$MAC <- sdat$minor_AF * (2 * sdat$n_complete_samples)

#QC for MAC (>= 20)
sdat<- sdat[MAC >= 20]

#select columns
sdat <- sdat %>% select(chr, pos, rsid, everything())

#Export formatting
print("------------------------")
print("Done. Exporting...")

fwrite(sdat,
       file      = opt$output,   
       sep       = "\t",           # tab delimiter
       quote     = FALSE,          # don’t quote strings
       showProgress = TRUE)

print("Data exported.")