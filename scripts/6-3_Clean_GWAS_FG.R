### Clean GWAS summary statistics (FinnGen)

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
suppressMessages(library("magrittr"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("GenomicRanges"))


option_list = list(
  
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input GWAS summary file", metavar="character"),
  make_option(c("-b", "--blocks"), type="character", default=NULL, 
              help="LD block file", metavar="character"),
  make_option(c("-N", "--sample_size"), type="character", default=0, 
              help="GWAS sample size", metavar="numeric"),
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
if (is.null(opt$sample_size)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (GWAS sample size).n", 
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

sdat <- fread(opt$file,sep="\t")
names(sdat)[1] <- "chr"

print("Data loaded.")

# Load annotation
fg_vars <- fread(file="/stor/scratch/Kirkpatrick/Jared/signbias/finngen/FG_variants.txt",
                 sep="\t")
names(fg_vars)[8] <- "rsids"

#remove multiallelic vars
sdat <- sdat[ , if (.N == 1) .SD, by = .(chr,pos)]
sdat <- sdat[ , if (.N == 1) .SD, by = rsids]
fg_vars <- fg_vars[ , if (.N == 1) .SD, by = .(chr,pos)]
fg_vars <- fg_vars[ , if (.N == 1) .SD, by = rsids]

#combine with sumstats
sdat <- left_join(sdat, fg_vars, join_by(rsids)) 
rm(fg_vars)
sdat<- sdat %>% select(chr=chr.x, pos=pos.x,ref=ref.x, alt=alt.x, rsids, 
                       nearest_genes, pval, mlogp, beta, sebeta, af_alt,
                       everything(), -ref.y, -`#variant`, -chr.y, -pos.y, -alt.y)

# SNPS only check
valid <- c("A","C","G","T")
sdat <- sdat[
  nchar(ref) == 1 & nchar(alt) == 1 &
    ref %in% valid  & alt %in% valid
]

#edit columns
sdat$pos <- as.numeric(sdat$pos)
sdat$af_alt <- as.numeric(sdat$af_alt)
sdat$beta  <- as.numeric(sdat$beta)
sdat$sebeta <- as.numeric(sdat$sebeta)
sdat$pval <- as.numeric(sdat$pval)
sdat$mlogp <- as.numeric(sdat$mlogp)
sdat$AC <- as.numeric(sdat$AC)
sdat$INFO <- as.numeric(sdat$INFO)

## load variants and match them
vars_matched <- fread(file="/stor/scratch/Kirkpatrick/Jared/signbias/data/aou_matched_rsids.tsv",
                      sep="\t")
sdat<- sdat[rsids %in% vars_matched$rsid]

# hg19 positions
sdat[, c("chr37", "pos37") := tstrsplit(b37_coord, ":", keep = 1:2)]
sdat[, pos37 := as.integer(pos37)]  

#minor allele freqs
sdat$af_ref <- 1-sdat$af_alt

sdat[, `:=`(
  minor_allele = ifelse(af_alt <= 0.5, alt, ref),
  minor_AF         = pmin(af_alt, (1-af_alt))   
)]

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
  filter(!is.na(sebeta)) %>% #remove NA SEs
    
  # polarize effects w.r.t to minor alleles (alt is effect allele in UKB)
  mutate(alt_is_minor = af_alt <= 0.5,
         m_Beta = ifelse(alt_is_minor, beta, -beta),
         MAF = minor_AF)

# assign ancestral states
print("------------------------")
print("Assigning ancestral states...")

adir<-"/stor/work/Kirkpatrick/scratch/Jared/signbias/data/ancestral_seqs/homo_sapiens_ancestor_"
anc <- rep(NA_character_, nrow(sdat))

for (i in 1:22) {
  idx <- which(sdat$chr == i)
  
  if (length(idx) == 0L) next               
  

  fasta_file <- sprintf("%s%d.fa", adir, 1)
  aseq       <- readDNAStringSet(fasta_file)[[1]]
  
  positions  <- as.integer(sdat$pos37[idx]) 
  
  good <- !is.na(positions)                 
  
  if (any(good)) {
    views            <- Views(aseq, start = positions[good], width = 1)
    anc[idx[good]]   <- toupper(as.character(views))
  }
  
  message("chr ", i, ": filled ", sum(good), " / ", length(idx), " rows")
}

sdat[, anc_allele := anc]

rm(aseq)
rm(views)

# Polarize wrt to derived alleles
sdat<- sdat%>%
  mutate(
    alt_is_anc = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'), alt == anc_allele, NA),
    d_beta = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'),
                    ifelse(alt_is_anc, -beta, beta), NA),
    DAF = ifelse(anc_allele %in% c('A', 'T', 'G', 'C'),
                 ifelse(alt_is_anc, 1 - af_alt, af_alt), NA)
)

# Assign SNPs to LD blocks
print("------------------------")
print("Assigning LD blocks...")
setDT(ld)

sdat <- sdat %>%
  left_join(
    ld,
    join_by(chr, between(pos, start, end))   
  ) %>%
  select(-start, -end)  

sdat <- sdat[!is.na(block_id)]

#remove unnecessary columns and QC for info 
print("------------------------")
print("Variant QC....")

sdat <- subset(sdat, INFO >= 0.8)

#QC for MAC (>= 20)
Ns <- as.numeric(opt$sample_size)
sdat$minor_AC <- sdat$MAF * (2*Ns)
sdat<- sdat[minor_AC >= 20]

#arrange columns
sdat <- sdat %>% select(chr, pos, rsids, ref, alt, minor_AC,	MAF, DAF, m_Beta, d_beta, 
                        beta,	se=sebeta, minor_allele, alt_is_minor, anc_allele, pos37, block_id,
                        pval)

#Export formatting
print("------------------------")
print("Done. Exporting...")

fwrite(sdat,
       file      = opt$output,   
       sep       = "\t",           # tab delimiter
       quote     = FALSE,          # don’t quote strings
       showProgress = TRUE)

print("Data exported.")