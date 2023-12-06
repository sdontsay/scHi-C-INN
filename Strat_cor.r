library(hicrep)
library(tseries)
library(RJSONIO)
library(stringr)
library(parallel)
library(combinat)
library(optparse)

# Define command-line options
option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, help="Input directory"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory"),
  make_option(c("-s", "--stage"), type="character", default=NULL, help="Stage parameter"),
  make_option(c("-r", "--resol"), type="integer", default=1000000, help="Resolution for get.scc"),
  make_option(c("-H", "--h"), type="integer", default=1, help="h for get.scc"),
  make_option(c("-l", "--lbr"), type="integer", default=2000000, help="Lower bound for get.scc"),
  make_option(c("-u", "--ubr"), type="integer", default=5000000, help="Upper bound for get.scc"),
  make_option(c("-m", "--mcore"), type="integer", default=30, help="Number of cores"),
  make_option(c("-g", "--genome"), type="character", default="hg19", help="Genome type (hg19, hg38, mm9, mm10)")
)

# Parse options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Function to get chromosomes based on genome type
get_chromosomes <- function(genome_type) {
  if (genome_type %in% c('hg19', 'hg38')) {
    return(paste("chr", c(1:22, "X"), sep=""))
  } else if (genome_type %in% c('mm9', 'mm10')) {
    return(paste("chr", c(1:19, "X", "Y"), sep=""))
  } else {
    stop("Unsupported genome type. Choose from 'hg19', 'hg38', 'mm9', 'mm10'.")
  }
}

# Use the parsed options
subdir <- opt$input_dir
outdir <- opt$output_dir
stage = opt$stage
resol = opt$resol
h = opt$h
lbr = opt$lbr
ubr = opt$ubr
mcore = opt$mcore
genome_type = opt$genome

if (!file.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Get the chromosomes based on genome type
chrom <- get_chromosomes(genome_type)

get_cells = function(x) {
  as.integer(str_split(str_split(x, '_')[[1]][2], 'cell')[[1]][2])
}

make_comb = function(x, y) {
  return(paste(x,'-',y, sep=''))
}

make_path = function(x) {
  return(paste(subdir, i, sprintf('/%s_cell%s_%s.txt', i, x, stage), sep='/'))
}

# In the cal_hic function, use the parameters from command-line
cal_hic = function(x, y) {
  mt1 <- ldf[[x]]
  mt2 <- ldf[[y]]
  scc.out = get.scc(mt1, mt2, resol = resol, h = h, lbr = lbr, ubr = ubr)
  return(scc.out$scc)
}

whole_cor <- data.frame()
for (i in chrom) {
  cellnames <- list.files(paste(subdir,i,sep='/'), pattern="*.txt", full.names=FALSE)
  cellnames = as.array(cellnames)
  cell_list = apply(cellnames, 1, get_cells) # get all cells in a directory
  cell_list = sort(cell_list)
  cell_path <- sapply(cell_list, make_path)
  ldf <- mclapply(cell_path, read.matrix, mc.cores = mcore)
  res <- combn(cell_list, 2)
  comb = mcmapply(make_comb, res[1,], res[2,], mc.cores = mcore)
  hic = mcmapply(cal_hic, res[1,], res[2,], mc.cores = mcore)
  write.table(hic, paste(outdir, '/', i, '.txt', sep=''),
              sep=',', col.names = 'weighted_pearson', row.names = comb)
  whole_cor = append(whole_cor,mean(hic, na.rm=TRUE))
}
total_avg <- mean(unlist(whole_cor, recursive = TRUE), na.rm=TRUE)
whole_cor = append(whole_cor,total_avg)
whole_cor = unlist(whole_cor, recursive = TRUE)
write.table(whole_cor, paste(outdir, 'sum_avg_whole.txt', sep='/'), sep=',',
            col.names = 'weighted_pearson')