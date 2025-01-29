library(magrittr)
require(dplyr)
require(matrixStats)
require(readr)
meta <- read.csv('data/novaseq_metadata/tagmentation_layout.csv') |> dplyr::select(-c('X'))

novaseq <- read.csv('results/combined_novaseq_counts.csv.gz') |> dplyr::select(-c('X'))

novaseq.joined <- dplyr::left_join(novaseq, meta, by=c('Plate'='tag_plate', 'Well'='location')) |> 
  dplyr::relocate(c('value', 'PID', 'timepoint', 'PID_timepoint'), .after='Well')

ont.files <- list.files('results/mutations-called-v2', full.names = T)
ont.counts <- lapply(ont.files, \(.f) {
  read.csv(.f) |> dplyr::mutate(file=.f) |> dplyr::relocate(file)
}) |> dplyr::bind_rows() |>
  dplyr::mutate(dataset_name=paste0("HSPC_pool", 
                            ifelse(grepl('05', file), '1',
                            ifelse(grepl('06', file), '2',
                            ifelse(grepl('07', file), '3_repeat',
                            ifelse(grepl('08', file), '4', NA)))))) |>
  dplyr::relocate(dataset_name) |>
  dplyr::rename(barcode=X)

hsc.obs <- read.csv('data/HSC_obs.csv.gz') |>
  dplyr::select(-c('X')) |> dplyr::mutate(barcode=gsub('-1', '', barcode))

mutations <- gsub('_(mt|wt|un|coverage)$', '', names(ont.counts)[grepl('^chr', names(ont.counts))]) |> unique()

# Get total coverage, and remove entries with < 10 UMI
umi.thr <- 10
coverage.matrix <- as.matrix(ont.counts[,paste0(mutations,'_coverage')])
coverage.matrix[coverage.matrix < umi.thr] <- NA

mt.matrix <- as.matrix(ont.counts[,paste0(mutations,'_mt')])
wt.matrix <- as.matrix(ont.counts[,paste0(mutations,'_wt')])

# Calcualte VAF and replace infs with NA
vaf.matrix <- mt.matrix/coverage.matrix
vaf.matrix[!is.finite(vaf.matrix)] <- NA

wt.vaf.matrix <- wt.matrix/coverage.matrix
wt.vaf.matrix[!is.finite(vaf.matrix)] <- NA

# Needed to convert back to a human redable mutation
moi.df <- read.csv('data/mutations-of-interest.csv')
moi <- gsub('_.*', '', moi.df$mut2_hg19)
names(moi) <- moi.df$mut2__hg38
mut.str <- sub('\\.', '>', sub('\\.', ':', mutations))

moi.full <- moi.df$mut2_hg19
names(moi.full) <- moi.df$mut2__hg38

mt.calls <- apply(vaf.matrix >= 0.2, 1, \(.r) {
  .rx <- .r
  .rx[is.na(.r)] <- F
  if(sum(.rx>0)){
    paste0(paste0('MT_', moi[mut.str][.r[!is.na(.r)]] |> unique()), collapse='; ')
  }else{
    paste0()
  }
  })

wt.calls <- apply(wt.vaf.matrix>=0.8, 1, \(.r) {
  .rx <- .r
  .rx[is.na(.r)] <- F
  if(sum(.rx>0)){
    paste0(paste0('WT_', moi[mut.str][.r[!is.na(.r)]] |> unique()), collapse='; ')
  }else{
    paste0()
  }
  })


un.calls <- apply((wt.matrix<0.8)&(vaf.matrix<0.2), 1, \(.r) {
  .rx <- .r
  .rx[is.na(.r)] <- F
  if(sum(.rx>0)){
    paste0(paste0('UN_', moi[mut.str][.r[!is.na(.r)]]), collapse='; ')
  }else{
    paste0()
  }
  })

call.str <- sapply(1:length(mt.calls), \(.i){
  if(length(wt.calls[[.i]])){
    paste(wt.calls[[.i]], mt.calls[[.i]], un.calls[[.i]], sep='; ')
  }else if (length(mt.calls[[.i]])){
    paste(mt.calls[[.i]], un.calls[[.i]], sep='; ')
  }else if (length(un.calls[[.i]])){
    paste(mt.calls[[.i]], un.calls[[.i]], sep='; ')
  }else{
    paste0()
  }
})

call.str.df <- ont.counts |> dplyr::select('dataset_name', 'barcode') |>
  dplyr::mutate(MT_STATUS=call.str)

colnames(vaf.matrix) <- moi.full[mut.str]
mt.vaf.df <- data.frame(vaf.matrix) |> 
  dplyr::bind_cols(ont.counts |> dplyr::select('dataset_name', 'barcode')) |>
  dplyr::relocate(c('dataset_name', 'barcode'))

colnames(wt.vaf.matrix) <- moi.full[mut.str]
wt.vaf.df <- data.frame(wt.vaf.matrix) |> 
  dplyr::bind_cols(ont.counts |> dplyr::select('dataset_name', 'barcode')) |>
  dplyr::relocate(c('dataset_name', 'barcode'))

hsc.obs.plus <- dplyr::left_join(hsc.obs, call.str.df, by=c('dataset_name', 'barcode'))

# readr::write_csv(hsc.obs.plus, 'results/HSC_obs_with-mt-status.csv.gz')
# readr::write_csv(mt.vaf.df, 'results/mt-vafs.csv.gz')
# readr::write_csv(wt.vaf.df, 'results/wt-vafs.csv.gz')

srsf2.names <- names(mt.vaf.df)[grep('SRSF2', names(mt.vaf.df))]
# SRSF2 check
srsf2.mt <- mt.vaf.df |> dplyr::select(all_of(c('dataset_name', 'barcode', srsf2.names))) |>
  dplyr::left_join(hsc.obs, by=c('dataset_name', 'barcode'))

srsf2.by.pat <- srsf2.mt |> dplyr::group_by(dataset_name, timepoint, patient) |>
  dplyr::summarize_at(vars(all_of(srsf2.names)), ~sum(.x>0.2, na.rm=T))

ont.counts |> dplyr::filter(barcode %in% (hsc.obs |> dplyr::filter(patient=='61292002', dataset_name=='HSPC_pool2'))$barcode, dataset_name=='HSPC_pool2') |> View()
