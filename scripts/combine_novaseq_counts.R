#!/bin/Rscript

require(dplyr)
require(tibble)

count.files <- list.files('results/novaseq_counts')

count.df <- lapply(count.files, \(.f){
    if(dim(read.table(file.path('results/novaseq_counts', .f)))[1]>1){
        read.table(file.path('results/novaseq_counts', .f), skip=1) |> 
            tibble::column_to_rownames('V1') |> t() |>
            data.frame() |>
            dplyr::mutate(Plate=paste0(stringr::str_split(.f, '_')[[1]][1:2], collapse='_'),
                          Well=stringr::str_split(.f, '_')[[1]][3],
                          Cell=paste0(Plate, '_', Well)
                          ) |>
            dplyr::relocate(Plate, Well) |>
            tibble::rownames_to_column() |> dplyr::select(-c('rowname')) |>
            tibble::column_to_rownames('Cell')
    }else{
            data.frame(Plate=paste0(stringr::str_split(.f, '_')[[1]][1:2], collapse='_'),
                       Well=stringr::str_split(.f, '_')[[1]][3]) |>
            dplyr::mutate(
                       Cell=paste0(Plate, '_', Well)
                    ) |>
            tibble::column_to_rownames('Cell')
    }
}) |> dplyr::bind_rows()

rep.list <- as.list(rep(0, dim(count.df)[2]))
names(rep.list) <- names(count.df)
count.df <- count.df |> tidyr::replace_na(rep.list)
write.csv(count.df, 'results/combined_novaseq_counts.csv')
