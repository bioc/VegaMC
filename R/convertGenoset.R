convertGenoset <- function(genosetElem, isbaf)
{


    if (isbaf == TRUE)
        {
            lrr <- assayData(genosetElem) $lrr
            baf <- assayData(genosetElem) $baf
    } else
        {
            lrr <- assayData(genosetElem) $cn
        }

    samples <- sampleNames(genosetElem)

    nsamples <- length(samples)
    chr_ind <- chrIndices(genosetElem)
    nprobes <- chr_ind[nrow(chr_ind), 2]


    if (isbaf == TRUE)
        {
            numcols <- nsamples * 2 + 3
    } else
        {
            numcols <- nsamples + 3
        }
    dataset_matrix <- matrix(0, nprobes, numcols)

    name_cols <- c()
    sample_index <- 1
    for (i in 4:numcols)
        {
            if (isbaf == TRUE)
                {
                    if (i %% 2 == 1)
                        {
                            name_cols <- c(name_cols,
                                           paste(samples[sample_index],
                                                 ".B Allele Freq", sep = ""))
                            dataset_matrix[,
                                           i] <- as.numeric(baf[,
                                                                sample_index])
                            sample_index <- sample_index + 1
                    } else
                        {
                            name_cols <- c(name_cols,
                                           paste(samples[sample_index],
                                                 ".Log R Ratio", sep = ""))
                            dataset_matrix[,
                                           i] <- as.numeric(lrr[,
                                                                sample_index])
                        }
            } else
                {
                    name_cols <- c(name_cols,
                                   paste(samples[sample_index],
                                         ".Log R Ratio", sep = ""))
                    dataset_matrix[, i] <- as.numeric(lrr[, sample_index])
                    sample_index <- sample_index + 1
                }

        }
    colnames(dataset_matrix) <- c("Name", "Chr", "Position", name_cols)

    dataset_matrix[, 1] <- 1:nprobes

    for (i in 1:nrow(chr_ind))
        {
            chr <- rownames(chr_ind)[i]
            chr <- substr(chr, 4, nchar(chr))
            start <- as.numeric(chr_ind[i, 1])
            end <- as.numeric(chr_ind[i, 2])

          dataset_matrix[start:end, 2] <- chr
        }

    dataset_matrix[, 3] <- pos(genosetElem)
    dataset <- paste(tempdir(), "dataset.txt", sep = "")
    write.table(dataset_matrix, dataset, sep = "\t", col.names =
                TRUE, row.names = FALSE, eol = "\n", quote = FALSE)

    return(dataset)

}
