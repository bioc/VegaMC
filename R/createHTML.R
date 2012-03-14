createHTML <-
function(output_file_name, segmentation, pval_threshold, loss_threshold,
        gain_threshold, loh_percentage, loh_threshold, baf, n_samples,
        n_probes, chromosomes, beta, min_region_bp_size,
        bs, getGenes){
## This function generates the html file reporting the list
## (with the respective details) of each detected region

    ## check if BAFs are observed
    if(baf==TRUE){
        print_baf <- "YES"
    }else{
        print_baf <- "NO"
    }

    n_reg <- nrow(segmentation)
    loss_pos <- which(as.numeric(segmentation[,6])<=pval_threshold)
    n_loss <- length(loss_pos)
    gain_pos <- which(as.numeric(segmentation[,7])<=pval_threshold)
    n_gain <- length(gain_pos)
    n_loh = "No Available"
    if(baf==TRUE){
        loh_pos <- which(as.numeric(segmentation[,8])<=pval_threshold)
        n_loh <- length(loh_pos)
    }

    ## Create the Header
    inFile <- system.file("template/HeaderRegions.html", package="VegaMC")
    inFile <- file(inFile, "r")
    outFile <- file(output_file_name, "w")
    html_genes=""
    if(getGenes == TRUE){
        html_genes <- paste(substr(output_file_name, 1,
                nchar(output_file_name)-5), "Genes.html", sep="")
    }
    values <- c(TITLE=output_file_name, DATE=date(), GENEFILE=html_genes)
    copySubstitute(inFile, outFile, values)
    close(inFile)


    ## Create the Summary
    inFile <- system.file("template/InputParameters.html",
                package="VegaMC")
    inFile <- file(inFile, "r")
    minBP = "No"
    if(is.na(min_region_bp_size)==FALSE){
            minBP=paste(format(min_region_bp_size,
                        big.mark="."), "bp")
    }
    values <- c(NUMSAMPLE=format(n_samples, big.mark="."))
    values <- c(values, NUMCHR=length(chromosomes))
    values <- c(values, NUMPROBE=format(n_probes, big.mark="."),
                        BAF=print_baf)
    values <- c(values, BETA=format(beta, decimal.mark=","))
    values <- c(values, MINBP=minBP)
    values <- c(values, LOSST=format(loss_threshold, decimal.mark=","))
    values <- c(values, GAINT=format(gain_threshold, decimal.mark=","))
    values <- c(values, LOHT=format(loh_threshold, decimal.mark=","))
    values <- c(values, LOHP=paste(format(as.integer(loh_percentage*100),
                        decimal.mark=","),"%"))
    values <- c(values, BS=format(bs, big.mark="."))
    values <- c(values, PVAL=format(pval_threshold, decimal.mark=","))
    values <- c(values, NUMREG=format(n_reg, big.mark="."))
    values <- c(values, NUMDEL=format(n_loss, big.mark="."))
    values <- c(values, NUMAMP=format(n_gain, big.mark="."))
    if(baf==TRUE){
        values <- c(values, NUMLOH=format(n_loh, big.mark="."))
    }else{
        values <- c(values, NUMLOH=n_loh)
    }
    copySubstitute(inFile, outFile, values)
    close(inFile)


    ## Create the Header for the Table of all Regions
    inFile <- system.file("template/HeaderTableAllRegions.html",
                package="VegaMC")
    inFile <- file(inFile, "r")
    values <- c(values, VAL="NO")
    copySubstitute(inFile, outFile, values)
    close(inFile)

    ## Insert a new Line for each region
    segmentation[,2:4] <- format(as.integer(segmentation[,2:4]),
                    big.mark=".")
    segmentation[,5:8] <- format(as.numeric(segmentation[,5:8]),
                    decimal.mark=",")
    for(i in 1:n_reg){
        inFile <- system.file("template/AllRegionTableRow.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c()
        values <- c(CHR=as.character(segmentation[i,1]))
        values <- c(values, BPSTART=as.character(segmentation[i,2]))
        values <- c(values, BPEND=as.character(segmentation[i,3]))
        values <- c(values, BPSIZE=as.character(segmentation[i,4]))   
        values <- c(values, L2MEAN=as.character(segmentation[i,5]))   
        values <- c(values, LOSSPVAL=as.character(segmentation[i,6]))
        values <- c(values, LOSSPRC=as.character(segmentation[i,9]))
        values <- c(values, GAINPVAL=as.character(segmentation[i,7]))
        values <- c(values, GAINPRC=as.character(segmentation[i,10]))
        values <- c(values, LOHPVAL=as.character(segmentation[i,8]))   
        values <- c(values, LOHPRC=as.character(segmentation[i,11]))
        copySubstitute(inFile, outFile, values)
        close(inFile)
    }

    ## Create the Footer for the Table of all Regions
    inFile <- system.file("template/FooterTable.html", package="VegaMC")
    inFile <- file(inFile, "r")
    values <- c(values, VAL="NO")
    copySubstitute(inFile, outFile, values)
    close(inFile)

    if(length(loss_pos)>0){
        ## Create the Header for the Table of Deleted Regions
        inFile <- system.file("template/HeaderTableAberrantRegions.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, TABNAME="Deleted", ABERR="Loss")
        copySubstitute(inFile, outFile, values)
        close(inFile)
        ## Insert a new Line for each region
        for(i in loss_pos){
            inFile <- system.file(
                    "template/AberrantRegionTableRow.html",
                    package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c()
            values <- c(CHR=as.character(segmentation[i,1]))
            values <- c(values,
                    BPSTART=as.character(segmentation[i,2]))
            values <- c(values,
                    BPEND=as.character(segmentation[i,3]))
            values <- c(values,
                    BPSIZE=as.character(segmentation[i,4]))   
            values <- c(values,
                    L2MEAN=as.character(segmentation[i,5]))   
            values <- c(values,
                    ABERRPVAL=as.character(segmentation[i,6]))
            values <- c(values,
                    ABERRPRC=as.character(segmentation[i,9]))
            copySubstitute(inFile, outFile, values)
            close(inFile)
        }
        ## Create the Footer for the Table of Deleted Regions
        inFile <- system.file("template/FooterTable.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, VAL="NO")
        copySubstitute(inFile, outFile, values)
        close(inFile)
    }

    if(length(gain_pos)>0){
        ## Create the Header for the Table of Amplified Regions
        inFile <- system.file("template/HeaderTableAberrantRegions.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, TABNAME="Amplified", ABERR="Gain")
        copySubstitute(inFile, outFile, values)
        close(inFile)
        ## Insert a new Line for each region
        for(i in gain_pos){
            inFile <- system.file(
                    "template/AberrantRegionTableRow.html",
                    package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c()
            values <- c(CHR=as.character(segmentation[i,1]))
            values <- c(values,
                    BPSTART=as.character(segmentation[i,2]))
            values <- c(values,
                    BPEND=as.character(segmentation[i,3]))
            values <- c(values,
                    BPSIZE=as.character(segmentation[i,4]))   
            values <- c(values,
                    L2MEAN=as.character(segmentation[i,5]))   
            values <- c(values,
                    ABERRPVAL=as.character(segmentation[i,7]))
            values <- c(values,
                    ABERRPRC=as.character(segmentation[i,10]))
            copySubstitute(inFile, outFile, values)
            close(inFile)
        }
        ## Create the Footer for the Table of Amplified Regions
        inFile <- system.file("template/FooterTable.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, VAL="NO")
        copySubstitute(inFile, outFile, values)
        close(inFile)
    }

    if(baf==TRUE){
        if(length(loh_pos)>0){
            inFile <- system.file(
                    "template/HeaderTableAberrantRegions.html",
                    package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c(values, TABNAME="LOH", ABERR="LOH")
            copySubstitute(inFile, outFile, values)
            close(inFile)
            for(i in loh_pos){
                inFile <- system.file(
                    "template/AberrantRegionTableRow.html",
                    package="VegaMC")
                inFile <- file(inFile, "r")
                values <- c()
                values <- c(CHR=as.character(segmentation[i,1]))
                values <- c(values,
                    BPSTART=as.character(segmentation[i,2]))
                values <- c(values,
                    BPEND=as.character(segmentation[i,3]))
                values <- c(values,
                    BPSIZE=as.character(segmentation[i,4]))   
                values <- c(values,
                    L2MEAN=as.character(segmentation[i,5]))   
                values <- c(values,
                    ABERRPVAL=as.character(segmentation[i,7]))
                values <- c(values,
                    ABERRPRC=as.character(segmentation[i,10]))
                copySubstitute(inFile, outFile, values)
                close(inFile)
            }
            ## Create the Footer for the Table of Amplified Regions
            inFile <- system.file("template/FooterTable.html",
                        package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c(values, VAL="NO")
            copySubstitute(inFile, outFile, values)
            close(inFile)
        }
    }
    ## Create the Footer
    inFile <- system.file("template/Footer.html", package="VegaMC")
    inFile <- file(inFile, "r")
    values <- c(BIOINFO="http://bioinformatics.biogem.it")
    copySubstitute(inFile, outFile, values)
    close(inFile)
    close(outFile)
}

