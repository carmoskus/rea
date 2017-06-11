# Main interface for running rea

doc = 'Usage: rea.R <counts> <cols> [ -exp -p <np> ]

options:
 -p <np>  Number of principle components
 -exp     Calculate lambda on expressed genes'

suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

#print(opts)

# Set R options
options(digits=4)

# Set constants
quantiles  = c(0.5, 0.25, 0.1, 0.05, 0.01)
df = 1

# Set parameters
counts.file = opts$counts
cols.file = opts$cols

# Load data
counts = as.matrix(read.table(counts.file, header=TRUE, sep="\t", row.names=1))
col.info = read.table(cols.file, header=TRUE, row.names=1, sep="\t")

# Start DESeq2
suppressPackageStartupMessages(library("DESeq2"))

analyze = function (np) {
    print(paste0("Analyze ", np))
    
    se = SummarizedExperiment(counts)
    coldata = colData(se)
    coldata$Group = col.info$group

    # Add covariates
    form = "Group"
    if (np > 0) {
        for (i in 1:arg2) {
            coldata[,paste0("Covar", i)] = sva[,i]
            form = paste0(form, " + Covar", i)
        }
    }
    
    print(form)
    
    dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = reformulate(form))
    dds = DESeq(dds, minReplicatesForReplace=Inf)
    
    print(resultsNames(dds))
    
    # Output results
    res = results(dds, c('Group', 'a', 'b'), cooksCutoff=Inf)
    res = res[order(res$pvalue),]
    #write.csv(as.data.frame(res), file=paste0(subdir, name, "_res.csv"))

    print(res)
    cat("\n")

    all.padj = p.adjust(res$pvalue, method="BH")
    cat(paste0("0.05\t", sum(all.padj <= 0.05, na.rm=TRUE), "\t", sum(res$padj <= 0.05, na.rm=TRUE), "\n"))

    # Compute inflation factors
    exp.p = quantile(res$pvalue[!is.na(res$padj)], quantiles, na.rm=TRUE)
    exp.chi = qchisq(1-exp.p, df)
    exp.lambda = exp.chi / qchisq(1-quantiles, df)

    all.p = quantile(res$pvalue, quantiles, na.rm=TRUE)
    all.chi = qchisq(1-all.p, df)
    all.lambda = all.chi / qchisq(1-quantiles, df)

    cat(paste(quantiles, collapse="\t"), "\n")
    cat(paste(exp.lambda, collapse="\t"), "\n")
    cat(paste(all.lambda, collapse="\t"), "\n")
    
    # Output metadata
    #write.table(sizeFactors(dds), file=paste0(subdir, name, "_sizes.txt"), sep="\t")

    #mc = as.data.frame(mcols(dds))
    #rownames(mc) = rownames(dds)
    #write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")
}

if (is.null(opts$p)) {
    analyze(0)
} else {
    analyze(as.integer(opts$p))
}
