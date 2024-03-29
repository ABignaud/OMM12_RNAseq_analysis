#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#----------------------------------
#- Load required packages
#- (all of them are available in bioconda/R)
#----------------------------------
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("DESeq2")
library("ashr")
library("rhdf5")
library("tximport")
library("ReportingTools")
library("hwriter")
library("ini")
library("ggrepel")
library("png")
library("ComplexHeatmap")
library("circlize")

#------------
#- Parse args
#------------


# inputdir = args$deseq2$INPUTDIR
metadata <- args[1]
design <- "false"
# condition <- "condition"
# treatment <- "citro"
# control <- "control"
pval <- args[4]
fc <- args[3]
# tx2gene.file = args$deseq2$TX2GENE


#--------------
#- Get metadata
#--------------
colData <- read.table(metadata, sep = ";", header = T)
samples <- colData[, 1]
colnames <- colnames(colData)
colData <- as.data.frame(colData[, -c(1)])
colData$condition <- factor(colData$condition)
colData$replicate <- factor(colData$replicate)
rownames(colData) <- samples
colnames(colData) <- colnames[-1]


#-------------------------------
#- Specify design and contrast
#-------------------------------
default <- "no"
if (design == "false") {
    design <- paste0(colnames(colData)[1:ncol(colData)], collapse = " + ")
    default <- "yes"
} else {
    if (!(condition %in% colnames(colData))) {
        stop("ERROR: condition '", condition, "' is not a column of the metadata file, please specify a valid condition")
    }
    if (!(treatment %in% colData[colnames(colData) == condition][, 1])) {
        stop("ERROR: treatment '", treatment, "' is not a level in condition ", condition, ". Please specify a valid treatment")
    }
    if (!(control %in% colData[colnames(colData) == condition][, 1])) {
        stop("ERROR: control '", control, "' is not a level in condition ", condition, ". Please specify a valid control")
    }
}

#-----------------------------------------------------------------------
#- Get counts and create dds
#- Make sure that first column has the sample names, that these samples
#-  are all in countData, and that they are in the same order
#-----------------------------------------------------------------------
counts <- "featureCounts"
countData <- read.table(paste0(args[2]), sep = "\t", header = T, check.names = FALSE)

if (colnames(countData)[1] == "ENSEMBL_ID") {
    geneID <- countData$ENSEMBL_ID
    countData <- select(countData, -ENSEMBL_ID)
    rownames(countData) <- geneID
} else if (colnames(countData)[1] == "Geneid") {
    geneID <- countData$Geneid
    countData <- select(countData, -Geneid)
    rownames(countData) <- geneID
} else {
    stop(
        "ERROR: unrecognized input file format. Please check that the file has the same format as the merged counts file outputted by the nextflow_rnaseq or nfcore_rnaseq pipelines, i.e., ",
        "\"ENSEMBL_ID\" column followed by sample counts, or \"Geneid\" and \"gene_name\" columns followed by sample counts."
    )
}

if (length(rownames(colData)[!rownames(colData) %in% colnames(countData)]) > 0) {
    stop(
        "ERROR: the following samples are not in the featureCounts matrix: ", paste(rownames(colData)[!rownames(colData) %in% colnames(countData)], collapse = ", "),
        ". Please make sure that the first column of your metadata file has the sample IDs."
    )
} else {
    countData <- countData[, rownames(colData)]
    if (!all(rownames(colData) == colnames(countData))) {
        stop("ERROR: Something is wrong, this should never happen [rownames(colData) ne colnames(countData)??]")
    }
}
dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = eval(parse(text = paste0("~ ", design)))
)

#-------------------------------------------
#- remove genes with <= 60 counts in all samples
#-------------------------------------------
dds <- dds[rowSums(counts(dds)) > 100, ]

#---------
#- Run DGE
#----------
dds <- DESeq(dds)

#---------
#- PCA
#----------
if (nrow(colData) < 30) {
    transformation <- rlog(dds, blind = FALSE)
} else {
    transformation <- vst(dds, blind = FALSE)
}
pcaData <- plotPCA(transformation, intgroup = colnames(colData), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaplot_name <- "PCAplot.png"
png(file = pcaplot_name)
if (ncol(colData) > 1) {
    ggplot(pcaData, aes(PC1, PC2, color = colData[, 1], shape = colData[, 2])) +
        geom_point(size = 3) +
        scale_shape_manual(values = 1:nlevels(colData[, 2])) +
        ggtitle("PCA plot") +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank())
} else {
    ggplot(pcaData, aes(PC1, PC2, color = colData[, 1])) +
        geom_point(size = 3) +
        ggtitle("PCA plot") +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank())
}
dev.off()

#-----------------------------------------------------------------------
#- Modify one of the functions in ReportingTools so that the norm count
#- 	plots per gene are in different folders for each contrast
#-	so that plots for a gene that come up in 2 analyses are not overwritten
#- modified from https://rdrr.io/bioc/ReportingTools/src/R/addReportColumns-methods.R
#-----------------------------------------------------------------------
setMethod("modifyReportDF",
    signature = signature(
        object = "DESeqResults"
    ),
    definition = function(df, htmlRep, object, DataSet, factor,
                          make.plots = TRUE, contrast, ...) {
        if ("EntrezId" %in% colnames(df)) {
            df <- entrezGene.link(df)
        }
        if (make.plots) {
            dots <- list(...)
            par.settings <- list()
            if ("par.settings" %in% names(dots)) {
                par.settings <- dots$par.settings
            }

            figure.dirname <- paste("figures", htmlRep$shortName, "/", contrast, sep = "")
            figure.directory <- file.path(
                dirname(path(htmlRep)),
                figure.dirname
            )
            dir.create(figure.directory, recursive = TRUE)

            df <- ReportingTools:::eSetPlot(df, DataSet, factor, figure.directory,
                figure.dirname,
                par.settings = par.settings,
                ylab.type = "Normalized Counts"
            )
            df
        }
        df
    }
)

if (default == "no") {
    #-------------------------------------------------------------
    #- Get results using  design and contrast defined by the user
    #-------------------------------------------------------------
    contrast <- c(condition, treatment, control)

    #- alpha is fdr threshold for summary display only
    res <- results(dds, contrast = contrast, alpha = 0.05)
    resNorm <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

    #----------------------------
    #- File with all the results
    #----------------------------
    write.table(resNorm, file = paste0(condition, "_", treatment, "_vs_", control, "_results.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

    #-----------------------
    #- Heatmap plot
    #-----------------------

    heatmap_name <- paste0(condition, "_", treatment, "_vs_", control, "_heatmap.png")

    de <- rownames(resNorm[resNorm$padj < 0.05 & !is.na(resNorm$padj) & abs(resNorm$log2FoldChange) > fc, ])

    if (length(de) < 5) {
        message <- paste0("Only ", length(de), " genes with Padj<0.05 and Log2FC>", fc, ".\nPlease rerun the pipeline with a lower FC threshold to generate a heatmap")

        png(file = heatmap_name, width = 600, height = 600, res = 200)
        plot(0, type = "n", axes = FALSE, ann = FALSE)
        mtext(message, side = 3)
        dev.off()
    } else {
        de_mat <- assay(transformation)[de, ]

        png(file = heatmap_name, width = 1500, height = 2000, res = 200)

        ha <- HeatmapAnnotation(df = select(colData, !!condition))
        p <- Heatmap(t(scale(t(de_mat))),
            name = "Normalized counts (scaled)",
            row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 4),
            heatmap_legend_param = list(legend_direction = "horizontal"),
            show_column_names = TRUE,
            top_annotation = ha,
            show_row_names = TRUE,
            cluster_columns = TRUE,
            column_names_side = "top"
        )
        print(p)
        dev.off()
    }

    #-----------------------
    #- MA plot
    #-----------------------
    maplot_name <- paste0(condition, "_", treatment, "_vs_", control, "_MAplot.png")
    png(file = maplot_name)
    maplot <- plotMA(resNorm, ylim = c(-5, 5))
    dev.off()

    #------------------------
    #- Volcano plot
    #------------------------

    vplot_name <- paste0(condition, "_", treatment, "_vs_", control, "_VolcanoPlot.png")

    ha <- data.frame(gene = rownames(resNorm), logFC = resNorm$log2FoldChange, Pval = resNorm$pvalue)
    ha <- within(ha, {
        Col <- "Other"
    })
    ha[abs(ha$logFC) > fc, "Col"] <- paste0("|logFC| > ", fc)
    ha[!is.na(ha$Pval) & ha$Pval < pval, "Col"] <- paste0("Pval < ", pval)
    ha[!is.na(ha$Pval) & ha$Pval < pval & abs(ha$logFC) > fc, "Col"] <- paste0("Pval < ", pval, " & |logFC| > ", fc)
    ha$Col <- factor(ha$Col, levels = c(paste0("Pval < ", pval), paste0("|logFC| > ", fc), paste0("Pval < ", pval, " & |logFC| > ", fc), "Other"), labels = c(paste0("Pval < ", pval), paste0("|logFC| > ", fc), paste0("Pval < ", pval, " & |logFC| > ", fc), "Other"))

    png(file = vplot_name)
    p <- ggplot(ha, aes(logFC, -log10(Pval))) +
        geom_point(aes(col = Col), alpha = 0.5) +
        ggtitle(paste0(condition, " - ", treatment, " vs ", control)) +
        geom_hline(yintercept = -log10(pval), linetype = 2, alpha = 0.5) +
        geom_vline(xintercept = fc, linetype = 2, alpha = 0.5) +
        geom_vline(xintercept = -fc, linetype = 2, alpha = 0.5) +
        scale_color_manual(values = c("indianred1", "gold2", "cornflowerblue", "gray47")) +
        theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank())

    #- add gene labels for gene switch pval anf fc below theresholds
    ha2 <- ha %>%
        filter(Pval < pval & (logFC >= fc | logFC <= -fc)) %>%
        select(gene, logFC, Pval, Col)

    p <- p + geom_text_repel(data = ha2, aes(label = gene), box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), force = 10, size = 3, segment.size = 0.25, segment.alpha = 0.5)
    print(p)
    dev.off()

    #-----------------------
    #- Report
    #-----------------------

    allplot_name <- paste0(condition, "_", treatment, "_vs_", control, "_AllPlot.png")
    img1 <- readPNG(maplot_name)
    img2 <- readPNG(vplot_name)
    png(file = allplot_name, width = 1200, height = 600)
    par(mai = rep(0, 4)) # no margins
    layout(matrix(1:2, ncol = 2, byrow = TRUE))
    for (i in 1:2) {
        plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = "i", yaxs = "i")
        rasterImage(eval(parse(text = paste0("img", i))), 0, 0, 1, 1)
    }
    dev.off()

    reportdir <- "deseq2_report"
    title <- paste0("RNA-seq DGE analysis using DESeq2 with user-defined design and ", counts, " counts")
    des2Report <- HTMLReport(shortName = "DESeq2_nextflow_pipeline_results", title = title, reportDirectory = reportdir)

    himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", pcaplot_name))
    publish(hwrite(himg, br = TRUE, center = T), des2Report)
    publish(paste0("<center><h3><u>Contrast:</u>  ", condition, " - ", treatment, " vs ", control), des2Report)
    publish("<h4>MA/Volcano plots", des2Report)
    himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", allplot_name))
    publish(hwrite(himg, br = TRUE, center = T), des2Report)
    publish("<h4>Top 100 differentially expressed genes", des2Report)
    publish(resNorm, des2Report, contrast = paste0(condition, "_", treatment, "_vs_", control), pvalueCutoff = 1, n = 100, DataSet = dds, factor = colData(dds)[[condition]])
    publish(paste0("<h4>Heatmap DE genes (FDR 0.05, |log2FC| ", fc, ")"), des2Report)
    himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", heatmap_name))
    publish(hwrite(himg, br = TRUE, center = T), des2Report)

    finish(des2Report)
} else {
    #-------------------------------------------------------------
    #- Get results using default design and all possible contrasts
    #-------------------------------------------------------------
    reportdirALL <- "deseq2_report"
    title <- paste0("RNA-seq DGE analysis using DESeq2 with default (no design specified) mode and ", counts, " counts")
    des2ReportALL <- HTMLReport(shortName = "DESeq2_nextflow_pipeline_results", title = title, reportDirectory = reportdirALL)

    himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", pcaplot_name))
    publish(hwrite(himg, br = TRUE, center = TRUE), des2ReportALL)
    n <- 1
    for (i in 1:ncol(colData)) {
        pairs <- combn(unique(colData[, i]), 2)
        for (j in 1:ncol(pairs)) {
            #- alpha is fdr threshold for summary display only
            res <- results(dds, contrast = c(colnames(colData)[i], as.character(pairs[, j])), alpha = 0.05)
            resNorm <- lfcShrink(dds, contrast = c(colnames(colData)[i], as.character(pairs[, j])), res = res, type = "ashr")

            #-----------------------
            #- File with all results
            #-----------------------
            write.table(resNorm, file = paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"), "_results.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

            #-----------------------
            #- Heatmap plot
            #-----------------------

            heatmap_name <- paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"), "_heatmap.png")

            de <- rownames(resNorm[resNorm$padj < 0.05 & !is.na(resNorm$padj) & abs(resNorm$log2FoldChange) > fc, ])

            if (length(de) < 5) {
                message <- paste0("Only ", length(de), " genes with Padj<0.05 and Log2FC>", fc, ". Please rerun the pipeline with a lower FC threshold to generate a heatmap")

                png(file = heatmap_name, width = 600, height = 600, res = 200)
                plot(0, type = "n", axes = FALSE, ann = FALSE)
                mtext(message, side = 3)
                dev.off()
            } else {
                de_mat <- assay(transformation)[de, ]

                png(file = heatmap_name, width = 1500, height = 2000, res = 200)

                ha <- HeatmapAnnotation(df = select(colData, !!colnames(colData)[i]))
                p <- Heatmap(t(scale(t(de_mat))),
                    name = "Normalized counts (scaled)",
                    row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 4),
                    heatmap_legend_param = list(legend_direction = "horizontal"),
                    show_column_names = TRUE,
                    top_annotation = ha,
                    show_row_names = TRUE,
                    cluster_columns = TRUE,
                    column_names_side = "top"
                )
                print(p)
                dev.off()
            }

            #-----------------------
            #- MA plot
            #-----------------------
            maplot_name <- paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"), "_MAplot.png")
            png(file = maplot_name)
            maplot <- plotMA(resNorm, ylim = c(-5, 5))
            dev.off()

            #------------------------
            #- Volcano plot
            #------------------------

            vplot_name <- paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"), "_VolcanoPlot.png")

            ha <- data.frame(gene = rownames(resNorm), logFC = resNorm$log2FoldChange, Pval = resNorm$pvalue)
            ha <- within(ha, {
                Col <- "Other"
            })
            ha[abs(ha$logFC) > fc, "Col"] <- paste0("|logFC| > ", fc)
            ha[!is.na(ha$Pval) & ha$Pval < pval, "Col"] <- paste0("Pval < ", pval)
            ha[!is.na(ha$Pval) & ha$Pval < pval & abs(ha$logFC) > fc, "Col"] <- paste0("Pval < ", pval, " & |logFC| > ", fc)
            ha$Col <- factor(ha$Col, levels = c(paste0("Pval < ", pval), paste0("|logFC| > ", fc), paste0("Pval < ", pval, " & |logFC| > ", fc), "Other"), labels = c(paste0("Pval < ", pval), paste0("|logFC| > ", fc), paste0("Pval < ", pval, " & |logFC| > ", fc), "Other"))

            png(file = vplot_name)
            p <- ggplot(ha, aes(logFC, -log10(Pval))) +
                geom_point(aes(col = Col), alpha = 0.5) +
                ggtitle(paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"))) +
                geom_hline(yintercept = -log10(pval), linetype = 2, alpha = 0.5) +
                geom_vline(xintercept = fc, linetype = 2, alpha = 0.5) +
                geom_vline(xintercept = -fc, linetype = 2, alpha = 0.5) +
                scale_color_manual(values = c("indianred1", "gold2", "cornflowerblue", "gray47")) +
                theme_classic() +
                theme(legend.position = "bottom", legend.title = element_blank())

            #- add gene labels for genes with pval anf fc below theresholds
            ha2 <- ha %>%
                filter(Pval < pval & (logFC >= fc | logFC <= -fc)) %>%
                select(gene, logFC, Pval, Col)

            p <- p + geom_text_repel(data = ha2, aes(label = gene), box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), force = 10, size = 3, segment.size = 0.25, segment.alpha = 0.5)
            print(p)
            dev.off()

            #-----------------------
            #- Report
            #-----------------------

            allplot_name <- paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_"), "_AllPlot.png")
            img1 <- readPNG(maplot_name)
            img2 <- readPNG(vplot_name)
            png(file = allplot_name, width = 1200, height = 600)
            par(mai = rep(0, 4)) # no margins
            layout(matrix(1:2, ncol = 2, byrow = TRUE))
            for (l in 1:2) {
                plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = "i", yaxs = "i")
                rasterImage(eval(parse(text = paste0("img", l))), 0, 0, 1, 1)
            }
            dev.off()

            publish(paste0("<center><h3><u>Contrast ", n, ":</u>  ", colnames(colData)[i], " - ", paste0(as.character(pairs[, j]), collapse = " vs ")), des2ReportALL)
            publish("<h4>MA/Volcano plots", des2ReportALL)
            himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", allplot_name))
            publish(hwrite(himg, br = TRUE, center = TRUE), des2ReportALL)
            publish("<h4>Top 100 differentially expressed genes", des2ReportALL)
            publish(resNorm, des2ReportALL, contrast = paste0(colnames(colData)[i], "_", paste0(as.character(pairs[, j]), collapse = "_vs_")), pvalueCutoff = 1, n = 100, DataSet = dds, factor = colData(dds)[[i]])
            publish(paste0("<h4>Heatmap DE genes (FDR 0.05, |log2FC| ", fc, ")"), des2ReportALL)
            himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/", heatmap_name))
            publish(hwrite(himg, br = TRUE, center = TRUE), des2ReportALL)
            n <- n + 1
        }
    }
    finish(des2ReportALL)
}
