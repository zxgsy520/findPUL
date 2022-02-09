#!/usr/bin/env Rscript

library(grid)
library(circlize)
library(ComplexHeatmap)
options(bitmapType='cairo')

add_help_args <- function(args){

    if(length(args) != 2) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot a cluster heatmap.\n")
        cat("Example:cluster_heatmap.R cazy_otu.tsv out.heatmap\n")
        cat("Input file format:\n")
        cat("\tsample1\tsample2\tsample3\n")
        cat("OTU1\t0.276046\t-1.339209\t1.063164")
        quit()
    }

}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

file <- args[1]
prefix <- args[2]

data1 <- read.table(file, header=T) #row.names=1, stringsAsFactors=F #读取数据
data <- as.matrix(data1)	#将数据转化为矩阵
rows <- nrow(data)	#计算行数
columns <- length(data[1,])	#计算列数

if(rows<=20) {
    heights <- 6
}else if(rows<=50) {
    heights <- 14
}else {
    heights <- 18
}

if(columns<=4) {
    widths <- 4
}else if(columns<=8) {
    widths <- 8
}else {
    widths <- 10
}

pdf(paste(prefix, ".pdf", sep=""), width=widths, height=heights)
a <- dev.cur()
png(paste(prefix, ".png", sep=""), width=200, height=600)

format <- colorRamp2(seq(min(data), max(data), length=3), c("blue", "#EEEEEE", "red"), space="RGB")
dev.control("enable")
Heatmap(data, name=" ", col=format, show_column_dend=FALSE)
dev.copy(which=a)
dev.off()
dev.off()
