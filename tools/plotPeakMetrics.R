library(readr)
library(tidyverse)
library(ggplot2)

jaccard_file = snakemake@input[["jaccard"]]
intervalFiles = snakemake@input[["interval"]]
jaccardJPG = snakemake@output[["jaccard"]]
intervalJPG = snakemake@output[["interval"]]


jaccard = read_tsv(jaccard_file)
jaccard = mutate(jaccard, Group = if_else(Comparison=="betweenReps", "repA.vs.repB", if_else(Comparison=="IDRprs", "pr1.vs.pr2", "Peaks.vs.PR")))
jaccard_plot = ggplot(jaccard, aes(Group, Tissue, fill = jaccard)) + geom_tile() + scale_fill_gradient(low = "white", high = "red")


alpha = 0.05
ratios=setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Tissue", "Comparison", "Ratio"))
for (i in intervalFiles)
{
  match = str_match(basename(i), "(.+?)_(.+?).txt")
  tissue = match[2]
  comp = match[3]
  file = read_tsv(i, col_names = FALSE)
  ratio = nrow(filter(file, X7 < alpha)) / nrow(file)
  ratios = add_row(ratios, Tissue = tissue, Comparison = comp, Ratio = ratio)
}
ratios = mutate(ratios, Group = if_else(Comparison=="betweenReps", "repA.vs.repB", if_else(Comparison=="IDRprs", "pr1.vs.pr2", "Peaks.vs.PR")))
interval_plot = ggplot(ratios, aes(Group, Tissue, fill = Ratio)) + geom_tile() + scale_fill_gradient(low = "white", high = "red")

ggsave(jaccardJPG, plot = jaccard_plot, dpi = 300)
ggsave(intervalJPG, plot = interval_plot, dpi = 300)
