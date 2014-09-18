library(ggplot2)
library(reshape2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")
gp = ggplot(df, aes(x=Position, y=Nucleotide)) 
gp + geom_tile(aes(fill=Weight)) + scale_fill_gradient2("Weight", low="#f1a340", mid="#f7f7f7", high="#998ec3", midpoint=0.5) + theme_bw()
ggsave(output.pdf, width=10,height=2)