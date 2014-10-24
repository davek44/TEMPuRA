library(ggplot2)
library(reshape2)
library(infotheo)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

# Read in the data and split into appropriate data frames
df = read.table(df.file, header=T, quote="\"",na.strings = ".")
drop = c('Score')
df_sequences = df[,!(names(df) %in% drop)]
df_scores = df$Score

# Discretize Scores
df_discrete_scores = discretize(df_scores)

# Calculate Mutual Information for each column
sequence_length = ncol(df_sequences)
mutual_information <- matrix(NA,sequence_length,1)
i <- 1
while (i<=sequence_length){
	mutual_information[i,] <- 1000*(mutinformation(df_sequences[,i], df_discrete_scores))
	i <- i + 1
}

# Plot the Mutual Information as a bar plot
mutual_information <- data.frame(mutual_information)
mutual_information$position <- rownames(mutual_information)
mutual_information$position <- as.character(mutual_information$position)
mutual_information$position <- factor(mutual_information$position, levels=unique(mutual_information$position))
gp <- ggplot(mutual_information, aes(x=position, y=mutual_information)) + 
	geom_bar(stat = "identity", fill="red") +
	theme_bw() +
	xlab('Position') +
	ylab('Millibits') 

ggsave(output.pdf)