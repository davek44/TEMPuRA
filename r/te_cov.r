library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
te = ca[2]
output.pdf = ca[3]

df <- read.table(df.file, header=T)

# keep data in the order given (otherwise it's sorts them)
df$data = factor(df$data, levels=unique(df$data))

gp <- ggplot(df, aes(x=indexes, y=coverage, color=data)) +
    geom_line(size=1, alpha=0.75) +
    scale_x_continuous(paste(te,"index")) +
    scale_y_continuous("") +
    theme_bw() +
    theme(text=element_text(size=20)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

color.count = length(unique(df$data))
if(color.count == 2) {
    gp <- gp + scale_color_manual("", values=c("tomato3","steelblue"))
} else if(color.count == 4) {
    gp <- gp + scale_color_manual("", values=c("#F33D6E", "#FF7340", "#3F8FD2", "#4D54D8"))
} else if(color.count == 6) {
    gp <- gp + scale_color_manual("", values=c("#FF4940", "#FF8040", "#F13C73", "#4867D6", "#7746D7", "#35C0CD"))
} else {
    gp <- gp + scale_color_brewer("", palette="Set1")
}

ggsave(output.pdf)
