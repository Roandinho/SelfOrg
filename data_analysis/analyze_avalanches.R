library(ggplot2)
library(dplyr)
library(poweRlaw)
theme_set(theme_bw(14))


## Function to extract parameters from the file name

extract_parms <- function(file){
  parms <- strsplit(file, split = "_")[[1]][-1]
  
  perturb <- strsplit(parms[grep("P", parms)], "-")[[1]][-1]
  perturb.bool <- as.logical(perturb[1])
  perturb.n <- as.numeric(perturb[2])
  
  hierachy <- strsplit(parms[grep("H", parms)], "-")[[1]][-1]
  hierachy.bool <- as.logical(hierachy[1])
  hierachy.n <- as.numeric(hierachy[2])
  
  neighbour.level <- as.numeric(strsplit(parms[grep("N", parms)], "-")[[1]][-1])
  
  update.interaction <- as.logical(strsplit(parms[grep("I", parms)], "-")[[1]][-1])
  
  
  steps <- as.numeric(
    strsplit(parms[grep("S",parms)], "-")[[1]][-1]
  )
  
  
  M <- as.numeric(sub(".csv", "", strsplit(parms[grep("M",parms)], "-")[[1]][-1]))
  return(list(perturb.bool, perturb.n, hierachy.bool, hierachy.n,
              neighbour.level, update.interaction, steps, M))
}



## Read all the files in the entropy folder

txt_files_entropy <- list.files("./data/entropy/")
txt_files_entropy <- paste("./data/entropy/", txt_files_entropy, sep = "")


data_list <- lapply(txt_files_entropy, function(file){
  print(file)
  # file <- txt_files_entropy[1]
  parms.ls <- extract_parms(file)
  data.raw <- read.csv(file, sep = ",", header = FALSE) 
  
  data <- as.data.frame(t(data.raw)) %>% rename(entropy = V1) %>% 
    mutate(time = seq(1, length(entropy))) 
  
  data$perturb <- parms.ls[[1]]
  data$perturb_n <- parms.ls[[2]]
  data$hierarchy <- parms.ls[[3]]
  data$hierarchy_n <- parms.ls[[4]]
  data$neighbours_level <- parms.ls[[5]]
  data$update <- parms.ls[[6]]
  data$steps <- parms.ls[[7]]
  data$M <- parms.ls[[8]]
  return(data)
})

# this dataframe contains all the entropy values for each file in the entropy folder
entropy <- do.call(rbind, data_list)



## Read all the files in the avalanches folder
txt_files_avalanches <- list.files("./data/avalanches/")
txt_files_avalanches <- paste("./data/avalanches/", txt_files_avalanches, sep = "")


data_list <- lapply(txt_files_avalanches, function(file){
  print(file)
  # file <- txt_files_entropy[2]
  #   file <- "./data/entropy/output_entropy_PF0_HT50_N10_T5000.csv"
  parms.ls <- extract_parms(file)
  data.raw <- read.csv(file, sep = ",", header = FALSE) 
  
  data <- as.data.frame(t(data.raw)) %>% rename(output = V1)  %>% 
    # this column contains the timestep at which the perturbation took place
    # it just sum the values of the avalanche length
    mutate(
      ava_length = as.numeric(gsub("[[:punct:]]|[[:punct:]] \\d+[[:punct:]]", "", as.character(output))),
      affected_cells = as.numeric(gsub("[[:punct:]]\\d+[[:punct:]]|[[:punct:]]", "", as.character(output))),
      time = cumsum(ava_length)) 
  
  data$perturb <- parms.ls[[1]]
  data$perturb_n <- parms.ls[[2]]
  data$hierarchy <- parms.ls[[3]]
  data$hierarchy_n <- parms.ls[[4]]
  data$neighbours_level <- parms.ls[[5]]
  data$update <- parms.ls[[6]]
  data$steps <- parms.ls[[7]]
  data$M <- parms.ls[[8]]
  return(data)
})


# This object contains the length and affected number of cells in each run
avalanches <- do.call(rbind, data_list) 

## The following code will plot the distribution of avalanche length and affected cells for only one file


# Different ways of ploting avalanche length

## histogram


## Points of the histrogram in log scale: http://stackoverflow.com/questions/11532559/convert-bar-into-points-in-hist-function
x <- avalanches$ava_length
.hist <- hist(x, breaks=500, main ="Avalanche length")
# str(.hist)
df.plot <- data.frame(mids = .hist$mids, counts = .hist$counts)
ggplot(df.plot, aes(mids, counts)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Avalanche length")



## log-log probability plot: http://stackoverflow.com/questions/14736038/log-log-probability-chart-in-r
data_pl = displ$new(avalanches$ava_length)
est <- estimate_xmin(data_pl)
data_pl$xmin <- est
plot(data_pl)
lines(data_pl, col=2)

## Uncomment these lines to perform bootstrap analysis and calculate p-value
#bs <- bootstrap_p(data_pl)
# bs$p


## Now do the same plots for the number of affected cells

x <- avalanches$affected_cells
.hist <- hist(x, breaks=500, main ="Affected cells")

df.plot <- data.frame(mids = .hist$mids, counts = .hist$counts)
ggplot(df.plot, aes(mids, counts)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Affected cells")



## log-log probability plot: http://stackoverflow.com/questions/14736038/log-log-probability-chart-in-r
data_pl = displ$new(avalanches$affected_cells)
est <- estimate_xmin(data_pl)
data_pl$xmin <- est
plot(data_pl)
lines(data_pl, col=2)

## Uncomment these lines to perform bootstrap analysis and calculate p-value
#bs <- bootstrap_p(data_pl)
# bs$p




