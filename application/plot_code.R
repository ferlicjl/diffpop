library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)


# Set working directory to that which contains our label output files
inDir = "C:/DFCI/Jeremy/flex/diffpop/application/example output/"
setwd(inDir)

# Generate a list of all filenames
lblfiles = list.files(inDir, pattern="^out.*_label.csv$", full.names=F)

# Read in all of the label files into one dataframe
myMergedData <-
  do.call(rbind,
          lapply(lblfiles, read.csv))

# Label each row with a corresponding file id
nfiles = length(lblfiles)
myMergedData$id = rep(1:nfiles, each = 801)

# Remove any data point where the LT label has died out
#        (necessary to divide by it in next step)
myMergedData = myMergedData[myMergedData$LT != 0.0,]
myMergedData[,3:9] = myMergedData[,3:9] / myMergedData[,2]

# Melt the data set in order to plot using ggplot2
myMergedData = melt(myMergedData, id.vars = c("time", "id"))

# Rename some columns 
names(myMergedData) = c("time", "id", "pop", "value")

# Experimental results directory
exp_dir = "C:/DFCI/Jeremy/flex/diffpop/application/experimental data/"

# Plotting function
plot_dat = function(pop, df){
  exp_pts = read.csv(paste(exp_dir, tolower(pop), "_fit.csv", sep = ""), header = F)
  names(exp_pts) = c("time", "value")
  exp_pts$id = 1
  exp_pts$color = "blue"
  
  if(file.exists(paste(exp_dir, tolower(pop), "_pred.csv", sep = ""))){
    pred_pts = read.csv(paste(exp_dir, tolower(pop), "_pred.csv", sep = ""), header = F)
    names(pred_pts) = c("time", "value")
    pred_pts$id = 1
    pred_pts$color = "green"
    
    exp_pts = rbind(exp_pts, pred_pts)
  }
  
  p <- ggplot(df[df$pop == pop,], aes(x = time, y = value, group = id))
  p = p +  stat_summary(aes(group = 1), geom = "ribbon", 
                        fun.ymin = function(x) quantile(x, 0.25), 
                        fun.ymax = function(x) quantile(x, 0.75), 
                        col = "grey80", alpha = 0.3) +
    stat_summary(aes(group = 1), geom = "line", fun.y = median, col = "red", size = 2) +
    geom_point(data = exp_pts, col = exp_pts$color, alpha = 1.0) +
    xlim(0, max(exp_pts$time + 50)) +
    ggtitle(pop) +
    ylab("")
  
  return(p)
}

# Apply plotting function to each population
plot.list = lapply(c("ST", "MPP", "CLP", "CMP", "GMP", "MEP", "proB"), 
                   plot_dat, df = myMergedData)

# Arrange the plots in a nice grid
grid.arrange(grobs = plot.list, 
             left = textGrob("Percent of Labeled Cells relative to LT-HSC",
                             rot = 90, vjust = 1))