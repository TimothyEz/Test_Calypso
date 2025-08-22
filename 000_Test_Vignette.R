# For some reason, after closure it's no longer a git repository on calypso drive.
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
## install timeOmics
BiocManager::install('timeOmics')

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("lmms", quietly = TRUE)) install.packages("lmms")

# load libraries 
library(timeOmics)
library(tidyverse)

# load data
data("timeOmics.simdata")
sim.data <- timeOmics.simdata$sim # Attempt to check edit 

dim(sim.data) 

# ---- Data Processing ----
# Platform-specific: Normalize based on data type

# Time-specific
remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X),
                      function(x) abs(sd(x)/mean(x))))
  return(X[, cv > cutoff])
}

data.filtered <- remove.low.cv(sim.data, 0.5)

# ---- Time Modelling ----
devtools::install_github("cran/lmms")
library(lmms)

# numeric vector containing the sample time point information 
time <- timeOmics.simdata$time
head(time)

# example of lmms
lmms.output <- lmms::lmmSpline(data = sim.data, time = time, sampleID = rownames(data.filtered), deri = FALSE, # had to change data.filtered to sim.data and basis to cubic to get plotLong to work
                               basis = "cubic", numCores = 4, timePredict = 1:9, 
                               keepModels = TRUE)

modelled.data <- t(slot(lmms.output, 'predSpline'))


### Let's Plot the modeled Profiles ###

# gather data
data.gathered <- modelled.data %>% as.data.frame() %>%
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to = "feature", values_to = 'value', -time)

# plot profiles
ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() + ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")

### Remove straight lines ###

# To remove straight line models i.e. inter-individual variation is too high (using 2-phase test procedure)
filter.res <- lmms.filter.lines(data = data.filtered, 
                                lmms.obj = lmms.output, time = time)
profile.filtered <- filter.res$filtered


# ----Single-Omic longitudinal clustering ----
# run pca 
pca.res <- pca(X = profile.filtered, ncomp = 5, scale = FALSE, center = FALSE) # double check this

# tuning ncomp 
pca.ncomp <- getNcomp(pca.res, max.ncomp = 5, X = profile.filtered,
                      scale = FALSE, center = FALSE)
pca.ncomp$choice.ncomp
plot(pca.ncomp)

# final model 
pca.res <- pca(X = profile.filtered, ncomp = 2, scale = FALSE, center = FALSE)

# extract cluster 
pca.cluster <- getCluster(pca.res)

plotIndiv(pca.res)
plotVar(pca.res)
plotLoadings(pca.res)
plotLong(pca.res, scale = FALSE, center = FALSE, 
         title = "PCA longitudinal clustering")


tune.spca.res <- tuneCluster.spca(X = profile.filtered, ncomp = 2, 
                                  test.keepX = c(2:10))
# selected features in each component
tune.spca.res$choice.keepX

plot(tune.spca.res)

# final model
spca.res <- spca(X = profile.filtered, ncomp = 2, 
                 keepX = tune.spca.res$choice.keepX, scale = FALSE)
plotLong(spca.res)
