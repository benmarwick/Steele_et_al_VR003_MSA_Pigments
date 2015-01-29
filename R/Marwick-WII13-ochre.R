# rm(list=ls()) # clear other data objects

# This file contains the functions called by the markdown file

############################################################
#' Read in the pXRF raw data
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' df <- read_in()
#' }

read_in <- function() {
# get all CSV files in all directories in this folder
files <- dir("data/pXRF", recursive=TRUE, full.names=TRUE, pattern="\\.(csv|CSV)$")
# originally from my C:\Users\marwick\Dropbox\Marwick WI13

# read all CSV files into R
library(plyr)
data <- llply(files, function(fn)  read.csv(fn,stringsAsFactor = FALSE), .progress = 'text')

# make a data frame where each column is one analysis and the colname is the filename

# create empty dataframe to hold the data
df <- data.frame(matrix(nrow=length(data), ncol = nrow(data[[1]])))

# loop through each CSV file and get the counts and put
# them into a column of a single data frame

for (i in 1:length(files)) {
df[i,] <- (data[[i]]$Bruker.AXS.Handheld.Inc..S1PXRF.Spectrum)
print(i) # show progress
}

# put the CSV file names as the column names in the new dataframe
colnames(df) <- rownames(data[[1]])
# and as first row to use later
df$filename <- files

# from line 21 - 2068 are the channels, separate them out
channels <- df[,c((which(names(df) %in% 0):ncol(df)))]
rownames(channels) <- df$filename

# first replace spaces with periods in colnames
names(df) <- gsub(" ", ".", names(df))
# extract sample LOc from filename
df$Sample.Loc <- sapply(strsplit(basename(files), "-"), "[", 2)
library(stringr)
df$Run <- as.numeric(str_extract(sapply(strsplit(basename(files), "-"), "[", 4), "[0-9]+"))

return(list(df = df, files = files))
}

############################################################
#' Get pXRF channels and compute AUC for those elements
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' aucAllSamples <- channels_and_AUC(df$df, df$files)
#' }

channels_and_AUC <- function(df, files) {
# get channels for elements, then auc for those elements
# for all samples

# look at specific channels for elements of interest
channelRanges <- read.csv("data/pXRF-channels.csv")
#channelRanges <- read.csv("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\Marwick WI13\\scripts\\pXRF-channels.csv")
# extract channels to ID elements
channelRangesStartEnd <- channelRanges[,c(1,4:5)]
# round to whole numbers
channelRangesStartEnd$Chan.Start <- round(channelRangesStartEnd$Chan.Start)
channelRangesStartEnd$Chan.End <- round(channelRangesStartEnd$Chan.End)

# subset only named channels from raw data
elements <- lapply(1:nrow(channelRangesStartEnd), function(x) df[,c((21+channelRangesStartEnd$Chan.Start[x]):(21+channelRangesStartEnd$Chan.End[x]))] )
elementsDF <- as.data.frame(do.call(cbind, elements))
# or just get one channel
# CaKa1 <- df[,c((21+channelRangesStartEnd$Chan.Start[1]):(21+channelRangesStartEnd$Chan.End[1]))]

# get list of elements
elementList <- as.character(channelRangesStartEnd[,1])

# auc for all channels of all samples
library(Bolstad2)
library(plyr)
# auc for all channels of all samples
aucAllSamples <- ldply(1:nrow(df), function(i) (sapply(1:nrow(channelRangesStartEnd),
                                   function(n) sintegral((channelRangesStartEnd$Chan.Start[n]:channelRangesStartEnd$Chan.End[n]),
                                  as.integer(df[,c((21+channelRangesStartEnd$Chan.Start[n]):(21+channelRangesStartEnd$Chan.End[n]))][i,]))$int)),
                       .progress = 'text')
colnames(aucAllSamples) <-  gsub("\\s","", elementList)
# reattach sample Loc and run number
# extract sample LOc from filename
aucAllSamples$Sample.Loc <- sapply(strsplit(basename(files), "-"), "[", 2)
# make sample names correct
aucAllSamples$Sample.Loc <- gsub("KKH1", "KKH", aucAllSamples$Sample.Loc)
aucAllSamples$Sample.Loc <- gsub("KKH2", "VR3", aucAllSamples$Sample.Loc)
aucAllSamples$Sample.Loc <- gsub("KKH3", "PL8", aucAllSamples$Sample.Loc)
aucAllSamples$Sample.Loc <- gsub("KKH4", "KFR", aucAllSamples$Sample.Loc)
# get run number
library(stringr)
aucAllSamples$Run <- as.numeric(str_extract(sapply(strsplit(basename(files), "-"), "[", 4), "[0-9]+"))


# calculate net peak area by subtracting area of background, see
# here for general concept
# http://www.osti.gov/energycitations/servlets/purl/6569198-OjE7Or/6569198.pdf
# first calculate continuum background peak for each element's channels
o <- 5 # offset from channel start and end
background <- ldply(1:nrow(df), function(i) (sapply(1:nrow(channelRangesStartEnd),
                                                    function(n) sintegral(c(channelRangesStartEnd$Chan.Start[n], channelRangesStartEnd$Chan.End[n]),
                                                                          c(as.integer(df[,c(21+channelRangesStartEnd$Chan.Start[n])- o][i]),
                                                                            as.integer(df[,c(21+channelRangesStartEnd$Chan.Start[n]+ o)][i])))$int)),
                    .progress = 'text')
colnames(background) <- gsub("\\s","", elementList)
# Now subtract background from area under curve
NPA <- aucAllSamples[,-c(ncol(aucAllSamples), (ncol(aucAllSamples)-1), (ncol(aucAllSamples)-2))] - background[,-ncol(background)]
rownames(NPA) <- basename(files)
head(NPA)
# good now convert all -ve values to zero
NPA1 <- data.frame(apply(NPA, 2, function(x) ifelse(x<0, 0, x)))
aucAllSamples <- cbind(NPA1, Sample.Loc=aucAllSamples$Sample.Loc, Run=aucAllSamples$Run)
# we seem to have two different names for the brick standard, let's fix that
aucAllSamples$Sample.Loc <- as.character(gsub("brick679", "brickn679", aucAllSamples$Sample.Loc))
unique(aucAllSamples$Sample.Loc)

return(list(aucAllSamples = aucAllSamples, elementList = elementList))
}


############################################################
#' Relative Standard Error (RSE)
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#'  my_RSE <- RSE(x)
#'  }

# RSE function
RSE <- function(x) {
  100*( (sd(x, na.rm = TRUE)/sqrt(length(x)))/mean(x, na.rm = TRUE) )
}

############################################################
#' Subset 5 runs from the 10 from each specimen
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#'  combys9 <- subset_5_runs(aucAllSamples$aucAllSamples, aucAllSamples$elementList)
#'  }

subset_5_runs <- function(aucAllSamples, elementList) {

### now subset from the ten runs to see if a more  ###################
#patterened result is obtained...

# generate all combinations of 5 runs from  ten runs
su <- 5 # number of runs to sample from ten runs
idx <- combn(unique(aucAllSamples$Run), su)

# check that it works on all ten runs per ID
library(plyr)

# make a list of dfs for all samples for each combination of five runs
# to prepare to calculate RSEs
combys1 <- lapply(1:ncol(idx), function(i) aucAllSamples[aucAllSamples$Run %in% idx[,i],] )

# prepare to calculate new colums with RSE and means
# make col names from RSEs
RSEs <- sapply(1:(ncol(aucAllSamples)-3), function(j) paste0("RSE",names(aucAllSamples[j])))
# make functions for RSEs
RSExs <- sapply(1:(ncol(aucAllSamples)-3), function(j) paste0("RSE(",names(aucAllSamples[j]),")"))
# combine col names and function
doRSE <- paste0(sapply(1:length(RSEs), function(x) paste0(RSEs[x],"=",RSExs[x])), collapse=",", sep="")

# make col names for means
meanss <- sapply(1:(ncol(aucAllSamples)-3), function(j) paste0("mean",names(aucAllSamples[j])))
# make functions for means
meanxs <- sapply(1:(ncol(aucAllSamples)-3), function(j) paste0("mean(",names(aucAllSamples[j]),")"))
# combine col names and function
domeans <- paste0(sapply(1:length(meanss), function(x) paste0(meanss[x],"=",meanxs[x])), collapse=",", sep="")

# combine RSE and means
func <- paste(doRSE, domeans, collapse = ",", sep =",")

# calculate RSEs and means for all combinates of 5 runs, @mnel's solution, adapted, this is quite slow

combys2 <- lapply(combys1,
                  function(x) do.call(ddply,c(.data = quote(x),
                     .variables = quote(.(Sample.Loc)),
                     .fun = quote(summarise),
                      eval(parse(text = sprintf('.(%s)', func ))))))
head(combys2)


# put all the data into one dataframe with col indicating origin
library(plyr)
names(combys2) <- paste0('df',seq_along(combys2))
combys3 <- ldply(combys2, rbind) # takes a moment....
head(combys3)

# find min RSE for each element for each Sample.Loc and get mean
combys4 <- ddply(combys3, .(Sample.Loc), numcolwise(min) )
combys5 <- cbind(combys4[,1], combys4[,RSEs])
head(combys5)

# convert to long table
elemsnospace <- gsub("\\s","", elementList)
combys6 <- lapply(elemsnospace[-length(elemsnospace)], function(i)
  data.frame(RSE=eval(parse(text=paste0("combys5$RSE",i))),
             mean=eval(parse(text=paste0("combys4$mean",
                                         i,"[match(combys5$RSE",
                                         i,", ", "combys4$RSE",i,")","]")))))
names(combys6) <- elemsnospace[-length(elemsnospace)]

# exclude elements where RSE could not be computed
have_RSE <- lapply(combys6, function(i) sum(i$RSE, na.rm = TRUE)) != 0
combys6a <- combys6[have_RSE]

combys7 <- ldply(combys6a, rbind)
head(combys7)
combys7$Sample.Loc <- rep(unique(aucAllSamples$Sample.Loc), length(combys6a))
# two elements seem to be missing...
# remove RSE, so we just have mean, element and sample Location
combys8 <- combys7[,-2]
library(reshape2)
combys9 <- dcast(combys8, Sample.Loc ~ .id, value.var='mean')
head(combys9)

# change col names to be just Element names
names(combys9)[2:ncol(combys9)] <- gsub("mean", "", names(combys9[2:ncol(combys9)]) )
names(combys9)[2:ncol(combys9)] <- gsub("Ka1", "", names(combys9[2:ncol(combys9)]) )
# check
names(combys9)
return(combys9)
}


############################################################
#' Compute correlations with Fe
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' combys9 <- corr_with_Fe(combys9)
#' }

corr_with_Fe <- function(combys9) {
# Look to see which elements are correlated with Fe,
# following Popelka-Filcoff et al. 2007, 2008

# now look for correlations with F
library("dplyr")
library("reshape2")
library("Hmisc")

# exlcude standards because they are not ochre
combys9a <- combys9 %>%
            filter(!Sample.Loc %in% c("brickn679", "n97bflint", "n98bclay"))
unique(combys9a$Sample.Loc)

corrs <- combys9a %>%
  select(-Sample.Loc) %>%
  as.matrix() %>%
  rcorr()

# r-values
r_vals <- corrs$r %>%
  melt() %>%
  arrange(-abs(value)) %>%
  filter(Var1 == "Fe") %>%
  mutate(elem  = paste0(Var1, Var2))

# p-values
p_vals <- corrs$P %>%
  melt() %>%
  arrange(abs(value)) %>%
  filter(Var1 == "Fe") %>%
  mutate(elem  = paste0(Var1, Var2))

# combine p and r values and get elements with
# significant positive correlation with Fe concentration
pos_corr_with_Fe <- merge(r_vals, p_vals, by = 'elem') %>%
  arrange(-value.x) %>%
  filter(value.y < 0.1 & value.x > 0)

# elements remaining are...
remaining <- as.character(pos_corr_with_Fe$Var2.x)

# so just keep these for the rest of the analysis...
combys9 <- combys9[, c('Sample.Loc', 'Fe', remaining)]
head(combys9)

return(list(combys9 = combys9, remaining = remaining))
}



############################################################
#' Calibration
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#'  combys12 <- calibration(combys9$combys9, combys9$remaining)
#' }

calibration <- function(combys9, remaining) {

# Now we have mean counts per element with background removed for each speciment
# we can now use standards data to convert from AUC to absolute measurement

# Make regression from standards to compute calibrated values

brick <- read.csv("data/pXRF/standards/SRM_679_Brick_Clay.txt")
# assume 1 kg of brick to make all measurements consistent
# convert  mg/kg to %mass (mg/kg = ppm)
options(scipen=999)
brick$perc_mass <- ifelse(brick$unit == "%mass",
                          brick$Quantity ,
                          brick$Quantity / 10000 ) # ten thou
brick$oxide <- brick$perc_mass * brick$oxide.conversion
# check to see how that sums
sum(brick$oxide)

# brick gives us Cr, Co, Ti and Fe, use other standards to get other Ni, Pb, Kr, V. We can get Ni and Pb from NIST-97b flint clay, according to Liu, T-B., Maynard, J.B., and Alten, J., 2006, Superheavy S isotopes from glacier-associated sediments of the Neoproterozoic of south China: Oceanic anoxia or sulfate limitation?, in Kesler, S.E., and Ohmoto, H., eds., Evolution of Early Earth's Atmosphere, Hydrosphere and Biosphere--Constraints from Ore Deposits: Boulder CO, Geological Society of America, Memoir 198, p. 205-222.
# But it seems the Ni and Pb don't have a +ve cor with Fe, so they don't make the cut

flint <- read.csv("data/pXRF/standards/97a.txt")
# assume 1 kg of sample to make all measurements consistent
# convert  mg/kg to %mass (mg/kg = ppm)
options(scipen=999)
flint$perc_mass <- ifelse(flint$unit != "ppm",
                          flint$Quantity ,
                          flint$Quantity / 10000 ) # ten thou
# no oxides
flint$oxide <- NA

plastic <- read.csv("data/pXRF/standards/97b.txt")
# assume 1 kg of sample to make all measurements consistent
# convert  mg/kg to %mass (mg/kg = ppm)
options(scipen=999)
plastic$perc_mass <- ifelse(plastic$unit != "ppm",
                            plastic$Quantity ,
                            plastic$Quantity / 10000 ) # ten thou
# no oxides
plastic$oxide <- NA

# combine flint, plastic and brick dataframes and keep names of dataframe
standard <- do.call(rbind, list(plastic=plastic, flint=flint, brick=brick))
standard$ID <- rownames(standard)

# get XRF data on standards

# get XRF counts of standards
standard_XRF <- combys9[combys9$Sample.Loc %in% c('brickn679', 'n97bflint', 'n98bclay'), ]

# get NIST %mass for the elements that we have XRF counts for
library("dplyr")
standard_percmass <- standard %>%
                        filter(Element1 %in% c('Fe', remaining)) %>%
                        select(perc_mass, Element1, ID)

standard_percmass_m <- melt(standard_percmass)
standard_XRF_m <- melt(standard_XRF)
# make the standard labels identical
standard_percmass_m$ID <- gsub("\\.[0-9]+", "", standard_percmass_m$ID)
standard_percmass_m$Sample.Loc <- with(standard_percmass_m,
                               ifelse(ID == 'brick', 'brickn679',
                                      ifelse(ID == "plastic", 'n98bclay',
                                             ifelse(ID == "flint", "n97bflint", "NA"))))
standard_percmass_m$element <- standard_percmass_m$Element1
standard_XRF_m$element <- standard_XRF_m$variable

# join up the XRF counts and the NIST values for the standards
standard_merge <- merge(standard_percmass_m, standard_XRF_m, by = c('Sample.Loc', 'element'))

# plot
library("ggplot2")
ggplot(standard_merge, aes(value.x, value.y, colour = Sample.Loc, group = 1)) +
  geom_point(size = 4) +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~element, scales = c("free"))

# flint in the Fe plot looks odd, exclude it...
standard_merge <- standard_merge %>%
                      filter(Sample.Loc != "n97bflint" | element != 'Fe')

# compute linear regression equations,
# then plug-in measured values from arch/geo samples

models <- dlply(standard_merge, "element", function(df)
  lm(value.x ~ value.y, data = df))

# Apply coef to each model and return a data frame
models_coef <- ldply(models, coef)

# Print the summary of each model
# l_ply(models, summary, .print = TRUE)

# Calibrate the elements in the standards
# get the elements that we've measured
standard_meas <- intersect(standard$Element1, names(combys9))
# subset measurement data for calibation
combys10 <- combys9[, standard_meas]

# apply calibration to measured data
calib_tmp <- vector("list", nrow(models_coef))
for(i in seq(nrow(models_coef))){
  elem_tmp <-combys10[, names(combys10) == models_coef$element[i]]
  calib_tmp[[i]] <- models_coef$`(Intercept)`[i] +
                  (models_coef$value.y[i] * elem_tmp)
}

# Co is a special case... deal with it - warning hard-coding here...
Co <- standard_merge %>%
        filter(element == "Co")
Co_multi <- Co$value.y * Co$value.x
# test to get certified value
Co_multi / Co$value.y

# apply to output from loop
calib_tmp[[length(calib_tmp)]] <- Co_multi / combys10$Co

# convert list to data frame
names(calib_tmp) <- standard_meas
combys11 <- do.call(cbind.data.frame, calib_tmp)
combys11$Sample.Loc <- combys9$Sample.Loc

# remove standards from data
combys11 <- combys11 %>%
  filter(!Sample.Loc %in% c("brickn679", "n97bflint", "n98bclay"))
unique(combys11$Sample.Loc)

# combys11 <- data.frame(Reduce(rbind, apply(combys10, 1, function(i) i * standard_multi)))

# express elments as log10 of ratio to FeOx in each sample
combys11a <- combys11[, names(combys11) != 'Sample.Loc']
combys12 <- log10(combys11a[, names(combys11a) != 'Fe'] / combys11a$Fe)
combys12$Sample.Loc <- combys9$Sample.Loc[!combys9$Sample.Loc %in% c("brickn679", "n97bflint", "n98bclay")]
rownames(combys12) <- combys12$Sample.Loc
head(combys12)

return(combys12)
}


############################################################
#' plots
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' plots(combys12)
#' }

plots <- function(combys12) {
# with just two elements, it's simplest just to plot them together
# without any multivariate stats...

# remove periods and numbers from sample names
things_to_remove <- c("[0-9]+", "\\.")
combys12$plotsymbol <- gsub(paste(things_to_remove, collapse = "|"), "", combys12$Sample.Loc)

# ggplot(combys12, aes(Ti, Co)) +
#   geom_text(aes(label = Sample.Loc)) +
#   theme_minimal() +
#   xlab(expression(log[10]*"[Ti/Fe]")) +
#   ylab(expression(log[10]*"[Co/Fe]"))
#
#
# ggplot(combys12, aes(Ti, Co)) +
#   geom_text(aes(label = plotsymbol)) +
#   theme_minimal() +
#   xlab(expression(log[10]*"[Ti/Fe]")) +
#   ylab(expression(log[10]*"[Co/Fe]"))


# quick PCA
combys13 <- combys12 %>% select(-Sample.Loc, -plotsymbol)
# remove any samples with non-numerics, esp Inf
# row.has.Inf <- apply(combys13, 1, function(x) sum(x) != Inf )
# combys13 <- combys13[row.has.Inf, ]
pca1 <- prcomp(combys13, scale=TRUE)
# remove VR3.02 an an outlier
# combys13 <- combys13[rownames(combys13) != 'VR3.02', ]
pca1 <- prcomp(combys13, scale=TRUE)

# biplot(pca1)

# More fancy
library(ggbiplot)
g <- ggbiplot(pca1, obs.scale = 1, var.scale = 1,
              groups =  combys12$plotsymbol, ellipse = TRUE,
              circle = TRUE, alpha = 0)
g <- g + scale_color_discrete(name = '') +
         geom_text(label = combys12$plotsymbol)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top') +
        xlim(c(-4,4)) +
        theme_minimal()
print(g)
ggsave("figures/pca_biplot.png",width = par("din")[1]*1.1, height = par("din")[2]/2*1.1)


# plot a cluster analysis
hc <- hclust(dist(scale(na.omit(combys13))))
png("figures/dendroplot.png", width=1233, height=484, res = 100)
plot(hc, hang = -1, cex = 0.75, main = "", xlab = "", sub = "")
dev.off()




## Use the silhouette widths for assessing the best number of clusters
# library(fpc)
# library("cluster")
# asw <- numeric(20)
# for (k in 2:20){
#   asw[k] <- pam(combys13, k) $ silinfo $ avg.width
# }
# k.best <- which.max(asw)
# cat("silhouette-optimal number of clusters:", k.best, "\n")
#
# # K-Means Cluster Analysis
# fit <- kmeans(combys13, k.best) #  cluster solution
# # get cluster means
# # aggregate(combys13,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# combys13_clus <- data.frame(combys13, fit$cluster)
# # Cluster Plot against 1st 2 principal components
# # vary parameters for most readable graph
# library(cluster)
# clusplot(combys13_clus, fit$cluster, color=FALSE, shade=FALSE,
#          col.clus = "black", labels=2, lines=0, main = "")

}

# consider to combine elemental and mag sus data to cluster...

############################################################
#' ochre mass
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' ochre_data <- get_ochre_data()
#' }

get_ochre_data <- function() {
  ochre <- read.csv("data/ochre.csv", stringsAsFactors = FALSE)[1:79,]
  ochre <- ochre[ochre$weight != "NA" & ochre$site != "NA" & ochre$site != "", ]
}


  ############################################################
  #' ochre mass
  #'
  #'
  #'
  #' @export
  #' @examples
  #' \dontrun{
  #' ochre_mass(ochre_data)
  #' }

ochre_mass <- function(ochre_data) {

ggplot(ochre_data, aes(x=site, y=weight, colour=context)) +
  geom_boxplot() +
  theme_minimal()
}


############################################################
#' Bayesian t-test
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' HDI_mean_diff <- b_t_t(ochre_data)
#' }

b_t_t <- function(ochre_data) {

library("BayesianFirstAid")
# compute bayesian t-test
weight_by_context <- bayes.t.test(weight ~ context, data = ochre_data)
# show plot
plot(weight_by_context)
# get values of 95% HDI of difference in means to use in text
HDI_mean_diff <- (weight_by_context$stats)[5,5:6]

return(HDI_mean_diff)
}

############################################################
#' Ochre surface
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' surface <- surface(ochre_data)
#' }

surface <- function(ochre_data) {

surface <- ochre_data[ochre_data$rounded != "" & ochre_data$rounded != "NA", ]
surface <- data.frame(prop.table(table(surface$rounded)) * 100)

return(surface)

}


############################################################
#' Mag sus plot
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' mag_sus_plot(ochre_data)
#' }

mag_sus_plot <- function(ochre_data) {


library("scales")
library("ggplot2")
ochre_data$FD <- as.numeric(as.character(ochre_data[,grepl("frequency.dependent", names(ochre_data))]))
ggplot(ochre_data, aes(x=LF.mass.specific.susceptibility, y=FD, colour=site)) +
  geom_point(size = 4) +
  scale_y_continuous(trans=log2_trans()) +
  scale_x_continuous(trans=log2_trans()) +
  theme_minimal() +
  ylab("Frequency dependancy (%)") +
  xlab(expression(Low~frequency~mass~susceptibility~~10^-6~m^3~kg^-1))

}



############################################################
#' Freq dep plot
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' freq_dep_plot(ochre_data)
#' }

freq_dep_plot <- function(ochre_data) {

  ochre_data$FD <- as.numeric(as.character(ochre_data$Percentage.frequency.dependent.))
library("scales")
library("ggplot2")
ggplot(ochre_data, aes(x=site, y=LF.mass.specific.susceptibility, colour=context)) +
  geom_boxplot() +
  scale_y_continuous(trans=log2_trans()) +
  theme_minimal() +
  ylab(expression(Low~frequency~mass~susceptibility~~10^-6~m^3~kg^-1))
}

############################################################
#' Mag sus corr
#'
#'
#'
#' @export
#' @examples
#' \dontrun{
#' HDI_corr <- mag_sus_corr(ochre_data)
#' }

mag_sus_corr <- function(ochre_data) {

library("BayesianFirstAid")

ochre_data$FD <- as.numeric(as.character(ochre_data$Percentage.frequency.dependent.))

mag_sus_corr <- with(ochre_data, bayes.cor.test(FD, LF.mass.specific.susceptibility))
HDI_corr <- mag_sus_corr$stats[1,5:6]

return(HDI_corr)
}

# ############################################################
# #' Table summarising ochre data
# #'
# #'
# #'
# #' @export
# #' @examples
# #' \dontrun{
# #' the_ochre_table <- ochre_table(ochre_data)
# #' }
#
# ochre_table <- function(ochre_data) {
# # table summarising ochre data
# the_ochre_table <- ochre_data %>%
#                     select(Loc, sq, layer.bag, weight, length, width, thickness,
#                            LF.mass.specific.susceptibility,
#                            Percentage.frequency.dependent.susceptibility..Îºfd..or.Ï.fd..)
#
# # %FD should be numeric also
# the_ochre_table$Percentage.frequency.dependent.susceptibility..Îºfd..or.Ï.fd.. <- as.numeric(the_ochre_table$Percentage.frequency.dependent.susceptibility..Îºfd..or.Ï.fd..)
#
# # do some rounding (http://stackoverflow.com/a/21328269/1036500)
# roundIfNumeric <- function(x, n=1)if(is.numeric(x)) round(x, n) else x
#
# the_ochre_table <- as.data.frame(
#                     lapply(the_ochre_table, roundIfNumeric, 2)
#                     )
# # add Fe, Co and Ti %mass columns
#
# ### run all the code...
# }

