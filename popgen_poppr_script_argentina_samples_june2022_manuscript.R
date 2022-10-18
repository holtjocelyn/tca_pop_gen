#Code following Gunwald lab through Github: https://grunwaldlab.github.io/Population_Genetics_in_R/First_Steps.html

#make sure the following packages are installed. https://grunwaldlab.github.io/Population_Genetics_in_R/Getting_ready_to_use_R.html
#install.packages("poppr")
#install.packages("adegenet")
#install.packages("ape")
#install.packages("ggplot2")
#install.packages("mmod")
#install.packages("magrittr")
#install.packages("dplyr")
#install.packages("treemap")
#install.packages("vegan")

#for access to R cheat sheets go here: https://www.rstudio.com/resources/cheatsheets/
  
#load the libraries to be able to analyze data. note that if you load other libraries they may overide content from other packages
library("poppr") #will override the pegas package, need to run pegas separtely for refined anlalysis. automatically loads adegenet and ade4 package 
library("ade4")
library("adegenet") #requires ade4 package
library("ape")
library("ggplot2")
library("mmod")
library("magrittr") #allows stringing together commands left to right using pipes (%>%)
library("dplyr")
library("treemap")


#Get Session Info 
sessionInfo()

#avoid loading these libraries simultaneously or they will override content
#library("pegas") #loading these other libraries simultaneously will mask features from ape and ade4
#library("vegan") #no noticiable issues
#library("ade4") #loading this package will mask 
 
#get more info on other data file formats that poppr takes
??read.genalex


#set your working directory
setwd("~/Desktop/tca_data/popprargentina")


####load info of packages and versions used for workflow
#https://popgen.nescent.org/StartSNP.html
options(width = 100)
devtools::session_info()


########LOAD FILE AND ASSIGN STRATA########## https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html
#upload the desired file
#example of reading in a structure file and automatically converting it to a genind object
#read.structure(file, n.ind = NULL, n.loc = NULL, onerowperind = NULL, col.lab = NULL, col.pop = NULL, col.others = NULL, row.marknames = NULL, NA.char = "-9", pop = NULL, sep = NULL, ask = TRUE, quiet = FALSE)

#duplicate this file and save formatted as a .tsv. then open it in numbers or excel to view the columns. samples are in column 1, and location info is in column 2
#col.lab: an integer giving the index of the column containing labels of genotypes. '0' if absent.
#col.pop: an integer giving the index of the column containing population to which genotypes belong. '0' if absent.
#if get a dimnames error, go back into .stru file with a text edit and make sure there is not an extra space before the first marker name. also make sure that the heading of the file output is deleted. also make sure to delete the last blank row of the file
#if the file cannot open or be found, double check in the file properties that the extension is .stru rather than .stru.tsv. manually edit to .stru if this is the case

#make sure have the correct number of loci, this is the number with lesser columns from the structure output
obj <- read.structure("~/Desktop/tca_data/popprargentina/tca_argentina_samples_structure_format_edited.stru", n.ind = 10, n.loc = 26657, onerowperind = FALSE, col.lab = 0, col.pop = 2, col.others = NULL, row.marknames = 1, NA.char = "-9", pop = NULL, sep = "\t", ask = TRUE, quiet = TRUE)

#which other optional columns should be read? press press enter


#press enter, so do not make any other optional columns to read

#Warning message: In df2genind(X = X, pop = pop, ploidy = 2, sep = sep, ncode = ncode) : duplicate labels detected for some individuals; using generic labels

#call the object of the formatted genind object from the structure file
obj

# make a file with sample names and population names. OR hot tip use the same population map from STACKS and STRUCTURE. Just make sure that there are no extra spaces between collection or location names

#upload this .txt file to R called: population_map_tca_us_numbered_semicomplex_unknowns_removed_corrected_headers_june2022.txt
strata.txt <- file.choose()

#now convert this text file into a table that the program can read
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information
strata(obj) <- strata

#check to see that this information was incorporated into the gen file
View(strata(obj))



#######################
#Once do this step SKIP redoing. THIS STEP TAKES A LONG TIME for 24k+ loci
#check out the distribution of genotpyes. genotype accumulation curve. 
#http://grunwaldlab.github.io/Population_Genetics_in_R/First_Steps.html
#check out the distribution of genotpyes. 
objgc <- genotype_curve(obj, sample = 100, maxloci = 0L, quiet = FALSE, thresh = 1, plot = TRUE, drop = TRUE, dropna = FALSE)  #??genotype_ curve

#view the genotype curve. 
objgc
#####################



#######
#determine allele frequencies and any missing data, make a locus table. https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html
#data("obj")
objlt <- locus_table(obj)
#error message:  In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
#VECTOR_ELT() can only be applied to a 'list', not a 'character'


#show the alleles at each locus and eveness
objlt

#now plot the missing data. ******Takes ~ 5 minutes to run
info_table(obj, plot = TRUE) # type = missing. Reports back no missing data found, likely because it is masked with zeros

#view(info_table(obj, plot = TRUE))

#plot genotypic diversity.*****Takes ~ 10 minutes to run. https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
poppr(obj)

#calculate genotype diversity or MLG for each population. https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
p.tab <- mlg.table(obj)

#call genotypic diversity among populations
p.tab

#calculate diversity statistics on a matric with multilocus genotypes per pop
diversity_stats(p.tab, H = TRUE, G = TRUE, lambda = TRUE, E5 = TRUE)




#######get some basic diversity statistics from data
#https://popgen.nescent.org/StartSNP.html
library("adegenet")

div <- summary(obj)

names(div)

div



#do some basic stats
library("hierfstat")
basicstat <- basic.stats(obj, diploid = TRUE, digits = 2)

names(basicstat)

basicstat

#get confidence limit intervals for Fis values
boot.ppfis(obj)

#call the calculation
boot.ppfis(dat = obj)


#make a pca of the data
objpca <- indpca(obj)

#now plot the graphic of the pca
plot(objpca, cex = 0.7)



###########Find potentially uninformative loci
nLoc(obj)

iobj <- informloci(obj)

nLoc(iobj)

#all sites are polymorphic, no sites removed. no loci found with MAF < 0.01 (minor allele frequency)
####################



######CALCULATE Gst########### https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html
#Calculate Gst: https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
#install.packages("mmod")
library("poppr")
library("ape")
library("magrittr")
library("mmod")

#calculate population genetic structure using another method with the Hedrick's statistic. **Takes ~ 2 mins to run
objhen1 <- Gst_Hedrick(obj)

#FST was originally formulated to measure genetic distance using biallelic markers (Wright 1969), but was generalized for multiple alleles (Nei 1973), referred to as GST (although confusingly, the terms GST and FST are often used interchangeably
objhen1



######CALCULATE GENETIC DISTANCE####### https://grunwaldlab.github.io/poppr/reference/genetic_distance.html
objdist <- provesti.dist(obj)

#call the genetic distance estimates
objdist

#HOW IS THIS INTERPRETED? THE SCALE has what appears to be 10 clusters, which is difficult to visualize on the NJ Tree because the labels are close together.



#Population differentiation measured by FST was originally proposed by Sewall Wright (Wright, 1949). This was later extended to a method based on diversity by Masatoshi Nei (Nei, 1973). 
#As researchers applied these metrics to microsatellites, genetic markers with a large number of alleles, it became clear that Nei’s measure would not correctly range from zero to one, so Philip Hedrick proposed a correction (Hedrick, 2005). 
#More recently, Lou Jost proposed another alternative (Jost, 2008). 
#You can tell a topic is popular when so many variants of it are generated. And there are more variants than mentioned here. 
#Discussion as to which measure may be appropriate for your data was posted to the Molecular Ecologist blog titled should I use FST, GST or D?.
#https://www.molecularecologist.com/2011/03/02/should-i-use-fst-gst-or-d-2/



########NEIGHBOR JOINING ANALYSIS############ *go back and include the campanotus outgroups then rerun the tca with concensus sequences of snps 
#develop a neighbor joining tree of distances. https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
library("poppr")
library("ape") #to visualize the neighbor joining function
library("magrittr")

objtree <- objdist %>%
  nj() %>% #calculate neighbor-joining tree
  ladderize() #organize branches by clade
plot(objtree)
add.scale.bar(length = 0.05) #add a scale bar showing a 5% difference


#Add bootstrap values. *****Takes ~ 4 hours to run #rerun the neighborjoining tree when can do this step
set.seed(999)
aboot(obj, dist = provesti.dist, sample = 200, tree = "nj", cuttoff = 50, quiet = TRUE)
#aboot(obj5miss, dist = provesti.dist, sample = 200, tree = "nj", cuttoff = 50, quiet = TRUE)
#aboot(obj100miss, dist = provesti.dist, sample = 200, tree = "nj", cuttoff = 50, quiet = TRUE)

###########



#######PERFORM CLUSTER ANALYSIS#####Makes a nice plot.
#more info on msn: https://grunwaldlab.github.io/poppr/reference/bruvo.msn.html

objmxclust <- find.clusters(obj) #change this for whether analyze data obj5miss or obj100miss

#choose the number of PCs to retain = 9

#choose the number of clusters = when select 2, see that cluster 2 and 3 have the same value, but if select 3 clusters then it wont give a k for it , lowest value. After run K = 3 low value with biological meaning, both low points in the dataset

#see cluster assignments
objmxclust

#list numeric statistics of the model
stat(objmxclust)


library("shiny")

#create a minimum spanning network. this should open an interactive GUI. https://grunwaldlab.github.io/Population_Genetics_in_R/Minimum_Spanning_Networks.html
imsn()



###code for ismn gui#########

#obj_sub <- popsub(obj, exclude = character(0))
#obj_dist <- diss.dist(obj_sub, percent = FALSE, mat = FALSE)
#min_span_net <- poppr.msn(obj_sub, obj_dist, showplot = FALSE, include.ties = TRUE)

#set.seed(69)
obj_sub <- popsub(obj, exclude = character(0))
obj_dist <- diss.dist(obj_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(obj_sub, obj_dist, showplot = FALSE, include.ties = TRUE)

#set.seed(69)
#plot_poppr_msn(obj,
               min_span_net,
               inds = "ALL",
               mlg = FALSE,
               gadj = 3,
               nodescale = 10,
               palette = rainbow,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)

########end code for ismn#######

obj_sub <- popsub(obj, exclude = character(0))
obj_dist <- diss.dist(obj_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(obj_sub, obj_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
plot_poppr_msn(obj,
               min_span_net,
               inds = "ALL",
               mlg = FALSE,
               gadj = 3,
               nodescale = 45,
               palette = rainbow,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)

#########################



########RUN DAPC ON DATA#####https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
#for scattered DAPC use this code: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
#https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/inst/doc/adegenet-dapc.pdf?revision=887&root=adegenet&pathrev=887


#use data: obj

#load library
library("adegenet")

#how to create a dapc in adegenet
#https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf


#load the data
###read structure data as a genind object
objgenind <- import2genind("~/Desktop/tca_data/popprargentina/tca_argentina_samples_structure_format_edited.stru", n.ind = 10, n.loc = 26657, onerowperind = FALSE, col.lab = 0, col.pop = 2, col.others = NULL, row.marknames = 1, NA.char = "-9", pop = NULL, sep = "\t", ask = TRUE, quiet = TRUE)

#import the mapping file to go with this dataset
#upload this .txt file to R called: population_map_tca_us_numbered_semicomplex_unknowns_removed_corrected_headers_june2022.txt
strata.txt <- file.choose()

#now convert this text file into a table that the program can read
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information
strata(objgenind) <- strata

#check to see that this information was incorporated into the gen file
View(strata(objgenind))

#set the population to state or location
pop(objgenind) <- strata$location

#check to make sure this worked by calling this argument
pop(objgenind)

#Ideally, the optimal clustering solution should correspond to the lowest BIC. In practice, the ’best’ BIC is often indicated by an elbow in the curve of BIC values as a function of k.
find.clusters(objgenind) #use this curve to pick the number on the x axis where get somewhat of an asymptote, and then plug this in for the PC for the analysis below
#select 9 pcs

#other cluster selection method
#find.clusters(objgenind, clust = NULL, n.pca = NULL, n.clust = NULL, method = c("kmeans", "ward"), stat =10 c("BIC", "AIC", "WSS"), choose.n.clust = TRUE, citerion = c("diffNgroup", "min", "goesup", "smoothNgoesup", "goodfit"), max.n.clust = round(nrow(objgenind@tab)/10), n.iter = 1e5, n.start = 10, scale = FALSE, truenames = TRUE)
#the lowest BIC value is 2, but there are 3 PCs in a graph
#these values are slightly different when calculated off the genind vs. the stru file

#get approximately 2 PCs; this should be equal to or less than the actual pops have. if go too high in number than the discriminant analysis eigenvalues will not be as clear and the data points will not be visible on the scatter plot
dapc.objgenind <- dapc(objgenind, var.contrib = TRUE, scale = FALSE, n.pca = 9, n.da = NULL)

#dapc.objgenind <- dapc(objgenind, pop = strata$location, var.contrib = TRUE, scale = FALSE, n.pca = 96, n.da = NULL) #n.da means how many discriminant analysis components, if null you will manually type this in and press enter

3 #choose the number of discriminant functions to retain = 2 shows a deltak of 2 and 3 being equal. entering 3 does not provide a value.

#?pch to find out more about the different symbols associated with the pch number code

#original scatter plot with labels covering points
#scatter(dapc.objgenind)

#input a color palate to use instead of the default colors, which could overlap and be difficult to tell apart
myCol <- c("darkblue", "purple", "green") #, "orange", "brown", "black") #"black", "grey", "red", "yellow2"

#This works. Get nice spread out sample bubbles that are semi-transparent
scatter(dapc.objgenind, bg = "white", pch = 20, solid = 0.4, cell = 0, cstar = 0, cex = 2,  clab = 0, col = myCol, legend = TRUE, scree.da = FALSE) #clab = 0 is what gets rid of the sample labels

#Info about each component of the code above for the scatter
#leave pch as 20 as this keeps the data points are circles
#cex = # determines the size of each sample point in a graph, so the lower the number the smaller the data point size
#putting cell = 0 gets rid of ellipse
#putting cstar = 0 gets rid of the lines between the samples in the dapc
#scree.da = FALSE , removes the dapc eigenvalues in the lower right corner, TRUE places the DAPC eigenvalues in the figure
###########################




###CALCULATE Fst and Fis VALUES######## https://search.r-project.org/CRAN/refmans/genepop/html/genedivFis.html
## gene diversity and Fis

#load libraries
library("genepop")
#library("pegas")


#set your working directory
setwd("~/Desktop/tca_data/popprargentina")

#instruction link to genepop program
#https://cran.r-project.org/web/packages/genepop/genepop.pdf
#duplicate the genepop.gen file. now save the duplicate as .tsv extension to keep it tab delimited. open this in a text editor such as BBEDIT and delete any pop header from the identifying sample of that population that appears to group together. For example, delete the "pop" from the Georgia population since it looks like it groups closely with Florida - based upon the PCA and DAPC. Save and check to make sure the file size is the same as the .gen file. Now duplicate and save the copy as a .txt file. This will give you the new values. 

##Run Fis values for all collections
#??genedivFis {genepop}
objFis <- genedivFis("tca_argentina_samples_snps_genepop_format_edited.txt", sizes = FALSE, outputFile = "tca_genepop_fis_output.txt", dataType = "Diploid", verbose = interactive())

objFis

## run Fst analysis using genepop file from STACKS (can also create input file using PGDSpider)
#??Fst{genepop}

## run a global Fst analysis using genepop - one runs comparisons for all populations
tca_snps_genepop_global_Fst <- Fst("tca_argentina_samples_snps_genepop_format_edited.txt", sizes = FALSE, pairs = FALSE, outputFile = "tca_fst_global_output.txt", dataType = "Diploid", verbose = interactive())

## run a pairwise Fst for all populations using genepop analysis
tca_snps_genepop_pairwise_Fst <- Fst("tca_argentina_samples_snps_genepop_format_edited.txt", sizes = FALSE, pairs = TRUE, outputFile = "tca_fst_pairwise_output.txt", dataType = "Diploid", verbose = interactive())


###test of genic and genotypic differentiation
#test_diff(
#inputFile,
#genic = TRUE,
#pairs = FALSE,
#outputFile = "",
#settingsFile = "",
#dememorization = 10000,
#batches = 100,
#iterations = 5000,
#verbose = interactive()
#)

###test of genic and genotypic differentiation; will take ~ 7 hours to run
test_diff("tca_argentina_samples_snps_genepop_format_edited.txt", genic = TRUE, pairs = FALSE, outputFile = "tca_argentina_genotypic_differentiation_test.txt", settingsFile = "", dememorization = 10000, batches = 100, iterations = 5000, verbose = interactive())





#******Test with current parameters takes about 7 hours to run **********
## run a pairwise differentiation test using genepop. Can create input file from STACKS or by using PGDSpider)
#reset parameters to run dememorization 1,000, batches = 100, iterations = 5,000
#??differentiation{genepop}
tca_snps_genepop_test_diff <- test_diff("tca_us_samples_snps_genepop_format.txt", genic = TRUE, pairs = TRUE, outputFile = "tca_snps_final_differentiation_test.txt", dememorization = 10000, batches = 100, iterations = 5000, verbose = interactive())



####RUN AMOVA######## https://grunwaldlab.github.io/poppr/reference/poppr.amova.html
#https://grunwaldlab.github.io/poppr/reference/poppr-package.html
#example dataset: data("Aeut")
#The implementation of AMOVA in poppr requires two very basic components: (1) A distance matrix derived from the data and (2) a separate table used to partition the data into different stratifications.
#https://grunwaldlab.github.io/poppr/reference/poppr-package.html

#load r packages
library("adegenet")
library("poppr")
#library("ade4")
#library("pegas")
#library("vcfR")

#set working directory
setwd("~/Desktop/tca_data/popprargentina")

#make sure to leave the default heading in the file, otherwise it won't load properly
obj2 <- read.genepop("tca_argentina_samples_snps_genepop_format_edited.gen", ncode = 2L, quiet = FALSE)

#obj2 <- read.genepop("tca_argentina_samples_snps_genepop_format_5pops_edited.gen", ncode = 2L, quiet = FALSE)

#Warning message:
#In df2genind(X = X, pop = as.character(pop), ploidy = 2, ncode = ncode,  :
#               Markers with no scored alleles have been removed

#call the object of the formatted object. if structure of genepop it is a genind object. if vcf it is class vcfR. this will also give percent missing data. ~17.58%
obj2

#upload this .txt file to R called: population_map_tca_us_numbered_semicomplex_unknowns_removed_corrected_headers_june2022.txt
#for examining only five pops use file: population_map_tca_us_numbered_semicomplex_unknowns_removed_corrected_headers_5pops_june2022.txt
strata.txt <- file.choose()
#select the updated pop map that represents combined pops that had low fst values and grouped close together on DAPC and were subsequently updated in the genepop file


#now convert this text file into a table that the program can read
strata <- read.table(file = strata.txt, header = TRUE, fill = TRUE)

#view the new table that was created from the text file which is now called strata
View(strata)

#edit the strata in your gen file to now have that information
strata(obj2) <- strata

#check to see that this information was incorporated into the gen file
View(strata(obj2))





#for interpreting AMOVA data: https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html

#for the simplified equation below for calculating amova: https://grunwaldlab.github.io/poppr/reference/poppr.amova.html


#for analysis, run separately with pegas and then again with ade4
#run with pegas
pamova <- poppr.amova(
  obj2,
  hier = ~location,
  clonecorrect = FALSE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = c("pegas"),
  nperm = 999
)
#make sure to run the model with pegas soley so dont get potential issues with ade4, can also try to run with ade4 by itself


###NOW GET THE RESULTS FROM THE ANOVA ANALYSIS######
pamova


#to run amova with default settings in poppr
#https://groups.google.com/g/poppr/c/NSag-55d6bs
#poppr.amova(obj2, ~location)


#interpreting amova results
#https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html




#for analysis, run separately with ade4
pamova <- poppr.amova(
  obj2,
  hier = ~location,
  clonecorrect = FALSE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = c("ade4"),
  nperm = 999
)
#make sure to run the model with pegas soley so dont get potential issues with ade4, can also try to run with ade4 by itself

###NOW GET THE RESULTS FROM THE ANOVA ANALYSIS######
pamova



#load package ade4
library("ade4")

library("pegas")

set.seed(1999)
objsig <- randtest(pamova, nrepet = 999) #set to 999 once get to work

objsig

plot(objsig)

write.table(pamova$componentsofcovariance, sep = ",", file = "~/Desktop/tca_data/popprus/tca_pamova_ade4.csv")


#Get Session Info 
sessionInfo()
