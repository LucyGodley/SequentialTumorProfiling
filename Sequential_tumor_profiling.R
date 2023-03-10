#Analysis of data for the manuscript titled "Using sequential tumor molecular profiling to identify germline variants" U of C OncoPlus somatic sequencing panel

#load the required analysis packages
library(ggplot2)
library(cli)
library(QuantPsyc)
library(pastecs)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(clinfun)
library(lubridate)
library(forcats)
library(RColorBrewer)
library(psych)
library(cvequality)
library(ggbeeswarm)
library(ROCit)
library(pROC)
library(VennDiagram)
library(stringr)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

#Load in Required raw data

#read in the clinically annotated reference genome from ClinVar downloaded 09/1/21
clinvar <- read.delim(file = "variant_summary.txt")

#subset the reference to only include those mapped to GRch37
clinvar <- subset(clinvar, clinvar$Assembly == "GRCh37")

#read in the gnomAD annotation downloaded and compiled on 12/6/19
gnomad <- read.csv(file = "gnomAD_OncoPlus.csv", header = TRUE)

#read in the COSMIC annotation downloaded and compiled on 10/9/19
cosmic <- read.csv(file = "Cosmic_OncoPlus.csv", header = TRUE)

#Spreadsheet with gene lists to isolate later, i.e. the ChIP genes and known Germline genes
genelist <- read.csv("GeneLists.csv", header = TRUE, na.strings=c("","NA"))

#load in the Variant Caller File and removes the row names column
cancer <- read.csv(file = "Compiled_Germline_Panel2.csv", na.strings = "999", header = TRUE)
cancer <- subset(cancer, select = -X)
cancer$DOB <- as.Date(cancer$DOB)

#removes duplicate rows in case the same month is added 2x
cancer <- distinct(cancer)

#fixes the C.Nom column so that they are all consistent with the proper annotation
cancer$C.Nom <- gsub("&gt;", ">", cancer$C.Nom)
cancer$C.Nom <- as.factor(cancer$C.Nom)

#ensures that the allele frequency is numeric
cancer$AF <- as.numeric(as.character(cancer$AF))
cancer$AF <- cancer$AF %>% replace_na(0)

#changes the Request.Date column to a date variable type
cancer$Request.Date <- as.Date(cancer$Request.Date , "%m/%d/%y")
class(cancer$Request.Date)
min(cancer$Request.Date)
max(cancer$Request.Date)

#EXCLUDES DATES FROM DOWNSTREAM ANALYSIS (use this for paper figures)
cancer <- cancer %>% dplyr::filter(Request.Date <= "2021-12-31")

#creates the data normalization function called "std" for standard date
std <- function(x) {  
  ((min(x) - x)) * -1 
}

#uses the dplyer package to perform the function on each patient and put it into a column called time_point
cancer <- cancer %>% group_by(MRN) %>% dplyr::mutate(time_point = std(Request.Date))

#Ensures the time_point date difference is measured in an interger representing the number of days
cancer$time_point <- as.integer(cancer$time_point)

#creates an age column with the patient age on the date of sequencing
cancer$age <- as.numeric(cancer$Request.Date - cancer$DOB)

#add a column to the cancer dataframe with the "Name" column to match the clinvar dataframe so that we can annotate it with ClinVar and COSMIC
cancer <- cancer %>% separate(C.Nom, c("G.Nom", "S.Nom"), sep = ":")
cancer$Name <- paste(cancer$G.Nom, cancer$Gene, sep = "(")
cancer$Name <- paste(cancer$Name, cancer$S.Nom, sep = "):")
cancer$Name <- paste(cancer$Name, cancer$P.Nom, sep = " (")
cancer$del <- ")"
cancer$Name <- paste(cancer$Name, cancer$del, sep = "")
cancer$del <- NULL

#add a column with the variance between mutations: use var(AF) if you want variance, and sd(AF)/mean(AF) if you want the coefficient of variation
cancer <- cancer %>% group_by(paste(Gene,P.Nom,MRN)) %>% dplyr::mutate(variance = sd(AF)/mean(AF))
cancer <- cancer[ , !(names(cancer) %in% c("paste(Gene, P.Nom, MRN)"))]

#makes a column with the unique nucleotide name
cancer$gene_snv <- paste(cancer$Gene, " ", cancer$S.Nom)

#makes a column with the unique protein mutation name
cancer$mutation <- paste(cancer$Gene, " ", cancer$P.Nom)

#combine the cancer dataframe with the ClinVar dataframe by the new column that we previously created called Name
cancer <- merge(cancer, clinvar[,c("Name","ClinicalSignificance","OriginSimple","PhenotypeList")], by = "Name", all.x = TRUE, all.y = FALSE)

#puts a value in the cells that are not in ClinVar
cancer$ClinicalSignificance <- cancer$ClinicalSignificance %>% fct_explicit_na(na_level = "not in ClinVar")
cancer$OriginSimple <- cancer$OriginSimple %>% fct_explicit_na(na_level = "not in ClinVar")

#makes columns in the gnomAD dataframe so this can be added to annotate our data also
gnomad$gene_snv <- paste(gnomad$Gene, " ", gnomad$S.Nom)
gnomad$mutation <- paste(gnomad$Gene, " ", gnomad$P.Nom)

#Annotate the OncoPlus data with the gnomAD dataframe
cancer <- merge(cancer, gnomad[,c("rsID","Chromosome","Position","Reference","Alternate","Consequence","Annotation","Allele.Frequency","gene_snv")], by = "gene_snv", all.x = TRUE, all.y = FALSE)

#makes a column in the cosmic dataframe that we can use to combine with the cancer dataframe to annotate it
cosmic$gene_snv <- paste(cosmic$GENE_NAME, " ", cosmic$MUTATION_CDS)

#This function actually works for the whole dataset. There just cannot be na in the variance column, so we make it 1 when there are not recurrent variants
germline <- function(CV) {
  if_else(CV <= 0.03 & max(CV) > 0.12, grade <- 1, 
  if_else(CV < 0.03 & max(CV) < 0.12, grade <- 2,
  if_else(CV > 0.03 & CV <= 0.05 & max(CV) > 0.12, grade <- 2,
  if_else(CV > 0.03 & CV <= 0.05 & max(CV) < 0.12, grade <- 3,
  if_else(CV > 0.05 & CV <= 0.12 & max(CV) > 0.12, grade <- 3,
  if_else(CV > 0.05 & CV <= 0.12 & max(CV) < 0.12, grade <- 4, grade <- 5))))))
}

germline <- function(CV) {
  dplyr::if_else(CV <= 0.03 & max(CV) > 0.1, grade <- 1, 
  dplyr::if_else(CV < 0.03 & max(CV) < 0.1, grade <- 2,
  dplyr::if_else(CV > 0.03 & CV <= 0.05 & max(CV) > 0.1, grade <- 2,
  dplyr::if_else(CV > 0.03 & CV <= 0.05 & max(CV) < 0.1, grade <- 3,
  dplyr::if_else(CV > 0.05 & CV <= 0.1 & max(CV) > 0.1, grade <- 3,
  dplyr::if_else(CV > 0.05 & CV <= 0.2 & max(CV) < 0.2, grade <- 4, grade <- 5))))))
}

#testing the function out ###### (pt w/o CV > 0.012), ###### (pt w/ CV > 0.12)
cancer$variance <- cancer$variance %>% replace_na(1)
cancer <- cancer %>% dplyr::group_by(MRN) %>% dplyr::mutate(grade = germline(variance))

#Adds a column with the inheritance pattern for each gene
cancer <- merge(cancer, genelist[,c("Germline_Blood","Inheritance_Blood")], by.x = "Gene", by.y = "Germline_Blood", all.x = TRUE, all.y = FALSE)
cancer$Inheritance_Blood <- dplyr::if_else(is.na(cancer$Inheritance_Blood), "Unknown", cancer$Inheritance_Blood)

#now create a new dataframe to work from to do the germline analysis, either use the entire dataset or exclude dates
germ <- cancer
noquote(unique(germ$Gene))
#creates a data table of patients who had the same gene pop up >= 2 times
germ2 <- germ %>% group_by(Patient, MRN, Gene) %>% dplyr::summarize(count = n())
germ2 <- subset(germ2, germ2$count >= 2)
length(unique(germ2$MRN))

#creates a column in the germ data frame indicating if the patient has recurring mutations in the same gene
germ$reocc <- germ$MRN %in% germ2$MRN

#makes a new dataframe with data from with only patients with reoccuring mutations
germ3 <- subset(germ, germ$reocc == "TRUE")
length(unique(germ3$MRN))
length(unique(paste(germ3$MRN, germ3$Request.Date)))

#creates a data table of patients who had the same protein change >= 2x
germ4 <- germ3 %>% group_by(Patient, MRN, P.Nom) %>% dplyr::summarize(count = n())
germ4 <- subset(germ4, germ4$count >= 2)

#creates a column in the germ3 data frame indicating if the patient has recurring mutations in the same gene
germ3$reocc <- (germ3$P.Nom %in% germ4$P.Nom) & (germ3$MRN %in% germ4$MRN)

#remakes the dataframe with data from with only patients with reoccuring mutations
germ4 <- subset(germ3, germ3$reocc == "TRUE")

#re-orders the data set by MRN, Gene, and then reverse test date
germ4 <- germ4[order(germ4$Patient, germ4$Gene, germ4$P.Nom, rev(germ4$Request.Date)),]

#tells how many patients there are with reoccuring variants in same gene
length(unique(germ4$MRN))

#creates a spreadsheet with the patient information and the number of times the SNP appeared all in the correct order
germ5 <- germ4 %>% group_by(Patient, Name) %>% dplyr::summarize(count = n())
germ5 <- merge(germ4, germ5, by=c("Patient","Name"))
germ5 <- subset(germ5, germ5$count > 1)
germ5 <- germ5[order(germ4$Patient, germ4$Gene, germ4$P.Nom, rev(germ4$Request.Date)),]

#subsets AF >= 0.3, may add/delete "& germ5$AF <= 0.6" if needed
germ6 <- subset(germ5, germ5$AF >= 0.3)

#subsets the mutations that have population frequency less than 1x10^-4 or are NA
germ7 <- subset(germ6, Allele.Frequency < 0.0001 | is.na(Allele.Frequency))

#creates a function that removes rows that have NA in a specified column
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#uses the function created to remove NA in the Patient column
germ5 <- completeFun(germ5, "Patient")

#Figure 1A: Normalized patient with germline CSF3R variant (patient #######) - exported SVG 650x400
ggplot(data = subset(cancer, cancer$MRN  == "#######"), aes(x = time_point/30, y = AF, col = paste(Gene," ",S.Nom," ",P.Nom))) + ylim(0,1) + geom_line() + geom_point(size = 3) + xlab("Months") + ylab("Allele Frequency") + ggtitle("Patient variants over time") + theme_classic() + scale_color_brewer(palette="Set1")

#Figure 1B: Coefficient of Variation formulas and these data
as.data.frame(cancer %>% dplyr::filter(MRN  == "#######", Gene == "CSF3R") %>% dplyr::select(AF, variance))
as.data.frame(cancer %>% dplyr::filter(MRN  == "#######", Gene == "RUNX1") %>% dplyr::select(AF, variance))
median(as.data.frame(cancer %>% dplyr::filter(MRN  == "#######", Gene == "CSF3R") %>% dplyr::select(AF))$AF)

#Figure 1C: Plot of descending order of the counts. Numbers used are below.
#Complete Dataset
#number of patients
length(unique(germ$MRN))
#number of OncoPlus Assays
length(unique(paste(germ$Patient, germ$Request.Date)))
#number of mutations that were found
length(germ$Name)
#number of unique mutations that were found
length(unique(germ$Name))

#Patients with Recurrent Variants
#number of patients
length(unique(germ5$MRN))
#total number of variants
length(germ5$Name)
#number of unique variants
length(unique(germ5$Name))
#number of genes 
length(unique(germ5$Gene))

#number of patients with CV < 0.108 and population frequency <= 0.0001
length(unique(germ7$MRN))
length(unique(germ7$MRN))/length(unique(cancer$MRN))
#total number of variants
length(germ7$gene_snv)
length(germ7$gene_snv)/length(germ$gene_snv)

#number of unique variants
length(unique(germ7$gene_snv))
length(unique(germ7$gene_snv))/length(unique(germ$gene_snv))
#number of genes
length(unique(germ7$Gene))
noquote(unique(germ7$Gene))

#Figure 1D: number of germline mutations via bar graph by gene
#load in the edited data
germ9 <- read.csv(file = "validation_cohort_final.csv", header = TRUE)

#add additional column if you decided to change the grading algorithm
germ9 <- merge(germ9, germ[,c("MRN","Gene","S.Nom","grade")], by = paste(c("MRN","Gene","S.Nom")), all.x = TRUE, all.y = FALSE)

#renames the columns
germ9 <- germ9 %>% dplyr::rename(grade = grade.x, grade.new = grade.y)

#makes sure we can remove the patients that did not consent for further research
germ9$consented.to.research <- as.factor(germ9$consented.to.research)
class(germ9$consented.to.research)
unique(germ9$consented.to.research)

#Here is the graph for Figure 1D
germ9 %>% 
  dplyr::filter(consented.to.research != "No", etiology == "germline") %>%
  dplyr::group_by(Gene, path.level) %>% 
  dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(freq = sum(count)) %>%
  ggplot(aes(x = reorder(Gene,(-freq)), y = count, fill=path.level)) +
  geom_bar(position="stack", stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values = c("#4DAF4A","#E41A1C","#377EB8", "#E3C535"))

#Here are the numbers if needed
germ9 %>% filter(consented.to.research != "No", etiology == "germline") %>% group_by(Gene, path.level) %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))  %>% as.data.frame() %>% arrange(desc(freq) ) 

germ9 %>% dplyr::filter(consented.to.research != "No", etiology == "germline") %>% dplyr::group_by(path.level) %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))  %>% as.data.frame() %>% arrange(desc(freq) ) 

#total number of patients with potential germline tests
length(unique(germ9$MRN))
#total number of patients consented to research
length(unique(filter(germ9, consented.to.research != "No")$MRN))
#total number of patients with material for germline confirmation
length(unique(filter(germ9, consented.to.research != "No", etiology != "not tested")$MRN))

#Numbers for total patients consented to research
length(unique(filter(germ9, consented.to.research != "No")$MRN))

#Here is the graph for Figure 1E
germ9 %>% 
  dplyr::filter(consented.to.research != "No", etiology != "not tested") %>%
  dplyr::group_by(grade.new, etiology) %>% 
  dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(freq = sum(count)) %>%
  ggplot(aes(x = grade.new, y = count, fill=etiology)) + 
  geom_bar(position="fill", stat = 'identity') + theme_classic() + scale_fill_brewer(palette="Set1")

#Numbers for each
germ9 %>% dplyr::filter(consented.to.research != "No", etiology != "not tested") %>% dplyr::group_by(grade.new, etiology) %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))

#Numbers for grade of each variants
germ9 %>% dplyr::filter(consented.to.research != "No", etiology != "not tested") %>% dplyr::group_by(etiology) %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))

#Number of patients with recurrent variant
length(unique(filter(germ9, consented.to.research != "No")$MRN))

#Number of recurrent variants
length((germ9 %>% dplyr::filter(consented.to.research != "No", etiology == "germline") %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE))$S.Nom)
length((germ9 %>% dplyr::filter(consented.to.research != "No", etiology != "germline") %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE))$S.Nom)

#number of patients with material available for additional testing
length(unique(filter(germ9, consented.to.research != "No", etiology != "not tested")$MRN))

#shows the names of the pathogenic variants that were confirmed
germ9 %>% dplyr::filter(consented.to.research != "No", path.level == "P", !is.na(varsome)) %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE)

#shows which were germline but grade 5 (basically false negatives)
germ9[which(germ9$grade.new == 5 & germ9$etiology == "germline"),] %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE)

#shows which were somatic but grade 1 (basically false positives)
germ9[which(germ9$grade == 1 & germ9$etiology == "somatic"),] %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE)

#plot the mutations for a single patient based on MRN
ggplot(data = subset(cancer, cancer$MRN  == "#######"), aes(x = time_point/30, y = AF, col = paste(mutation," ",variance))) + ylim(0,1) + geom_line() + geom_point(size = 3) + xlab("Months") + ylab("Allele Frequency") + ggtitle("Patient variants over time") + theme_classic() + scale_color_brewer(palette="Set1")

#NOW WILL DO THE ROC CURVES
#removes the confirmed variants that were not tested or reported on the germline panels
germ10 <- germ9 %>% dplyr::filter(etiology != "not tested", consented.to.research != "No")

#fixes the ID column to only include those consented to research
germ10 <- germ10 %>% dplyr::group_by(MRN) %>% dplyr::mutate(ID = cur_group_id())

#places a 1 in a column if germline and 0 if not
germ10$cohort.status <- dplyr::if_else(germ10$etiology == "germline", 1, 0)

#removes multiple entries from each patient
germ10 <- germ10 %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE)
germ10 %>% dplyr::group_by(etiology) %>% dplyr::summarise(count = n())


germ9 %>% dplyr::filter(consented.to.research != "No", path.level == "P", etiology == "germline") %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::summarise(count = n())

#Runs the ROC curve (Supplemental Figure 2B)
roc <- rocit(score = germ10$grade.new, class = germ10$cohort.status, method = "bin")
summary(roc)
print(ciAUC(roc, level = 0.95))
plot(roc, YIndex = T, values = T, col = c(2,4), level = 0.95)
plot(ciROC(roc, level = 0.95))
ksplot(roc)
message("KS Stat (binormal) : ", ksplot(roc)$`KS stat`)
message("KS Stat (binormal) cutoff : ", ksplot(roc)$`KS Cutoff`)
measureit(score = germ10$grade.new, class = germ10$cohort.status, measure = c("ACC", "SENS","SPEC","PPV","NPV","pDLR"))

#histogram of dates so you can see ranges of confirmed patients
germ9$confirmation.date <- as.Date(germ9$confirmation.date , "%m/%d/%Y")
ggplot(germ9 %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE), aes(x=confirmation.date)) + theme_bw() + geom_histogram(binwidth=15, fill="blue",color="black") + scale_x_date(date_breaks = "years", date_labels = "%d-%m-%Y")

germ9[which(germ9$Gene == "TP53"),] %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::filter(etiology != "somatic")

#shows which were germline but grade 5 (basically false negatives) for Supplemental Table 3
germ10[which(germ10$grade.new == 5 & germ10$etiology == "germline"),] %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::select(ID, MRN, Gene, S.Nom, grade)

germ10[which(germ10$MRN == "######"),] %>% dplyr::distinct(MRN, Gene, S.Nom, .keep_all = TRUE) %>% dplyr::select(ID, MRN, Gene, S.Nom, grade.new)

###NOW BEGIN ANALYSIS OF UNFILTERED ONCOPLUS DATA###

#Read in all the Unfiltered OncoPlus data
oncoplus <- fread("ClinVar_190624_211231.txt", fill= TRUE, showProgress = TRUE)
names(oncoplus)
oncoplus_update <- fread("TMB_190624_211231.txt", sep = "\t", header = TRUE)
names(oncoplus_update)
oncoplus_update2 <- fread("Clinical_170701_190621.txt", fill= TRUE, showProgress = TRUE)
clindata <- read.csv("CGL_MRN_190624_211231.csv", header = TRUE)
names(clindata) 
clindata_update <- read.csv("CGL_MRN_170701_190621.csv", header = TRUE)

#Loads the epic dataset with patient overall survival pulled 1/08/22 (Oncoplus_survival5)
osdata <- read.csv("Oncoplus_survival5.csv", header = TRUE)

#combine the above OncoPlus and clindata files into one dataframe with all the sequencing data
oncoplus <- rbind(oncoplus, oncoplus_update, oncoplus_update2)
clindata <- rbind(clindata, clindata_update)
clindata <- dplyr::distinct(clindata)

#Adds a column with the CGL number, take this out if needed
oncoplus$CGL <- gsub("S.*", "", oncoplus$SAMPLE)

#changes the Request.Date column to a date variable type
oncoplus <- transform(oncoplus, DATE = as.Date(as.character(DATE), "%y%m%d"))

#uses the dplyer package to perform the function on each patient and put it into a column called time_point
oncoplus <- oncoplus %>% dplyr::group_by(CGL) %>% dplyr::mutate(time_point = std(DATE))

#Ensures the time_point date difference is measured in an interger representing the number of days
oncoplus$time_point <- as.integer(oncoplus$time_point)

#some of the dates are 2-5 days apart and performed on the same sample. Below removes duplicate rows based on columns within the same sample and gNomen
oncoplus <- oncoplus %>% distinct(SAMPLE, gNomen, .keep_all = TRUE)
oncoplus$SAMPLE <- gsub("T.*", "", oncoplus$SAMPLE)

#add a column with the variance between mutations: use var(AF) if you want variance, and sd(AF)/mean(AF) if you want the coefficient of variation
oncoplus$AF <- as.numeric(oncoplus$AF)
which(is.na(oncoplus$AF))
oncoplus <- oncoplus %>% dplyr::filter(!is.na(AF))
oncoplus <- oncoplus %>% dplyr::group_by(paste(CGL, gNomen)) %>% dplyr::mutate(variance = sd(AF)/mean(AF))

#sets the variance to 1 for any variant that is not repeated in the data
oncoplus$variance <- oncoplus$variance %>% replace_na(1)

#runs the automatic grading function
oncoplus <- oncoplus %>% dplyr::group_by(CGL) %>% dplyr::mutate(grade = germline(variance))

#adds a column indicating if the gene is a common ChIP gene and is in BOTH papers, can edit to change to one paper if interested
oncoplus$ChIP <- oncoplus$gene %in% genelist$ChIP_NEJM & oncoplus$gene %in% genelist$ChIP_NatMed

#Adds a Column with the number of variants in each sample
oncoplus <- oncoplus %>% dplyr::group_by(SAMPLE) %>% dplyr::mutate(number_variants = length(unique(gNomen)))

#puts the DOB column in the correct format and creates a new column with the age in years
clindata$DOB <- as.Date(clindata$DOB, "%m/%d/%Y")
clindata$Order.Date <- mdy(clindata$Order.Date)
clindata$Sequence.Date <- mdy(clindata$Sequence.Date)

#Creates a function to overcome the date characterization in R
DOB <- function(x, year=2018){
  m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900+m, 2000+m)
  x
}

#normalizes the date to overcome the strange date characterization in R
clindata$DOB <- DOB(clindata$DOB)

#creates a column with patient age (in days)
clindata$age <- as.numeric(clindata$Order.Date - clindata$DOB)

#creates a new Sample column from the order number column to the CGL Sample Number so that the dataframes can be combine
clindata$CGL.y <- gsub(".T.*", "", clindata$Order)
clindata$Sample <- gsub(" ", "", paste(clindata$CGL.y, clindata$S.))

#combines the clinical data with the sequencing data by CGL number
oncoplus <- merge(oncoplus, clindata, by.x = "SAMPLE", by.y = "Sample", all.x = TRUE, all.y = FALSE)

#fills in the variants without a primary tumor site designation with "Not Designated"
oncoplus$OncoTree.Primary.Tumor.Site <- sub("^$", "Not Designated", oncoplus$OncoTree.Primary.Tumor.Site)
oncoplus$OncoTree.Primary.Tumor.Site <- oncoplus$OncoTree.Primary.Tumor.Site %>% replace_na("Not Designated")
oncoplus$OncoTree.Diagnosis <- sub("^$", "Not Designated", oncoplus$OncoTree.Diagnosis)
oncoplus$OncoTree.Diagnosis <- oncoplus$OncoTree.Diagnosis %>% replace_na("Not Designated")

#puts a value in the cells that are not in ClinVar
class(oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- as.character(oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*Benign*", "Benign", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*benign*", "Likely benign", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*Other*", "Benign", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*Pathogenic*", "Pathogenic", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub("*Pathogenic*.", "Pathogenic", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*pathogenic*", "Likely pathogenic", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*pathogenicity", "VUS", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*VUS", "VUS", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*Not*", "Not ", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- gsub(".*Risk*", "Risk", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- sub("^$", "Not in ClinVar", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs <- sub(".*provided*", "Not in ClinVar", oncoplus$clinVarClinSignifs)
oncoplus$clinVarClinSignifs[grepl("pathogenic", oncoplus$clinVarClinSignifs)] <- "Likely pathogenic"
oncoplus$clinVarClinSignifs[grepl("Pathogenic", oncoplus$clinVarClinSignifs)] <- "Pathogenic"
oncoplus$clinVarClinSignifs[grepl("Drug response", oncoplus$clinVarClinSignifs)] <- "Drug response"
oncoplus$clinVarClinSignifs[grepl("Risk factor", oncoplus$clinVarClinSignifs)] <- "Risk factor"
oncoplus$clinVarClinSignifs[grepl("VUS", oncoplus$clinVarClinSignifs)] <- "VUS"
oncoplus$clinVarClinSignifs[grepl("Not in ClinVar", oncoplus$clinVarClinSignifs)] <- "Not in ClinVar"
oncoplus$clinVarClinSignifs <- as.factor(oncoplus$clinVarClinSignifs)

#fills in the empty rows if the variant is non-coding
oncoplus$codingEffect <- sub("^$", "non-coding", oncoplus$codingEffect)

#adds a true/false column if the primary site is blood/bone marrow or lymph
oncoplus$in.blood <- if_else(oncoplus$OncoTree.Primary.Tumor.Site == "Blood/Bone Marrow" | oncoplus$OncoTree.Primary.Tumor.Site == "Lymph", TRUE, FALSE)

#adds a column with the inheritance according to Feurstein SK, Trottier AM, Estrada-Merly N, et al. Germline predisposition variants occur in myelodysplastic syndrome patients of all ages [published online ahead of print, 2022 Aug 19]. Blood. 2022;blood.2022015790
oncoplus <- merge(oncoplus, genelist[,c("Germline_Blood","Inheritance_Blood")], by.x = "gene", by.y = "Germline_Blood", all.x = TRUE, all.y = FALSE)

#Removes synonymous variants (although they may be germline)
oncoplus <- oncoplus %>% filter(codingEffect != "synonymous")

#creates a dataframe of patients who have >= 2 samples in the dataset
reocc <- oncoplus %>% dplyr::group_by(CGL, gNomen) %>% dplyr::summarize(count = n())
reocc <- reocc %>% dplyr::group_by(CGL) %>% dplyr::summarize(count = n())
reocc <- subset(reocc, reocc$count >= 2)

#creates a column in the oncoPlus dataframe indicating if the patient has repeated tests
oncoplus$reocc <- oncoplus$CGL %in% reocc$CGL

#adds column if the MRN is in our OncoPlus dataset
oncoplus$in.osdata <- oncoplus$MRN %in% osdata$MRN

#removes any row that is not in blood/bone marrow or lymph node
oncoplus_heme <- oncoplus %>% dplyr::filter(in.blood == TRUE)

#makes a new dataframe with data from with only patients with reoccurring oncoplus testing
oncoplus2 <- subset(oncoplus_heme, oncoplus$reocc == "TRUE")

#creates a dataframe of patients who have >= 2 variants that repeated in the dataset 
reocc <- oncoplus2 %>% dplyr::group_by(CGL, gNomen) %>% dplyr::summarize(count = n())
reocc <- subset(reocc, reocc$count >= 2)

#creates a column in the germ data frame indicating if the patient has recurring mutations in the same gene
oncoplus2$reocc <- oncoplus2$CGL %in% reocc$CGL
oncoplus2 <- subset(oncoplus2, oncoplus2$reocc == "TRUE")

#adds a column with the number of unique mutations per CGL number
oncoplus2 <- merge(oncoplus2, reocc, by=c("CGL","gNomen"))

#makes a new dataframe from the one with reoccuring mutations based on specific cutoffs below
oncoplus3 <- oncoplus2 %>% dplyr::filter(grade == 1 | grade == 2)
oncoplus3 <- oncoplus3 %>% dplyr::filter(gnomadAltFreq_all < 0.0001 | is.na(gnomadAltFreq_all))

#Prints summary numbers for the analysis
length(unique(oncoplus_heme$CGL)) #number of total patients
length(unique(oncoplus2$CGL)) #number of patients with multiple samples
length(unique(oncoplus_heme$SAMPLE)) #number of samples
length(unique(oncoplus_heme$OncoTree.Diagnosis)) #number of diagnosis
length(unique(oncoplus_heme$gene)) #number of unique genes
length(unique(oncoplus_heme$gNomen)) #number of unique variants

length(unique(oncoplus$gene))

#tells us the number of patients (CGL numbers) and samples
length(unique(oncoplus3$CGL))
length(unique(oncoplus3$CGL))/length(unique(oncoplus_heme$CGL))
length(unique(oncoplus3$SAMPLE))
length(unique(oncoplus3$SAMPLE))/length(unique(oncoplus_heme$SAMPLE))
length(unique(oncoplus3$gene))
length(unique(oncoplus3$gene))/length(unique(oncoplus_heme$gene))
length(unique(oncoplus3$gNomen))
length(unique(oncoplus3$gNomen))/length(unique(oncoplus_heme$gNomen))

#Uses stringr package to pull out the MRN into its own column
osdata$MRN <- stringr::str_extract(string = osdata$Patient.Name..MRN., pattern = "(?<=\\().*(?=\\))")
osdata$MRN <- as.integer(osdata$MRN)

#makes the date columns actually dates using the package lubridate
osdata <- transform(osdata, OrderDate = as.Date(as.character(Order.Date), "%m/%d/%Y"), BirthDate = as.Date(as.character(Date.of.Birth), "%m/%d/%Y"), DeathDate = as.Date(as.character(DeathDate), "%m/%d/%Y"))

#Adds a column for if a patient is alive or not, 1=alive, 2=dead
osdata$alive <- ifelse(is.na(osdata$DeathDate), 1, 2)

#Adds a column with patients age when the data was pulled or age when they died
osdata$age <- ifelse(is.na(osdata$DeathDate), interval(osdata$BirthDate, "2021-06-25") %/% months(1), interval(osdata$BirthDate, osdata$DeathDate) %/% months(1))

#adds column if the MRN is in our OncoPlus dataset
osdata$in.oncoplus <- osdata$MRN %in% oncoplus$MRN

#removes the patients that we don't have oncoplus data for
osdata <- subset(osdata, osdata$in.oncoplus == TRUE)
osdata <- subset(osdata, select = -in.oncoplus)

#adds column if the MRN is in our oncoplus3 dataset or has a potential germline variant
osdata$in.germline <- osdata$MRN %in% oncoplus3$MRN

#adds column if the MRN is pathogenic and in our oncoplus3 dataset
osdata$in.path.germline <- osdata$MRN %in% (oncoplus3 %>% dplyr::filter(clinVarClinSignifs == "Pathogenic" | clinVarClinSignifs == "Likely pathogenic"))$MRN

#adds column if the MRN associated with a particular tumor site in the oncoplus data
osdata$in.blood <-  osdata$MRN %in% (as.data.frame(oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(OncoTree.Primary.Tumor.Site == "Blood/Bone Marrow" | OncoTree.Primary.Tumor.Site == "Lymph"))$MRN)

#adds a column if there is a patient has a variant in a specific gene as specified: these two may be interesting (GRIN2A, KMT2A) "ATRX"   "ERBB2"  "MLH1"   "ARID1A" "ATM" "ERBB4"  "RUNX1" "SF3B1"  "FANCA"  "KDM6A"  "NCOA2" 
osdata$in.gene <- dplyr::if_else(osdata$MRN %in% (as.data.frame(oncoplus3 %>% group_by(CGL) %>% dplyr::filter(gene == "DDX41"))$MRN), "DDX41", "")

osdata$in.gene <- dplyr::if_else(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(gene == "CHEK2"))$MRN), paste(osdata$in.gene, "CHEK2"), osdata$in.gene)

osdata$in.gene <- dplyr::if_else(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(gene == "RUNX1"))$MRN), paste(osdata$in.gene, "RUNX1"), osdata$in.gene)

#adds a column with the germline grade of the variant identified
osdata$in.grade <- ifelse(osdata$MRN %in% (oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 1))$MRN, 1, "")
osdata$in.grade <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 2))$MRN), paste(osdata$in.grade, 2), osdata$in.grade)
osdata$in.grade <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 3))$MRN), paste(osdata$in.grade, 3), osdata$in.grade)
osdata$in.grade <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 4))$MRN), paste(osdata$in.grade, 4), osdata$in.grade)
osdata$in.grade <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 5))$MRN), paste(osdata$in.grade, 5), osdata$in.grade)

osdata$in.grade <- ifelse(is.na(osdata$in.grade), ifelse(osdata$MRN %in% (oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 2))$MRN, 2, NA), osdata$in.grade)
osdata$in.grade <- ifelse(is.na(osdata$in.grade), ifelse(osdata$MRN %in% (oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 3))$MRN, 3, NA), osdata$in.grade)
osdata$in.grade <- ifelse(is.na(osdata$in.grade), ifelse(osdata$MRN %in% (oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 4))$MRN, 4, NA), osdata$in.grade)
osdata$in.grade <- ifelse(is.na(osdata$in.grade), ifelse(osdata$MRN %in% (oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(grade == 5))$MRN, 5, NA), osdata$in.grade)

#creates a column separating which grade of germline variants are in each patient for figures | grepl("3",osdata$in.grade)
osdata$in.grade1 <- ifelse(grepl("1",osdata$in.grade) | grepl("2",osdata$in.grade),"Grade_12","Grade_345")

#adds column with true/false, because the grade scale computes the HR using the germline as the control
osdata$in.grade2 <- ifelse(grepl("5",osdata$in.grade1) ,FALSE,TRUE)

#adds a column if there is a clinical significant field in the filtered patient
osdata$in.pathogenic <-  osdata$MRN %in% (as.data.frame(oncoplus %>% dplyr::group_by(CGL) %>% dplyr::filter(str_detect(clinVarClinSignifs, "Pathogenic")))$MRN)

#merges the osdata with the clindata
osdata <- merge(osdata, clindata[,c("MRN","Sex","OncoTree.Specimen.Site","OncoTree.Primary.Tumor.Site","OncoTree.Diagnosis","Percent.Tumor")], by.x = "MRN", by.y = "MRN", all.x = TRUE, all.y = FALSE)

#adds a column with the inheritance if there is a known variant
osdata$inheritance <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(!is.na(Inheritance_Blood), Inheritance_Blood == "AR"))$MRN), paste(osdata$inheritance, "AR"), NA)
osdata$inheritance <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(!is.na(Inheritance_Blood), Inheritance_Blood == "AD/AR"))$MRN), paste(osdata$inheritance, "AD"), osdata$inheritance)
osdata$inheritance <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(!is.na(Inheritance_Blood), Inheritance_Blood == "AD"))$MRN), paste(osdata$inheritance, "AD"), osdata$inheritance)
osdata$inheritance <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus3 %>% dplyr::group_by(CGL) %>% dplyr::filter(!is.na(Inheritance_Blood), Inheritance_Blood == "X-linked"))$MRN), paste(osdata$inheritance, "X-linked"), osdata$inheritance)

osdata[which(osdata$in.grade3 == "Grade_12"),c("MRN", "in.blood")]

#plot the mutations for a single patient based on MRN normalized date on x-axis ###### 
ggplot(data = subset(oncoplus, oncoplus$MRN  == "######"), aes(x = time_point/30, y = AF, col = paste(gene," ",gNomen," ",variance))) + ylim(0,1) + geom_line(show.legend = FALSE) + geom_point(size = 3, show.legend = FALSE) + xlab("Months") + ylab("Allele Frequency") + ggtitle("Patient variants over time") + theme_classic() + scale_color_brewer(palette="Set1")

oncoplus[which(oncoplus$MRN == "######"),c("gene","gNomen","grade","Inheritance_Blood","time_point")] %>% dplyr::distinct(gNomen, .keep_all = TRUE)

#creates a column with the inheritance and grade together for the figure)
osdata$in.grade3 <- dplyr::if_else(osdata$in.grade1 == "Grade_12" & grepl("AD",osdata$inheritance),"Grade_12_AD", osdata$in.grade1)

#removed all but the earliest date
osdata %<>% dplyr::group_by(MRN) %>% dplyr::filter(OrderDate == min(OrderDate))
osdata <- osdata[!duplicated(osdata$MRN), ]

#creates a landmark column of time difference from earliest test to when the status of each patient was collected
osdata$landmark <- ifelse(is.na(osdata$DeathDate), interval(osdata$OrderDate, "2021-12-31") %/% months(1), interval(osdata$OrderDate, osdata$DeathDate) %/% months(1))

#histogram of dates so you can see ranges of tests
ggplot(osdata %>% filter(in.blood == TRUE), aes(x=OrderDate)) + theme_bw() + geom_histogram(binwidth=15, fill="blue",color="black") + scale_x_date(date_breaks = "years", date_labels = "%d-%m-%Y")

#Makes a more comprehensive figure of landmark analysis from time of first OncoPlus test for presentations (Figure 2B)
ggsurvplot(survfit(Surv(landmark, alive) ~ in.grade3, data= osdata %>% dplyr::filter(in.blood == TRUE)),
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Grade 3/4/5", "AD + Grade 1/2", "Grade 1/2"),  
           xlab="Time (months)",
           palette=c("blue", "red", "dark green"), 
           title="Time-Dependent Covariate", 
           risk.table.height=.15)

coxph(Surv(landmark, alive) ~ in.grade3, data = osdata %>% dplyr::filter(in.blood == TRUE)) %>% gtsummary::tbl_regression(exp = TRUE)

#input the year of survival that you want to document
summary(survfit(Surv(landmark, alive) ~ in.grade3, data = osdata %>% dplyr::filter(in.blood == TRUE)), times = 24)

#plots the data as a cummulative survival plot (Figure 2D)
ggsurvplot(survfit(Surv(age/12, alive) ~ in.grade3, data= osdata %>% dplyr::filter(in.blood == TRUE)), conf.int = FALSE, pval = TRUE,
           palette = c("blue", "red", "dark green"),
           risk.table = TRUE, risk.table.col = "strata",
           fun = "cumhaz")

#input the year of survival that you want to document
summary(survfit(Surv(age/12, alive) ~ in.grade3, data = osdata %>% dplyr::filter(in.blood == TRUE)), times = 75)

#changes the levels of factor for in.grade3 so that the reference is the patients without grade 1/2 variants
osdata$in.grade3 <- factor(osdata$in.grade3, levels = c("Grade_345","Grade_12_AD","Grade_12"))

#Gives us the Hazard Ratio, CI, and P-value
coxph(Surv(age/12, alive) ~ in.grade3, data = osdata %>% dplyr::filter(in.blood == TRUE)) %>% gtsummary::tbl_regression(exp = TRUE)

#tells us the number of patients with survival data available
length((osdata %>% dplyr::filter(in.blood == TRUE) %>% dplyr::distinct(MRN, .keep_all = TRUE))$MRN)

#tells us the number of patients with survival data available in the germline data
length((osdata %>% dplyr::filter(in.blood == TRUE, in.grade2 == TRUE) %>% dplyr::distinct(MRN, .keep_all = TRUE))$MRN)

#tells us the number of patients with survival data available in the germline data and that have a AD or AR/AD variant
length((osdata %>% dplyr::filter(in.blood == TRUE, in.grade3 == "Grade_12_AD") %>% dplyr::distinct(MRN, .keep_all = TRUE))$MRN)

#tells us the number of grade 1 variants with survival data available
length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, grade == 1) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)
length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, grade == 1) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)/length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)

#tells us the number of grade 2 variants with survival data available
length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, grade == 2) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)
length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, grade == 2) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)/length((oncoplus %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE) %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)

#provides a list of AD gene variants that are included in the kaplan-meier curve
oncoplus4 <- oncoplus3 %>% dplyr::filter(in.osdata == TRUE, !is.na(Inheritance_Blood), Inheritance_Blood != "AR", grade == 1 | grade == 2) %>% dplyr::distinct(MRN, gNomen, .keep_all = TRUE) %>% select(CGL, MRN, gene, cNomen, pNomen, grade, Inheritance_Blood)
oncoplus4 <- oncoplus4 %>% group_by(MRN) %>% dplyr::mutate(ID = cur_group_id())

table(osdata$in.grade)
table(osdata$in.grade3)
table(oncoplus4$Inheritance_Blood)
length(unique(oncoplus4$MRN))
oncoplus4$MRN <- as.integer(oncoplus4$MRN)

#tells the patients that are not shared between datasets -> should be NONE
dplyr::setdiff(as.data.frame(osdata %>% dplyr::filter(in.grade3 == "Grade_12_AD"))$MRN, oncoplus4$MRN)

#tells the number of grade 1 and 2 variants in the pipeline caller
length((oncoplus4 %>% dplyr::distinct(CGL, cNomen, .keep_all = TRUE))$pNomen) #number of AD variants
length((oncoplus3 %>% dplyr::filter(in.osdata == TRUE) %>% dplyr::distinct(CGL, cNomen, .keep_all = TRUE))$pNomen) #denominator of called variants
length((oncoplus4 %>% dplyr::distinct(gene, .keep_all = TRUE))$gene)
length((oncoplus %>% dplyr::distinct(gene, .keep_all = TRUE))$gene)

#Creates Supplemental Table 3 from the probable variants in the Kaplan-Meier
oncoplus4 <- oncoplus4 %>% dplyr::arrange(ID)
write.csv(oncoplus4 %>% dplyr::select(ID, gene, cNomen, pNomen, grade, Inheritance_Blood), "Supplemental Table 3.csv")

#makes a dataframe with potential germline variants that are NOT autosomal dominant
oncoplus4_noAD <- oncoplus3 %>% dplyr::filter(in.osdata == TRUE, Inheritance_Blood == "AR" | is.na(Inheritance_Blood), grade == 1 | grade == 2) %>% dplyr::distinct(MRN, gNomen, .keep_all = TRUE) %>% select(CGL, MRN, gene, cNomen, pNomen, grade, Inheritance_Blood)
oncoplus4_noAD <- oncoplus4_noAD %>% group_by(MRN) %>% dplyr::mutate(ID = cur_group_id())
oncoplus4_noAD <- oncoplus4_noAD %>% dplyr::arrange(ID)

write.csv(oncoplus4_noAD %>% dplyr::select(ID, gene, cNomen, pNomen, grade, Inheritance_Blood), "Supplemental Table 3 noAD.csv")

#tells the patients that are not shared between datasets -> should be NONE
dplyr::setdiff(as.data.frame(osdata %>% dplyr::filter(in.grade3 == "Grade_12"))$MRN, oncoplus4_noAD$MRN)

#read in the manually annotated table to add back to Oncoplus4 dataframe
oncoplus4_anno <- read.csv(file = "Supp.Table.3.Annotated.csv", header = TRUE)

#add the manual annotation to OncoPlus4 making Oncoplus5
oncoplus5 <- merge(oncoplus4, oncoplus4_anno[,c("ID","gene","pNomen","ClinVar","varsome","path.level")], by=c("ID","gene","pNomen"), all.x = TRUE, all.y = FALSE)

#
#shows a table of the number of annotated variants
oncoplus5 %>% dplyr::group_by(path.level) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))  %>% as.data.frame() %>% arrange(desc(freq)) 

#lists the genes that had grade 1 or 2 germline variants
oncoplus3 %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, grade == 1 | grade == 2) %>% dplyr::distinct(gene) %>% unlist(recursive = FALSE, use.names = FALSE) %>% print(quote=FALSE, row.names = FALSE) 

#makes another column with the manual annotation of pathogenic potential germline variants in the os dataset
osdata$in.path.germline <-  ifelse(osdata$MRN %in% (oncoplus5 %>% dplyr::group_by(CGL) %>% dplyr::filter(path.level == "P"))$MRN, "P", "")
osdata$in.path.germline <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus5 %>% dplyr::group_by(CGL) %>% dplyr::filter(path.level == "VUS"))$MRN), paste(osdata$in.path.germline, "VUS"), osdata$in.path.germline)
osdata$in.path.germline <- ifelse(osdata$MRN %in% (as.data.frame(oncoplus5 %>% dplyr::group_by(CGL) %>% dplyr::filter(path.level == "B"))$MRN), paste(osdata$in.path.germline, "B"), osdata$in.path.germline)
osdata$in.path.germline <- ifelse(grepl("P",osdata$in.path.germline) ,"P",osdata$in.path.germline)
osdata$in.path.germline <- ifelse(grepl("VUS",osdata$in.path.germline) ,"VUS",osdata$in.path.germline)

#tells us the number of pathogenic variants with survival data available
length((oncoplus3 %>% dplyr::filter(in.blood == TRUE, in.osdata == TRUE, clinVarClinSignifs == "Pathogenic" | clinVarClinSignifs == "Likely pathogenic") %>% dplyr::distinct(cNomen, .keep_all = TRUE))$cNomen)

#plots a quick kaplan-meier with identifying if a patient has a gene associated with OS
ggsurvplot(survfit(Surv(landmark, alive) ~  in.gene + in.grade1, data= osdata %>% dplyr::filter(in.blood == TRUE)),
           conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           palette="Set1",
           xlab="Time (months)",
           title="Time-Dependent Covariate", 
           risk.table.height=.15)

summary(survfit(Surv(landmark, alive) ~ in.gene + in.grade1, data = subset(osdata, osdata$in.blood == TRUE)), times = 18)

#plots a quick kaplan-meier with pathogenic potential germline variant as surrogate for OS
ggsurvplot(survfit(Surv(landmark, alive) ~  in.path.germline + in.grade3, data= osdata %>% dplyr::filter(in.blood == TRUE, in.grade3 == "Grade_12_AD")),
           conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           palette="Set1",
           xlab="Time (months)",
           title="Time-Dependent Covariate", 
           risk.table.height=.15)


coxph(Surv(landmark, alive) ~ in.path.germline + in.grade3, data = osdata %>% dplyr::filter(in.blood == TRUE)) %>% gtsummary::tbl_regression(exp = TRUE)

#Supplemental Figure 1A and 1B
CSF3R <- subset(cancer, cancer$Gene == "CSF3R")

#loads the dataset with the CSF3R confirmed germline status
CSF3R_germ <- read.csv(file = "CSF3R_germline.csv", header = TRUE)

#creates a data table of patients who had the same gene pop up >= 2 times
CSF3R_iv <- CSF3R %>% dplyr::group_by(Patient, MRN, S.Nom, P.Nom, mutation, variance) %>% dplyr::summarize(count = n())
CSF3R_iv <- CSF3R_iv %>% dplyr::filter(count >= 2)
CSF3R_iv <- merge(x = CSF3R_iv, y = CSF3R_germ[ , c("S.Nom", "cohort.status")], by.x = "S.Nom", by.y = "S.Nom", all.x=TRUE)

#creates a column in the CSF3R dataframe indicating if the patient has recurring mutations in the same gene
CSF3R$reocc <- CSF3R$MRN %in% CSF3R_iv$MRN

#creates a column in the CSF3R dataframe indicating if the mutation has been confirmed as germline or not
CSF3R <- merge(x = CSF3R, y = CSF3R_germ[ , c("P.Nom.R", "cohort.status")], by.x = "P.Nom", by.y = "P.Nom.R", all.x=TRUE)

#add a column with the average AF of mutations
CSF3R <- CSF3R %>% group_by(gene_snv) %>% mutate(mean_AF = max(AF))

#box plot of germline, somatic, and unconfirmed variants: Supplemental Figure 2A
ggplot(CSF3R %>% dplyr::filter(reocc == TRUE & !is.na(cohort.status)) %>% dplyr::distinct(MRN, gene_snv, .keep_all = TRUE) %>% dplyr::select(gene_snv, mean_AF, variance, cohort.status), aes(cohort.status, variance, fill = cohort.status)) + geom_boxplot(alpha = 0.5) + geom_quasirandom(alpha = 0.5) + xlab("Etiology") + ylab("Variance") + scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Numbers for each
CSF3R %>% dplyr::filter(reocc == TRUE & !is.na(cohort.status)) %>% dplyr::group_by(cohort.status) %>% dplyr::distinct(MRN, gene_snv, .keep_all = TRUE) %>% dplyr::summarise(count = n()) %>% dplyr::mutate(freq = sum(count))

#Provides the P-value
kruskal.test(variance ~ cohort.status, data = CSF3R %>% dplyr::filter(reocc == TRUE & !is.na(cohort.status)) %>% dplyr::group_by(cohort.status) %>% dplyr::distinct(MRN, gene_snv, .keep_all = TRUE))

#Now generate a ROC curve with sensitivity and specificity using variance, you can do it with 3 different stats empirical (emp), binormal (bin), or non-parametric (non)
CSF3R_iv$cohort.status2 <- ifelse(CSF3R_iv$cohort.status == "Germline", 0, 1)
measureit(score = CSF3R_iv$variance, class = CSF3R_iv$cohort.status2, measure = c("ACC", "SENS","SPEC","PPV","NPV","pDLR"))
roc <- rocit(score = CSF3R_iv$variance, class = CSF3R_iv$cohort.status2, method = "bin")
summary(roc)
print(ciAUC(roc, level = 0.95))

#Supplemental Figure 2B
plot(roc, YIndex = F, values = F, col = c(2,4))

#makes a temporary dataframe with the gene counts and a column with true/false if it is in the list of known germline predisposition genes
temp <- oncoplus3 %>% dplyr::filter(in.osdata == TRUE) %>% dplyr::distinct(CGL, cNomen, .keep_all = TRUE) %>% dplyr::group_by(gene) %>% dplyr::summarize(count = n())
temp$known <- dplyr::if_else(temp$gene %in% genelist$Germline_Blood, "Germline", NA_character_)
temp$known <- dplyr::if_else(temp$gene %in% genelist$Somatic_Blood, paste(temp$known, "Somatic"), temp$known)
temp$known <- dplyr::if_else(is.na(temp$known), "Not Associated", temp$known)
temp <- merge(temp, genelist[,c("Germline_Blood","Inheritance_Blood")], by.x = "gene", by.y = "Germline_Blood", all.x = TRUE, all.y = FALSE)

#write.csv(temp, "Grade1_2_variants.csv")

#Graphs the number of germline mutations via bar graph by gene highlighting the known germline genes (not shown are genes with less than 5 variants): Figure 2A
ggplot(temp %>% dplyr::filter(count > 3), aes(x= reorder(gene, -count), y=count, fill = known)) + geom_bar(stat = "identity", show.legend = TRUE) + xlab("Gene") + ylab("Number of Mutations") + scale_fill_brewer(palette="Set1") + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8))

as.data.frame(temp[order(temp$count, decreasing = TRUE),])
describeBy(temp, group = temp$known)
sum((temp %>% filter(grepl("Germline",temp$known)))$count)
sum(temp$count)
sum((temp %>% filter(grepl("Germline",temp$known)))$count)/sum(temp$count)
sum(grepl("Germline",temp$known))
sum(grepl("Germline",temp$known))/length(temp$known)
#tells us the expected percentage based on number of genes sequenced and number of germline mutations
length(unique(genelist$Germline_Blood))/length(unique(oncoplus$gene))

#Graphs only germline genes via bar graph by gene highlighting the inheritance (Supplemental Figure 3A)
ggplot(temp %>% dplyr::filter(!is.na(Inheritance_Blood)), aes(x= reorder(gene, -count), y=count, fill = Inheritance_Blood)) + geom_bar(stat = "identity", show.legend = TRUE) + xlab("Gene") + ylab("Number of Mutations") + scale_fill_brewer(palette="Paired") + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8))



#makes a temporary dataframe with the gene counts and a column with true/false if it is in the list of ChIP genes: Supplemental Figure 3A
temp <- oncoplus3 %>% dplyr::filter(in.osdata == TRUE) %>% dplyr::distinct(CGL, cNomen, .keep_all = TRUE) %>% dplyr::group_by(gene) %>% dplyr::summarize(count = n())
temp$chip <- temp$gene %in% genelist$ChIP_NEJM | temp$gene %in% genelist$ChIP_NatMed

#Graphs the number of germline mutations via bar graph by gene highlighting the known CHIP genes
ggplot(temp %>% filter(count > 3), aes(x= reorder(gene, -count), y=count, fill = chip)) + geom_bar(stat = "identity", show.legend = TRUE) + xlab("Gene") + ylab("Number of Mutations") + scale_fill_manual(values = c("grey40","blue")) + theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8))

as.data.frame(temp[order(temp$count, decreasing = TRUE),])
sum((temp %>% filter(chip == TRUE))$count)
sum((temp %>% filter(chip == TRUE))$count)/sum(temp$count)
sum(temp$chip)
sum(temp$chip)/length(temp$chip)
#tells us the expected percentage based on number of genes sequenced and number of ChIP mutations
length(unique(unlist(genelist$ChIP_NEJM, genelist$ChIP_NatMed)))/length(unique(oncoplus$gene))

#Exports Supplemental Table 2 from out data
names(germ10)
germ_manuscript <- germ10 %>% dplyr::select(ID, Gene, S.Nom, P.Nom, count, variance, grade.new, etiology, ClinVar, varsome, path.level)
germ_manuscript <- germ_manuscript %>% dplyr::arrange(ID)
write.csv(germ_manuscript, "Supplemental Table 2.csv")
