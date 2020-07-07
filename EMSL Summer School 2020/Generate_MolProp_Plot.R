### Script to generate molecular property plots; should be able to be converted to a function easily
# This script must be run after FTMS_Analysis.R/FREDA or else it will not work correctly
# by RED

# Options
property = "Mass" # Specify which property distribution you want to track through time
                  # Valid options are "Mass", "NOSC", "AI", "AI_Mod", DBE", "DBE_O", "DBE_AI", "GFE"

# Load necessary packages
require(reshape2); require(ggplot2); require(easycsv)


# ################## #
#### Load in data ####
# ################## #

# Set directory
setwd(easycsv::choose_dir())

# Load files
data = read.csv(list.files(pattern = "Processed.*Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "Processed.*Mol.csv"), row.names = 1)

# Ordering data
data = data[,order(colnames(data))]


# ################## #
#### Error Checks ####
# ################## #

if(!identical(row.names(data), row.names(mol))){
  stop("Your data file and mol file do not match. Please address this.")
} # If the row names for your data and mol files do not match, the incorrect files are loading

if(length(which(colnames(mol) %in% "bs1_class")) == 0){
  stop("It seems like FTMS_Analysis.R or the FREDA app has not been run.")
} # This script requires outputs from the ftmsRanalysis or FREDA


# ################### #
#### Cleaning data ####
# ################### #

# Setting data to presence/absence
data[data > 0] =  1 # Given that intensities are not necessarily related to 
# concentration, we set values to presence/absence

# Removing data without molecular formula
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Adding a m/z column
mol$Mass = as.numeric(row.names(mol))


# ################### #
#### Data analysis ####
# ################### #

# Plotting the distribution of one property
sing.prop = NULL

for(i in 1:ncol(data)){
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = data.frame(Samples = colnames(data)[i], Property = mol[row.names(temp), property]) # Subsetting the mol file

  sing.prop = rbind(sing.prop, temp)
} # Creating an object for the single propertys

p = ggplot(sing.prop, aes(x = Samples, y = Property)) +
  geom_boxplot() + xlab(NULL) + ylab(property) +
  ggtitle(paste0(property, " Distribution Plot")) +
  theme_bw() + theme(axis.text = element_text(color = "black"),
                     axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5))

ggsave(paste0(property, " Distribution Plot.pdf"), plot = p, width = 4, height = 3)

# Creating empty object to store average molecular property values
mol.info = data.frame(NOSC = rep(NA, ncol(data)), AI = NA, DBE = NA, N = NA, S = NA, P = NA, Peaks = NA, 
                      row.names = colnames(data), stringsAsFactors = F) # Creating an object for mean properties

for(i in 1:ncol(data)){
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = mol[row.names(temp),] # Subsetting the mol file
  
  mol.info$NOSC[i] = mean(temp$NOSC, na.rm = T)
  mol.info$AI[i] = mean(temp$AI, na.rm = T)
  mol.info$DBE[i] = mean(temp$DBE, na.rm = T)
  mol.info$N[i] = mean(temp$N, na.rm = T)
  mol.info$S[i] = mean(temp$S, na.rm = T)
  mol.info$P[i] = mean(temp$P, na.rm = T)
  mol.info$Peaks[i] = length(temp$P)
} # This for-loop determines the average molecular property values by sample

mol.info = melt(as.matrix(mol.info)) # Converting data to long format for better plotting in ggplot

# Plotting bar plots
p = ggplot(mol.info, aes(x = Var1, y = value)) +
  geom_bar(stat = "identity", aes(fill = Var2)) + facet_grid(Var2~., scales = "free_y") +
  xlab(NULL) + ggtitle("Molecular Characterisitics by Sample") +
  theme_bw() + theme(axis.text = element_text(color = "black"),
                     axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
                     legend.position = "none")

ggsave("Molecular Characteristics by Sample.pdf", plot = p, width = 5, height = 9)
