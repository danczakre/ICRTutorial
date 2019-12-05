### Script to generate Van Krevelen plots; should be able to be converted to a function easily
# This script must be run after FTMS_Analysis.R/FREDA or else it will not work correctly
# by RED

# Options
freq_thresh = 0 # Adjusts the amount of data which is plotted based upon occurence frequency
                # 0 will plot all data; 100 will plot data that is present in 100% of the samples
samp.num = 1 # This option will plot the Van Krevelen for the specified column/sample

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


# ###################### #
#### Plotting results ####
# ###################### #

# Adjusting data based upon frequency threshold
freq = (rowSums(data)/ncol(data))*100 # Finding frequnecies
freq.data = data[which(freq > freq_thresh),] # Creating data/mol objects based upon frequency
freq.mol = mol[which(freq > freq_thresh),]

# Plotting all data
p = ggplot(data = freq.mol, aes(x = OtoC_ratio, y = HtoC_ratio)) +
  geom_point(aes(color = bs2_class)) + xlab("O:C") + ylab("H:C") + 
  ggtitle(paste0("All Data - ", freq_thresh, "% threshold")) +
  labs(color = "Bailey et al, 2017 - Boundary Set") + xlim(c(0, max(mol$OtoC_ratio))) + 
  ylim(c(0, max(mol$HtoC_ratio))) + theme_bw() + theme(axis.text = element_text(color = "black"))

ggsave(paste0("All Data - ", freq_thresh, " threshold.pdf"), plot = p, width = 9, height = 7)

# Plotting one sample based upon selection
if(samp.num > ncol(freq.data)){
  stop("You selected an incorrect sample number (i.e., too high).")
} # Checks to make sure that the sample number is valid

samp.mol = freq.mol[which(freq.data[,samp.num] > 0),]

p = ggplot(data = samp.mol, aes(x = OtoC_ratio, y = HtoC_ratio)) +
  geom_point(aes(color = bs2_class)) + xlab("O:C") + ylab("H:C") + 
  ggtitle(paste0(colnames(freq.data)[samp.num], " - ", freq_thresh, "% threshold")) +
  labs(color = "Bailey et al, 2017 - Boundary Set") + xlim(c(0, max(mol$OtoC_ratio))) + 
  ylim(c(0, max(mol$HtoC_ratio))) + theme_bw() + theme(axis.text = element_text(color = "black"))

ggsave(paste0(colnames(freq.data)[samp.num], " - ", freq_thresh, " threshold.pdf"), 
       plot = p, width = 9, height = 7)
