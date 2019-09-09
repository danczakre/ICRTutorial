### Identifying biochemical transformations associated with compounds detected by FTICR
# RED 2019

require(dplyr); require(tidyr); require(ggplot2)

options(digits=10) # Sig figs in mass resolution data

Sample_Name = "HJ_Andrews_Tutorial"

#######################
### Loading in data ###
#######################

# Loading in ICR data
setwd(easycsv::choose_dir())
data = read.csv("Processed_HJ_Andrews_Tutorial_Data.csv", row.names = 1) # Keeping data and mol-data seperate to ensure they are unaltered
mol = read.csv("Processed_HJ_Andrews_Tutorial_Mol.csv", row.names = 1)

colnames(data) = paste("Sample_", colnames(data), sep="")

# Checking row names consistency
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect in your row names")
}

# Loading in transformations
trans.full =  read.csv("Transformation_List_Amino_Acid.csv")
trans.full$Name = as.character(trans.full$Name)

# ############# #
#### Errors ####
# ############ #

# Checking row names consistency between molecular info and data
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect: the mol. info and peak counts don't match")
}

if(length(which(mol$C13 == 1)) > 0){
  stop("Isotopic signatures weren't removed")
}

if(length(grep("QC_SRFAII", colnames(data))) > 0){
  stop("Suwannee River standards are still in the data")
}

if(length(grep("rep1|rep2", colnames(data))) > 0){
  stop("Technical replicates are still in the data")
}

if(max(data) > 1){
  print("Data was not presence/absence, so we are converting it")
  data[data > 1] = 1
}

# Creating output directories
if(!dir.exists("Transformation Peak Comparisons")){
  dir.create("Transformation Peak Comparisons")
}

if(!dir.exists("Transformations per Peak")){
  dir.create("Transformations per Peak")
}


###########################################
### Running through the transformations ###
###########################################

# pull out just the sample names
samples.to.process = colnames(data)

# error term
error.term = 0.000010

# matrix to hold total number of transformations for each sample
tot.trans = numeric()

# matrix to hold transformation profiles
profiles.of.trans = trans.full
head(profiles.of.trans)

for (current.sample in samples.to.process) {
  
  print(date())
  
  one.sample.matrix = cbind(as.numeric(as.character(row.names(data))), data[,which(colnames(data) == current.sample), drop = FALSE]) # "drop = FALSE" ensures that the row and column names remain associated with the data
  colnames(one.sample.matrix) = c("peak", colnames(one.sample.matrix[2]))
  # print(head(one.sample.matrix))
  
  Sample_Peak_Mat <- one.sample.matrix %>% gather("sample", "value", -1) %>% filter(value > 0) %>% select(sample, peak)
  Distance_Results <- Sample_Peak_Mat %>% left_join(Sample_Peak_Mat, by = "sample") %>% filter(peak.x > peak.y) %>% mutate(Dist = peak.x - peak.y) %>% select(sample, Dist,peak.x,peak.y)
  Distance_Results$Dist.plus = Distance_Results$Dist + error.term
  Distance_Results$Dist.minus = Distance_Results$Dist - error.term
  Distance_Results$Trans.name = -999
  head(Distance_Results)
  
  dist.unique = unique(Distance_Results[,'sample']) #unique samples
  
  date()
  
  #counter = 1
  
  for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
    
    mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
    if (length(mass.diff) > 1) { break() }
    Distance_Results$Trans.name[ which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)  ] = current.trans
    #print(c(counter,current.trans,mass.diff,length(mass.diff)))
    
    #counter = counter + 1
    
  }
  
  date()
  
  Distance_Results = Distance_Results[-which(Distance_Results$Trans.name == -999),]
  head(Distance_Results)
  
  # Creating directory if it doesn't exist, prior to writing the output file
  if(length(grep(Sample_Name,list.dirs("Transformation Peak Comparisons", recursive = F))) == 0){
    dir.create(paste("Transformation Peak Comparisons/", Sample_Name, sep=""))
    print("Directory created")
  }

  write.csv(Distance_Results,paste("Transformation Peak Comparisons/",Sample_Name,"/Peak.2.Peak_",dist.unique,".csv",sep=""),quote = F,row.names = F)

  # Alternative .csv writing
  # write.csv(Distance_Results,paste("Transformation Peak Comparisons/", "Peak.2.Peak_",dist.unique,".csv",sep=""),quote = F,row.names = F)
  
  # sum up the number of transformations and update the matrix
  tot.trans = rbind(tot.trans,c(dist.unique,nrow(Distance_Results)))
  
  # generate transformation profile for the sample
  trans.profile = as.data.frame(tapply(X = Distance_Results$Trans.name,INDEX = Distance_Results$Trans.name,FUN = 'length')); head(trans.profile)
  colnames(trans.profile) = dist.unique
  head(trans.profile)
  
  # update the profile matrix
  profiles.of.trans = merge(x = profiles.of.trans,y = trans.profile,by.x = "Name",by.y = 0,all.x = T)
  profiles.of.trans[is.na(profiles.of.trans[,dist.unique]),dist.unique] = 0
  head(profiles.of.trans)
  str(profiles.of.trans)
  
  # find the number of transformations each peak was associated with
  peak.stack = as.data.frame(c(Distance_Results$peak.x,Distance_Results$peak.y)); head(peak.stack)
  peak.profile = as.data.frame(tapply(X = peak.stack[,1],INDEX = peak.stack[,1],FUN = 'length' )); dim(peak.profile)
  colnames(peak.profile) = 'num.trans.involved.in'
  peak.profile$sample = dist.unique
  peak.profile$peak = row.names(peak.profile)
  head(peak.profile);
  
  # Creating directory if it doesn't exist, prior to writing the output file
  if(length(grep(Sample_Name,list.dirs("Transformations per Peak", recursive = F))) == 0){
    dir.create(paste("Transformations per Peak/", Sample_Name, sep=""))
    print("Directory created")
  }

  # Writing data to newly created directory
  write.csv(peak.profile,paste("Transformations per Peak/",Sample_Name,"/Num.Peak.Trans_",dist.unique,".csv",sep=""),quote = F,row.names = F)
  
  # Alternative .csv writing
  # write.csv(peak.profile,paste("Transformations per Peak/", "Num.Peak.Trans_",dist.unique,".csv",sep=""),quote = F,row.names = F)
  
  print(dist.unique)
  print(date())
  
}

# format the total transformations matrix and write it out
tot.trans = as.data.frame(tot.trans)
colnames(tot.trans) = c('sample','total.transformations')
tot.trans$sample = as.character(tot.trans$sample)
tot.trans$total.transformations = as.numeric(as.character(tot.trans$total.transformations))
str(tot.trans)
write.csv(tot.trans,paste(Sample_Name,"_Total_Transformations.csv", sep=""),quote = F,row.names = F)

# write out the trans profiles across samples
write.csv(profiles.of.trans,paste(Sample_Name, "_Trans_Profiles.csv", sep=""),quote = F,row.names = F)

# ######################## #
### Plotting the results ###
# ######################## #

# Preprocessing data
row.names(profiles.of.trans) = profiles.of.trans$Name # Fixing row names
trans.pro = profiles.of.trans[,-1:-2] # Removing extra columns
trans.pro = apply(trans.pro, 2, function(x) x/sum(x))
trans.pro = as.data.frame(t(trans.pro))

# Creating factors
factors = data.frame(Sample_Name = row.names(trans.pro), Sample_Type = "Pore Water", stringsAsFactors = F)
factors$Sample_Type[grep("SW48", factors$Sample_Name)] = "Surface Water"

# Melting data
trans.pro$Type = factors$Sample_Type
trans.melt = melt(trans.pro, id.vars = c("Type"))

pdf("AA Transformations by Water Type.pdf")
print(
  ggplot(trans.melt, aes(x = Type, y = value, group = Type))+
    geom_boxplot(aes(color = Type))+
    facet_wrap(.~variable, ncol = 3, scales = "free_y")+
    ylab("Proportion of AA Transformations (%)")+
    xlab(NULL)+
    scale_color_manual(values = c("dodgerblue4", "firebrick4"))+ # Remove this line if you want to run more samples
    ggtitle("AA Proportion by water type")+
    theme_bw()+
    theme(text = element_text(size = 14),
          axis.text.x = element_text(colour = "black", angle = 12.5, hjust = 1, vjust = 1),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
)
dev.off()


trans.melt$variable = factor(trans.melt$variable, 
                             levels = unique(trans.melt$variable)[order(colMeans(trans.pro[,-22]), decreasing = T)])

pdf("Comparison of AA Transformations.pdf")
print(
  ggplot(trans.melt, aes(x = variable, y = value, group = variable))+
    geom_boxplot()+
    ylab("Proportion of AA Transformations (%)")+
    xlab(NULL)+
    ggtitle("AA Proportions")+
    theme_bw()+
    theme(text = element_text(size = 14),
          axis.text.x = element_text(colour = "black", angle = 37.5, hjust = 1, vjust = 1),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
)
dev.off()
