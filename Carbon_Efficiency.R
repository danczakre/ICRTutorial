# Compute the lambda for chemical compositions
# By Hyun-Seob Song et al.
# Edited by RED

require(reshape2); require(ggplot2); require(easycsv)

# Sample Name
sample_name = "HJ_Andrews_Tutorial"

# read in molecular file
setwd(easycsv::choose_dir())
data = read.csv(list.files(pattern = "Processed.*Data.csv"), row.names = 1) # Keeping data and mol-data seperate to ensure they are unaltered
mol = read.csv(list.files(pattern = "Processed.*Mol.csv"), row.names = 1)

df = mol # Set df equal to mol

# ################ #
### Calculations ###
# ################ #

# extract molecular formulae
CHEMICAL_ELEMENTS = c("C","H","N","O","S","P")
molecularFormula <- unique(df$MolForm[!is.na(df$MolForm)])
length(which(df$C > 1))

# extract numerical formulae
numericalFormula <- array(0, dim=c(length(molecularFormula), length(CHEMICAL_ELEMENTS)))
for (k in 1:length(molecularFormula)){
  formula <- molecularFormula[k]
  ge <- gregexpr("[A-Z]\\d*", formula, perl=TRUE)
  s_index <- ge[[1]]
  s_len <- attr(s_index, "match.length")
  for (i in 1:length(s_len)){
    token <- substr(formula, s_index[i], s_index[i] + s_len[i] - 1)
    element <- substr(token, 1, 1)
    if (grepl(element, "CHNOSP")) {
      idx = which(CHEMICAL_ELEMENTS %in% element)
      if (numericalFormula[k, idx] > 0) next  # for C13
      if (s_len[i] == 1) {
        numericalFormula[k, idx] = 1
      } else {
        numElement <- try(strtoi(substr(formula, s_index[i] + 1, s_index[i] + s_len[i] - 1)))
        if (class(numElement)=="integer"){
          numericalFormula[k, idx] = numElement
        } else {
          print(paste("[ERROR] an unknown chemical element found:", token, "in", formula))
        }
      }
    } else {
      print(paste("[ERROR] an unknown chemical element found:", element, "in", formula))
    }
  }
}

######################## compute lambda ########################
getThermoStoich <- function(chemForm) {
  a <- chemForm[1]
  b <- chemForm[2]
  c <- chemForm[3]
  d <- chemForm[4]
  e <- chemForm[5]
  f <- chemForm[6]
  z <- 0 #chemForm[7]
  
  # Step 1a) stoichD: stoichiometries for an electron donor
  ySource <- -1
  yH2o <- -(3*a+4*e-d)
  yHco3 <- a
  yNh4 <- c
  yHpo4 <- e
  yHs <- f
  yH <- 5*a+b-4*c-2*d+7*e-f
  yE <- -z+4*a+b-3*c-2*d+5*e-2*f
  stoichD <- c(ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE)
  stoichD[c(9,10)] <- 0 # add additional components: e-acceptor and biomass
  
  # Step 1b) stoichA: stoichiometries for an electron acceptor (i.e., oxygen)
  stoichA <- rep(0, 10)
  stoichA[9] <- -1  # oxygen
  stoichA[7] <- -4  #  h+
  stoichA[8] <- -4  #  e-
  stoichA[2] <- 2  #  h2o
  
  # Step 1c) stoichCat: stoichiometries for catabolic reaciton 
  yEd <- stoichD[8]
  yEa <- stoichA[8]
  stoichCat <- stoichD-(yEd/yEa)*stoichA
  
  # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton 
  #          (N source = NH4+)
  chemFormBiom <- c(1, 1.8, 0.2, 0.5, 0, 0, 0)  # C H_1.8 N_0.2 O_0.5
  aB <- chemFormBiom[1]
  bB <- chemFormBiom[2]
  cB <- chemFormBiom[3]
  dB <- chemFormBiom[4]
  eB <- chemFormBiom[5]
  fB <- chemFormBiom[6]
  zB <- chemFormBiom[7]
  
  ySource <- -1
  yH2o <- -(3*aB+4*eB-dB)
  yHco3 <- aB
  yNh4 <- cB
  yHpo4 <- eB
  yHs <- fB
  yH <- 5*aB+bB-4*cB-2*dB+7*eB-fB
  yE <- -zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB
  stoichAnStarB <- c(ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE)
  stoichAnStarB[c(9,10)] <- 0  # add additional components: e-acceptor and biomass
  stoichAnStarB <- -stoichAnStarB
  stoichAnStarB[10] <- stoichAnStarB[1]
  stoichAnStarB[1] <- 0
  
  # Step 2b) "overall" anabolic reaction
  eA4Anabolic=c( # electron acceptor for anabolic reaction
    'O2',    # Kleerebezem and Van Loosdrecht (2010)
    'HCO3-' # % McCarty (year?)
  )
  
  for (i in 1:length(eA4Anabolic)) {
    eA4Ana <- eA4Anabolic[i]
    if (eA4Ana == 'O2') {
      stoichAnStar_O2 <- stoichAnStarB+(1/a)*stoichD
      yEana <- stoichAnStar_O2[8]
      if (yEana > 0)
        stoichAn_O2 <- stoichAnStar_O2-yEana/yEa*stoichA
      else if (yEana < 0)
        stoichAn_O2 <- stoichAnStar_O2-yEana/yEd*stoichD
      else
        stoichAn_O2 <- stoichAnStar_O2
    } else if (eA4Ana == 'HCO3-') {
      yEd <- stoichD[8] # TODO
      yEa <- stoichAnStarB[8]
      stoichAn_HCO3 <- stoichD-(yEd/yEa)*stoichAnStarB
      stoichAn_HCO3 <- stoichAn_HCO3/stoichAn_HCO3[10]  # TODO: normalize?
    }
  }
  
  # Step 3: get lambda
  
  # - estimate delGcox0 using LaRowe and Van Cappellen (2011)
  ne <- -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
  nosc <- -ne/a+4  # nominal oxidataion state of carbon 
  delGcox0PerE <- 60.3-28.5*nosc  # kJ/C-mol
  delGcox0 <- delGcox0PerE*a*abs(stoichD[1])  # kJ/rxn
  
  # - estimate delGf0 for electron donor
  delGf0_D_zero <- 0
  delGf0_zero <- c(delGf0_D_zero, -237.2, -586.8, -79.3, -1096.1, 12.1, 0, 0, 16.4, -67)
  delGcox0_zero <- drop(delGf0_zero %*% stoichD)
  delGf0_D_est <- (delGcox0-delGcox0_zero)/stoichD[1]
  # - finally, delGf0
  delGf0 <- delGf0_zero
  delGf0[1] <- delGf0_D_est
  
  # - standard delG at pH=0
  delGcat0 <- drop(delGf0 %*% stoichCat)
  delGan0_O2 <- drop(delGf0 %*% stoichAn_O2)
  delGan0_HCO3 <- drop(delGf0 %*% stoichAn_HCO3)
  
  # - stadard delG at pH=7
  R <- 0.008314  # kJ/(K.mol)
  T <- 298  # K
  iProton <- 7  # [eD,h2o,hco3-,nh4+,hpo4^2-,hs-,h+,e-,eA,biom]
  delGcox <- delGcox0+R*T*stoichD[iProton]*log(1e-7)
  delGcat <- delGcat0+R*T*stoichCat[iProton]*log(1e-7)
  delGan_O2 <- delGan0_O2+R*T*stoichAn_O2[iProton]*log(1e-7)
  delGan_HCO3 <- delGan0_HCO3+R*T*stoichAn_HCO3[iProton]*log(1e-7)
  
  # The Thermodynamic Electron Equivalents Model (TEEM)
  # --------
  eta <- 0.43
  delGsyn <- 200  # kJ/(mol.X)
  if (delGan_O2 < 0)
    m_O2 <- 1
  else
    m_O2 <- -1
  
  if (delGan_HCO3 < 0)
    m_HCO3 <- 1
  else
    m_HCO3 <- -1
  
  lambda_O2 <- (delGan_O2*eta^m_O2+delGsyn)/(-delGcat*eta)
  lambda_HCO3 <- (delGan_HCO3*eta^m_HCO3+delGsyn)/(-delGcat*eta)
  
  if (lambda_O2 > 0)
    stoichMet_O2 <- lambda_O2*stoichCat+stoichAn_O2
  else
    stoichMet_O2 <- stoichAn_O2
  
  if (lambda_HCO3 > 0)
    stoichMet_HCO3 <- lambda_HCO3*stoichCat+stoichAn_HCO3
  else
    stoichMet_HCO3 <- stoichAn_HCO3
  
  delGdis_O2 <- drop(delGf0 %*% stoichMet_O2) + R*T*stoichMet_O2[iProton]*log(1e-7)
  delGdis_HCO3 <- drop(delGf0 %*% stoichMet_HCO3) + R*T*stoichMet_HCO3[iProton]*log(1e-7)
  # delGdis <- 200+18*(6-a)^1.8 + exp(((-0.2+nosc)^2)^0.16*(3.6+0.4*a))
  
  # list(delGcox0PerE = delGcox0PerE,
  #      delGcox0 = delGcox0,
  #      delGcox = delGcox,
  #      delGcat0 = delGcat0,
  #      delGcat = delGcat,
  #      delGan0_O2 = delGan0_O2,
  #      delGan0_HCO3 = delGan0_HCO3,
  #      delGan_O2 = delGan_O2,
  #      delGan_HCO3 = delGan_HCO3,
  #      delGdis_O2 = delGdis_O2,
  #      delGdis_HCO3 = delGdis_HCO3,
  #      lambda_O2 = lambda_O2,
  #      lambda_HCO3 = lambda_HCO3,
  #      stoich.D = stoichD,
  #      stoich.A = stoichA,
  #      stoich.Cat = stoichCat,
  #      stoich.An_O2 = stoichAn_O2,
  #      stoich.An_HCO3 = stoichAn_HCO3,
  #      stoich.Met_O2 = stoichMet_O2,
  #      stoich.Met_HCO3 = stoichMet_HCO3)
  c(delGcox0PerE,delGcox0,delGcox,delGcat0,delGcat,delGan0_O2,delGan0_HCO3,
    delGan_O2,delGan_HCO3,delGdis_O2,delGdis_HCO3,lambda_O2,lambda_HCO3,
    stoichD,stoichA,stoichCat,stoichAn_O2,stoichAn_HCO3,
    stoichMet_O2,stoichMet_HCO3)
}
######################## compute lambda ########################

getLambda <- function(formulaMatrix) {
  nrows = nrow(formulaMatrix)
  lambda_rst <- array(0, dim=c(nrows, 83))
  for(i in 1:nrows) {
    lambda_rst[i,] <- getThermoStoich(formulaMatrix[i,])
  }
  lambda_rst
}

# out <- getLambda(molecularFormula, numericalFormula)
out <- getLambda(numericalFormula)

# build data frame
df <- as.data.frame(out)
# build col names
names <- rep("", 83)
names[1:13] <- c("delGcox0PerE","delGcox0","delGcox","delGcat0","delGcat","delGan0_O2","delGan0_HCO3",
                "delGan_O2","delGan_HCO3","delGdis_O2","delGdis_HCO3","lambda_O2","lambda_HCO3")
stoich_colnames <- c("donor","h2o","hco3","nh4","hpo4","hs","h","e","acceptor","biom")
stoich_types <- c("stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3",
                  "stoichMet_O2","stoichMet_HCO3")
for (i in 1:length(stoich_types)) {
  names[((i-1)*10+14):(i*10+13)] <- array(sapply(stoich_types[i], paste, stoich_colnames, sep="_"))
}
colnames(df) <- names
df['formula'] <- molecularFormula

# write.csv(df, file = paste0(sample_name, "_Lambda_Ycsi.csv"), row.names = F)

# ############################ #
#### Plotting Y_csi results ####
# ########################### #

# Merging Y_csi into mol data
df = data.frame(formula = df$formula, lamO2 = df$lambda_O2, lamHCO3 = df$lambda_HCO3,
                 Y_met = df$stoichMet_O2_donor)

mol$Mass = row.names(mol)
mol = merge(x = mol, y = df, by.x = "MolForm", by.y = "formula", all = T, )
row.names(mol) = mol$Mass
mol = mol[order(as.numeric(mol$Mass)),]
mol$Y_met = mol$Y_met*-1

# Creating a Y_csi-by-sample object
Y_met = as.data.frame(matrix(data = 0, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data)))

for(i in 1:ncol(data)){
  temp = data[,i]
  Y_met[which(temp > 0),i] = mol$Y_met[which(temp > 0)]
}

# Cleaning data
Y_met[is.na(Y_met)] = 0
Y_met = Y_met[-which(rowSums(Y_met) == 0),]
Y_met = Y_met[,-grep("PP48.000012", colnames(Y_met))]

# Adding factors
Y_melt = as.data.frame(t(Y_met))
Y_melt$Type = "Pore Water"
Y_melt$Type[grep("SW48", row.names(Y_melt))] = "Surface Water"
Y_melt$Sample = row.names(Y_melt)

# Plotting totality of data
Y_melt = melt(Y_melt, id.vars = c("Type", "Sample"))
Y_melt = Y_melt[-which(Y_melt$value == 0),]

ggplot(Y_melt, aes(x = Sample, y = value))+
  geom_boxplot()+
  facet_grid(.~Type, scales = c("free_x"))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = 37.5, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

# Plotting average data
Y_met_av = data.frame(Sample = colnames(Y_met), Type = "Pore Water", 
                      value = NA, row.names = colnames(Y_met),
                      stringsAsFactors = F)
Y_met_av$Type[grep("SW48", Y_met_av$Sample)] = "Surface Water"
Y_met_av$Sample = gsub(pattern = "HJ.*PP48_000|HJ.*SW48_000", replacement = "", Y_met_av$Sample)

for(i in 1:ncol(Y_met)){
  temp = Y_met[,i]
  temp = temp[-which(temp == 0)]
  Y_met_av$value[i] = mean(temp)
}

ggplot(Y_met_av, aes(x = Sample, y = value, group = Type))+
  geom_point(aes(color = Type))+
  geom_line(aes(color = Type))+
  xlab("Time Point")+
  ylab(expression("Y"[csi]))+
  scale_color_manual(values = c("dodgerblue4", "firebrick4"))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())
