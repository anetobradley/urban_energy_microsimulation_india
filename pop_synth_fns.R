# mod_1_pop_synth
# IHEUS Model Module 1
# This module contains functions to prepare data for population synthesis, uses IPF 
# to synthesise a population for a given city using Indian Census Data and performs 
# internal validation of the resulting dataset.
# A P Neto-Bradley 2020

set.seed(2019)

library(dplyr)
library(ipfp)

#### PREPARE CONSTRAINTS TABLE ####
# Import census constraint data for Bangalore
const_table_from_census <- function(city){
  census_abstract <- read.csv(paste0("Census_Abstract_",city,".csv"))
  census_ward_socio_ec <- read.csv(paste0("Census_HH_Q_",city,".csv"), stringsAsFactors = F)
  
  # Filter out relevant constraints
  constraints <- left_join(census_abstract[,c(5,10,11,17,20,26)],census_ward_socio_ec[,c(8,23:64,109:118,120,123,126,127,129,136)], by ="Ward")
  
  # Tidy header rows and labels - Table values as percentage
  constraints$P_SC <- round((constraints$P_SC/constraints$TOT_P)*100, digits = 1)
  constraints$P_ST <- round((constraints$P_ST/constraints$TOT_P)*100, digits = 1)
  constraints$P_ILL <- round((constraints$P_ILL/constraints$TOT_P)*100, digits = 1)
  
  # Make sure constraint names match 
  colnames(constraints) <- c("Ward","No_HH","Total_Pop",
                             "SC","ST", # Relates to IHDS ID13 (SC - 4, ST - 5)
                             "Illiterate", # Relates to IHDS HHEDUC
                             "Rf_thatch","Rf_plastic","Rf_Handtile","Rf_Machtile","Rf_Brick","Rf_Stone","Rf_Metal","Rf_Concrete","Rf_Other", # Relates to IHDS HQ5
                             "Wl_thatch","Wl_plastic","Wl_Mud","Wl_Wood","Wl_Unpacked_stone","Wl_Stone","Wl_Metal","Wl_Brick","Wl_Concrete","Wl_Other", # Relates to IHDS HQ4
                             "Fl_Mud","Fl_Wood","Fl_Brick","Fl_Stone","Fl_Concrete","Fl_Tiles","Fl_Other", # Relates to IHDS HQ6
                             "No_excl_rm","1_rm","2_rm","3_rm","4_rm","5_rm","6_rm", # Relates to IHDS SA1
                             "HH_of_1","HH_of_2","HH_of_3","HH_of_4","HH_of_5","HH_of_6","HH_greater_than_6", # Relates to IHDS NPERSONS
                             "Own_house","Rent_house", # Relates to IHDS CG1
                             "Cook_Firewood","Cook_Crop","Cook_Dung","Cook_Coal","Cook_Kero","Cook_LPG","Cook_Elec","Cook_Biogas","Cook_Other","No_cooking", # Relates to IHDS FU6, FU7, FU8, FU9, FU10, FU11, FU12
                             "Cooking_in","Cooking_out","No_kitchen", # Relates to IHDS SA2
                             "Banking_access", # Relates to IHDS DB1A & DB9E
                             "Own_TV", # Relates to IHDS CGTV
                             "Own_Moped") # Relates to IHDS CG8
  
  # Note that the IHDS has only one joint category for building material - tiles and stone types
  constraints[,9] <- as.numeric(as.character(constraints[,9])) + as.numeric(as.character(constraints[,10])) 
  constraints[,20] <- as.numeric(as.character(constraints[,20])) + as.numeric(as.character(constraints[,21])) 
  
  # Note that the IHDS has a limited range of cooking fuel categories so we will remove 
  # non-corresponded variables
  constraints <- constraints[,-c(10,21,55:57)]
  
  return(constraints)
  
}

#### PREPARE SEED DATA ####
seed_from_ihds_dist <- function(distwise,dist_code){
  #### Subset IHDS ##########
  # Import and isolate the Bangalore Urban District data from IHDS
  load("36151-0002-Data.Rda")
  
  # Filter for Bangalore Urban District - Code 2920 
  if(distwise == T){
    ihds_filter <- filter(da36151.0002, DISTRICT == dist_code)
  }
  else{
    ihds_filter <- da36151.0002 %>%
      filter(STATEID == dist_code) %>%
      filter(URBAN2011 == "(1) urban 1")
  }
  
  # Filter for relevant variables
  seed <- ihds_filter[,c("IDHH",
                         "ID13","HHEDUC",
                         "HQ5","HQ4","HQ6",
                         "SA1","NPERSONS",
                         "CG1",
                         "FU6","FU7","FU8","FU9","FU10","FU11","FU12",
                         "SA2",
                         "DB1A","DB9E",
                         "CGTV",
                         "CG8",
                         "CG18",
                         "WA2A",
                         "HQ1",
                         "SA4")]
  
  # Before we can use this factor dummy variable tranformation we need to replace NA's with values.
  # In this case we will default to assigning NA to the other category
  #### Handle Missing Data ####
  seed[,c(2)][which(is.na(seed[,c(2)]))] <- levels(seed[,c(2)])[6]
  seed[,c(4)][which(is.na(seed[,c(4)]))] <- levels(seed[,c(4)])[10]
  seed[,c(5)][which(is.na(seed[,c(5)]))] <- levels(seed[,c(5)])[9]
  seed[,c(6)][which(is.na(seed[,c(6)]))] <- levels(seed[,c(6)])[7]
  seed[,c(9)][which(is.na(seed[,c(9)]))] <- levels(seed[,c(9)])[4]
  
  seed[,"SA1"][which(is.na(seed[,"SA1"]))] <- 1
  
  if(mean(is.na(seed[,c(17)])>0)){
    seed[,c(17)] <- as.character(seed[,c(17)])
    seed[,c(17)][which(is.na(seed[,c(17)]))] <- "(4) No cooking"
    seed[,c(17)] <- as.factor(seed[,c(17)])
  }
  else{
    seed[,c(17)] <- factor(seed[,c(17)], levels = c(levels(seed[,c(17)]), "(4) No cooking"))
  }
 
  #### Create Dummy Variables from Factors ####
  
  # Format variables as dummies to match the constraints
  seed_caste <- model.matrix(~seed[,c(2)]-1)
  seed_materials <- model.matrix(~.+0, data=seed[,c(4:6)], contrasts.arg = lapply(seed[,c(4:6)], contrasts, contrasts=F))
  seed_owners <- model.matrix(~seed[,c(9)]-1)
  seed_cook_loc <- model.matrix(~seed[,c(17)]-1)
  
  # We need to bin room no. and household size
  seed_rooms <- model.matrix(~cut(seed[,"SA1"], breaks = c(0:6,Inf), labels=c(as.character(0:6)))-1)
  seed_persons <- model.matrix(~cut(seed[,"NPERSONS"], breaks = c(0:6,Inf), labels=c(as.character(1:6),"6+"))-1)
  
  # We also need to sort the fuel categories
  #### Engineer Cooking Dummy Vars from Energy Qs ####
  cook_firewood <- (1*(seed$FU7 == "(2) Mainly cooking 2")) + (1*(seed$FU7 == "(5) Combination 5"))
  cook_crop <- (1*(seed$FU8 == "(2) Mainly cooking 2")) + (1*(seed$FU8 == "(5) Combination 5"))
  cook_dung <- (1*(seed$FU9 == "(2) Mainly cooking 2")) + (1*(seed$FU9 == "(5) Combination 5"))
  cook_coal<- (1*(seed$FU12 == "(2) Mainly cooking 2")) + (1*(seed$FU12 == "(5) Combination 5"))
  cook_kero <- (1*(seed$FU10 == "(2) Mainly cooking 2")) + (1*(seed$FU10 == "(5) Combination 5"))
  cook_LPG <- (1*(seed$FU11 == "(2) Mainly cooking 2")) + (1*(seed$FU11 == "(5) Combination 5"))
  
  bio_stove <- (1*(seed$FU6 != "(4) Other/Not biomass (Kerosene, LPG etc.) 4"))
  
  cook_firewood <- cook_firewood*bio_stove
  cook_firewood[is.na(cook_firewood)] <- 0
  cook_LPG <- cook_LPG*(1-cook_kero)
  
  no_cooking <- (1*(is.na(bio_stove)==TRUE))
  
  seed_cooking <- data.frame(cook_firewood,cook_crop,cook_dung,cook_coal,cook_kero,cook_LPG,no_cooking)
  # The following are not included in the IHDS breakdown
  #cook_elec
  #cook_biogas
  #cook_other
  
  #### Extract Remaining Vars ####
  # Illiteracy (no schooling) needs to be deduced
  seed_illit <- (1*(seed$HHEDUC == "(00) none 0"))
  
  # Banking is asked in two parts in IHDS need to combine
  seed_banking <- (1*(seed$DB1A=="(1) Yes 1"))+(1*(seed$DB9E=="(1) Yes 1"))
  seed_banking[which(seed_banking==2)] <- 1
  seed_banking[which(is.na(seed_banking))] <- 0
  
  # Appliances
  seed_TV <- (1*(seed$CGTV=="(1) Yes 1"))
  seed_TV[which(is.na(seed_TV))] <- 0
  seed_Moped <- (1*(seed$CG8=="(1) Yes 1"))
  seed_Moped[which(is.na(seed_Moped))] <- 0
  
  
  #### Compile all seed data ####
  seed_proc <- data.frame(seed[,"IDHH"],
                          seed_caste,
                          seed_illit,
                          seed_materials,
                          seed_rooms,seed_persons,
                          seed_owners,
                          seed_cooking,
                          "triv_blank"=0,# Triv Patch
                          "triv_two" =0,# Triv Patch
                          seed_cook_loc,
                          seed_banking,
                          seed_TV,
                          seed_Moped)
  
  # We need to remove or combine categories/sub-categories not in census constraints
  # First will need to combine sub categories to match the category in the census
  seed_proc[,11] <- seed_proc[,11] + seed_proc[,16]
  seed_proc[,14] <- seed_proc[,14] + seed_proc[,17]
  seed_proc[,50] <- seed_proc[,50] + seed_proc[,51]
  seed_proc[,61] <- seed_proc[,61] + seed_proc[,62]
  
  seed_out <- seed_proc[complete.cases(seed_proc),-c(2:4,7,16,17,51:52,62)]

  # Need to rename the colnames for matching purposes.
  colnames(seed_out) <- c("HHID",
                           "SC","ST", # Relates to IHDS ID13 (SC - 4, ST - 5)
                           "Illiterate", # Relates to IHDS HHEDUC
                           "Rf_thatch","Rf_Handtile","Rf_Stone","Rf_plastic","Rf_Metal","Rf_Concrete","Rf_Brick","Rf_Other", # Relates to IHDS HQ5
                           "Wl_thatch","Wl_Mud","Wl_plastic","Wl_Wood","Wl_Brick","Wl_Metal","Wl_Unpacked_stone","Wl_Concrete","Wl_Other", # Relates to IHDS HQ4
                           "Fl_Mud","Fl_Wood","Fl_Brick","Fl_Stone","Fl_Concrete","Fl_Tiles","Fl_Other", # Relates to IHDS HQ6
                           "No_excl_rm","1_rm","2_rm","3_rm","4_rm","5_rm","6_rm", # Relates to IHDS SA1
                           "HH_of_1","HH_of_2","HH_of_3","HH_of_4","HH_of_5","HH_of_6","HH_greater_than_6", # Relates to IHDS NPERSONS
                           "Own_house","Rent_house", # Relates to IHDS CG1
                           "Cook_Firewood","Cook_Crop","Cook_Dung","Cook_Coal","Cook_Kero","Cook_LPG","No_cooking", # Relates to IHDS FU6, FU7, FU8, FU9, FU10, FU11, FU12
                           "Cooking_out","Cooking_in","No_kitchen", # Relates to IHDS SA2
                           "Banking_access", # Relates to IHDS DB1A & DB9E
                           "Own_TV", # Relates to IHDS CGTV
                           "Own_Moped") # Relates to IHDS CG8
  
  # Finally reorder columns to match constraints
  seed_out <- seed_out[, c(1:4,5,8,6,11,7,9,10,12,13,15,14,16,19,18,17,20,21,22:51,53,52,54:57)]
  return(seed_out) 
}

#### USE ipfp TO RUN IPF ####
ipf_popgen <- function(constraints,seed){
  cons <- (constraints[-c(200:209),"No_HH"]/10000)*apply(constraints[-c(200:209),], 2, as.numeric)+1
  ind_cat <- seed[,colnames(constraints[,c(4:59)])]
  ind_agg <- colSums(ind_cat)
  ind_cat <- ind_cat[,names(which(ind_agg!=0))]
  n_ind <- nrow(seed)
  x0=rep(1, n_ind)
  
  # Apply ipfp over all wards
  weights <- apply(cons[,names(which(ind_agg!=0))], MARGIN = 1, FUN = function(x) ipfp(x, t(ind_cat),x0, maxit=20))
  
  #### INTEGERISATION ####
  # From Lovelace and Ballas (2013)
  int_trs <- function(x){
    # For generalisation purpose, x becomes a vector
    xv <- as.vector(x) # allows trs to work on matrices
    xint <- floor(xv) # integer part of the weight
    r <- xv - xint # decimal part of the weight
    def <- round(sum(r)) # the deficit population
    # the weights be 'topped up' (+ 1 applied)
    topup <- sample(length(x), size = def, prob = r)
    xint[topup] <- xint[topup] + 1
    dim(xint) <- dim(x)
    dimnames(xint) <- dimnames(x)
    xint
  }
  
  #### EXPANSION ####
  int_expand_vector <- function(x){
    index <- 1:length(x)
    rep(index, round(x))
  }
  
  #### GENERATE INDIVIDUALS ####
  test_synth <- NULL
  weights_int <- unlist(apply(weights, MARGIN = 2, 
                              FUN = function(x) int_expand_vector(int_trs(x))))
  test_synth <- data.frame(seed[weights_int,1:57], 
                           zone = rep(1:nrow(cons), round(colSums(weights))))
  
  #### INTERNAL VALIDATION ####
  # Test against the census we started from on a city scale.
  print((sum(sqrt(cons[,2]))))
  # Plot synthetic totals against actual census totals
  plot(colSums(test_synth[2:57])/nrow(test_synth), colSums(cons[,4:59])/(sum(sqrt(cons[,2])/10)), col="blue")
  lines(c(0,1),c(0,100))
  text(colSums(test_synth[2:57])/nrow(test_synth), colSums(cons[,4:59])/(sum(sqrt(cons[,2])/10)), names(colSums(test_synth[2:57])), pos=4, col="red")
 
  print(cor(colMeans((cons[,4:59])), colSums(test_synth[2:57])))
  
  return(test_synth)
}
