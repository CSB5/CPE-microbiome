#####################
#### Supp Fig 3 ####
#####################
# %%
meta_diversity_jan22 <- read.delim("../data/Supplementary_file_1.txt", stringsAsFactors = TRUE)
# Remove "1 month post clearance" designation as that is not used in the plots
meta_diversity_jan22_no1 = meta_diversity_jan22[(meta_diversity_jan22$Colonization.status != "CPE_negative_1" &!is.na(meta_diversity_jan22$Colonization.status)), ]
meta_diversity_jan22_no1$Colonization.status <- factor(meta_diversity_jan22_no1$Colonization.status,
        levels = c("CPE_positive", "CPE_negative_0", "CPE_negative_2", "Family_member"))

# Remove rows with "Data not collected" designation for certain comparisons
temp.df = meta_diversity_jan22[meta_diversity_jan22$Antibiotics.since.last.visit != "Data not collected", ]
# Remove "1 month post clearance" designation as that is not used in the plots
temp.df.no1 = temp.df[(temp.df$Colonization.status != "CPE_negative_1" & !is.na(temp.df$Colonization.status)),  ]
temp.df.no1$Colonization.status <- factor(temp.df.no1$Colonization.status,
        levels = c("CPE_positive", "CPE_negative_0", "CPE_negative_2", "Family_member"))

# %% #### Effect of antibiotics ####
# Find effect
library('lme4')
mem_col_rand_ab = lmer(Genus.level.Shannon.diversity ~ Colonization.status + (1 | Antibiotics.since.last.visit), data = temp.df)
summary(mem_col_rand_ab)
ranef(mem_col_rand_ab) # No: 0.1867723, Yes: -0.1867723

# Remove effect
temp.df.no1$Shan_ab = temp.df.no1$Genus.level.Shannon.diversity
for (row in 1:nrow(temp.df.no1)){
  if (temp.df.no1[row,"Antibiotics.since.last.visit"] == "No"){
    temp.df.no1[row,"Shan_ab"] = temp.df.no1[row,"Shan_ab"] - 0.1867723
  }
  else if (temp.df.no1[row,"Antibiotics.since.last.visit"] == "Yes"){
    temp.df.no1[row,"Shan_ab"] = temp.df.no1[row,"Shan_ab"] + 0.1867723
  }
  else {
    print("error")
  }
}

# %% # Plot
boxplot(Shan_ab~Colonization.status, data=temp.df.no1, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))
axis(1, at = seq(0, 4, by = 1))

# %% #### Effect of hospitalization ####
# Find effect
mem_col_rand_hosp = lmer(Genus.level.Shannon.diversity ~ Colonization.status + (1 | Hospitalization.since.last.visit), data = temp.df)
summary(mem_col_rand_hosp)
ranef(mem_col_rand_hosp) # No: 0.17796921, Yes: 0.04535177, Still.in.hosp: -0.22332098

# Remove effect
temp.df.no1$Shan_hosp = temp.df.no1$Genus.level.Shannon.diversity
for (row in 1:nrow(temp.df.no1)){
  if (temp.df.no1[row,"Hospitalization.since.last.visit"] == "No"){
    temp.df.no1[row,"Shan_hosp"] = temp.df.no1[row,"Shan_hosp"] - 0.17796921
  }
  else if (temp.df.no1[row,"Hospitalization.since.last.visit"] == "Yes"){
    temp.df.no1[row,"Shan_hosp"] = temp.df.no1[row,"Shan_hosp"] - 0.04535177
  }
  else if (temp.df.no1[row,"Hospitalization.since.last.visit"] == "Still in hospital"){
    temp.df.no1[row,"Shan_hosp"] = temp.df.no1[row,"Shan_hosp"] + 0.22332098
  }
  else {
    print("error")
  }
}

# Plot
boxplot(Shan_hosp~Colonization.status, data=temp.df.no1, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))
axis(1, at = seq(0, 4, by = 1))

#### Effect of individuals ####

# Find effect
mem_col_rand_indiv = lmer(Genus.level.Shannon.diversity ~ Colonization.status + (1 | Indiv.Code), data = meta_diversity_jan22)
summary(mem_col_rand_indiv)
ranef(mem_col_rand_indiv)
c = ranef(mem_col_rand_indiv)$Indiv.Code

# Remove effect
meta_diversity_jan22_no1$Shan_indiv = meta_diversity_jan22_no1$Genus.level.Shannon.diversity
for (row in 1:nrow(meta_diversity_jan22_no1)){
  famTerm = meta_diversity_jan22_no1[row, "Indiv.Code"]
  intTerm = c$`(Intercept)`[row.names(c) == famTerm]
  meta_diversity_jan22_no1[row,"Shan_indiv"] = meta_diversity_jan22_no1[row,"Shan_indiv"] - intTerm
}

# Plot
boxplot(Shan_indiv~Colonization.status, data=meta_diversity_jan22_no1, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))
axis(1, at = seq(0, 4, by = 1))

#### Effect of gender ####

# Find effect
mem_col_rand_gender = lmer(Genus.level.Shannon.diversity ~ Colonization.status + (1 | Gender), data = meta_diversity_jan22)
summary(mem_col_rand_gender)
ranef(mem_col_rand_gender) #F emale: -0.08492628, Male: 0.08492628

# Remove effect
meta_diversity_jan22_no1$Shan_gender = meta_diversity_jan22_no1$Genus.level.Shannon.diversity
for (row in 1:nrow(meta_diversity_jan22_no1)){
  if (meta_diversity_jan22_no1[row,"Gender"] == "Female"){
    meta_diversity_jan22_no1[row,"Shan_gender"] = meta_diversity_jan22_no1[row,"Shan_gender"] + 0.08492628
  }
  else if (meta_diversity_jan22_no1[row,"Gender"] == "Male"){
    meta_diversity_jan22_no1[row,"Shan_gender"] = meta_diversity_jan22_no1[row,"Shan_gender"] - 0.08492628
  }
  else {
    print("error")
  }
}

# Plot
boxplot(Shan_gender~Colonization.status, data=meta_diversity_jan22_no1, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))
axis(1, at = seq(0, 4, by = 1))

#### Effect of ethnicity ####

# Find effect
mem_col_rand_race = lmer(Genus.level.Shannon.diversity ~ Colonization.status + (1 | Race), data = meta_diversity_jan22)
summary(mem_col_rand_race)
ranef(mem_col_rand_race) # Chinese: 0.08851847, Indian: -0.11229462, Malay: 0.10066944, Others: -0.07689329

# Remove effect
meta_diversity_jan22_no1$Shan_race = meta_diversity_jan22_no1$Genus.level.Shannon.diversity
for (row in 1:nrow(meta_diversity_jan22_no1)){
  if (meta_diversity_jan22_no1[row,"Race"] == "Chinese"){
    meta_diversity_jan22_no1[row,"Shan_race"] = meta_diversity_jan22_no1[row,"Shan_race"] - 0.08851847
  }
  else if (meta_diversity_jan22_no1[row,"Race"] == "Indian"){
    meta_diversity_jan22_no1[row,"Shan_race"] = meta_diversity_jan22_no1[row,"Shan_race"] + 0.11229462
  }
  else if (meta_diversity_jan22_no1[row,"Race"] == "Malay"){
    meta_diversity_jan22_no1[row,"Shan_race"] = meta_diversity_jan22_no1[row,"Shan_race"] - 0.10066944
  }
  else if (meta_diversity_jan22_no1[row,"Race"] == "Others"){
    meta_diversity_jan22_no1[row,"Shan_race"] = meta_diversity_jan22_no1[row,"Shan_race"] + 0.07689329
  }
  else {
    print("error")
  }
}

# Plot
boxplot(Shan_race~Colonization.status, data=meta_diversity_jan22_no1, ylab="", xlab="",
        col = c("indianred", "dodgerblue4", "cornflowerblue", "lightgreen"), xaxt="n",yaxt="n",
        ylim=c(-0.2,3.6))
axis(2, at = seq(0, 3, by = 1))
axis(1, at = seq(0, 4, by = 1))
