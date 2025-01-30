####Snakeskin Biochemistry Manuscript Code 2024####

#' `Project:`         * Snake skin lipids and fungal-bacterial interactions influence the growth of Ophidiomyces ophidiicola  *
#' `Investigator:`    * Kaitlyn M. Murphy *
#' `Date:`            * August 19th, 2024 *

#Set working directory
setwd("H:/Shared drives/Microbiome Ecology Lab Shared Drive/Kaitlyn Murphy PC backup/SFD_Analyses/SFD_Manuscript")
getwd()

#Install these packages
install.packages('nlme')
install.packages('ggplot2')
install.packages('tidyverse')
install.packages('data.table')
install.packages('readxl')
install.packages('ggpubr')
install.packages('growthcurver')
install.packages('rstatix')
install.packages('RColorBrewer')
install.packages('lsmeans')
install.packages('lme4')
install.packages('lmerTest')
install.packages('ggtext')
install.packages('gt')
install.packages('multcomp')
install.packages('emmeans')
install.packages('chromote')
install.packages('gridExtra')
install.packages("viridis")
install.packages("MuMIn")
install.packages("ggsignif")
install.packages("superb")
install.packages("dplyr")
install.packages("sjPlot")
install.packages("webshot")
install.packages("sjPlot")
install.packages("scales")
install.packages("stringr")
install_phantomjs()

#Call on these libraries
library(nlme)
library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(growthcurver)
library(rstatix)
library(RColorBrewer)
library(lsmeans)
library(lme4)
library(lmerTest)
library(ggtext)
library(gt)
library(multcomp)
library(emmeans)
library(chromote)
library(gridExtra)
library(viridis)  
library(MuMIn)
library(ggsignif)
library(superb)
library(sjPlot)
library(webshot)
library(sjPlot)
library(scales)
library(stringr)

#Call on datasets
#NOTE: There are multiple used in this script

#Call on first dataset ("Exp1_Snakeskin")
datum_snakeskin=read.csv(file.choose())
head(datum_snakeskin)

#Call on second dataset ("Exp2_MinMedia")
datum_mm=read.csv(file.choose())
head(datum_mm)

#Call on third dataset ("Exp3_DBFI")
datum_bfi=read.csv(file.choose())
head(datum_bfi)

#Call on fourth dataset ("Exp4_SpentMedia_BacGrowth")
datum_sm=read.csv(file.choose())
head(datum_sm)

#Call on  fifth dataset ("Exp4_SpentMedia_FunGrowth")
datum_smj=read.csv(file.choose())
head(datum_smj)



#Link to color palette used in this manuscript
#https://hauselin.github.io/colorpalettejs/

####Exp 1: Snakeskin lipid extraction and O. ophiodiicola growth####

#Remove contaminated plates and lizard from dataset
datum_snakeskin <- datum_snakeskin[-c(1, 5, 14, 27, 34, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48), ]

#Change columns from continuous to categories
datum_snakeskin$treatment <- as.factor(datum_snakeskin$ID)
datum_snakeskin$TREAT <- as.factor(datum_snakeskin$TREAT)
datum_snakeskin$REPLICATE <- as.numeric(datum_snakeskin$REPLICATE)
datum_snakeskin$DATE <- as.factor(datum_snakeskin$DATE)
datum_snakeskin$TRIAL <- as.factor(datum_snakeskin$TRIAL)

######Analyses#####

#Preliminary plots
plot(AREAMM2~TREAT, data=datum_snakeskin)
plot(AREAMM2~AREAPXL, data=datum_snakeskin)

#Test of treatment
resultsskin = lmer(AREAMM2~relevel(TREAT, ref="Control") + (1|TRIAL:REPLICATE), data=datum_snakeskin)
summary(resultsskin)

######Plot#####

datum_snakeskin$TREAT_2 = factor(datum_snakeskin$TREAT,c("Control", "Solvent Control", "P. platyrhinos", "T. sirtalis", "D. punctatus", "Ch. bottae", "Cr. oreganus", "P. catenifer"))

tiff("Figure1A.tiff", width = 7, height = 4, units = 'in', res = 300)

ggplot(datum_snakeskin, aes(x = TREAT_2 , y = AREAMM2)) + 
  geom_boxplot(color="black", fill="#848076", alpha=0.2) +
  geom_point() +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, colour = "black"),
        axis.line = element_line(linewidth=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(face="italic", colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Species", y=(expression("Area in mm" ^ 2))) +
  scale_x_discrete(guide = guide_axis(angle = 45)) 

dev.off()

####Exp 2: O. ophiodiicola growth on minimal media####

#Clean up dataset
min.media.2 <- datum_mm %>% 
  unite(col='plate_id', c('MEDIUM', 'CONC', 'REPLICATE', 'TRIAL'), sep = '_', remove = FALSE) %>%
  mutate_at('AVERAGE', as.numeric) 

min.media.3 <- min.media.2 %>% 
  mutate(MEDIUM = factor(MEDIUM, levels = c('Keratin', 'Cholesterol', 'Oleic Acid', 'Squalene'))) %>% 
  mutate(CONC = factor(CONC, 
                       levels = c('1', '0.1', '0.01', '0.001', '0.0001')))

media_palette <- c("#002051", "#7f7c75", "#bbaf71", "#3c4d6e")

######Analyses#####

####  only examine growth at day 21

day.21.data <- datum_mm %>% 
  select(-c(DIAMETER..mm...RED, DIAMETER..mm...GREEN)) %>% 
  mutate_at('AVERAGE', as.numeric) %>% 
  mutate(AVERAGE = AVERAGE - 7) %>% 
  filter(DAY == 21) %>% 
  mutate(MEDIUM = factor(MEDIUM, levels = c('Control - Plate', 'Control - Glass', 'Keratin', 'Cholesterol', 'Oleic Acid', 'Squalene'))) %>% 
  mutate(CONC = factor(CONC, 
                       levels = c('Control', '1', '0.1', '0.01', '0.001', '0.0001'))) %>% 
  filter(!is.na(CONC)) %>% 
  filter(!is.na(MEDIUM)) %>% 
  mutate(MEDIUM.2 = factor(case_when(MEDIUM == 'Control - Plate' ~ "Control",
                                     MEDIUM == 'Control - Glass' ~ "Control",
                                     TRUE ~ MEDIUM)))


#Check reference levels
levels(day.21.data$CONC)
day.21.data$CONC <- relevel(day.21.data$CONC, ref = "Control")
levels(day.21.data$CONC)
levels(day.21.data$MEDIUM)
day.21.data$MEDIUM <- relevel(day.21.data$MEDIUM, ref = "Control - Plate")
levels(day.21.data$MEDIUM)
day.21.data$MEDIUM.2 <- relevel(day.21.data$MEDIUM.2, ref = "Control")


####   look at controls across trials
day.21.data.control <- day.21.data %>% 
  filter(CONC == 'Control')

results_mm_control=lmer(AVERAGE ~ TRIAL + MEDIUM + (1|REPLICATE:MEDIUM), data = day.21.data.control, na.action=na.omit)
summary(results_mm_control) #near significant difference in trial (P = 0.059) with a difference of 2.3 mm

summary(glht(results_mm_control, linfct = mcp(MEDIUM = "Tukey")), test = adjusted("holm"))

########  run the model without controls since we are comparing growth among the treatments, not to the controls here

day.21.data.no.control <- day.21.data %>% 
  filter(!CONC == 'Control')

day.21.data.no.control$MEDIUM <- relevel(day.21.data.no.control$MEDIUM, ref = "Keratin") #define keratin as reference level
day.21.data.no.control$CONC <- relevel(day.21.data.no.control$CONC, ref = "1") #define 1% as reference level


results_mm.2 <- lm(AVERAGE ~ MEDIUM + CONC*MEDIUM, data = day.21.data.no.control, na.action=na.omit)
summary(results_mm.2) #lots of significant results including interactions

###  post-hoc test
results_mm.2_post.hoc <- emmeans(results_mm.2, ~ MEDIUM + CONC*MEDIUM)
contrast(results_mm.2_post.hoc, "pairwise", simple = "each", combine = TRUE) #set method to pairwise to include all potential comparisons

#Results table
tab_model(results_mm.2, file = "MinMedia.html") #shows multiple significant main effects but a whole lot of significant interactions


##################  Subset by each concentration individually, include all controls since they do not differ
#subset the data by medium
datum_mm_keratin <- day.21.data[ which( day.21.data$MEDIUM == "Keratin" | day.21.data$CONC == "Control") , ]
datum_mm_squalene <- day.21.data[ which( day.21.data$MEDIUM == "Squalene" | day.21.data$CONC == "Control"), ]
datum_mm_chol <- day.21.data[ which( day.21.data$MEDIUM == "Cholesterol" | day.21.data$CONC == "Control") , ]
datum_mm_oa <- day.21.data[ which( day.21.data$MEDIUM == "Oleic Acid" | day.21.data$CONC == "Control") , ]

#Analyze each media separately
results_mm_keratin =lm(AVERAGE ~ CONC, data = datum_mm_keratin, na.action=na.omit)
summary(results_mm_keratin)
results_mm_squalene =lm(AVERAGE ~ CONC, data = datum_mm_squalene, na.action=na.omit)
summary(results_mm_squalene)
results_mm_oa =lm(AVERAGE ~ CONC, data = datum_mm_oa, na.action=na.omit)
summary(results_mm_oa)
results_mm_chol =lm(AVERAGE ~ CONC, data = datum_mm_chol, na.action=na.omit)
summary(results_mm_chol)


####    examine the mixed media plates ####  

day.21.mixed.data <- min.media %>% 
  select(-c(DIAMETER..mm...RED, DIAMETER..mm...GREEN)) %>% 
  mutate_at('AVERAGE', as.numeric) %>% 
  mutate(AVERAGE = AVERAGE - 7) %>% 
  filter(DAY == 21) %>% 
  mutate(MEDIUM = factor(MEDIUM, levels = c('Keratin + Chol', 'Keratin + Squal', 'Keratin + Oleic', 
                                            'Keratin', 'Cholesterol', 'Oleic Acid', 'Squalene'))) %>% 
  mutate(CONC = factor(CONC, 
                       levels = c('1'))) %>% 
  filter(!is.na(CONC)) %>% 
  filter(!is.na(MEDIUM))


#Check reference levels
day.21.mixed.data$MEDIUM <- relevel(day.21.mixed.data$MEDIUM, ref = "Keratin")

###  run model
mix.media.lm =lm(AVERAGE ~ MEDIUM, data = day.21.mixed.data, na.action=na.omit)
summary(mix.media.lm)

###  post-hoc test
mixed.media_post.hoc <- emmeans(mix.media.lm, ~ MEDIUM)
contrast(mixed.media_post.hoc, "pairwise") #set method to pairwise to include all potential comparisons


######Plot#####

#Figure 2A - single carbon source dataset

min.media.3 %>% 
  mutate(AVERAGE = AVERAGE - 7) %>% 
  filter(!CONC == 'Control') %>% 
  filter(!CONC == '0.00001') %>% 
  filter(DAY == 21) %>% 
  filter(MEDIUM == 'Cholesterol' | MEDIUM == 'Keratin' | MEDIUM == 'Oleic Acid' | MEDIUM == 'Squalene') %>% 
  ggplot(aes(x = CONC, y = AVERAGE, fill = MEDIUM))+
  geom_boxplot(alpha = 0.7)+
  scale_fill_manual(values = media_palette) +
  scale_y_continuous(limits = c(0,65))+
  theme_classic()+
  labs(x = 'Concentration (%)', y = "Fungal Growth (mm)")+
  facet_grid(. ~ MEDIUM)+
  geom_hline(yintercept = 41.250000, linetype = 2) +
  annotate("rect", xmin = 0.4, xmax = 5.6, ymin = 40.34, ymax = 42.16,
           alpha = .1, fill = "black")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size = 15), strip.text.x = element_text(size = 14), 
        legend.position="none")

###  figure 2c

media_palette.2 <- c("#002051", "#7f7c75", "#7f7c75", "#bbaf71", "#bbaf71", "#3c4d6e" , "#3c4d6e")

datum_mm %>% 
  select(-c(DIAMETER..mm...RED, DIAMETER..mm...GREEN)) %>% 
  mutate_at('AVERAGE', as.numeric) %>% 
  mutate(AVERAGE = AVERAGE - 7) %>% 
  filter(DAY == 21) %>% 
  mutate(CONC = factor(CONC, levels = c('1'))) %>%
  mutate(MEDIUM = factor(case_when(MEDIUM == 'Keratin + Chol' ~ "Keratin +\nCholestrol",
                                   MEDIUM == 'Keratin + Oleic' ~ "Keratin +\nOleic Acic",
                                   MEDIUM == 'Keratin + Squal' ~ "Keratin +\nSqualene",
                                   TRUE ~ MEDIUM))) %>% 
  mutate(MEDIUM = factor(MEDIUM, levels = c('Keratin', 'Keratin +\nCholestrol', 'Cholesterol', 'Keratin +\nOleic Acic',
                                            'Oleic Acid', 'Keratin +\nSqualene', 'Squalene'))) %>% 
  filter(!is.na(CONC)) %>% 
  filter(!is.na(MEDIUM)) %>% 
  ggplot(aes(x = MEDIUM, y = AVERAGE, fill = MEDIUM))+
  geom_boxplot(alpha = 0.7)+
  scale_fill_manual(values = media_palette.2) +
  scale_y_continuous(limits = c(0,65))+
  theme_classic()+
  labs(x = 'Medium', y = "Fungal Growth (mm)")+
  geom_hline(yintercept = 41.250000, linetype = 2) +
  annotate("rect", xmin = 0.4, xmax = 7.6, ymin = 40.34, ymax = 42.16,
           alpha = .1, fill = "black")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size = 15), strip.text.x = element_text(size = 14), 
        legend.position="none")


####Exp 3: Direct fungal-bacterial interactions####

#####Bacterial growth analyses#####

######Analyses#####

#Set plotting theme
theme_set(theme_pubr() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text.y   = element_text(size=20),
                  axis.text.x   = element_text(size=20),
                  axis.title.y  = element_text(size=20, margin = margin(r = 15)),
                  axis.title.x  = element_text(size=20, margin = margin(t = 15)),
                  legend.text = element_text(size=15),
                  legend.title = element_text(size = 20),
                  plot.title = element_text(face= "italic",size=25, hjust = 0.5),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")))

#Change isolate numbers to names in dataset
datum_bfi <- datum_bfi %>% 
  mutate(taxa = if_else(taxa == 'BR1.10', 'Chryseobacterium sp.', taxa)) %>% #change the species id to post-sequencing genus name
  mutate(taxa = if_else(taxa == 'BR1.11', 'Lelliottia nimipressuralis', taxa)) %>%
  mutate(taxa = if_else(taxa == 'BR1.7', 'Mammaliicoccus sciuri', taxa)) %>%
  mutate(taxa = if_else(taxa == 'BR1.9', 'Acinetobacter guillouiae', taxa)) %>%
  mutate(taxa = if_else(taxa == 'EKS20.5', 'Stenotrophomonas sp.', taxa)) %>%
  mutate(taxa = if_else(taxa == 'TR087-7.2', 'Flavobacterium odoratimimum', taxa))

#Create a column that defines this interaction
datum_bfi$index2 <- as.numeric(interaction(datum_bfi$oo, datum_bfi$index))

#cursory examination of bacterial growth curves
datum_bfi %>%
  filter(taxa != "blank") %>%
  ggplot(aes(x = day, y = od620, color = as.factor(oo), lty = media)) + 
  geom_smooth(method = "loess") +
  geom_point() +
  facet_wrap(~taxa) + 
  scale_x_continuous(limits = c(1,7))

#constrain analysis to initial four days 
df.bac <- datum_bfi %>%
  filter(day <= 4)

# Averaging across blank wells should capture any variation in blank wells fairly well 
df.bac %>%
  filter(taxa == "blank") %>%
  ggplot(aes(day, od620)) + 
  geom_point() +
  geom_smooth(aes(color = media), method = "loess", show.legend = F) + 
  facet_wrap(~media, scales = "free") 

# average od620 measurements for all blank measurements at all time-points
df.bac.blankmean <- df.bac %>%
  filter(taxa == "blank") %>%
  group_by(index2, media, day, oo) %>% 
  summarise_at(vars(od620), list(od620.blank = mean)) %>%
  right_join(df.bac, by = c("index2", "media", "day", 'oo')) %>%
  ungroup() %>%
  mutate(od620 = ifelse(taxa == "blank", od620.blank, od620)) %>%
  dplyr::select(-od620.blank) 

#remove repeat blank columns 
df.bac.blank_singles <- df.bac.blankmean %>%
  filter(taxa == "blank") %>%
  group_by(index2, media, day, oo) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  rbind(filter(df.bac.blankmean, taxa != "blank"))

#clean up working environment
df.bac <- df.bac.blank_singles
functions <- as.vector(lsf.str())
gdata::keep(list = c('datum_bfi', 'df.bac', functions), sure = T)

#write loop to summarize curve for each well in each plate in the dataframe 
growthcurves <- list()

#summarize growth curve for each well in dataframe 
for(i in 1:length(unique(df.bac$index2))){
  growthcurves[[i]] <- df.bac %>%
    filter(index2 == unique(.$index2)[i]) %>%
    mutate(time = day) %>%
    dplyr::select(time, od620, well) %>%
    pivot_wider(names_from = well, values_from = od620) %>%
    #rename_with(~str_replace_all(., c('G7'='blank', 'G1'='blank'))) %>%
    SummarizeGrowthByPlate() %>%
    filter(sample != "")
  #join logistic growth curve models with experimental metadata
  growthcurves[[i]] <- df.bac %>%
    filter(index2 == unique(.$index2)[i]) %>%
    dplyr::select(well, row, column, media, taxa, bacterial_strain, index2, replicate, oo, -index, -index2) %>%
    dplyr::rename(sample = well) %>%
    group_by(sample) %>%
    sample_n(1) %>%
    ungroup() %>%
    right_join(growthcurves[[i]], by = "sample")
}

#generate dataframe with growth curve summary statistics 
growthcurves <- rbindlist(growthcurves) %>%
  arrange(taxa, media, oo, replicate) %>%
  mutate(oo = as.factor(oo))

#classify taxa as growers or non-growers 
growthcurves <- df.bac %>%
  group_by(taxa, media, oo) %>%
  summarize(min = min(od620),
            max = max(od620),
            growth.rate = max/min) %>%
  group_by(taxa) %>%
  filter(growth.rate == max(growth.rate)) %>%
  mutate(growth = ifelse(growth.rate > 2, "1", "0")) %>%
  dplyr::select(taxa, growth) %>%
  right_join(., growthcurves, by = "taxa")

bac.growthplots <- list()

# Find the unique values for each variable 
sapply(lapply(growthcurves, unique), length)

for(i in 1:length(unique(growthcurves$taxa))){
  
  #conduct anova test 
  res.aov.i <- growthcurves %>%
    filter(taxa == unique(.$taxa)[i]) %>%
    anova_test(auc_l~ oo * media)
  
  #report anova results
  print(unique(growthcurves$taxa)[i])
  print(res.aov.i)
  
  #conduct post-hoc test
  pwc.i <- growthcurves %>%
    filter(taxa == unique(.$taxa)[i]) %>%
    group_by(media) %>%
    emmeans_test(auc_l ~ oo, p.adjust.method = "none") %>%
    add_y_position() %>%
    add_x_position(x = "media", dodge = 1)
  
  p <- growthcurves %>%
    filter(taxa == unique(.$taxa)[i]) %>%
    mutate(media = recode(media, "keratin" = "Keratin", "m9" = "M9"),
           oo = factor(recode(oo, "1" = "Present", "0" = "Absent"), levels = c("Present", "Absent"))) %>%
    ggplot(aes(media, auc_l)) +
    geom_boxplot(aes(fill = oo), position = position_dodge(1)) + 
    labs(x = "Media", 
         y = "Bacterial growth area\nunder the curve", 
         title = unique(growthcurves$taxa)[i], 
         fill = expression(italic("O. ophidiicola"))) +
    theme(legend.title=element_blank(),
          plot.title = element_text(face="italic", size=18),
          axis.text.y   = element_text(size=14),
          axis.text.x   = element_text(size=14),
          axis.title.y  = element_text(size=14, margin = margin(r = 15)),
          axis.title.x  = element_text(size=14, margin = margin(t = 15)),) + 
    scale_fill_manual(values = c("#bbaf71", "#3c4d6e")) +
    ylim(0,3.5)
  
  #if(
    #all(filter(growthcurves, taxa == unique(growthcurves$taxa)[i])$growth == "1"))
    #{
    bac.growthplots[[i]] <- p + stat_pvalue_manual(pwc.i, label = "p.adj.signif", hide.ns = F) 
  #} else {
   # bac.growthplots[[i]] <- p 
  #}
}

######Plot######

#Full plots
tiff("Figure5C_Full.tiff", width = 28, height = 4, units = 'in', res = 300)

ggarrange(plotlist = bac.growthplots[c(1,2,3,4,5,6,7)], 
          common.legend = T, 
          nrow = 1, ncol = 7,
          legend = "right")

dev.off()

#####Fungal growth analysis####

######Analyses######

#create dataframe to study fungal growth growth 
df.fg <- datum_bfi %>%
  filter(oo == 1) %>%
  group_by(taxa, media, replicate) %>%
  sample_n(1) %>%
  arrange(taxa, media, replicate) %>%
  mutate(oo_growth = 0,
         day = 0,
         od620 = NA) %>%
  dplyr::select(-row, -column, -plate, -index, -index2) %>%
  rbind(., dplyr::select(datum_bfi %>% filter(oo == 1), names(.))) %>%
  ungroup()

#Check coverage (lack NA values) across timeseries
df.fg %>% 
  mutate(day = as.factor(day)) %>%
  group_by(day) %>% 
  dplyr::summarise(sumNA = sum(is.na(oo_growth)),
                   total = length(oo_growth)) %>%
  mutate(coverage = (1-sumNA/total)*100)

#cursory examination of fungal growth
df.fg %>%
  filter(day <= 3) %>%
  ggplot(aes(x = day, y = oo_growth, color = media, lty = media)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~taxa)

#pull taxa for multipanel figure
df.fg.tax <- df.fg %>%
  mutate(media = as.factor(recode(media, "keratin" = "Keratin", "m9" = "M9")),
         taxa = as.factor(recode(taxa, "blank" = "Blank"))) 

df.fg.tax <- df.fg.tax %>%
  filter(taxa %in% c(factor(df.fg.tax$taxa,c("Blank", "Acinetobacter guillouiae", "Chryseobacterium sp.", "Flavobacterium odoratimimum", "Lelliottia nimipressuralis", "Mammaliicoccus sciuri", "Stenotrophomonas sp."))))
 
#character string to specify taxa below
tax <- df.fg.tax %>%
  filter(taxa != "Blank") %>%
  dplyr::select(taxa) %>%
  unique() %>%
  pull()

oogrowth_plots <- list() 

#looping code to generate multipanel figure
for(i in 1:(length(unique(df.fg.tax$taxa))-1)){
  
  dat <- df.fg.tax %>%
    filter(taxa == "Blank" | taxa == tax[i])
  
  oogrowth_plots[[i]] <- dat %>%
    ggplot(aes(x = day, y = oo_growth, color = fct_inorder(taxa))) + 
    geom_point() + 
    geom_smooth(data = dat %>% filter(taxa == "Blank"), method = "lm", se = F, lty = "dashed") + 
    geom_smooth(data = dat %>% filter(taxa == tax[i]), method = "lm", se = T, alpha = 0.25) +
    labs(title = tax[i], 
         x = "Time (days)",
         y = "Fungal Growth Radius (cm)") + 
    facet_wrap(~media) + 
    scale_linetype_manual(guide = "none") + 
    theme(legend.title=element_blank(),
          plot.title = element_text(face="italic", size=18),
          axis.title.y  = element_text(size=20, margin = margin(t = 5))) + 
    scale_color_manual(values = c("#bbaf71","gray"))
  
}

#effect of blanks 
fg.media.lm <- lm(oo_growth ~ day*media, data = filter(df.fg.tax, taxa == "Blank"))
summary(fg.media.lm)

df.fg.tax %>%
  filter(taxa == "Blank") %>%
  mutate(fit = predict(fg.media.lm, se.fit = T)$fit,
         se.upr = fit + predict(fg.media.lm, se.fit = T)$se.fit,
         se.lwr = fit - predict(fg.media.lm, se.fit = T)$se.fit) %>%
  ggplot(aes(x = day, y = fit)) + 
  geom_ribbon(aes(ymin = se.lwr, ymax = se.upr, group = media), alpha = 0.25) + 
  geom_line(aes(color = media), linewidth = 1) + 
  geom_point(aes(y = oo_growth, color = media)) + 
  labs(title = "Fungal Growth without Bacteria",
       subtitle = "Interaction term not significant - p = 0.23") + 
  theme(plot.title = element_text(size=25, hjust = 0.0))

df.fg.tax$media <- as.factor(df.fg.tax$media)

#Main analysis comparison 
fg.lm <- lm(oo_growth ~ day*relevel(media, ref = "M9")*relevel(taxa, ref = "Blank"), data = df.fg.tax)
summary(fg.lm)

fg.lm.pr <- as.data.frame(emtrends(fg.lm, pairwise ~ media*taxa, var = "day")$contrasts) %>%
  mutate(significance = ifelse(p.value < 0.05, "significant", "non-significant")) %>%
  arrange(significance) %>%
  mutate(contrasts = gsub("[()]", "", contrast)) %>%
  separate(col = contrasts, sep = " - ", into = c("con.1", "con.2")) %>%
  dplyr::select(-contrast) %>%
  separate(col = con.1, sep = " ", into = c("media.1", "taxa.1")) %>%
  separate(col = con.2, sep = " ", into = c("media.2", "taxa.2")) %>%
  dplyr::select(taxa.1, media.1, taxa.2, media.2, p.value, significance) %>% 
  arrange(taxa.1, media.1)
fg.lm.pr

#In M9 trendlines shouldn't differ 
fg.lm.pr.m9 <- fg.lm.pr %>%
  filter(media.1 == "M9" & media.2 == "M9") %>%
  filter(taxa.1 == "Blank" | taxa.2 == "Blank")
fg.lm.pr.m9

#Expected difference in slope between blank and bacteria presence growth on keratin media
fg.lm.pr.keratin <- fg.lm.pr %>%
  filter(media.1 == "Keratin" & media.2 == "Keratin") %>%
  filter(taxa.1 == "Blank" | taxa.2 == "Blank")
fg.lm.pr.keratin

#comparison of keratin and M9 slopes within a bacterial treatment  
fg.lm.pr.within.taxa <- fg.lm.pr %>%
  filter(media.1 != media.2) %>%
  filter(taxa.1 == taxa.2)
fg.lm.pr.within.taxa

#Examing the effect on oo_mass (non-significant)
fg.lm.mass <- lm(oo_mass ~ day*media*taxa, data = df.fg.tax)
summary(fg.lm.mass)

fg.lm.mass.pr <- as.data.frame(emtrends(fg.lm.mass, pairwise ~ media*taxa, var = "day")$contrasts) %>%
  mutate(significance = ifelse(p.value < 0.05, "significant", "non-significant")) %>%
  arrange(significance) %>%
  mutate(contrasts = gsub("[()]", "", contrast)) %>%
  separate(col = contrasts, sep = " - ", into = c("con.1", "con.2")) %>%
  dplyr::select(-contrast) %>%
  separate(col = con.1, sep = " ", into = c("media.1", "taxa.1")) %>%
  separate(col = con.2, sep = " ", into = c("media.2", "taxa.2")) %>%
  dplyr::select(taxa.1, media.1, taxa.2, media.2, p.value, significance) %>% 
  arrange(taxa.1, media.1)
fg.lm.mass.pr

#In M9 trendlines shouldn't differ 
fg.lm.mass.pr.m9 <- fg.lm.mass.pr %>%
  filter(media.1 == "M9" & media.2 == "M9") %>%
  filter(taxa.1 == "Blank" | taxa.2 == "Blank")
fg.lm.mass.pr.m9

#Expected difference in slope between blank and bacteria presence growth on keratin media
fg.lm.mass.pr.keratin <- fg.lm.mass.pr %>%
  filter(media.1 == "Keratin" & media.2 == "Keratin") %>%
  filter(taxa.1 == "Blank" | taxa.2 == "Blank")
fg.lm.mass.pr.keratin

#comparison of keratin and M9 slopes within a bacterial treatment  
fg.lm.mass.pr.within.taxa <- fg.lm.mass.pr %>%
  filter(media.1 != media.2) %>%
  filter(taxa.1 == taxa.2)
fg.lm.mass.pr.within.taxa

######Plot######

#select relevant plots and arrange into a multipanel figure 

tiff("Figure5B_Full.tiff", width = 30, height = 6, units = 'in', res = 300)

ggarrange(plotlist = oogrowth_plots,
          common.legend = F,
          nrow = 1, ncol = 6,
          legend = "bottom")

dev.off()

####Exp 4: Spent media (indirect fungal-bacterial interactions)####

#####Bacterial growth on fungal spent media#####

#Set plotting theme
theme_set(theme_pubr() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text.y   = element_text(size=20),
                  axis.text.x   = element_text(size=20),
                  axis.title.y  = element_text(size=20, margin = margin(r = 15)),
                  axis.title.x  = element_text(size=20, margin = margin(t = 15)),
                  legend.text = element_text(size=15),
                  legend.title = element_text(size = 20),
                  plot.title = element_text(size=25, hjust = 0.5),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")))

#Clean up dataset and group by GENUS
datum_sm <- datum_sm %>% 
  replace_na(list(sp_id = "blank")) %>% #rename the species id for blank columns
  mutate(sp_id = if_else(sp_id == 'BR1.10', 'Bacteroidia/Chryseobacterium/BR1.10', sp_id)) %>% #change the species id to post-sequencing genus name
  mutate(sp_id = if_else(sp_id == 'BR1.1', 'Bacteroidia/Chryseobacterium/BR1.1', sp_id)) %>% 
  mutate(sp_id = if_else(sp_id == 'BR1.11', 'Gammaproteobacteria/Lelliottia nimipressuralis/BR1.11', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.12', 'Gammaproteobacteria/Acinetobacter guillouiae/BR1.12', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.2', 'Gammaproteobacteria/Acinetobacter dispersus/BR1.2', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.4', 'Bacteroidia/Empedobacter tilapiae/BR1.4', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.5', 'Bacteroidia/Flavobacterium odoratum/BR1.5', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.6', 'Bacilli/Mammaliicoccus sciuri/BR1.6', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.7', 'Bacilli/Mammaliicoccus sciuri/BR1.7', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.8', 'Gammaproteobacteria/Acinetobacter guillouiae/BR1.8', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.9', 'Gammaproteobacteria/Acinetobacter guillouiae/BR1.9', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.1', 'Alphaproteobacteria/Sphingomonas taxi/EKS17.1', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.15', 'Bacilli/Pristimantibacillus/EKS17.15', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.19', 'Actinomycetia/Microbacterium/EKS17.19', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.8', 'Gammaproteobacteria/Stenotrophomonas/EKS17.8', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.29', 'Actinomycetia/Streptomyces/EKS17.29', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.1', 'Actinomycetia/Dermatophilaceae/EKS20.1', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.10', 'Gammaproteobacteria/Stenotrophomonas/EKS20.10', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.3', 'Bacilli/Lysinibacillus/EKS20.3', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.4', 'Actinomycetia/Rhodococcus erythropolis/EKS20.4', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.5', 'Gammaproteobacteria/Stenotrophomonas/EKS20.5', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.8', 'Actinomycetia/Paenarthrobacter/EKS20.8', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-5.1', 'Gammaproteobacteria/Aeromonas hydrophila/TR087-5.1', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-5.3', 'Gammaproteobacteria/Morganella morganii/TR087-5.3', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-6.1', 'Gammaproteobacteria/Serratia liquefaciens/TR087-6.1', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.2', 'Bacteroidia/Flavobacterium odoratimimum/TR087-7.2', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.4', 'Gammaproteobacteria/Morganella morganii/TR087-7.4', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.6', 'Gammaproteobacteria/Morganella morganii/TR087-7.6', sp_id)) %>%
  unite(col='sample_name', c('row', 'column','media', 'sp_id', 'experiment'), sep = '_', remove = FALSE) %>% #create a single column that can create individual ID for all wells
  group_by(sample_name) %>% #group by the new plate ID
  filter(hour > 16) %>% #keep only those that were measured for all 48 hours
  ungroup() #ungroup the dataframe

#Clean up dataset and group by CLASS
datum_sm <- datum_sm %>% 
  replace_na(list(sp_id = "blank")) %>% #rename the species id for blank columns
  mutate(sp_id = if_else(sp_id == 'BR1.10', 'Bacteroidia', sp_id)) %>% #change the species id to post-sequencing genus name
  mutate(sp_id = if_else(sp_id == 'BR1.1', 'Bacteroidia', sp_id)) %>% 
  mutate(sp_id = if_else(sp_id == 'BR1.11', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.12', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.2', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.4', 'Bacteroidia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.5', 'Bacteroidia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.6', 'Bacilli', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.7', 'Bacilli', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.8', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'BR1.9', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.1', 'Alphaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.15', 'Bacilli', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.19', 'Actinomycetia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.8', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS17.29', 'Actinomycetia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.1', 'Actinomycetia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.10', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.3', 'Bacilli', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.4', 'Actinomycetia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.5', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'EKS20.8', 'Actinomycetia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-5.1', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-5.3', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-6.1', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.2', 'Bacteroidia', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.4', 'Gammaproteobacteria', sp_id)) %>%
  mutate(sp_id = if_else(sp_id == 'TR087-7.6', 'Gammaproteobacteria', sp_id)) %>%
  unite(col='sample_name', c('row', 'column','media', 'sp_id', 'experiment'), sep = '_', remove = FALSE) %>% #create a single column that can create individual ID for all wells
  group_by(sample_name) %>% #group by the new plate ID
  filter(hour > 16) %>% #keep only those that were measured for all 48 hours
  ungroup() #ungroup the dataframe

######Analyses######

#cursory examination of bacterial growth curves in blank wells
datum_sm %>%
  filter(taxa != "blank") %>%
  ggplot(aes(x = hour, y = od620, color = as.factor(media))) + 
  geom_smooth(method = "loess") +
  geom_point() +
  facet_wrap(~sp_id)

# Averaging across blank wells should capture any variation in blank wells fairly well 
datum_sm %>%
  filter(taxa == "blank") %>%
  ggplot(aes(hour, od620)) + 
  geom_point() +
  geom_smooth(aes(color = media), method = "loess", show.legend = F) + 
  facet_wrap(~media, scales = "free") 

#average od620 measurements for all blank measurements at all time-points
datum_sm.blank <- datum_sm %>%
  filter(taxa == "blank") %>%
  group_by(experiment, media, hour) %>% 
  summarise_at(vars(od620), list(od620.blank = mean)) %>%
  right_join(datum_sm, by = c("experiment", "media", "hour")) %>%
  ungroup() %>%
  mutate(od620 = ifelse(taxa == "blank", od620.blank, od620)) %>%
  dplyr::select(-od620.blank) 

#remove repeat blank columns 
datum_sm.blank_singles <- datum_sm.blank %>%
  filter(taxa == "blank") %>%
  group_by(experiment, media, hour) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  rbind(filter(datum_sm.blank, taxa != "blank"))

#clean up working environment
datum_sm.bac <- datum_sm.blank_singles
functions <- as.vector(lsf.str())
gdata::keep(list = c('datum_sm', 'datum_sm.bac', functions), sure = T)

#write loop to summarize curve for each well in each plate in the dataframe 
growthcurves1 <- list()

# generate variable that will allow loop to run on a plate/media basis
datum_sm$index <- interaction(datum_sm$testid, datum_sm$media)

#summarize growth curve for each well in dataframe 
for(i in 1:length(unique(datum_sm$index))){
  growthcurves1[[i]] <-
    filter(datum_sm, index == unique(datum_sm$index)[i]) %>%
    mutate(time = hour) %>%
    dplyr::select(time, od620, well) %>%
    pivot_wider(names_from = well, values_from = od620) %>%
    #rename_with(~str_replace_all(., c('D4'='blank', 'D7'='blank', 'D10'='blank', 'D1'='blank', 'E1'='blank'))) %>%
    SummarizeGrowthByPlate() %>%
    filter(sample != "")
  #join logistic growth curve models with experimental metadata
  growthcurves1[[i]] <- 
    filter(datum_sm, index == unique(datum_sm$index)[i]) %>%
    dplyr::select(well, row, column, media, taxa, sp_id, bac_sp, index, testid) %>%
    rename(sample = well) %>%
    group_by(sample) %>%
    sample_n(1) %>%
    ungroup() %>%
    right_join(growthcurves1[[i]], by = "sample")
}

#generate dataframe with growth curve summary statistics 
growthcurves1 <- rbindlist(growthcurves1) %>%
  rename(well = sample) %>%
  arrange(testid, desc(sp_id), media) %>%
  as.data.frame()

#Clean up working environment
functions <- as.vector(lsf.str())
gdata::keep(list = c('growthcurves1', functions), sure = T)

# Find the unique values for each variable 
sapply(lapply(growthcurves1, unique), length)

#Clean up dataset for plotting
df.gc.tax <- growthcurves1 %>%
  mutate(media = fct_recode(as.factor(media), "m9" = "M9", "keratin" = "1%-Keratin", "5d" = "5d-Keratin", "9d" = "9d-Keratin", '13d' = "13d-Keratin")) %>% 
  mutate(media = factor(media, levels = c("m9", "keratin", "5d", "9d", "13d")))

#Check reference levels
levels(df.gc.tax$media)

#Linear model 
df.gc.lm <- lm(auc_l ~ media*sp_id, data = df.gc.tax)
summary(df.gc.lm)

#Run loop for anova and plots
plots <- list()

for(i in 1:length(unique(df.gc.tax$sp_id))){
  
  res.aov.i <- df.gc.tax %>%
    filter(sp_id == unique(df.gc.tax$sp_id)[i]) %>%
    anova_test(auc_l ~ media)
  
  pwc.i <- df.gc.tax %>%
    filter(sp_id == unique(df.gc.tax$sp_id)[i]) %>%
    emmeans_test(auc_l ~ media, p.adjust.method = "none") %>%
    add_xy_position(x = "media", scales = "free")
  
  plot.i <- df.gc.tax %>%
    filter(sp_id == unique(df.gc.tax$sp_id)[i]) %>%
    ggbarplot(x = "media", y = "auc_l", fill = "media", add = "mean_se") +
    scale_fill_viridis(discrete= TRUE, option="cividis") +
    stat_pvalue_manual(pwc.i, hide.ns = T) +
    labs(
      x = "Media", 
      y = "Area under the curve", 
      title = str_wrap(unique(df.gc.tax$sp_id)[i],width=50)
    ) + 
    theme(legend.position = "none",
          axis.title.y  = element_text(size=15, margin = margin(r = 15)),
          axis.title.x  = element_text(size=15, margin = margin(t = 15))) + 
    coord_cartesian(ylim = c(0,39))
  
  plots[[i]] <- plot.i
  print(unique(df.gc.tax$sp_id)[i])
  print(res.aov.i)
}

######Plot######

#select relevant plots and arrange into a multipanel figure 

tiff("BaconFunSpentMedia.tiff", width = 20, height = 20, units = 'in', res = 300)

ggarrange(plotlist = plots[c(1:19)],
          ncol = 5, 
          nrow = 4) 

dev.off()

#Change columns from continuous to categories
datum_sm$treatment <- as.factor(datum_sm$taxa)

results_sm=lmer(od620 ~ taxa + hour + (1|testid), data = datum_sm,na.action=na.omit)
summary(results_sm)

tiff("Figure4A.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_sm, aes(x = taxa , y = od620, color= media)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#3c4d6e", "#7f7c75","#bbaf71", "grey", "black")) +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(face="italic", colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Species", y=(expression("OD620"))) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(color='Media type') 

dev.off()''

#Subset data by species

datum_sm1 <- subset(datum_sm, bac_sp=="Acinetobacter")
datum_sm2 <- subset(datum_sm, bac_sp=="Chryseobacterium")
datum_sm3 <- subset(datum_sm, bac_sp=="Stenotrophomonas")
datum_sm4 <- subset(datum_sm, taxa=="blank")

datum_sm5 <- datum_sm[ which( datum_sm$media == "5d-Keratin") , ]

#Analysis of bacterial growth curves

#Rename day column to time for GrowthCurver purposes
datum_sm5$time <- datum_sm5$hour

#Rearrange dataset to analyze
summarybac <- datum_sm5 %>% 
  group_by(taxa, time) %>%
  summarise(od620 = mean(od620))

summarybac <- pivot_wider(summarybac, names_from = taxa, values_from = od620)

gc_fit <- SummarizeGrowthByPlate(summarybac)
plot(gc_fit)
head(gc_fit)
output_file_name <- "test.txt"
write.table(gc_fit, file = test, 
            quote = FALSE, sep = "\t", row.names = FALSE)

datum_sm1 <- datum_sm %>%
  group_by(hour, media) %>% 
  summarise(
    sd = sd(od620, na.rm=TRUE),
    od620 = mean(od620),
    .groups="keep"
  )

datum_sm2 <- datum_sm %>%
  group_by(hour, media) %>% 
  summarise(
    sd = sd(od620, na.rm=TRUE),
    od620 = mean(od620),
    .groups="keep"
  )

datum_sm3 <- datum_sm %>%
  group_by(hour, media) %>% 
  summarise(
    sd = sd(od620, na.rm=TRUE),
    od620 = mean(od620),
    .groups="keep"
  )

datum_sm4 <- datum_sm %>%
  group_by(hour, media) %>% 
  summarise(
    sd = sd(od620, na.rm=TRUE),
    od620 = mean(od620),
    .groups="keep"
  )

tiff("FigureAcin.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_sm1, aes(x = as.numeric(hour) , y = od620, color= media)) + 
  geom_line() +
  scale_color_manual(values=c("#3c4d6e", "#7f7c75","#bbaf71", "grey", "black")) +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Hour", y=("OD620")) +
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(color='Media type') +
  ggtitle("Acinetobacter")

dev.off()

tiff("FigureChrys.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_sm2, aes(x = as.numeric(hour) , y = od620, color= media)) + 
  geom_line() +
  scale_color_manual(values=c("#3c4d6e", "#7f7c75","#bbaf71", "grey", "black")) +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Hour", y=("OD620")) +
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(color='Media type') +
  ggtitle("Chryseobacterium")

dev.off()

tiff("FigureSteno.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_sm3, aes(x = as.numeric(hour) , y = od620, color= media)) + 
  geom_line() +
  scale_color_manual(values=c("#3c4d6e", "#7f7c75","#bbaf71", "grey", "black")) +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Hour", y=("OD620")) +
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(color='Media type') +
  ggtitle("Stenotrophomonas")

dev.off()

tiff("FigureControl.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_sm4, aes(x = as.numeric(hour) , y = od620, color= media)) + 
  geom_line() +
  scale_color_manual(values=c("#3c4d6e", "#7f7c75","#bbaf71", "grey", "black")) +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  labs(x="Hour", y=("OD620")) +
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(color='Media type') +
  ggtitle("Control")

dev.off()

#####Fungal growth on bacterial spent media####

#Clean up dataset
datum_smj1 <- datum_smj %>% 
  unite(col='sample_name', c('bacteria', 'medium_treatment', 'base_medium', 'replicate'), sep = '_', remove = FALSE) %>% #create a single column that can create individual ID for all plates
  mutate(unique.ID = if_else(unique.ID == 'Control(K)', 'First\npassage control', unique.ID)) %>% #change the control names
  mutate(unique.ID = if_else(unique.ID == 'Control2(K)', 'Second\npassage control', unique.ID)) %>% #change the control names
  filter(!is.na(diameter..mm.)) %>% #remove all NA rows
  group_by(sample_name) %>% #group by the new plate ID
  #filter(n() == 7) %>% #keep only those that were measured for all 7 days
  ungroup() #ungroup the dataframe

#Check coverage (lack NA values) across timeseries
datum_smj1 %>% 
  mutate(day = as.factor(day)) %>%
  group_by(day) %>% 
  dplyr::summarise(sumNA = sum(is.na(diameter..mm.)),
                   total = length(diameter..mm.)) %>%
  mutate(coverage = (1-sumNA/total)*100)

#cursory examination of fungal growth
datum_smj1 %>%
  ggplot(aes(x = day, y = diameter..mm., color = bacteria, lty = bacteria)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~base_medium)

######Analyses#####

#Subset to the final day of growth
datum_smj2 <- datum_smj1[ which(datum_smj1$day == "14") ,]

#Convert columns
datum_smj2$bacteria <- as.factor(datum_smj2$bacteria)
datum_smj2$base_medium <- as.factor(datum_smj2$base_medium)
datum_smj2$unique.ID <- as.factor(datum_smj2$unique.ID)
datum_smj2$medium_treatment <- as.factor(datum_smj2$medium_treatment)
datum_smj2$day <- as.numeric(datum_smj2$day)

#Subset to single an combined media types (separated by controls)
datum_smj2_single <- datum_smj2[ which( datum_smj2$medium_treatment == "keratin" | datum_smj2$medium_treatment == "keratinContol") , ]
datum_smj2_mixed <- datum_smj2[ which( datum_smj2$medium_treatment == "keratin A into B" | datum_smj2$medium_treatment == "Keratin B into A" | datum_smj2$medium_treatment == "Keratin A and B" |datum_smj2$medium_treatment == "Keratin Control2") , ]

#Check reference levels
levels(datum_smj2_single$bacteria)
levels(datum_smj2_mixed$medium_treatment)

#Analysis of fungal growth at 14 days in single media types
results_sm=lmer(diameter..mm. ~ relevel(base_medium, ref="M9") + relevel(bacteria, ref="None") + (1|replicate), data = datum_smj2_single, na.action=na.omit)
summary(results_sm)

#Analysis of fungal growth at 14 days in mixed media types
results_sm1=lmer(diameter..mm. ~ relevel(base_medium, ref="M9") + relevel(medium_treatment, ref="Keratin Control2") + (1|replicate), data = datum_smj2_mixed, na.action=na.omit)
summary(results_sm1)


######Plot######

#Subset to just keratin for plotting
datum_smj3_single <- subset(datum_smj2_single, base_medium=="keratin")
datum_smj3_mixed <- subset(datum_smj2_mixed, base_medium=="keratin")

#Reorder treatment groups
datum_smj3_single$bacteria_2 = factor(datum_smj3_single$bacteria,c("Chryseobacterium", "Stenotrophomonas", "None"))
datum_smj3_single$unique.ID_2 = factor(datum_smj3_single$unique.ID,c("A(K)", "B(K)", "First\npassage control"))
datum_smj3_mixed$bacteria_2 = factor(datum_smj3_mixed$bacteria,c("Chryseobacterium", "Stenotrophomonas", "StenotrophomonasChryseobacterium", "None"))
datum_smj3_mixed$unique.ID_2 = factor(datum_smj3_mixed$unique.ID,c("AintoB(K)", "BintoA(K)", "AplusB(K)", "Second\npassage control"))

#Adding p-values to figure
stat.test <- datum_smj3_single %>%
  t_test(diameter..mm. ~ unique.ID_2) %>%
  add_xy_position(x = "unique.ID_2")

stat.test1 <- datum_smj3_mixed %>%
  t_test(diameter..mm. ~ unique.ID_2) %>%
  add_xy_position(x = "unique.ID_2")

#Figure 6A and B

Plot1.J<- 
ggplot(datum_smj3_single, aes(x = unique.ID_2 , y = diameter..mm.)) + 
  geom_boxplot(fill=c("#bbaf71", "#3c4d6e", "#002051")) +
  geom_jitter(color="black") +
  theme_bw() +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")) +
  labs(x="Species", y=("Fungal growth (mm)"), title="Single bacterium media") +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test, hide.ns = T) +
  ylim(32,53)

Plot2.J<- 
  ggplot(datum_smj3_mixed, aes(x = unique.ID_2 , y = diameter..mm.)) + 
  geom_boxplot(fill=c("#bbaf71", "#3c4d6e", "#7f7c75", "#002051")) +
  geom_jitter(color="black") +
  theme_bw() +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")) +
  labs(x="Species", y=("Fungal growth (mm)"), title="Combined media") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test1, hide.ns = T)

tiff("FigureFunMedia.jpeg", width = 8, height = 4, units = 'in', res = 300)

ggarrange(Plot1.J, Plot2.J,
          font.label = list(size = 18),
          nrow = 1, ncol = 2)

dev.off()

datum_smj2 <- subset(datum_smj, base_medium=="M9")

tiff("FigureJohnM9.tiff", width = 8, height = 4, units = 'in', res = 300)

ggplot(datum_smj2, aes(x = unique.ID , y = diameter..mm., fill= bacteria)) + 
  geom_boxplot() +
  geom_jitter(color="black") +
  scale_fill_manual(values=c("#bbaf71", "#7f7c75", "#3c4d6e", "#002051")) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")) +
  labs(x="Species", y=("Area mm^2")) +
  labs(fill='Bacterial species') +
  ggtitle("Oo on bacterial spent media (M9) - Day 14")

dev.off()

ggplot(datum_sm2, aes(x = day, y = diameter..mm.)) + 
  geom_point(aes(color = bacteria)) +
  geom_smooth(aes(color = bacteria), method = 'lm', show.legend = F, alpha = 0.1) + 
  geom_smooth(aes(color = bacteria), method = 'lm', alpha = 0) + 
  scale_color_manual(name= "Bacteria", values = c("#3c4d6e", "#7f7c75","#bbaf71", "grey")) + 
  labs(title = "Oo on bacterial spent media",
       x = "Day", 
       y = "Growth (mm)") + 
  guides(shape = "none") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme_classic(12.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")) 

#select relevant plots and arrange into a multipanel figure 

tiff("Figure5B.tiff", width = 12, height = 8, units = 'in', res = 300)

datum_smj1 %>%
  ggplot(aes(x = day, y = diameter..mm., color = bacteria)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Fungal growth on bacterial spent media", 
       x = "Time (days)",
       y = "Fungal Growth Radius (cm)") + 
  facet_wrap(~base_medium) + 
  scale_linetype_manual(guide = "none") + 
  theme(legend.title=element_blank()) 
#scale_color_manual(values = c(brewer.pal(n=5,"Greys")[3], c("#bbaf71")))

dev.off()