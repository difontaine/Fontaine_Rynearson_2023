#### Script for processing data for Chapter 1 manuscript
## Field Sample Processing ###
#D. Fontaine April 13th, 2023


#Field Sample Analysis for 18S 6-year time series from Narragansett Bay
library(phyloseq)
library(plyr)
library(reshape2)
library(tidyverse)
library(dplyr)
library(data.table)
library(scales)
library(maptools)
library(raster)
library(prettymapr)
library(ggsci)
library(gridExtra)
library(ggtext)
library(car)
library(lubridate)
library(caret)
library(lubridate)
library(vegan)
library(cluster)
library(factoextra)
library(ggdendro)
library(grid)
library(ggpubr)

##Theme for plotting
y <- theme(legend.key = element_rect(fill = "white") ,
           legend.background = element_blank(), #To delete any of these, just do element_b 
           legend.text=element_text(family = "Times"), legend.title=element_text(family = "Times", face="bold"), panel.background = element_rect(fill = "white", colour="black"), panel.grid.major = element_line(colour = "white", linewidth=.5),
           panel.grid.minor = element_line(colour = "white", linewidth=.5),
           legend.position = "right",
           strip.background = element_rect(fill = "white", color = "black", size = .5), 
           plot.title=element_text(family = "Times"),
           axis.text = element_text(family = "Times"),
           axis.title = element_text(family = "Times", face = "bold"))

### SUPPLEMENTAL FIGURE 1: ENVI PLOTS ###
envi_data <- read_csv("Envdata.csv")
#get month info
envi_data$month <- month(envi_data$Date)
envi_data$year <- year(envi_data$Date)

#group by month and calculate average temp, DIN and sal
env_data_summarized <- envi_data%>%
  group_by(month) %>%
  summarise(avg_temp = mean(Temp), sd_temp = sd(Temp), avg_DIN = mean(DIN_μM), sd_DIN = sd(DIN_μM), avg_sal = mean(Salinity), sd_sal = sd(Salinity))

#plot Temp
temp_plot <-  ggplot(env_data_summarized, aes(as.factor(month), avg_temp))+
  geom_line(size = 1, color='black', group = 1)+
  geom_point(aes(x=month, y=avg_temp), size = 2.5, data= env_data_summarized)+
  y+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size =15))+
  geom_errorbar(aes(ymin=avg_temp-sd_temp, ymax=avg_temp+sd_temp), width = 0.2)+
  labs(x = "Month", y = "Average monthly temperature (˚C) (± s.d.)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5))+
  annotate(geom = "text", x = 1, y =  25, label = "a)", fontface =2, size = 6)
#ylim(c(-1,3.5))
temp_plot  
ggsave('Manuscript/Supplemental_Info/Figures/avg_temp.png', plot = temp_plot ,
       scale = 1, width = 5.5, height = 5 , units ="in",
       dpi = 300)

#env_data_summarized[,-1]

#plot DIN
DIN_plot <-  ggplot(env_data_summarized, aes(factor(month), avg_DIN))+
  geom_line(size = 1, color='black', group = 1)+
  geom_point(aes(x=, y=avg_DIN), size = 2.5, data= env_data_summarized)+
  y+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size =15))+
  geom_errorbar(aes(ymin=avg_DIN-sd_DIN, ymax=avg_DIN+sd_DIN), width = 0.2)+
  labs(x = "Month", y = "Average monthly DIN (˚µM) (± s.d.)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5))+
 annotate(geom = "text", x = 1, y =  15, label = "b)", fontface =2, size = 6)
#ylim(c(-1,3.5))
DIN_plot  
ggsave('Manuscript/Supplemental_Info/Figures/avg_DIN.png', plot = DIN_plot   ,
       scale = 1, width = 5.5, height = 5 , units ="in",
       dpi = 300)

### RAW DATA PROCESSING
field_count_tab <- read.table("zimmerman_prim_field_ASVs_4.14_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")
#read in sample info table -- need it to be a tsv for later on analysis otherwise phyloseq doesn't work
sample_info_tab_field <- read.table("sample_info_field.tsv", header=T, row.names=1,
                                    check.names=F, sep="\t") #Note to create a .tsv file, you need to save the excel spreadsheet as a regular .txt file and then change the extension of that file to ".tsv"

#read in taxonomy table
asv_tax_field <- read.table("all_diat_4.14_updatedNov21.tsv", header=T, row.names=1,
                            check.names=F, sep="\t")

#Make asv_tax a matrix for physeq object to work
asv_tax_field <- as.matrix(asv_tax_field)

#Make samplinfo phy object
sample_info_tab_phy <- sample_data(sample_info_tab_field)

# first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(field_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(asv_tax_field) #Have to make the taxa table a matrix in order for the physeq object thing to work
#Making a physeq object
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

#Subset physeq object for just the bacillariophyta class
all_diat <- subset_taxa(ASV_physeq, tax.Class =="Bacillariophyta")

all_diat_df <- as.data.frame(otu_table(all_diat))


#Count the number of sequences that were ID'ed as "Bacillariophyta
all_diat_df_counts <- as.data.frame(colSums(all_diat_df))

#get the min, max and mean for the all_diat_df_counts dataframe which gives the reads that were ID'ed as Diatoms and kept for downstream processing
min(all_diat_df_counts$`colSums(all_diat_df)`)
max(all_diat_df_counts$`colSums(all_diat_df)`)
mean(all_diat_df_counts$`colSums(all_diat_df)`)

#add sample name for csv
all_diat_df_counts<- rownames_to_column(all_diat_df_counts, var = "Sample")

#get sample info
sample_info_tab_field_counts <- rownames_to_column(sample_info_tab_field, var = "Sample")

all_diat_df_counts <- all_diat_df_counts %>%
  left_join(sample_info_tab_field_counts, by = "Sample")

#save csv to add to supplemental Table 1
write_csv(all_diat_df_counts, "diatom_reads.csv")

#Add ASV as column
all_diat_df <- rownames_to_column(all_diat_df, var = "ASV")

all_diat_df_taxa <- as.data.frame(tax_table(all_diat))
all_diat_df_taxa <- rownames_to_column(all_diat_df_taxa, var = "ASV")


all_diat_df <- all_diat_df %>%
  pivot_longer(cols = 2:83,
               names_to = "Sample",
               values_to = "count")

# count up the total number of reads per sample
diat_counts_total <- ddply(all_diat_df, .(Sample), summarise, total_reads = sum(count))

all_diat_df <- all_diat_df %>%
  left_join(diat_counts_total, by = 'Sample') %>%
  mutate(per_rel = (count/total_reads)*100)

#filter reads for cut-off threshold based on mock community analysis (drop everything that is below 0.075% relative abundance). See "read_threshold.R" script to see how this was done
all_diat_df <- all_diat_df %>%
  filter(per_rel >0.075)

#add taxonomy information to the ASV file
all_diat_df_with_taxa <- all_diat_df %>%
  left_join(all_diat_df_taxa, by = "ASV")

#it is only present if it has a % read rel abundance greater than 0.075%...this is similar to above but this gives each row a "presence" info based on the % relative abundance threshold
all_diat_df_with_taxa$presence <- sapply(all_diat_df_with_taxa$per_rel, function(x) {ifelse(any(x > 0.075), 1, 0)})


#need to statistically test the field sample replicates
field_sample_reps <- all_diat_df_with_taxa %>%
  filter(Sample %in% c("SF80b", "SF81b", "SF82b")) %>%
  select(ASV,Sample,per_rel)

field_sample_reps$presence <- 1

#do anova to test for differences in presence b/t field samples
fst <- aov(presence ~ Sample, data = field_sample_reps)
summary(fst)
TukeyHSD(fst, which = "Sample")


#Figuring out how many ASVs per taxa level 
all_diat_selected <- all_diat_df %>%
  select(ASV) %>%
  distinct()

all_diat_selected <- all_diat_selected %>% #so 658 distinct ASVs 
  inner_join(all_diat_df_taxa, by = "ASV")

all_diat_selected_sp <- all_diat_selected %>%
  separate(tax.Species, into = c('G', 'ext'), sep = "_") %>%
  filter(grepl('sp', ext)) #205 were not assigned to species but then minus 3 for those assigned to species with "sp" so, 202 total were only assigned to genus

all_diat_selected_nosp <- all_diat_selected %>%
  separate(tax.Species, into = c('G', 'ext'), sep = "_") %>%
  filter(!grepl('sp', ext)) %>% #344 + 3 species that have the "sp" in their name = 347
  drop_na()

all_diat_select_fam <- all_diat_selected %>%
  filter(is.na(tax.Genus)) #109 - 8 NA at order = 101 at fam

all_diat_select_order <- all_diat_selected %>%
  filter(is.na(tax.Family))

### FIGURE 1: MAP ###
#Load USA and individual state shapefile data
USA <- getData('GADM', country='USA', level=1)
RI <- (USA[USA$NAME_1=="Rhode Island",])
MA <- (USA[USA$NAME_1=="Massachusetts",])
CT <- (USA[USA$NAME_1=="Connecticut",])
NY <- (USA[USA$NAME_1=="New York",])
NJ <- (USA[USA$NAME_1=="New Jersey",])
PA <- (USA[USA$NAME_1=="Pennsylvania",])
MD <- (USA[USA$NAME_1=="Maryland",])

#Plot zoomed out
pdf('Manuscript/Figs/zoomed_out_map.pdf')
plot(RI, bg = "white", border = "black", xaxt = "n", yaxt = "n", col = "grey70", ylim = c(40, 42.5), xlim = c(-73,-70), cex.axis=1)
plot(MA, axes = T, col = "grey70", border = "black", add = T)
plot(CT, axes = T, col = "grey70", border = "black", add = T)
plot(NY, axes = T, col = "grey70", border = "black", add = T)
plot(NJ, axes = T, col = "grey70", border = "black", add = T)
plot(PA, axes = T, col = "grey70", border = "black", add = T)
plot(MD, axes = T, col = "grey70", border = "black", add = T)

#save zoomed out plot
dev.off()

#Plot zoomed in
pdf('Manuscript/Figs/zoomed_in_map.pdf')
plot(RI, bg = "white", border = "black", axes = T, col = "grey90", ylim = c(41.3, 41.9), xlim = c(-71,-71))
plot(MA, axes = T, col = "grey90", border = "black", add = T)
plot(CT, axes = T, col = "grey90", border = "black", add = T)
plot(NY, axes = T, col = "grey90", border = "black", add = T)

points(-71.4, 41.567, pch =21, bg = "black", cex = 1, lwd = 1)

#I didn't end up adding an arrow and scale bar so these may need to be adjusted for size, but they're included here if you want them
library(prettymapr)
addnortharrow(pos = "topright", padin = c(0.1,0.1), scale = 1, lwd = 1, border = "black", cols = c("white","black"), text.col = "black")
addscalebar(plotunit = NULL, plotepsg = NULL, widthhint = 0.25, unitcategory = "metric", htin = 0.1, padin = c(2,0.2), style = "bar", bar.cols = c("black","white"),
            lwd = 1, linecol = "black", tick.cex = 0.5, labelpadin = 0.15, label.cex = 1.2, label.col = "black")
dev.off()


#filter and only keep ASVs that have a percent rel abundance greater than 0.075 for the significance testing
mock_with_taxa$presence <- sapply(mock_with_taxa$per_rel, function(x) {ifelse(any(x > 0.075), 1, 0)})

#Add group info
mock_with_taxa_group1 <- mock_with_taxa %>%
  dplyr::filter(Sample %in% c("SF85b", "SF86b", "SF87b"))

mock_with_taxa_group1$MC <- "MC1"

mock_with_taxa_group2 <- mock_with_taxa %>%
  dplyr::filter(Sample %in% c("SF88b", "SF89b", "SF90b"))

mock_with_taxa_group2$MC <- "MC2"

mock_with_taxa_group3 <- mock_with_taxa %>%
  dplyr::filter(Sample %in% c("SF91b", "SF92b", "SF93b"))

mock_with_taxa_group3$MC <- "MC3"

mock_with_taxa_group4 <- mock_with_taxa %>%
  dplyr::filter(Sample %in% c("SF94b", "SF95b", "SF96b"))

mock_with_taxa_group4$MC <- "MC4"

mock_with_taxa_complete <- rbind(mock_with_taxa_group4, mock_with_taxa_group3,mock_with_taxa_group2,mock_with_taxa_group1)

mock_counts <- mock_with_taxa_complete %>%
  filter(Class == "Bacillariophyta") %>%
  group_by(Sample) %>%
  summarise(sum_counts = sum(per_rel))

#filter out ASVs with percent rel abundance of 0
mock_with_taxa_complete <- mock_with_taxa_complete %>%
  filter(presence == 1)

#read in sample info for mock samples
mock_samples <- read_csv("mock_samples/sample_info_mock.csv")

#generating plot of mock samples
mock_comp <- read_csv("mock_samples/mock_taxa_filtered.csv") #Read in file that has the observed and expected mock %s

mock_comp <- mock_comp %>%
  dplyr::select(MC, Replicate, Species, per_rel, Sample)

mock_exp <- mock_comp %>%
  filter(Replicate == "Exp") %>%
  filter(MC != "MC2") %>% #filter out MC2 since we're not using it 
  dplyr::rename(avg_rel = per_rel) %>%
  select(-Sample)

mock_exp$Replicate <- "Exp."
mock_exp$sd <- "NA"


#only plotting the average replicate information
mock_comp_noexp <- mock_comp %>% 
  filter(Replicate != "Exp") %>%
  dplyr::group_by(MC,Species,Replicate, Sample) %>%
  dplyr::summarise(total_rel = sum(per_rel)) %>%
  dplyr::rename(avg_rel = total_rel)

mock_1 <- mock_comp_noexp %>% 
  filter(MC== "MC1")

mock_1_aov <- aov(avg_rel ~ Sample, data = mock_1)
summary(mock_1_aov) #mock 1 is 0.98

mock_3 <- mock_comp_noexp %>% 
  filter(MC== "MC3")

mock_3_aov <- aov(avg_rel ~ Sample, data = mock_3)
summary(mock_3_aov) #mock 2 for paper (really mock 3) is 0.97 (only used 3 mock communities)

mock_4 <- mock_comp_noexp %>% 
  filter(MC== "MC4")

mock_4_aov <- aov(avg_rel ~ Sample, data = mock_4)
summary(mock_4_aov) #mock 3 in paper (really mock 4) is 0.97 (only used 3 mock communities))

mock_avg <- mock_comp_noexp %>%
  dplyr::group_by(MC,Species) %>%
  dplyr::summarise(sd = sd(avg_rel, na.rm = TRUE), avg_rel = mean(avg_rel))


mock_additional <- ddply(mock_avg, .(MC), summarise, avg_rel = 100 - sum(avg_rel))
mock_additional$Species <- '<0.075% of Reads'
mock_additional$sd <- 'NA'

total_mock <- rbind(mock_additional , mock_avg)

total_mock$Replicate <- "Obs."

total_mock <- rbind(total_mock, mock_exp)

#filter for MC1, MC3 and MC4
total_mock_supp_table <- total_mock %>%
  filter(MC %in% c("MC1", "MC3", "MC4")) %>%
  filter(Species != "<0.075% of Reads")

#save for supplemental table
write_csv(total_mock_supp_table, "mock_samples/mock_supp_table.csv" )


#do anova to test for difference between observed and expected
two_way_mock_1 <- aov(avg_rel ~ Replicate*Species, data = total_mock)
summary(two_way_mock_1)

#test that the percent relative abundances add up to 100 for the MCs 
test_sum <- total_mock %>%
  filter(Replicate != "Exp.") %>%
  dplyr::group_by(MC) %>%
  dplyr::summarise(sum = sum(avg_rel))


total_mock$MC <- factor(total_mock$MC, levels = c("MC1", "MC2", "MC3", "MC4"))

total_mock$Species <- with(total_mock, ifelse(Species == "<0.075% of Reads", "<0.075% of Reads", 
                                              paste0("*",Species,"*")))

#names for MC communities -- need it to be ordered for the paper (the Mc1 is 1, the MC3 is 2 and the MC4 is 3)
mock_names <- c(
  `MC1` = "MC 1",
  `MC3` = "MC 2",
  `MC4` = "MC 3"
)



## FIGURE 2 PLOT
mock_plot <- total_mock%>%
  filter(MC !="MC2")%>% #filter out MC2 because want to select MC1, 3 an 4 since these three samples have the most different proportions ofcells
  ggplot(aes(x=Replicate, y=avg_rel, fill=Species)) +
  geom_bar(width=1, stat="identity") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))+
 facet_wrap(~MC, nrow = 1, labeller = as_labeller(mock_names))+
  theme(axis.title.x = element_blank())+
  #theme(legend.text=element_text(size=10), legend.key.size = unit(.2, 'cm'), #change legend ke<- size
  #legend.key.height = unit(.4, 'cm'), #change legend key height
  #legend.key.width = unit(.4, 'cm'))+
  #scale_fill_manual(values = c("#E64B35B2" ,"#4DBBD5B2", "#00A087B2", "#3C5488B2" ,"#7E6148B2",  "#91D1C2B2", "#F39B7FB2" ,"#8491B4B2"
  #,"#DC0000B2" ), labels = c(expression(italic("D. brightwellii")), expression(italic("H. akashiwo")), expression(italic("S. marinoi")), expression(italic("S. menzelii")), expression(italic("T. anguste-lineata")), expression(italic("T. pseudonana")), expression(italic("T. rotula")),expression(italic("T. weissflogii")), "<0.075% of Reads"))+
  #scale_fill_manual(values = c("#8c510a",
  #"#bf812d",
  #"#dfc27d",
  #"#f6e8c3",
  #"#f5f5f5",
  #"#c7eae5",
#"#80cdc1",
#"#35978f",
#"#01665e"), labels = c(expression(italic("D. brightwellii")), expression(italic("H. akashiwo")), expression(italic("S. marinoi")), expression(italic("S. menzelii")), expression(italic("T. anguste-lineata")), expression(italic("T. pseudonana")), expression(italic("T. rotula")),expression(italic("T. weissflogii")), "<0.075% of Reads"))+
#theme(legend.text = element_markdown())+
scale_fill_manual(values = c('#328170',
                             '#8e6e91',
                             '#e89e8e',
                             '#abbef0',
                             '#91b886',
                             '#556894',
                             '#95554b',
                             '#634959',
                             '#8d7a53'), labels = c(expression(italic("D. brightwellii")), expression(italic("H. akashiwo")), expression(italic("S. marinoi")), expression(italic("S. menzelii")), expression(italic("T. anguste-lineata")), expression(italic("T. pseudonana")), expression(italic("T. rotula")),expression(italic("T. weissflogii")), "<0.075% of Reads"))+
  
  #scale_fill_manual(values = c('lightsteelblue3','royalblue3','darkseagreen4', 'darkseagreen1', 'lemonchiffon3', 'burlywood','mistyrose3', 'lightsalmon4','grey20'), labels = c(expression(italic("D. brightwellii")), expression(italic("H. akashiwo")), expression(italic("S. marinoi")), expression(italic("S. menzelii")), expression(italic("T. anguste-lineata")), expression(italic("T. pseudonana")), expression(italic("T. rotula")),expression(italic("T. weissflogii")), "<0.075% of Reads"))+
  #theme(legend.text = element_markdown())+
  # scale_fill_viridis(discrete = "TRUE")+
  labs(x="Sample", y="% Relative abundance")+
  #theme(legend.margin =margin(r=10,l=5,t=1,b=1))+
  y+
  theme(legend.text.align = 0)+
  theme(strip.text.x = element_text(size = 8, family = "Times"))+
  theme(legend.position = "bottom", legend.justification = "center")+
  theme(legend.text = element_text(size = 8))+
  theme(legend.title = element_blank())+
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.margin=margin(3,3,3,3),
        legend.box.margin=margin(-8,-8,-8,-8),
        legend.spacing.x = unit(0.1, "cm"))
mock_plot

ggsave("mock_comp_4.14_nooutline_averaged.pdf", plot = mock_plot, device = "pdf", width =3.5, height = 3, units = c("in"),
       dpi = 1200)
ggsave("mock_comp_4.14_nooutline_averaged.tiff", plot = mock_plot, device = "tiff", width =3.5, height = 3, units = c("in"),
       dpi = 300)

#### END CODE FOR FIGURE 2 #####


### FIGURE 3: TOTAL GENUS FREQUENCY ACROSS TIME SERIES ###

#The "all_diat_df_with_taxa" dataframe comes from above
all_diat_df_with_GenSp <- all_diat_df_with_taxa %>%
  drop_na() %>%
  filter(Sample != "SF81b") %>% #drop these sample replicates
  filter(Sample != "SF82b")

total_freq <- ddply(all_diat_df_with_GenSp , .(tax.Species), summarise, total_count = sum(presence))


total_freq$freq <- ((total_freq$total_count)/80)

total_genus_freq <- ddply(all_diat_df_with_GenSp, .(tax.Genus, Sample), summarise, total_count = sum(presence))


total_genus_freq$tax.Genus[total_genus_freq$tax.Genus == "Hemidiscaceae_X"] <- "Actinocyclus" #Hemidiscaceae is actually Actinocyclus, so fix this

total_genus_freq$total_count <- as.numeric(total_genus_freq$total_count)

total_genus_freq <- total_genus_freq %>%
  filter(total_count >= 1)

total_genus_freq$gen_presence <- 1

count_genus_formatch <- read.csv("genus_count_match.csv")

count_genus_formatch <-  count_genus_formatch %>%
  dplyr::select(Genus, Present_LM) %>%
  distinct()

total_genus_freq_all <- ddply(total_genus_freq, .(tax.Genus), summarise, gen_count = sum(gen_presence))


total_genus_freq_all$freq <- total_genus_freq_all$gen_count/80

total_genus_freq_all <- total_genus_freq_all %>%
  dplyr::rename(Genus = tax.Genus)%>%
  left_join(count_genus_formatch, by = "Genus")

total_genus_freq_all$Present_LM <- factor(total_genus_freq_all$Present_LM , levels = c("Yes", "Not from December 2008 – 2014", "Never"))

total_genus_freq_all[is.na(total_genus_freq_all)] = "No"

total_genus_freq_all$freq <- total_genus_freq_all$freq*100


gen_freq <- ggplot(total_genus_freq_all, aes(x = reorder(Genus, -freq), y = freq, fill = Present_LM))+
  geom_bar(stat = "identity", color = "black", width = 0.8, size = 0.2)+
  y+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10,family = "Times", face = "italic"),
        axis.title.x = element_text(size = 12, family = "Times", face = "bold"),
        axis.title.y = element_text(size = 12, family = "Times", face = "bold"),
        axis.text.y = element_text(size = 10,family = "Times"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 12))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  labs(x = "Genus", y = "Frequency of occurence (%)")+
  scale_fill_manual(values = c("black", "grey70", "white"), name = "Observed presence with LM", labels = c("Yes", "Not from 2008 - 2014 (study period)", "Never"))+
  theme(legend.position = c(0.8, 0.7))+
  scale_y_continuous(expand = c(0,1), limits = c(0,100))
gen_freq

ggsave("Manuscript/Figs/genus_freq.pdf", plot = gen_freq, device = "pdf", width = 7, height = 3.5, units = c("in"),
       dpi = 600)

#### END CODE FOR FIGURE 3 #####

### FIGURE 4: A) ASV RICHNESS AND B) CHLOROPHYLL WITH RICHNESS ###
### FIGURE 4A
###### FIGURE 3 AND STATISTICS FOR RICHNESS ###########
diat_count_genus <- read_csv("Unknown_diatoms_NB_Oct20.csv")
#Do a little filtering
diat_count_genus <- diat_count_genus %>%
  filter(DATE >= as.Date("2008-12-09") & DATE <= as.Date("2014-12-30")) %>%
  filter(COUNT_TYPE == 'S/R Surface') %>%
  dplyr::select(-Centric_unknown, -Unknown_diatoms, -Pennate_unknown, -Diatom_unknown)

#Make Dataframe longer
diat_count_genus_long <- pivot_longer(diat_count_genus, 
                                      cols = 5:152,
                                      names_to = "Genus",
                                      values_to = "count")

diat_count_genus_long$count[is.na(diat_count_genus_long$count)] <-  0

#Separate Genus column
diat_count_species_long <- separate(diat_count_genus_long, Genus, into = c("Genus", "Species"), sep = "_", extra = "drop")

#Determinig which species is present
diat_count_species_long$presence <- sapply(diat_count_species_long$count, function(x) {ifelse(any(x > 0), 1, 0)})

#species richness
diat_count_species_long_rich <- ddply(diat_count_species_long, .(DATE), summarise, count_richness = sum(presence))

diat_count_species_long_rich$DATE <- as_date(diat_count_species_long_rich$DATE)

#The "all_diat_df_with_taxa" dataframe comes from the "Chpater1_dataprocessing.R" file
meta_count <- ddply(all_diat_df_with_taxa, .(Sample), summarise, meta_richness = sum(presence))

meta_count <- meta_count %>%
  left_join(sample_info, by = "Sample") %>%
  filter(Sample != "SF81b") %>%
  filter(Sample != "SF82b")

#Join together metabarcoding genus data with count genus data
meta_count_joined <- meta_count %>%
  left_join(diat_count_species_long_rich, by = c("sample_contents" = "DATE")) 

#calculate averaged and range for LM richness
count_rich <- meta_count_joined_method %>%
  filter(rich_method == "count_richness")

#but first test fr outliers
Boxplot(richness ~ season, data = count_rich)

#remove outliers
count_rich <- count_rich[-c(61,80),]

mean(count_rich$richness, na.rm = TRUE)
min(count_rich$richness, na.rm = TRUE)
max(count_rich$richness, na.rm = TRUE)

kruskal.test(count_rich$richness~ count_rich$season) # yes seasonal difference in LM based richness

Boxplot(richness ~ season, data = count_rich) #count richness (LM) is highest in summer and lowest in spring

#Make longer so that can separate by richness method
meta_count_joined_method <- meta_count_joined %>%
  dplyr::select(sample_contents, season, meta_richness, count_richness) %>%
  pivot_longer(cols = 3:4,
               names_to = "rich_method",
               values_to = "richness"
  )

#Test to see if there is a difference in richness according to method and yes, there is a difference
wilcox.test(richness~rich_method, data = meta_count_joined_method) 


#Arrange sample contents column in order
meta_count_joined_method$sample_contents <- as.character(meta_count_joined_method$sample_contents)

meta_count_joined_method$rich_method <- factor(meta_count_joined_method$rich_method, levels = c("meta_richness", "count_richness"))

#test for outliers test
Boxplot(meta_richness ~ season, data = meta_count_joined)
#remove outliers
meta_count_joined <- meta_count_joined[-11,]

meta_count_joined <- meta_count_joined %>%
  filter(Sample != "SF69b")

meta_count_joined$season <- factor(meta_count_joined$season, levels = c("Winter", "Spring", "Summer", "Fall"))

#test for significance in richness by season and yes there is a difference. p <0.05
kruskal.test(meta_count_joined$meta_richness ~ meta_count_joined$season)

#plot species richness from metabarcoding
species_rich_plot <- ggplot(meta_count_joined, aes(x=as.factor(season), y=meta_richness))+
  #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  geom_boxplot(width = 0.3, outlier.shape = 1,position =  position_dodge(width = 0.1), outlier.size = 2)+
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.2) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.2) +
  stat_summary(fun=mean, geom="point", shape=16, size=2, color="black", fill="black")+
  labs(x = "Season", y = "ASV richness") +
  ylim(0,65)+
  y+
  #theme_classic() + # remove panel background and gridlines
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5))+
  theme(axis.title = element_text(size = 12, face = "bold", family = "Times"))+
  theme(axis.text = element_text(size = 10, family = "Times"),
        legend.text = element_text(size = 16, family = "Times"),
        legend.title = element_text(size = 18, face = "bold", family = "Times"))+
  annotate("text", x=0.6, y=65, label = "a)", family = "Times", fontface = 2, size = 4)+
  annotate(geom = "text", x = "Winter", y = 0, label = "n = 20", fontface = 1,size = 3, family = "Times")+
  annotate(geom = "text", x = "Spring", y = 0, label = "n = 19", fontface = 1,size = 3, family = "Times")+
  annotate(geom = "text", x = "Summer", y = 0, label = "n = 20", fontface = 1,size = 3, family = "Times")+
  annotate(geom = "text", x = "Fall", y = 0, label = "n = 21", fontface = 1,size = 3, family = "Times")
species_rich_plot

ggsave("Manuscript/Figs/species_rich_metaonly.tiff", plot = species_rich_plot, device = "tiff", width = 3.5, height = 3, units = c("in"),
       dpi = 600)


#SUPPLEMENTAL FIGURE 3: UNKNOWN DIATOMS#

#Making a plot to show % unknown diatoms in stacked bar chart form
diat$Chaetoceros_compressus <- as.numeric(diat$Chaetoceros_compressus)
diat$Guinardia_cylindrus <- as.numeric(diat$Guinardia_cylindrus)

diat[is.na(diat)] <- 0

diat_long_spp <- diat %>%
  pivot_longer(cols = 5:155,
               names_to = "Species",
               values_to = "Count") %>%
select(-Unknown_diatoms) %>%
  filter(grepl('spp', Species))


#calculating the number of taxa in the light microscopy dataset
total_diat<- diat %>%
  pivot_longer(cols = 5:155,
               names_to = "Species",
               values_to = "Count") %>%
  select(-Unknown_diatoms)  %>%
  filter(Count >0) %>%
  distinct(Species)



just_spp_count <- ddply(diat_long_spp, .(DATE, Total_abundance, year, month), summarise, total_unk_species = sum(Count))

just_spp_count$percent_unk_species <- (just_spp_count$total_unk_species/just_spp_count$Total_abundance)*100

diat_unk_genus <- read_csv("Unknown_diatoms_NB_genuslevel.csv")

diat_unk_genus <- diat_unk_genus %>%
  select(DATE, `COUNT TYPE`, `Total abundance`, unknown_diatom) %>%
  pivot_longer(cols = 4,
               names_to = "Species",
               values_to = "Count") 


diat_unk_genus <- diat_unk_genus %>%
  filter(DATE >= as.Date("2008-12-09") & DATE <= as.Date("2014-12-30")) %>%
  filter(DATE != "2013-12-30")%>%
  filter(`COUNT TYPE`== 'S/R Surface')

diat_unk_genus$percent_unk_genus <- (diat_unk_genus$Count/diat_unk_genus$`Total abundance`)*100

diat_unk_genus$year <- year(diat_unk_genus$DATE) #Create year from date
diat_unk_genus$year <- as.factor(diat_unk_genus$year) #make year factor

diat_unk_genus$month <- month(diat_unk_genus$DATE)
diat_unk_genus$month <- as.factor(diat_unk_genus$month) #make year factor


diat_unk_genus$year[diat_unk_genus$year == "2008"] <- "2009"

diat_unk_genus <- diat_unk_genus %>%
 select(percent_unk_genus, year, month)

just_spp_count <- just_spp_count %>%
  select(percent_unk_species, year, month)

joined_unk <- just_spp_count %>%
  left_join(diat_unk_genus, by = c("year","month"))

joined_unk_long <- joined_unk %>%
  select(year, month, percent_unk_genus, percent_unk_species) %>%
  pivot_longer(cols = 3:4,
               names_to = "Taxonomic Level",
               values_to = "Percent") 

joined_unk_long$`Taxonomic Level` <- as.factor(joined_unk_long$`Taxonomic Level`)

#figure out the yearly median of cells counted
median_genus <- joined_unk_long %>%
  filter(`Taxonomic Level` == "percent_unk_genus") %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(avg_genus = mean(Percent))

median(median_genus$avg_genus)

median_species <- joined_unk_long %>%
  filter(`Taxonomic Level` == "percent_unk_species") %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(avg_genus = mean(Percent, na.rm = TRUE))

median(median_species$avg_genus)

#plot boxplot
unk_joined_boxplot <- ggplot(joined_unk_long, aes(x=as.factor(year), y=Percent, fill = `Taxonomic Level`)) + 
  geom_boxplot(width = 0.6, outlier.shape = 1, outlier.size = 2)+
  #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = `Taxonomic Level`), width = 0.2) +
  #stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.2) +
  labs(x = "Year", y = "% Unidentified diatom cells per sample") +
  #theme_classic() + # remove panel background and gridlines
  scale_fill_manual(values = c("grey", "white"), labels = c("Genus","Species"))+
  theme(axis.title = element_text(size = 26))+
  theme(axis.text = element_text(size = 22))+
  y+
  theme(element_blank())+
  theme(legend.position = "none")+
  theme_classic() + # remove panel background and gridlines
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5))+
  theme(axis.title = element_text(size = 14, face = "bold", family = "Times"))+
  theme(axis.text = element_text(size = 12, family = "Times"),
        legend.text = element_text(size = 12, family = "Times"),
        legend.title = element_text(size = 14, face = "bold", family = "Times"))+
  annotate(geom = "text", x = "2009", y = -3, label = "n = 53", fontface = 1,size = 4, family = "Times")+
  annotate(geom = "text", x = "2010", y = -3, label = "n = 51", fontface = 1,size = 4, family = "Times")+
  annotate(geom = "text", x = "2011", y = -3, label = "n = 48", fontface = 1,size = 4, family = "Times")+
  annotate(geom = "text", x = "2012", y = -3, label = "n = 7", fontface = 1,size = 4, family = "Times")+
  annotate(geom = "text", x = "2013", y = -3, label = "n = 52", fontface = 1,size = 4, family = "Times")+
  annotate(geom = "text", x = "2014", y = -3, label = "n = 52", fontface = 1,size = 4, family = "Times")
unk_joined_boxplot 

ggsave("Manuscript/Supplemental_Info/Figures/unk_diatom_boxplot_year_joined.png", plot = unk_joined_boxplot , device = "png", width = 15, height = 12, units = c("cm"),
       dpi = 400)


### END OF FIGURE 4A PLOT ###

# Calculate the average richness per season
fall_rich <- meta_count_joined %>%
  filter(season == "Fall") %>%
  summarise(mean = mean(meta_richness))

spring_rich <- meta_count_joined %>%
  filter(season == "Spring") %>%
  summarise(mean = mean(meta_richness))

### FIGURE 4B: CHLOROPHYLL AND RICHNESS
# read in chl data
chl_nb <- read_csv("chldata_NB.csv")
chl_nb <- chl_nb %>%
  filter(Sample_Date >= "2008-12-09" & Sample_Date <= "2014-12-30")%>%
  dplyr::select(Sample_Date, surf_chla_corrected)

chl_nb$Sample_Date <- as_date(chl_nb$Sample_Date)
chl_nb$surf_chla_corrected <- as.numeric(chl_nb$surf_chla_corrected)

#join chl data with richness data
meta_count_chl <- meta_count %>%
  left_join(chl_nb, by = c("sample_contents" = "Sample_Date"))

#meta_count_chl$`surface chla all` <- as.numeric(meta_count_chl$`surface chla all`)

meta_count_chl$season <- factor(meta_count_chl$season, levels = c("Winter", "Spring", "Summer", "Fall"))

#figure out which model is best

#calculate quadratic model
meta_count_chl$chl2 <- meta_count_chl$surf_chla_corrected^2

meta_count_chl <- meta_count_chl %>%
  filter(!is.na(surf_chla_corrected)) %>%
  filter(Sample != "SF1b") %>% #remove these samples because they're outliers as determined with boxplot
  filter(Sample != "SF69b")


quadraticModel <- lm(meta_richness~ surf_chla_corrected + chl2, data=meta_count_chl)
summary(quadraticModel)

linearmodel <- lm(meta_richness~ surf_chla_corrected, data=meta_count_chl)
summary(linearmodel)

logmodel <- lm(meta_richness~ log10(surf_chla_corrected), data=meta_count_chl)
summary(logmodel)


#Actual plot of richness and chlorophyll 
chl_richness_all <- meta_count_chl %>%
  ggplot(aes(y = meta_richness, x = surf_chla_corrected))+
  geom_point(size = 2)+
  labs(x = "Chlorophyll (µg/L)", y = "ASV richness")+
  stat_smooth(method='lm', formula = y~poly(x,2), se = FALSE, color = "black")+
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "black", fill = "NA", size = 0.5))+
  scale_fill_manual(name = "Season", values = c("dodgerblue4", "darkseagreen","darkred", "darkorange"))+
  y+
  theme(axis.title = element_text(size = 12, face = "bold", family = "Times"))+
  theme(axis.text = element_text(size = 10, family = "Times"))+
  ylim(c(0,65))+
  theme(legend.position = c(0.87, 0.8))+
  # annotate(geom = "text", x = 25, y = 60, label = "y = -0.08x^2 + 1.6x + 32.4", fontface = 3,size = 3.5)+
  annotate("text", x=19.6, y=62, label = as.character(expression("y == -0.07*x^{2}")), parse = T, family = "Times", size = 3)+
  annotate("text", x=24, y=55, label = as.character(expression("R^{2} == 0.13")), parse = T, family = "Times", size = 3)+
  annotate("text", x=0, y=65, label = "b)", family = "Times", fontface = 2, size = 4)+
  annotate("text", x=26.5, y=62, label = "+ 1.6x + 32.4", family = "Times", size = 3)
chl_richness_all

ggsave("Manuscript/Figs/species_rich_chl.tiff", plot =chl_richness_all , device = "tiff", width = 3.5, height = 3, units = c("in"),
       dpi = 600)

### END OF CODE FOR FIGURE 4 ###

### FIGURE 5: TB-RDA
## Reading in environmental data for NBay time series for RDA with all ASVs
env <- read_csv("Envdata.csv")

#read in sample info
sample_info <- read_csv("sample_info_field.csv")

sample_info$sample_contents <-  as_date(sample_info$sample_contents)
sample_info <- sample_info %>% 
  filter(!is.na(sample_contents))
sample_info$month <- lubridate::month(sample_info$sample_contents)
#add season
sample_info$season <- with(sample_info,
                           ifelse(month %in% c(12,01,02), "Winter",
                                  ifelse(month %in% c(03,04,05), "Spring",
                                         ifelse(month %in% c(06,07,08), "Summer",
                                                ifelse(month %in% c(09,10, 11), "Fall", NA)))))

env_season <- env %>%
  left_join(sample_info, by = c("Date"= "sample_contents"))%>%
  dplyr::select(Date, season) 


env_season$season <- factor(env_season$season, levels = c("Winter", "Spring", "Summer", "Fall"))


env_season <- column_to_rownames(env_season, var = "Date")
season <- as.factor(env_season$season)


filt_combined_wide_all <- all_diat_df_with_taxa %>%
  filter(Sample != "SF81b")%>%
  filter(Sample != "SF82b") %>%
  dplyr::select(Sample, ASV, presence) %>%
  pivot_wider(id_cols = 1,
              names_from = ASV,
              values_from = presence)

filt_combined_wide_all <- filt_combined_wide_all%>%
  left_join(sample_info ,by = "Sample")

#function for normalizing env data
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#Normalize env data
env_norm <- as.data.frame(lapply(env[2:8], min_max_norm))
#Add date column back in 
env_norm$Date <- env$Date

#Join together matrix of community info with env info
filt_combined_envi_all <- filt_combined_wide_all %>%
  left_join(env_norm, by = c("sample_contents" = "Date")) #make sure date is in the format YYYY-MM-DD from the csv file

filt_combined_envi_all <- column_to_rownames(filt_combined_envi_all, var = "sample_contents")

species_all <- filt_combined_envi_all[,2:658] #Species go from 2 to 658

#make dataframe
species_all <- as.data.frame(species_all)
#change all "NA" to 0
species_all[is.na(species_all)] <-  0

species_counts <- as.data.frame(colSums(species_all))
hist(species_counts$`colSums(species_all)`)
species_counts <- species_counts %>%
  filter(`colSums(species_all)` ==1)

#transform species data with hellinger transformation
species_transformed_all <- decostand(species_all, method = "hellinger") #transform 
#get envi data
environment_all <- filt_combined_envi_all[,661:667]

#run the rda
dbRDA_all <- rda(species_transformed_all ~ Temp  + DIN_μM + DIP_μM + Si_μM +Salinity, environment_all, dist = "jaccard")

#look at screeplot
screeplot(dbRDA_all)
#Get R2 values
dbRDA_R2_all <- RsquareAdj(dbRDA_all)$r.squared

#get base plot for it
plot(dbRDA_all)

#Overall test of the significance of the analysis
anova(dbRDA_all, permutations = 999) 

#test axes for sig
anova(dbRDA_all, by = "axis", perm.max = 5000) 
#test for significance of env. variables
anova(dbRDA_all, by = "terms", permu = 5000) 

scores(dbRDA_all) #Getting the scores
scores_dbRDA_all <- scores(dbRDA_all)


site_scores_all <- data.frame(scores_dbRDA_all$sites)
site_scores_all <- rownames_to_column(site_scores_all, var = "Date")
species_scores_all <- scores_dbRDA_all$species #separating out species scores

scores_dbRDA_matrix_all <- as.data.frame(scores_dbRDA_all$species)

spec_names_all <- row.names(scores_dbRDA_matrix_all)

rda_with_labels <- orditorp(dbRDA_all, display = "species", labels = spec_names_all, choices = c(1,2), pch = 2, col = "red")
rda_with_labels


site_scores_all$Date <- as.Date(site_scores_all$Date)

env_season <- rownames_to_column(env_season, var = "Date")
env_season$Date <- as.Date(env_season$Date)

site_scores_all <- site_scores_all %>%
  left_join(env_season, by = "Date")

site_scores_all <- column_to_rownames(site_scores_all, var = "Date")

species_scores_forname_all <- as.data.frame(scores_dbRDA_all$species)

species_scores_forname_all$RDA1_abs <- abs(species_scores_forname_all$RDA1)
species_scores_forname_all$RDA2_abs <- abs(species_scores_forname_all$RDA2)

species_scores_forname_all <- species_scores_forname_all %>%
  mutate(test = RDA1_abs + RDA2_abs)

species_scores_select_all <- species_scores_forname_all %>%
  filter(test >0.14) 
species_scores_select_all <- species_scores_select_all[-16,]

#save species scores file
write.csv(species_scores_select_all, "selected_species_scores_withallASVs.csv")

#run the envi fit
env_forfit <- environment_all %>%
  dplyr::select(Temp, DIN_μM, DIP_μM, Si_μM, Salinity)


envfit_all <- envfit(dbRDA_all, env_forfit, permutations = 999)


#get envi scores
env_scores_all <- data.frame((envfit_all$vectors)$arrows, (envfit_all$vectors)$r, (envfit_all$vectors)$pvals)

env_scores_all <- rownames_to_column(env_scores_all, var = "Variable")

#change envi variable names
env_scores_all[1,1] <- "TEMP"
env_scores_all[2,1] <- "DIN"
env_scores_all[3,1] <- "DIP"
env_scores_all[4,1] <- "SI"
env_scores_all[5,1] <- "SAL"


#save envi scores
real_env_scores <- read_csv("real_envi_scores.csv")

env_scores_plot_all <- env_scores_all %>%
  dplyr::select(Variable, RDA1, RDA2)


# figuring out what species are significantly correlated in dbRDA
#envfit() takes the output of metaMDS() and the species matrix you created
fit <- envfit(dbRDA_all, species_all, perm = 999) 

# extract p-values for each species
fit_pvals <- fit$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")

# extract coordinates for species, only keep species with p-val = 0.001
fit_spp <- fit %>% #111 are significant 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fit_pvals, by = "species") %>% 
  filter(pvals < 0.01)


#Determine the number of clusters
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(species_transformed_all, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#Do the kmeans analsyis 
km_res_rda <- kmeans(species_transformed_all, 2, nstart = 1000)
print(km_res_rda)

#add cluster info back to DF
site_scores_cluster <- cbind(site_scores_all, cluster = km_res_rda$cluster)
fviz_cluster(km_res_rda, data = species_transformed_all ,
             #palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07", "green"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot"
)

fviz_nbclust(species_transformed_all, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

#rda with the segments
rda.plot_segment <- ggplot()+
  geom_segment(data = species_scores_select_all, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = env_scores_plot_all, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = -0.9, 0.22, label = "TEMP", fontface =2, family = "Times")+
  annotate(geom = "text", x = -0.88,y = -.55, label = "SI", fontface =2, family = "Times")+
  annotate(geom = "text", x = -0.6, y = -0.86, label = "DIP", fontface =2, family = "Times")+
  annotate(geom = "text", x = -0.01, y = -1.07, label = "DIN", fontface =2, family = "Times")+
  annotate(geom = "text", x = -0.25, y = -1.02, label = "SAL", fontface =2, family = "Times")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_rect(aes(xmin = -0.25, xmax = 0.3, ymin = -0.25, ymax = 0.25),
            fill = "transparent", color = "grey", size = 0.6) +
  annotate("Text", x = 0.15, y = -0.18, label = "Species", color="black", family="Times", fontface = 2, size=3)
rda.plot_segment

site_scores_cluster$cluster[site_scores_cluster$cluster== "1"] <- "Summer-Fall"
site_scores_cluster$cluster[site_scores_cluster$cluster== "2"] <- "Winter-Spring"

#the complete plot
complete_rda<- rda.plot_segment +
  geom_point(data = site_scores_all, aes(x=RDA1, y=RDA2, fill = season), pch =21, color = "black", size = 3)+
  #geom_hline(yintercept=0, linetype="dotted") +
  stat_ellipse(data=site_scores_cluster, aes(x =RDA1, y = RDA2, fill=cluster),
               geom = "polygon", alpha=0.2, level=0.95)+
  y+
  #geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [64%]", y = "tb-RDA2 [16%]")+
  scale_fill_manual(name = "Season", values = c("darkorange","darkseagreen" ,"darkred", "grey80", "dodgerblue4" ,"gray28"))+
  # scale_fill_manual(name = "Season", values = c("dodgerblue4","darkseagreen","darkred" ,"darkorange" ))+
  #scale_fill_manual(values = c("dodgerblue4","darkseagreen", "darkred" ,"darkorange"))+
  #scale_color_manual(values = c("gray27", "gray80"))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 14, family = "Times", face = "bold"))+
  theme(axis.text = element_text(size = 12, family = "Times"))+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")

# annotate(geom = "text", x = -1, y = .98, label = "a)", fontface = 2,size = 5)
complete_rda


ggsave("Manuscript/Figs/dbRDA_plot_allASVs_no_species_4.14.tiff", plot = complete_rda, device = "tiff", width = 3.5, height = 3.2, units = c("in"),
       dpi = 600)

#getting a no ckuster plot to obtain the legend
nocluster_plot <- rda.plot_segment+
  geom_point(data = site_scores_all, aes(x=RDA1, y=RDA2, fill = season), pch =21, color = "black", size = 3)+
  #geom_hline(yintercept=0, linetype="dotted") +
  y+
  #geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [64%]", y = "tb-RDA2 [16%]")+
  # scale_fill_manual(name = "Season", values = c("darkorange","darkseagreen" ,"darkred", "grey80", "dodgerblue4" ,"gray28"))+
  scale_fill_manual(name = "Season", values = c("dodgerblue4","darkseagreen","darkred" ,"darkorange" ))+
  #scale_fill_manual(values = c("dodgerblue4","darkseagreen", "darkred" ,"darkorange"))+
  #scale_color_manual(values = c("gray27", "gray80"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 10, family = "Times"))+
  theme(legend.title = element_blank())+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")

# annotate(geom = "text", x = -1, y = .98, label = "a)", fontface = 2,size = 5)
nocluster_plot

ggsave("Manuscript/Figs/noclusterplot_legend.tiff", plot = nocluster_plot, device = "tiff", width = 3.5, height = 3.2, units = c("in"),
       dpi = 600)

#### END OF CODE FOR FIGURE 5 ###

### SUPPLEMENTAL FIGURE 4: SPECIES VECTORS FOR RDA ###
#With species
species_rda <- ggplot()+
  geom_segment(data = species_scores_select_all, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  #geom_hline(yintercept=0, linetype="dotted") +
  y+
  #geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [64%]", y = "tb-RDA2 [16%]")+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1))+
  #stat_ellipse(data=site_scores_cluster, aes(x =RDA1, y = RDA2, fill=cluster),
  # geom = "polygon", alpha=0.2, level=0.95)+
  #scale_fill_manual(name = "Season", values = c("darkorange","darkseagreen" ,"darkred", "grey80", "dodgerblue4" ,"gray25"))
  #scale_fill_nejm()+
  annotate(geom = "text", x = 0.154, y = 0.037, label = "T. cur1", fontface = 3, size = 3.5)+
  annotate(geom = "text", x = -0.1, y = -0.07, label = "T. ten1", fontface = 3,size = 3.5)+  
  annotate(geom = "text", x = .08, y = -0.08, label = "T. a-l2",fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.208, y = 0.03, label = "T. a/n1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -.157, y = -0.055, label = "R. del1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.197, y = 0.057, label = "T. gui1", fontface = 3,size = 3.5)+
  #annotate(geom = "text", x = -0.8, y = -.02, label = "S. gp2", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.208, y = 0.006, label = "T. a/p.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -.15, y = -.078, label = "M. com1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.17, y = -0.04, label = "Cy. sp1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.11, y = -.146, label = "T. pun1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.22, y = -0.094, label = "T. min1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -.2, y = 0.07, label = "Ch. ten1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.05, y = -0.14, label = "T. e/m.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.07 , y = 0.093, label = "M. spi.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.132 , y = -0.07, label = "G. del1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.154 , y = -.05, label = "C. hys1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.175 , y = .032, label = "S. g/m/j2", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.215 , y = -.005, label = "S. g/m/j1.", fontface = 3,size = 3.5)+
  #annotate(geom = "text", x = 0.15 , y = -0.06, label = "T. ang.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.148 , y = -.025, label = "L. min1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.215, y = .022, label = "S. p/t1", fontface = 3,size = 3.5)+
  xlim(-0.23, 0.23)+
  annotate(geom = "text", x = -0.11, y = .068, label = "M. com2", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.11, y = .045, label = "T. ang.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.15, y = -.1, label = "ASV_173*", fontface = 1,size = 3.5)+
  annotate(geom = "text", x = 0.18, y = -.016, label = "R. d/s.", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.158, y = .02, label = "T. a-l1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.16, y = -.003, label = "Ch. sp1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.146, y = -.032, label = "L. min2", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.05, y = -.07, label = "Ch. dec1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = -0.16, y = -.01, label = "C. pel1", fontface = 3,size = 3.5)+
  annotate(geom = "text", x = 0.12, y = .01, label = "Cy. str1", fontface = 3,size = 3.5)+
  theme(axis.title = element_text(size = 14, family = "Times", face = "bold"))+
  theme(axis.text = element_text(size = 12, family = "Times"))
species_rda

ggsave("Manuscript/Figs/dbRDA_species_inset.png", plot = species_rda, device = "png", width = 15, height = 12, units = c("cm"),
       dpi = 400)

### END OF CODE FOR SUPPLEMENTAL 4 ###

### FIGURE 6: TIME SERIES OF RDA VALUES ### 
#Making a graph of CAP1 scores and date using site_scores DF
site_scores_date_all <- rownames_to_column(site_scores_all, var = "Date")
site_scores_date_all <- site_scores_date_all[order(site_scores_date_all$Date),] #order by date
site_scores_date_all$date_value <- seq(1,80)
site_scores_date_all$Date <- as_date(site_scores_date_all$Date)
site_scores_date_all <- site_scores_date_all %>%
  dplyr::select(-season) %>%
  left_join(sample_info, by = c("Date" = "sample_contents")) %>%
  dplyr::select(-Sample, -month)


#Need to tell ggplot how to shade
shade_winter <-  data.frame(x1=c("2008-12-01","2009-12-01", "2010-12-01", "2011-12-01", "2012-12-01", "2013-12-01", "2014-12-01"), x2=c("2009-02-28","2010-02-28","2011-02-28","2012-02-29","2013-02-28", "2014-02-28","2015-02-28"), y1=c(-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6), y2=c(0.6,0.6,0.6,0.6,0.6,0.6,0.6))
shade_spring <-  data.frame(x1=c("2009-02-28", "2010-02-28", "2011-02-28","2012-02-29","2013-02-28", "2014-02-28"), x2=c("2009-05-31","2010-05-31","2011-05-31","2012-05-31","2013-05-31","2014-05-31"), y1=c(-0.6,-0.6,-0.6,-0.6,-0.6, -0.6), y2=c(0.6,0.6,0.6,0.6,0.6,0.6))
shade_summer <-  data.frame(x1=c("2009-05-31","2010-05-31","2011-05-31","2012-05-31","2013-05-31","2014-05-31"), x2=c("2009-08-31","2010-08-31","2011-08-31","2012-08-31","2013-08-31","2014-08-31"), y1=c(-0.6,-0.6,-0.6,-0.6,-0.6, -0.6), y2=c(0.6,0.6,0.6,0.6,0.6,0.6))
shade_fall <-  data.frame(x1=c("2009-08-31","2010-08-31","2011-08-31","2012-08-31","2013-08-31","2014-08-31"), x2=c("2009-12-01", "2010-12-01", "2011-12-01", "2012-12-01", "2013-12-01", "2014-12-01"), y1=c(-0.6,-0.6,-0.6,-0.6,-0.6, -0.6), y2=c(0.6,0.6,0.6,0.6,0.6,0.6))

#for RDA 2 values (DIN plot)
shade_winter_cap2 <-  data.frame(x1=c("2008-12-01","2009-12-01", "2010-12-01", "2011-12-01", "2012-12-01", "2013-12-01", "2014-12-01"), x2=c("2009-02-28","2010-02-28","2011-02-28","2012-02-29","2013-02-28", "2014-02-28","2015-02-28"), y1=c(-1,-1,-1,-1,-1,-1,-1), y2=c(1,1,1,1,1,1,1))
shade_spring_cap2 <-   data.frame(x1=c("2009-02-28", "2010-02-28", "2011-02-28","2012-02-29","2013-02-28", "2014-02-28"), x2=c("2009-05-31","2010-05-31","2011-05-31","2012-05-31","2013-05-31","2014-05-31"), y1=c(-1,-1,-1,-1,-1, -1), y2=c(1,1,1,1,1,1))
shade_summer_cap2 <-  data.frame(x1=c("2009-05-31","2010-05-31","2011-05-31","2012-05-31","2013-05-31","2014-05-31"), x2=c("2009-08-31","2010-08-31","2011-08-31","2012-08-31","2013-08-31","2014-08-31"), y1=c(-1,-1,-1,-1,-1, -1), y2=c(1,1,1,1,1,1))
shade_fall_cap2 <-  data.frame(x1=c("2009-08-31","2010-08-31","2011-08-31","2012-08-31","2013-08-31","2014-08-31"), x2=c("2009-12-01", "2010-12-01", "2011-12-01", "2012-12-01", "2013-12-01", "2014-12-01"), y1=c(-1,-1,-1,-1,-1, -1), y2=c(1,1,1,1,1,1))


shade_winter$season <- "Winter"
shade_spring$season <- "Spring"
shade_summer$season <- "Summer"
shade_fall$season <- "Fall"

shade_winter_cap2$season <- "Winter"
shade_spring_cap2$season <- "Spring"
shade_summer_cap2$season <- "Summer"
shade_fall_cap2$season <- "Fall"

shade <- rbind(shade_winter, shade_spring, shade_summer, shade_fall)
shade$season <- factor(shade$season, levels = c("Winter", "Spring", "Summer", "Fall"))


shade_cap2 <- rbind(shade_winter_cap2, shade_spring_cap2, shade_summer_cap2, shade_fall_cap2)
shade_cap2$season <- factor(shade_cap2$season, levels = c("Winter", "Spring", "Summer", "Fall"))



#Adding temp data to site_scores date DF
site_scores_date_all <- site_scores_date_all %>%
  left_join(env, by = "Date") %>%
  dplyr::select(date_value, Date, RDA1, RDA2, DIN_μM, Temp, season)

#write_csv(site_scores_date, "site_scores_date.csv")

shade$x1 <- as_date(shade$x1 )
shade$x2 <- as_date(shade$x2)

shade_cap2$x1 <-  as_date(shade_cap2$x1 )
shade_cap2$x2 <-  as_date(shade_cap2$x2 )


RDA_timeseries <- ggplot(site_scores_date_all, aes(Date,RDA1))+
  geom_rect(data = shade,aes(x = NULL,y = NULL,xmin = x1,xmax = x2,ymin = y1,ymax = y2,fill = season),
            alpha = 0.2)+
  geom_line(size = 0.7 ,color='black')+
  geom_point(aes(x=Date, y=RDA1), size = 1, data= site_scores_date_all)+
  #theme_black()+
  y+
  scale_fill_manual(name = "Season", 
                    values = c('dodgerblue4',"darkseagreen","darkred","darkorange"))+
  labs(y = "tb-RDA1 [64%]", x = "Date")+
  #scale_x_continuous(expand = c(0,0),breaks = seq(1,80, by = 4))+
  scale_y_continuous(expand = c(0,0), breaks = seq(-.5,.5, by = 0.5))+
  scale_x_date(expand = expansion(0), breaks = date_breaks("1 year"), labels = date_format("%Y"))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text = element_text(size = 8, family = "Times"),
        axis.title.y = element_text(size = 10, family = "Times", face = "bold"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")
  
RDA_timeseries

ggsave("Manuscript/Figs/RDA_timeseries_allASVs_4.14.tiff", plot =   RDA_timeseries, device = "tiff", width = 3.6, height = 2, units = c("in"),
       dpi = 600)

#RDA time series with CAP2 values
RDA_timeseries_CAP2 <- ggplot(site_scores_date_all, aes(Date,RDA2))+
  geom_line(size = 0.7, color='black')+
  geom_point(aes(x=Date, y=RDA2), size = 1, data= site_scores_date_all)+
  geom_rect(data = shade_cap2,aes(x = NULL,y = NULL,xmin = x1,xmax = x2,ymin = y1,ymax = y2,fill = season),
            alpha = 0.2)+
  y+
  scale_fill_manual(name = "Season", 
                    values = c('dodgerblue4',"darkseagreen","darkred","darkorange"))+
  labs(y = "tb-RDA2 [16%]", x = "Time Series Day")+
  scale_y_continuous(expand = c(0,0), breaks = seq(-1,1, by = 0.5))+
  scale_x_date(expand = expansion(0), breaks = date_breaks("1 year"), labels = date_format("%Y"))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text = element_text(size = 8, family = "Times"),
        axis.title.y = element_text(size = 10, family = "Times", face = "bold"),
        axis.title.x = element_blank())+
  theme(legend.position = "none")
RDA_timeseries_CAP2

ggsave("Manuscript/Figs/RDA_timeseries_allASVs_CAP2.tiff", plot =   RDA_timeseries_CAP2, device = "tiff", width = 3.5, height = 2, units = c("in"),
       dpi = 600)

timeseries_merged <- ggarrange(RDA_timeseries, RDA_timeseries_CAP2, ncol = 1, nrow = 2, common.legend = FALSE)

ggsave("Manuscript/Figs/timeseries_merged.tiff", plot =  timeseries_merged, device = "tiff", width = 3.5, height = 4, units = c("in"),
       dpi = 600)



#then need to graph temp and DIN values separately
RDA_Temp <-  ggplot(site_scores_date_all, aes(Date,Temp))+
  geom_line(size = 0.7, color='black', linetype = 6)+
  #geom_point(aes(x=date_value, y=Temp), size = 2.5, data= site_scores_date)+
  #theme_black()+
  y+
  scale_fill_manual(name = "Season", 
                    values = c('dodgerblue4',"darkseagreen","darkred","darkorange"))+
  labs(y = "Temperature (˚C)", x = "Time Series Day")+
  scale_x_continuous(expand = c(0,0),breaks = seq(1,79, by = 4))+
  scale_y_continuous(position = "right")+
  theme(panel.border = element_rect(fill=NA, colour = "transparent", size=1),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 10, family = "Times", face = "bold"))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank())
RDA_Temp


ggsave("Manuscript/Figs/RDA_Temp.tiff", plot =   RDA_Temp, bg = "transparent", device = "tiff", width = 3.27, height = 1.74, units = c("in"),
       dpi = 600)

RDA_DIN <-  ggplot(site_scores_date_all, aes(Date,DIN_μM))+
  geom_line(size = 0.7, color='black', linetype = 6)+
  # geom_point(aes(x=date_value, y=Temp), size = 2.5, data= site_scores_date)+
  y+
  scale_fill_manual(name = "Season", 
                    values = c('dodgerblue4',"darkseagreen","darkred","darkorange"))+
  labs(y = "DIN (µM)", x = "Time Series Day")+
  scale_x_continuous(expand = c(0,0),breaks = seq(1,80, by = 4))+
  scale_y_continuous(position = "right")+
  theme(panel.border = element_rect(fill=NA, colour = "transparent", size=1),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 10, family = "Times", face = "bold"))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank())
RDA_DIN

ggsave("Manuscript/Figs/RDA_DIN.tiff", plot =   RDA_DIN, bg = "transparent", device = "tiff", width = 3.25, height = 1.73, units = c("in"),
       dpi = 600)

#then match up these plots in powerpoint and export as PDF

#correlations for temp and DIN with RDA1 and RDA2, respectively 
cor.test(site_scores_date_all$Temp, site_scores_date_all$RDA1, method = "pearson")
cor.test(site_scores_date_all$DIN_μM, site_scores_date_all$RDA2, method = "pearson")

### END OF CODE FOR FIGURE 6 ###

### FIGURE 7: GENUS FREQUENCY PER MONTH ###
temp_genus_freq <- all_diat_df_with_taxa %>%
  left_join(sample_info, by = "Sample") %>%
  filter(Sample!= "SF81b") %>%
  filter(Sample != "SF82b") %>%
  filter(!is.na(tax.Genus)) %>%
  dplyr::select(tax.Genus, sample_contents, presence, month, season) %>%
  dplyr::group_by(month, tax.Genus, sample_contents) %>%
  dplyr::summarise(sum_genus = sum(presence)) 

temp_genus_freq$sum_genus <- 1


temp_genus_freq$presence <- 1

month_freq_genus <- temp_genus_freq %>%
  dplyr::group_by(tax.Genus, month) %>%
  dplyr::summarise(sum_presence = sum(presence)) 


month_freq_genus_filtered <- month_freq_genus %>%
  dplyr::group_by(tax.Genus) %>%
  filter(any(sum_presence>=3))

month_freq_genus_filtered <- month_freq_genus_filtered %>%
  dplyr::group_by(month, tax.Genus) %>%
  left_join(sampling_effort_month, by = "month")


month_freq_genus_filtered$month_freq <- month_freq_genus_filtered$sum_presence/month_freq_genus_filtered$month_total*100

month_freq_genus_filtered$month_abb <- mymonths[month_freq_genus_filtered$month ]

month_freq_genus_filtered$month_abb <- factor(month_freq_genus_filtered$month_abb, levels = c("Jan","Feb","Mar", "Apr","May","Jun", "Jul","Aug","Sep", "Oct","Nov","Dec"))


heatmap_genus <- ggplot(data = month_freq_genus_filtered, aes(x = as.factor(month_abb), y = tax.Genus)) +
  geom_tile(aes(fill = month_freq), color = 'black') +
  scale_fill_viridis(limits = c(0,100), name = "Monthly Frequency \n of Occurrence")+
  #scale_fill_gradient(low="grey", high = "steelblue", name = "Presence", limits = c(0,1))+
  y+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 1, family = "Times"))+
  theme(axis.text.y = element_text(size = 8, face = "italic", family = "Times"),
        axis.title = element_text(size = 10, face = "bold", family = "Times"),
        legend.position = "top",
        legend.text = element_text(size = 7, family = "Times"),
        legend.title = element_text(size = 8, family = "Times"))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  labs(x = "Month", y = "Genus")+
  theme(legend.position = "top")+
  theme(legend.title.align = 0.5)
heatmap_genus

ggsave("genus_month_freq.pdf", plot = heatmap_genus , device = "pdf", width = 3.5, height = 5, units = c("in"),
       dpi = 300)

#for clustering dendrogram
temp_genus_freq_wide <- month_freq_genus_filtered%>%
  dplyr::select(tax.Genus,month,month_freq) %>%
  pivot_wider(names_from = month,
              values_from = month_freq)



# Run clustering
genus_heatmap_matrix <- as.matrix(temp_genus_freq_wide[, -1])
rownames(genus_heatmap_matrix) <- temp_genus_freq_wide$tax.Genus

genus_heatmap_matrix[is.na(genus_heatmap_matrix)] <- 0
genus_heatmap_dendro <- as.dendrogram(hclust(d = dist(x = genus_heatmap_matrix)))

#Plot dendro
#create dendro
dendro.plot <- ggdendrogram(data = genus_heatmap_dendro, rotate = TRUE)

# Preview the plot
print(dendro.plot)
dendro.plot <- dendro.plot + theme(axis.text.y = element_blank(), axis.text.x = element_blank())

gp1 <- ggplotGrob(dendro.plot)

#Order dendro
heatmap.order <- order.dendrogram(genus_heatmap_dendro)


# Order the levels according to their position in the cluster
month_freq_genus_filtered$tax.Genus <- factor(x = month_freq_genus_filtered$tax.Genus, levels = temp_genus_freq_wide$tax.Genus[heatmap.order],ordered = TRUE)

#new order manually
month_freq_genus_filtered$tax.Genus <- factor(month_freq_genus_filtered$tax.Genus, levels = c("Pseudo-nitzschia", "Corethron","Pleurosigma","Detonula","Odontella","Attheya","Licmophora", "Nitzschia",  "Tabularia",  "Brockmanniella", "Navicula", "Coscinodiscus", "Eucampia",  "Actinoptychus", "Actinocyclus", "Talaroneis", "Thalassionema", "Cylindrotheca", "Minutocellus", "Cerataulina", "Rhizosolenia", "Leptocylindrus", "Ditylum", "Guinardia","Skeletonema", "Chaetoceros", "Thalassiosira", "Minidiscus","Cyclotella" ))





#heatmap
heatmap.plot <- ggplot(data = month_freq_genus_filtered, aes(x = as.factor(month_abb), y = tax.Genus)) +
  geom_tile(aes(fill = month_freq), color = 'black') +
  scale_fill_viridis(limits = c(0,100), name = "Monthly Frequency \n of Occurrence")+
  #scale_fill_gradient(low="grey", high = "steelblue", name = "Presence", limits = c(0,1))+
  y+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 1, family = "Times"))+
  theme(axis.text.y = element_text(size = 8, face = "italic", family = "Times"),
        axis.title = element_text(size = 10, face = "bold", family = "Times"),
        legend.position = "top",
        legend.text = element_text(size = 7, family = "Times"),
        legend.title = element_text(size = 8, family = "Times"))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  labs(x = "Month", y = "Genus")+
  theme(legend.position = "top")+
  theme(legend.title.align = 0.5)
heatmap.plot

print(heatmap.plot)

gp2 <- ggplotGrob(heatmap.plot)

grid.arrange(gp2, gp1, ncol = 2)

#so this is the correct order for the genera, but need to manually switch the nodes in affinity designer

### END OF CODE FOR FIGURE 7 ### 

### FIGURE 8: TEMPORAL PATTERNS OF ASVS ###
########Niche exploration _time plot figure 6
minidiscus <- all_diat_df_with_taxa %>%
  filter(tax.Genus %in% c("Minidiscus", "Rhizosolenia", "Cyclotella", "Chaetoceros", "Leptocylindrus"))

#need to figure out the temperatures in which chaet sp1 and chaet tenuiss occur
temp_check <- minidiscus %>%
  filter(ASV %in% c("ASV_4", "ASV_16"))

temp_check <- temp_check %>%
  left_join(sample_info, by = "Sample") %>%
  filter(Sample!= "SF81b") %>%
  filter(Sample != "SF82b") %>%
  select(ASV, sample_contents, Sample, tax.Species, month, presence)

temp_check <- temp_check %>%
  left_join(env, by = c('sample_contents' = 'Date' ))

minidiscus <- minidiscus %>%
  left_join(sample_info, by = "Sample") %>%
  filter(Sample!= "SF81b") %>% #exclude field reps
  filter(Sample != "SF82b")

minidiscus_for_removal <- minidiscus %>%
  dplyr::select(ASV, Sample, tax.Species, month, per_rel)

minidiscus_for_removal$presence <- 1

#want to figure out what months species are in
minidiscus_for_removal <- minidiscus_for_removal %>%
  dplyr::select(ASV, tax.Species, month, presence) %>%
  unique()

minidiscus_for_removal <- ddply(minidiscus_for_removal, .(ASV), summarise, total_count = sum(presence))

species_removal <- minidiscus_for_removal %>%
  filter(total_count < 4)

minidiscus <- minidiscus %>%
  anti_join(species_removal, by = "ASV")


minidiscus_calcs <- ddply(minidiscus, .(ASV,tax.Genus, tax.Species, month), summarise, total_count = sum(presence))

month_totals <- ddply(minidiscus_calcs, .(month, tax.Genus), summarise, total_sp = sum(total_count))

minidiscus_calcs <- minidiscus_calcs %>%
  left_join(month_totals, by = c("month", "tax.Genus"))

minidiscus_calcs$asv_freq <- (minidiscus_calcs$total_count/minidiscus_calcs$total_sp)*100

minidiscus_calcs$month_abb <- mymonths[minidiscus_calcs$month]

minidiscus_calcs$month_abb <- factor(minidiscus_calcs$month_abb, levels = c("Jan","Feb","Mar", "Apr","May","Jun", "Jul","Aug","Sep", "Oct","Nov","Dec"))


write_csv(minidiscus_calcs, "ASV_species_names_freq_onlyASV.csv")

minidiscus_calcs_names <- read.csv("ASV_species_names_freq.csv") #this file include species names

minidiscus_calcs_complete <-minidiscus_calcs_names

mini <- read_csv("ASV_zero_freq.csv") #this file includes the months where the ASVs were at a frequency of 0

minidiscus_calcs_complete <- minidiscus_calcs_complete %>%
  dplyr::select(tax.Species, asv_freq, month, month_abb)

mini$month_abb <- mymonths[mini$month]

minidiscus_calcs_complete <- rbind(mini, minidiscus_calcs_complete)

minidiscus_calcs_complete <- separate(minidiscus_calcs_complete, tax.Species, sep = "_", into = c("tax.Genus", "ext"), remove= FALSE)

minidiscus_calcs_complete$month_abb <- factor(minidiscus_calcs_complete$month_abb, levels = c("Jan","Feb","Mar", "Apr","May","Jun", "Jul","Aug","Sep", "Oct","Nov","Dec"))

minidiscus_calcs_complete$tax.Genus <- factor(minidiscus_calcs_complete$tax.Genus, levels = c("Chaetoceros", "Minidiscus", "Cyclotella", "Rhizosolenia", "Leptocylindrus", "Ditylum"))


minidiscus_calcs_complete$tax.Species <- factor(minidiscus_calcs_complete$tax.Species)


### calculating the number of ASVs for minidiscus and cyclotella
mini_cyclo <- filt_species_long %>% filter(tax.Genus %in% c("Minidiscus", "Cyclotella", "Chaetoceros"))

mini_cyclo_sum <- mini_cyclo %>% 
  dplyr::group_by(tax.Genus, ASV) %>% 
  dplyr::summarise(sum_ASV = sum(presence)) 

mini_cyclo_sum$presence <- 1
mini_cyclo_sum <- mini_cyclo_sum %>%
  dplyr::group_by(tax.Genus) %>%
  dplyr::summarise(sum_asvs = sum(presence))

#Actual figure 8 plot
minidiscus_test <- minidiscus_calcs_complete %>%
  #drop_na() %>%
  #filter(tax.Genus == "Chaetoceros") %>%
  ggplot(aes(x =as.factor(month_abb), y = asv_freq, fill = tax.Species))+
  geom_bar(position="stack", stat="identity", width = 1, color = "black")+
  scale_fill_manual(values = c("grey95", "pink4", "mediumpurple4", "peachpuff4", "mistyrose","darkmagenta","palevioletred4","palevioletred","lightpink4","pink3","orchid1", "plum1","plum4","rosybrown1","thistle","thistle4","violetred4","lightcoral","purple4","maroon","mediumpurple","lavenderblush3","violet","lightcoral","pink","lavender","lightslateblue","darkgreen", "darkolivegreen3", "darkseagreen", "darkseagreen1", "forestgreen", "darkolivegreen", "lightgreen","wheat4", "wheat3", "lightgoldenrod3","tan4" , "darkblue","lightblue","lightblue4", "mediumaquamarine","steelblue","darkred", "orangered3", "darksalmon","rosybrown1"))+
  #scale_fill_manual(values = c("grey1", "grey4","grey7", "grey10", "grey11","grey14","grey17","grey20","grey23","grey26","grey29", "grey32","grey35","grey38", "grey41","grey44","violetred4","grey47","grey50","grey53","grey56","grey59","grey61","grey64","grey67","grey70","lightslateblue","grey73","darkgreen", "darkolivegreen3", "darkseagreen", "darkseagreen1", "forestgreen", "darkolivegreen", "lightgreen","wheat4", "wheat3", "lightgoldenrod3","tan4" , "darkblue","lightblue","lightblue4", "mediumaquamarine","steelblue","darkred", "orangered3", "darksalmon","rosybrown1"))+
  #scale_color_manual(values = ("hotpink", "pink4", "mediumpurple4", "peachpuff4", "mistyrose","darkmagenta","palevioletred4","palevioletred","lightpink4","pink3","orchid1", "plum1","plum4","rosybrown1","thistle","thistle4","violetred4","lightcoral","purple4","maroon","mediumpurple","lavenderblush3","violet","lightcoral","pink","lavender","lightslateblue","purple","darkgreen", "darkolivegreen3", "darkseagreen", "darkseagreen1", "forestgreen", "darkolivegreen", "lightgreen","wheat4", "wheat3", "lightgoldenrod3","tan4" , "darkblue","lightblue","lightblue4", "mediumaquamarine","steelblue","darkred", "orangered3", "darksalmon","rosybrown1"))+
  #scale_shape_manual(values = c('Actinocyclus_curvatulus'=15, 'Cyclotella_sp1'=0, 'Cyclotella_sp3'=16, 'Cyclotella_striata'=17, 'Minidiscus_comicus'=1, 'Minidiscus_spinulatus' = 2, 'Minidiscus_variabilis' = 3, 'Minutocellus_polymorphus' = 4, 'Navicula_lanceolata' = 5, 'Navicula_sp1' = 6, 'Cyclotella_choctawhatcheeana_litoralis' = 7))+
  #geom_point(aes(fill = tax.Species), pch = 21, color = "black")+
  #geom_line(aes(color = tax.Species), size = .75)+
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.5))+
  facet_wrap(~tax.Genus, ncol = 1)+
  y+
 theme(legend.position = "none")+
  theme(axis.title = element_text(size = 12, family = "Times", face = "bold"))+
  theme(axis.text.y = element_text(size = 10, family = "Times"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1, family = "Times"))+
  theme(strip.text = element_text(size = 12, family = "Times", face = "italic"))+
  theme(strip.text = element_text( face = "italic"))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  labs(y = "% ASV proportion", x = "Month")
#scale_x_discrete(breaks = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))
minidiscus_test


ggsave("Manuscript/Figs/ASV_freq_bars.pdf", plot = minidiscus_test, device = "pdf", width = 2.2, height = 6, units = c("in"),
       dpi = 600)

#comment out the 'legend.postion = "none"' line to get the full legend for the supplemental material
ggsave("ASV_freq_bars_full.pdf", plot = minidiscus_test, device = "pdf", width = 35, height = 20, units = c("cm"),
       dpi = 600)

### END OF CODE FOR FIGURE 8 ### 





