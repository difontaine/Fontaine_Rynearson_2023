### Script for determing spurious reads for NB Diatom Paper
# D.Fontaine

#Read in the count table
mock_count_tab <- read.table("mock_samples/zimmerman_prim_mock_ASVs_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")


#Only need certain samples for mocks
mock_count_tab <- mock_count_tab %>%
  dplyr::select("SF85b", "SF86b","SF87b","SF88b","SF89b","SF90b","SF91b","SF92b","SF93b","SF94b","SF95b","SF96b")
#read in sample info table -- need it to be a tsv for later on analysis otherwise phyloseq doesn't work
sample_info_tab_mock <- read.table("mock_samples/sample_info_mock.tsv", header=T, row.names=1,
                                   check.names=F, sep="\t") #Note to create a .tsv file, you need to save the excel spreadsheet as a regular .txt file and then change the extension of that file to ".tsv"

#read in taxonomy table
asv_tax_mock <- read.table("mock_samples/zimmerman_ASVs_taxonomy.tsv", header=T, row.names=1,
                           check.names=F, sep="\t")

#Make asv_tax a matrix for physeq object to work
asv_tax_mock<- as.matrix(asv_tax_mock)

#Make samplinfo phy object
sample_info_tab_mock <- sample_data(sample_info_tab_mock)

# first we need to create a phyloseq object using our un-transformed count table
count_tab_mock <- otu_table(mock_count_tab, taxa_are_rows=T)
tax_tab_mock <- tax_table(asv_tax_mock) #Have to make the taxa table a matrix in order for the physeq object thing to work
#Making a physeq object
ASV_physeq_mock <- phyloseq(count_tab_mock, tax_tab_mock, sample_info_tab_mock)

mock_otu <- as.data.frame(otu_table(ASV_physeq_mock))
mock_otu <- rownames_to_column(mock_otu, var = "ASV")

mock_taxa <- as.data.frame(tax_table(ASV_physeq_mock))

mock_otu <- mock_otu%>%
  pivot_longer(cols = 2:13,
               names_to = "Sample",
               values_to = "count")

mock_counts_total <- ddply(mock_otu, .(Sample), summarise, total_reads = sum(count))


mock_otu <- mock_otu %>%
  left_join(mock_counts_total, by = 'Sample') %>%
  mutate(per_rel = (count/total_reads)*100)

mock_taxa <- rownames_to_column(mock_taxa, "ASV")

mock_with_taxa <-  mock_otu%>%
  inner_join(mock_taxa, by = "ASV")

#### SUPPLEMENTAL FIGURE 2: SPURIOUS ASVS ######
mock_spurious <- mock_with_taxa %>%
  dplyr::filter(Genus != "Skeletonema") %>%
  dplyr::filter(Genus != "Heterosigma") %>%
  dplyr::filter(Genus != "Thalassiosira") %>%
  dplyr::filter(Genus != "Ditylum") %>%
  dplyr::filter(Genus != "Conticribra") #Conticribra_weisflogii =  Thalassiosira weisflogii (in mocks) and the only Conticribra found is Conticribra weisflogii so this syntax makes sense

mock_spurious$type <- "spurious"

mock_spurious <- mock_spurious %>%
  filter(per_rel !=0)

#PLot of spurious ASVs to help decide what the filtering cut-off will be
mock_spurious_plot_ <- ggplot(mock_spurious, aes(x=as.factor(type), y=per_rel)) + 
  geom_jitter(color="black", size=2, alpha=0.9) +
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  # stat_summary(fun = mean, geom = "point", shape = 23, size = 4)+
  theme_bw() +
  y+
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))+
  labs(x = "Spurious ASVs", y = "% Relative Abundance")+
  geom_hline(yintercept =  0.075, linetype = "dashed")+
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())
mock_spurious_plot_
ggsave("Manuscript/Supplemental_Info/mock_spurious.png", plot = mock_spurious_plot_, device = "png", width = 12, height = 15, units = c("cm"),
       dpi = 300)

#The determined threshold is 0.075%, which is 0.03% highers than the highest spurious read in the mock communities. Now, go back to main script and apply this threshold for data processing

#Counting the diatom percent reads from the mock_taxa_complete DF
moock_counts <- mock_with_taxa_complete %>%
  filter(ASV != "ASV_220") %>% #this is H. akashiwo and isn't a diatom
  group_by(Sample) %>%
  summarise(sum_reads = sum(per_rel))
