#########################################
# Characterizing the effects of calcium 
# and prebiotic fiber on human gut microbiota 
# composition and function using a randomized 
# crossover design- a feasibility study
#
# Lara Yoon 
#  lyoon6@ucla.edu
#
# Analyses
#
# Written with references to Leah Stiemsma
# Last updated: January 2021
# 
# Do not distribute without permission
###########################################

###############
##  LIBRARIES
###############
#install.packages("pacman")
library(pacman)
pacman::p_load("devtools","here", "BiocManager", #tools 
               "assertr", "tidyverse", "reshape2", "janitor", #data management
               "dada2", "phyloseq", "qiime2R", "ggplot2", "DESeq2", "GUniFrac", "vegan", "LDM", "Maaslin2", "lme4", "table1", "microbiome", #analysis
               "scales", "grid", "ggsci", "ggpubr", "viridis", "picante", "knitr", "gridExtra", "Glimma")  # data viz
p_loaded()

## If first time with Maaslin2 or Glimma, run below
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")
# BiocManager::install("Glimma")
# library(devtools)
# devtools::install_github("Maaslin2")
# BiocManager::version()


#################
## TAXA- DATA 
################
## Run this function when dealing with the old OTU tables, need to use this to distinguish your taxa into their appropriate ranks
parse_taxonomy_simple = function(char.vec) {
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
            "Species", "Strain")
  rv = strsplit(char.vec, ";")[[1]]
  rv = c(rv, rep(NA, length(ranks) - length(rv)))
  names(rv) = ranks
  rv
}


## Convert to Phyloseq 
CRC <- import_qiime(otufilename = "CRC_OTUtable.txt", mapfilename= "CRC_mapfile.txt", parseFunction = parse_taxonomy_simple)
CRC
    # otu_table()   OTU Table:         [ 1623 taxa and 60 samples ]
    # sample_data() Sample Data:       [ 60 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks ]

summarize_phyloseq(CRC)
sumCRC <- summarize_phyloseq(CRC)
lapply(sumCRC, function(x) write.table( data.frame(x), 'summarize_phyobj.csv'  , append= T, sep=',' ))


#################
## Sample- DATA 
################
view(sample_data(CRC))

## rename 'key' values 
sample_data(CRC)$key <- recode(sample_data(CRC)$key, "Inulin - 1"="Pre-Inulin", "Inulin - 2"="Post-Inulin", 
                       "Calcium - 1"="Pre-Calcium", "Calcium - 2"="Post-Calcium", 
                       "Inulin & Calcium - 1"="Pre-Combined", "Inulin & Calcium - 2"="Post-Combined" )
## convert 'key' from character to factor
sample_data(CRC)$key <- as_factor(sample_data(CRC)$key)
## reorder factor 
sample_data(CRC)$key <- factor(sample_data(CRC)$key, levels=c("Pre-Inulin","Post-Inulin", "Pre-Calcium" , "Post-Calcium", "Pre-Combined","Post-Combined", "Post Final Washout"))

## rename 'supplement' values 
sample_data(CRC)$supplement <- recode(sample_data(CRC)$supplement, "A"="Inulin", "B"="Calcium", "C"="Combined", "W"="Final Washout")
## convert 'supplement' from character to factor
sample_data(CRC)$supplement <- as_factor(sample_data(CRC)$supplement)


##Add sequence variable 
sample_data(CRC)$seq <- ifelse(sample_data(CRC)$ID %in% c('E3597','R7798'),'ABC', 
                               ifelse(sample_data(CRC)$ID %in% c('S6368'),'ACB',
                                      ifelse(sample_data(CRC)$ID %in% c('F4577'),'BAC',
                                             ifelse(sample_data(CRC)$ID %in% c('G5758', 'U3078'), 'BCA',
                                                    ifelse(sample_data(CRC)$ID %in% c('M9595','Q8747'),'CAB', 'CBA')))))

## Weeks variable 
sample_data(CRC)$Weeks <- factor(sample_data(CRC)$Weeks)
sample_data(CRC)$Weeks


#Baseline and Final variable 
sample_data(CRC)$basefinal <- ifelse(sample_data(CRC)$Weeks==1, "Baseline", 
                                     ifelse(sample_data(CRC)$Weeks==7, "Post Final Washout", NA))
view(sample_data(CRC))
str(sample_data(CRC))


## Create useful subsets 
# Before
CRC.before <-subset_samples(CRC, BA %in% (c("Before")))
CRC.before 
    # otu_table()   OTU Table:         [ 1623 taxa and 27 samples ]
    # sample_data() Sample Data:       [ 27 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks ]

# After 
CRC.after <-subset_samples(CRC, BA %in% (c("After")))
CRC.after
    # otu_table()   OTU Table:         [ 1623 taxa and 26 samples ]
    # sample_data() Sample Data:       [ 26 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks ]

# Inulin
CRC.inulin <-subset_samples(CRC, supplement %in% (c("Inulin")))
CRC.inulin
    # otu_table()   OTU Table:         [ 1623 taxa and 17 samples ]
    # sample_data() Sample Data:       [ 17 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks ]

# Calcium
CRC.calcium <-subset_samples(CRC, supplement %in% (c("Calcium")))
CRC.calcium
    # otu_table()   OTU Table:         [ 1623 taxa and 18 samples ]
    # sample_data() Sample Data:       [ 18 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks 

# Combined
CRC.comb<-subset_samples(CRC, supplement %in% (c("Combined")))
CRC.comb
    # otu_table()   OTU Table:         [ 1623 taxa and 18 samples ]
    # sample_data() Sample Data:       [ 18 samples by 13 sample variables ]
    # tax_table()   Taxonomy Table:    [ 1623 taxa by 8 taxonomic ranks ]





## get tax table 
taxnames <- as.data.frame(tax_table(CRC))
write.csv(taxnames,"taxnames.csv")


#################
##  Filtering
################
get_taxa_unique(CRC, "Phylum")
table.features_phyla <- as.data.frame(table(tax_table(CRC)[,"Phylum"], exclude=NULL))
view(table.features_phyla)

CRC <- subset_taxa(CRC, !is.na(Phylum) & !Phylum %in% c("", " p__"))
table.features_phyla2 <- as.data.frame(table(tax_table(CRC)[,"Phylum"], exclude=NULL))
view(table.features_phyla2)

prevdf = apply(X=otu_table(CRC), 
               MARGIN = ifelse(taxa_are_rows(CRC), yes=1, no=2), 
               FUN=function(x) {sum(x>0)})
# Add taxonomy and total read counts
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CRC),
                    tax_table(CRC))
view(prevdf)

prevdfSUMS <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdfSUMS
names(prevdfSUMS)[2] <- "meanprev"
names(prevdfSUMS)[3] <- "sumprev"
view(prevdfSUMS)


# Subset to remaining phyla 
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(CRC, "Phylum"))

prevfilterplot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(CRC),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") + 
  ggtitle("Taxa prevalence versus total count (filtering step) ")
# No natural separation noted. 
prevfilterplot

# Define prevalence threshold as 2% of total samples (.02*60=1.2)
prevthresh = 0.02 * nsamples(CRC) 
prevthresh  #1.2

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevthresh)]
CRC2 = prune_taxa(keepTaxa, CRC)
CRC
CRC2 #1110 taxa and 60 samples

view(sample_data(CRC2))

summarize_phyloseq(CRC2)

CRC2

getwd()

## get tax table 
taxnames <- as.data.frame(tax_table(CRC2))
write.csv(taxnames,"taxnames_filt.csv")

# get OTU talbe 
otutab <- as.data.frame(otu_table(CRC2))
write.csv(otutab,"otutab_filt.csv")


#################
##  Richness/Diversity 
################

##Richness/diversity plot
p.shan <- plot_richness(CRC, x="key", measures=c("Shannon"), color="BA") + 
  geom_boxplot() +
  geom_jitter(width=0.15) +
  theme_bw() +
  scale_colour_npg() + 
  #scale_fill_npg() + 
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x = element_text(hjust = 0)) + 
  theme(legend.position="none", strip.text.x = element_blank()) +
  xlab("Sample Collection Timepoint") + ylab("Shannon Diversity Index") + 
  coord_cartesian(ylim = c(2, 5)) + 
  annotate("text", x = 1, y = 5, label = "B")
p.shan
  ## Note Shannon index equally weights richness and evenness


# Subset to baseline and final
CRCsub <- subset_samples(CRC, sample_data(CRC)$basefinal %in% (c("Baseline", "Post Final Washout")))
show_col(pal_npg("nrc")(10))
cols <- c("Baseline"= "#3C5488FF", 
          "Post Final Washout"="#00A087FF")

p.shan2 <- plot_richness(CRCsub, x="basefinal", measures=c("Shannon"), color="basefinal") + 
  geom_boxplot() +
  geom_jitter(width=0.15) +
  scale_color_manual(values=cols)+
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none", strip.text.x = element_blank()) +
  xlab(" ") + ylab("Shannon Diversity Index")+ 
  coord_cartesian(ylim = c(2, 5))+ 
  annotate("text", x = 1, y = 5, label = "A")
p.shan2


# plot next ot each other 
library(cowplot)
plot_grid(p.shan2, p.shan, align = "v", ncol = 2, rel_widths = c(7/24, 17/24))
grid.arrange(p.shan2, p.shan, ncol=2)

##Richness/diversity stats
r.shan <- estimate_richness(CRC, measures = 'Shannon')
d = sample_data(CRC)
r.shan <- cbind(r.shan, d)
aov.key<- aov(Shannon~key, data=r.shan)
summary(aov.key)
capture.output(aov.key, file="richness_shannon_aov_REVISIONS.txt")
TukeyHSD(aov.key) ##No significant findings, but did export plot


## Comparison to baseline 
r.shan2 <- estimate_richness(CRCsub, measures = 'Shannon')
d2 = sample_data(CRCsub)
r.shan2 <- cbind(r.shan2, d2)
view(r.shan2)
r.shan2_wide <- dcast(r.shan2, ID ~ basefinal, value.var = "Shannon") 
r.shan2_wide$base <- r.shan2_wide$Baseline
r.shan2_wide$final <- r.shan2_wide$`Post Final Washout`
r.shan2_wide_na <- drop_na(r.shan2_wide)
wilcox.test(r.shan2_wide_na$base, r.shan2_wide_na$final, paired=TRUE, data=r.shan2_wide_na)
ttest_shan <- t.test(r.shan2_wide_na$base, r.shan2_wide_na$final, paired=TRUE, data=r.shan2_wide_na)
ttest_shan
capture.output(ttest_shan, file="richness_shannon_ttest_REVISIONS.txt")



##repeated measures anova:
#(aov(Shannon ~ key + Error(ID/key), data=r.shan)) ##no significance

##write out file to include in LR analysis
write.csv(r.shan, "richness_shannon.csv") 


### After 
r.shan.after <- estimate_richness(CRC.after, measures = 'Shannon')
d = sample_data(CRC.after)
r.shan.after <- cbind(r.shan.after, d)
aov.after <- aov(Shannon ~ key , data=r.shan.after)
summary(aov.after) # p= 0.136
TukeyHSD(aov.after) # no significance





### LMER & AOV tests 
# Weeks = period; Seq = sequence
# fit <- lmer(Shannon ~ key + Weeks + seq + (1|ID), data=r.shan.after)
# anova(fit)
# fit2 <- lmer(Shannon ~ key * Weeks + (1|ID), data=r.shan.after)
# anova(fit2)

## Note: we do not have enough degrees of freedom to include period and sequence effects in the model! 




#################
##  PCoA + PERMANOVA 
################

## Ordinate
ord <- ordinate(CRC, "PCoA", "bray")
## Plot
p.ord <- plot_ordination(CRC, ord, type="samples", color="key")  + 
  theme_bw() +
  geom_point(size=3) +
  scale_colour_npg() + 
  scale_fill_npg() + 
  labs(color = "Sample Collection Timepoint") +
  ggtitle("Beta Diversity, Before and After Each Intervention") + 
  theme(legend.position="right")
p.ord

##permanova for after samples and all the data
# calcluate bray curtis distance matrix
r.bray.dist <- phyloseq::distance(CRC, "bray") 
# make a data frame from sample data 
df <- data.frame(sample_data(CRC))
# Adonis test
r.adon <- adonis(r.bray.dist ~ supplement, data=df)
r.adon # Pr(>F)=0.42
write.csv(r.adon[["aov.tab"]], "adonis.csv") 
    # our adonis test is not significant 
    # so we cannot reject the null hypothesis that 
    # our three sites have the same centroid.


#Homogeneiy of dispersion test 
beta <- betadisper(r.bray.dist, df$supplement)
r.perm <- permutest(beta)
r.perm
write.csv(r.perm[["tab"]], "homogeneity_dispersion_permutest.csv") 
    # r betadisper results are not significant, 
    # meaning we cannot reject the null hypothesis 
    # that our groups have the same dispersions. 



## Permanova notes: 
# Ordinatio (PCoA) - shows whether samples overlap 
# Permanva statistically determine if the centers (centroids) of the cluster of samples

# ANOVA lets one tell if the mean value differs 
# among treatment groups, 
# so the PERMANOVA lets one determine if 
# centroids differ in ordinations

# Permanova not done on the actual ordination but is 
# done on the underlying distance matrices

#################
##  RELATIVE ABUNDANCE
################

##Preprocess your data for relative abundance analysis
top100 <- names(sort(taxa_sums(CRC), decreasing=TRUE))[1:100] #identify the top 100 OTUs
top100_trs <- transform_sample_counts(CRC, function(OTU) OTU/sum(OTU)) ##Transform to relative abundance
top100_prune <- prune_taxa(top100, top100_trs) ##prune according to these parameters

view(sample_data(top100_prune))
view(otu_table(top100_prune))
top100_prune
  # otu_table()   OTU Table:         [ 100 taxa and 60 samples ]
  # sample_data() Sample Data:       [ 60 samples by 15 sample variables ]
  # tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]

##right now, each sample type is represented by the number of samples sequenced, categorical merge to adjust this
merge = merge_samples(top100_prune, "key")
sample_data(merge)$key <- levels(sample_data(CRC)$key)

##transform to % relative abundance 
merge.100 = transform_sample_counts(merge, function(x) 100 * x/sum(x))

## convert 'key' from character to factor
sample_data(merge.100)$key <- as_factor(sample_data(merge.100)$key)
## reorder factor 
sample_data(merge.100)$key <- factor(sample_data(merge.100)$key, levels=c("Pre-Inulin","Post-Inulin", "Pre-Calcium" , "Post-Calcium", "Pre-Combined","Post-Combined", "Post Final Washout"))

## rename col names in tax table 
colnames(tax_table(merge.100)) <- c("Domain", "Phylum", "Class", "Order", "Family",  "Genus", "Species", "Strain")
phyloseq::tax_table(merge.100)[1:8,1:8]
tax_table(merge.100)[, colnames(tax_table(merge.100))] <- gsub(tax_table(merge.100)[, colnames(tax_table(merge.100))], pattern = "[a-z]__", replacement = "")

sample_data(merge.100)

##Plot

# DF for annotations 
nsize <- data.frame(
  timepoint = c("", "Pre-Inulin", "Post-Inulin", "Pre-Calcium", "Post-Calcium", "Pre-Combined", "Post-Combined", "Post Final Washout"), 
  nsize = c("Sample size", 9,8,9,9,9,9,7))
nsizet <- t(nsize)
view(nsizet)

##Phylum - bacteroidetes seems underrepresented given that these are gut samples (Not terribly unusual but could be an artifact of the sequencing)
p.phylum <- plot_bar(merge.100, "key", "Abundance", "Phylum") + 
  xlab("Sample Collection Timepoint") + 
  ylab("Percent Abundance") + 
  ggtitle("Phylum Abundance Before and After Each Intervention") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position='stack') +  
  theme_bw() +
  scale_colour_npg() + 
  scale_fill_npg() + 
  annotate("text", x = 1, y = -2, label = "n=9", size=4)+ 
  annotate("text", x = 2, y = -2, label = "n=8", size=4)+ 
  annotate("text", x = 3, y = -2, label = "n=9", size=4)+ 
  annotate("text", x = 4, y = -2, label = "n=9", size=4)+ 
  annotate("text", x = 5, y = -2, label = "n=9", size=4)+ 
  annotate("text", x = 6, y = -2, label = "n=9", size=4)+ 
  annotate("text", x = 7, y = -2, label = "n=7", size=4)
p.phylum
##bacteroidetes does seem to expand after A and B supplements, suggesting that the supplement promotes growth of this taxon.
##c supplement doesn't seem to do anything (comparing before and after)




##Order abundance with Phylum Wrap
p.order <- plot_bar(merge.100, "key", "Abundance", "Order") + 
  xlab("Sample Collection Timepoint") + 
  ylab("Percent Abundance") + 
  ggtitle("Order Abundance Before and After Each Intervention") + 
  geom_bar(aes(color=Order), stat="identity", position='stack') + 
  facet_wrap(~Phylum) +   
  scale_colour_npg() + 
  scale_fill_npg()+
  theme_bw()
p.order
##Seems like lachno is representing most of the firmicutes, shrinks when bifido expands?

##Genus, tried viewing a few different ways with facet_wraps
p.genus <- plot_bar(merge.100, "key", "Abundance", "Genus") + 
  xlab("Sample Collection Timepoint") + 
  ylab("Percent Abundance") + 
  ggtitle("Genera Abundance Before and After Each Intervention") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position='stack')

##Bifidobacterium changes quite a bit between supplements
p.phylum.2 <- plot_bar(merge.100, "key", "Abundance", "Phylum") + 
  xlab("Sample Collection Timepoint") + ylab("Percent Abundance") + 
  ggtitle("Genera Abundance Before and After Each Intervention") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position='stack') + 
  facet_wrap(~Genus)
p.phylum.2

##Family abundance with Phyla wrap (or opposite) - shows bifido and rumino swapping places before and after each supplement
p <- plot_bar(merge.100, "combined", "Abundance", "Phylum") + 
  xlab("Supplement Collection Points") + ylab("Percent Abundance") + 
  ggtitle("Family Abundance Before and After Each Supplement") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position='stack') + 
  facet_wrap(~Family)


### Tops 
top5phy <- sort(tapply(taxa_sums(merge.100), tax_table(merge.100)[, "Phylum"], sum), decreasing=TRUE)[1:5]
pphy    <- subset_taxa(merge.100, Phylum %in% names(top5phy))
p.phylum.top<-
  plot_bar(pphy, x="key", fill="key", facet_grid= ~ Phylum) +
  geom_bar(aes(color=key, fill=key), stat="identity", position="stack")+
  labs( title = "Phyla")+ 
  theme_classic()+
  xlab("Sample Collection Timepoint")+ 
  theme(legend.position = "none") + 
  scale_colour_npg() + 
  scale_fill_npg()
p.phylum.top 

#family- top5
top5fam <- sort(tapply(taxa_sums(merge.100), tax_table(merge.100)[, "Family"], sum), decreasing=TRUE)[1:5]
pfam    <- subset_taxa(merge.100, Family %in% names(top5fam))
subpo.filt.fam <-
  plot_bar(pfam, x="key", fill="key", 
           facet_grid= ~ Family, title = "Family") + 
  geom_bar(aes(color=key, fill=key), 
           stat="identity", position="stack") +
  theme(axis.text.x = element_text(angle = 90)) + 
  theme_classic()+
  xlab("Sample Collection Timepoint")+ 
  theme(legend.position = "none") + 
  scale_colour_npg() + 
  scale_fill_npg() 
subpo.filt.fam




#####
## Before/Final 
#####

top100_prune
# otu_table()   OTU Table:         [ 100 taxa and 60 samples ]
# sample_data() Sample Data:       [ 60 samples by 15 sample variables ]
# tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]

#subset samples 
# top100_prunesub <- subset_samples(top100_prune, basefinal %in% (c("Baseline", "Post Final Washout")))
# top100_prunesub <- subset_samples(top100_prune, Weeks %in% (c("1", "2", "4", "6", "7"))

# recode 
sample_data(top100_prune)$weekschr <- recode(sample_data(top100_prune)$Weeks, 
                                             "1"="Pre-Study (Baseline)", "2"="Post-Intervention 1", 
                                             "3"="Pre-Intervention 2", "4"="Post-Intervention 2", 
                                             "5"="Pre-Intervention 3", "6"="Post-Intervention 3", 
                                             "7"="Post- Final Washout")
view(sample_data(top100_prune))

##right now, each sample type is represented by the number of samples sequenced, categorical merge to adjust this
merge2 = merge_samples(top100_prune, "weekschr")
sample_data(merge2)$weekschr <- factor(sample_data(merge2)$weekschr)
levels(sample_data(merge2)$weekschr)
sample_data(merge2)$weektime <- recode(sample_data(merge2)$weekschr, 
                                       "1"="Pre-Study (Baseline)", "2"="Post-Intervention 1", 
                                       "3"="Pre-Intervention 2", "4"="Post-Intervention 2", 
                                       "5"="Pre-Intervention 3", "6"="Post-Intervention 3", 
                                       "7"="Post- Final Washout")

sample_data(merge2)


##transform to % relative abundance 
merge2.100 = transform_sample_counts(merge2, function(x) 100 * x/sum(x))

## convert 'key' from character to factor
sample_data(merge2.100)$weektime <- as_factor(sample_data(merge2.100)$weektime)
view(sample_data(merge2.100))

## rename col names in tax table 
colnames(tax_table(merge2.100)) <- c("Domain", "Phylum", "Class", "Order", "Family",  "Genus", "Species", "Strain")
phyloseq::tax_table(merge2.100)[1:8,1:8]
tax_table(merge2.100)[, colnames(tax_table(merge2.100))] <- gsub(tax_table(merge2.100)[, colnames(tax_table(merge2.100))], pattern = "[a-z]__", replacement = "")

merge2.100

p.phylum <- plot_bar(merge2.100, "weektime", "Abundance", "Phylum") + 
  xlab("Sample Collection Timepoint") + 
  ylab("Percent Abundance") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position='stack') +  
  theme_bw() +
  scale_colour_npg() + 
  scale_fill_npg()
p.phylum






#################
##  DESEQ 
## https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
################
## rename taxa names in tax table 
phyloseq::tax_table(CRC)[1:8,1:8]
tax_table(CRC)[, colnames(tax_table(CRC))] <- gsub(tax_table(CRC)[, colnames(tax_table(CRC))], pattern = "[a-z]__", replacement = "")


##DeSeq2 on ALL DATA- not useful bc multi-factorial
# diagdds = phyloseq_to_deseq2(CRC, ~ key) ##convert to deseq object & estimate disperions
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric") #fit & run your model (default multiple-inference correction is Benjamini-Hochberg)
# res = results(diagdds, cooksCutoff = FALSE) # save results 
# res = res[order(res$padj, na.last=NA),] # order by adjp value & remove entries with NA
# alpha = 0.01 # set alpha 
# sigtab = res[which(res$padj < alpha), ] # subset significant results 
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(CRC)[rownames(sigtab), ], "matrix")) # merge with tax table
# head(sigtab) # view significant results 
# posigtab <- sigtab[sigtab[, "log2FoldChange"] !=0, ]  # where log2FoldChange is not 0
# posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
# posigtab #clean table for legibility 


##Top 100 OTUs- need to specify contrast in results
top100 <- names(sort(taxa_sums(CRC), decreasing=TRUE))[1:100]
top100_prune <- prune_taxa(top100, CRC)
diagdds = phyloseq_to_deseq2(top100_prune, ~ key)
diagdds <- diagdds[ rowSums(counts(diagdds)) > 1, ]
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

alpha = 0.1

res.inulin = results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Inulin", "Pre-Inulin"))
res.calcium = results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Calcium", "Pre-Calcium"))
res.combined = results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Combined", "Pre-Combined"))

# Inulin results 
sigtab.in = res.inulin[which(res.inulin$padj < alpha), ]
sigtab.in = cbind(as(sigtab.in, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.in), ], "matrix"))
sigtab.in
posigtab.in <- sigtab.in[sigtab.in[, "log2FoldChange"] !=0, ] 
posigtab.in <- posigtab.in[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.in

# Calcium results 
sigtab.ca = res.calcium[which(res.calcium$padj < alpha), ]
sigtab.ca = cbind(as(sigtab.ca, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.ca), ], "matrix"))
sigtab.ca
posigtab.ca <- sigtab.ca[sigtab.ca[, "log2FoldChange"] !=0, ] 
posigtab.ca <- posigtab.ca[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.ca

# Combined results
sigtab.co = res.combined[which(res.combined$padj < alpha), ]
sigtab.co = cbind(as(sigtab.co, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.co), ], "matrix"))  # nothing significant
sigtab.co 
posigtab.co <- sigtab.co[sigtab.co[, "log2FoldChange"] !=0, ] 
posigtab.co <- posigtab.co[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.co

    # R will choose ref level for factors based on alphabetical order
    # for factors, comparison is last level of this variable over the first level
    # Results tables are generated using the functionresults, which extracts a results table with log2 fold changes,pvalues and adjustedpvalues.
    # log2FoldChange = the estimates are of the logarithmic fold changelog2(treated/untreated)


# After comparisons
res.Ca.In = results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Calcium", "Post-Inulin"))
res.Co.In= results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Combined", "Post-Inulin"))
res.Co.Ca = results(diagdds, cooksCutoff = FALSE, contrast=c("key", "Post-Combined", "Post-Calcium"))

alpha2=0.1

# Calcium vs Inulin results 
sigtab.1 = res.Ca.In[which(res.Ca.In$padj < alpha2), ]
sigtab.1 = cbind(as(sigtab.1, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.1), ], "matrix"))
sigtab.1
posigtab.1 <- sigtab.1[sigtab.1[, "log2FoldChange"] !=0, ] 
posigtab.1 <- posigtab.1[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.1

# Combined vs Inulin results 
sigtab.2 = res.Co.In[which(res.Co.In$padj < alpha2), ]
sigtab.2 = cbind(as(sigtab.2, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.2), ], "matrix"))
sigtab.2
posigtab.2 <- sigtab.2[sigtab.2[, "log2FoldChange"] !=0, ] 
posigtab.2 <- posigtab.2[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.2

# Combined vs Calcium results
sigtab.3 = res.Co.Ca[which(res.Co.Ca$padj < alpha2), ]
sigtab.3 = cbind(as(sigtab.3, "data.frame"), as(tax_table(top100_prune)[rownames(sigtab.3), ], "matrix"))
sigtab.3
posigtab.3 <- sigtab.3[sigtab.3[, "log2FoldChange"] !=0, ] 
posigtab.3 <- posigtab.3[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
posigtab.3


## Save to file 
write.csv(res.Ca.In@listData, "DeSeq2_PostCa-PostIn.csv")
write.csv(res.Co.In@listData, "DeSeq2_PostCo-PostIn.csv")
write.csv(res.Co.Ca@listData, "DeSeq2_PostCo-PostCa.csv")

write.csv(res.inulin@listData, "DeSeq2_Post-Pre_In.csv")
write.csv(res.calcium@listData, "DeSeq2_Post-Pre_Ca.csv")
write.csv(res.combined@listData, "DeSeq2_Post-Pre_Co.csv")




## Plotting 
#Note: Run and save in order 
# Calcium vs Inulin results 
x = tapply(sigtab.1$log2FoldChange, sigtab.1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.1$Phylum = factor(as.character(sigtab.1$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab.1$log2FoldChange, sigtab.1$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.1$Genus = factor(as.character(sigtab.1$Genus), levels=names(x))
p.log2fold.CaIn <- ggplot(sigtab.1, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=6) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  scale_colour_npg() + 
  scale_fill_npg()+
  ggtitle("Log2-Fold Change of OTUs, Post-Calcium vs Post-Inulin Timepoints ")
p.log2fold.CaIn

# Combined vs Inulin results 
x = tapply(sigtab.2$log2FoldChange, sigtab.2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.2$Phylum = factor(as.character(sigtab.2$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab.2$log2FoldChange, sigtab.2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.2$Genus = factor(as.character(sigtab.2$Genus), levels=names(x))
p.log2fold.CoIn <- ggplot(sigtab.2, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=6) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  scale_colour_npg() + 
  scale_fill_npg() + 
  ggtitle("Log2-Fold Change of OTUs, Post-Combined vs Post-Inulin Timepoints ")
p.log2fold.CoIn

# Combined vs Calcium results 
x = tapply(sigtab.3$log2FoldChange, sigtab.3$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.3$Phylum = factor(as.character(sigtab.3$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab.3$log2FoldChange, sigtab.3$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.3$Genus = factor(as.character(sigtab.3$Genus), levels=names(x))
p.log2fold.CoCa <- ggplot(sigtab.3, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=6) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  scale_colour_npg() + 
  scale_fill_npg()+
  ggtitle("Log2-Fold Change of OTUs, Post-Combined vs Post-Calcium Timepoints ")
p.log2fold.CoCa


########################## OLD CODE ######################## 
# ##Inulin-only comparison with top 100 abundant 
# i.100_prune <- prune_taxa(top100,CRC.inulin)
# 
# diagdds.i = phyloseq_to_deseq2(i.100_prune, ~ key)
# diagdds.i = DESeq(diagdds.i, test="Wald", fitType="parametric")
# 
# res.i = results(diagdds.i, cooksCutoff = FALSE)
# res.i
# alpha = 0.01
# sigtab.i = res.i[which(res.i$padj < alpha), ]
# sigtab.i = cbind(as(sigtab.i, "data.frame"), as(tax_table(i.100_prune)[rownames(sigtab.i), ], "matrix"))
# head(sigtab.i)
# sigtab.i
# posigtab.i <- sigtab.i[sigtab.i[, "log2FoldChange"] !=0, ] 
# posigtab.i <- posigtab.i[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
# posigtab.i
# 
# 
# ##Calcium-only comparison
# ca.100_prune <- prune_taxa(top100,CRC.calcium)
# 
# diagdds.ca = phyloseq_to_deseq2(ca.100_prune, ~ key)
# diagdds.ca = DESeq(diagdds.ca, test="Wald", fitType="parametric")
# 
# res.ca = results(diagdds.ca, cooksCutoff = FALSE)
# alpha = 0.01
# sigtab.ca = res.ca[which(res.ca$padj < alpha), ]
# sigtab.ca = cbind(as(sigtab.ca, "data.frame"), as(tax_table(ca.100_prune)[rownames(sigtab.ca), ], "matrix"))
# head(sigtab.ca)
# 
# posigtab.ca <- sigtab.ca[sigtab.ca[, "log2FoldChange"] !=0, ] 
# posigtab.ca <- posigtab.ca[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
# posigtab.ca
# 
# 
# ##Combined comparison
# co.100_prune <- prune_taxa(top100,CRC.comb)
# 
# diagdds.co = phyloseq_to_deseq2(co.100_prune, ~ key)
# diagdds.co = DESeq(diagdds.co, test="Wald", fitType="parametric")
# 
# res.co = results(diagdds.co, cooksCutoff = FALSE)
# alpha = 0.01
# sigtab.co = res.co[which(res.co$padj < alpha), ]
# sigtab.co = cbind(as(sigtab.co, "data.frame"), as(tax_table(co.100_prune)[rownames(sigtab.co), ], "matrix"))
# head(sigtab.co)
# 
# posigtab.co <- sigtab.co[sigtab.co[, "log2FoldChange"] !=0, ] 
# posigtab.co <- posigtab.co[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
# posigtab.co


################## EXTRA CODE ###############
##  (Leah)
#################

# ##Looked into mt test - need to review this more
# ##write out file
# top100 <- names(sort(taxa_sums(CRC), decreasing=TRUE))[1:100] #identify the top 100 OTUs
# top100_trs <- transform_sample_counts(CRC, function(OTU) OTU/sum(OTU)) ##Transform to relative abundance
# top100_prune <- prune_taxa(top100, top100_trs) ##prune according to these parameters
# df <- psmelt(top100_prune_trs) ##coerces to df but totally reorganized
# 
# ##coreced this way, better
# OTU1 = as(otu_table(top100_prune_trs), "matrix") #extract otu data from the phyloseq object
# if(taxa_are_rows(top100_prune_trs)){OTU1 <- t(OTU1)}
# # Coerce to data.frame
# OTUdf = as.data.frame(OTU1) ##missing taxa info
# 
# ##extract taxa too
# taxa = as(tax_table(top100_prune_trs), "matrix") 
# 
# ##Write out, have to write out as CSVs - tab delimited doesn't work well.
# write.csv(OTUdf, "~/Desktop/MIcrobiome files/top100OTUs.csv")
# write.csv(taxa, "~/Desktop/MIcrobiome files/top100OTUsTAXA.csv")



#################
##  SCFA 
################
library(readxl)

##SCFA data - for duplicates, took an average of the two readings
SCFAleah <- read.csv("CRC_SCFAdata.csv", header=T) ##had to reorg (see file)
SCFA <- SCFAleah

# # Raw SCFA values (de-id) 
# SCFAraw <- read_excel("SCFA_raw.xlsx")
# SCFAraw <- select(SCFAraw, -c(DataFile))
# # Code Key
# key <- read_excel("SFCA_IDkey_WCMC.xlsx")
# # Merge 
# SCFA <- full_join(SCFAraw, key, by=c("SampleName")) 
# # Long 
# names(SCFA)
# SCFAlong <- melt(SCFA, id.vars=c("SampleID", "Collection"), 
#                  measure.vars=c("Acetic acid", "Butyric acid","Formic acid", "Isovaleric acid","Propionic acid"))



## rename 'key' values 
SCFA$Supplement_labels <- recode(SCFA$Supplement_labels, "Inulin - 1"="Pre-Inulin", "Inulin - 2"="Post-Inulin", 
                               "Calcium - 1"="Pre-Calcium", "Calcium - 2"="Post-Calcium", 
                               "Inulin & Calcium - 1"="Pre-Combined", "Inulin & Calcium - 2"="Post-Combined", 
                               "Washout"="Post Final Washout")
## convert 'key' from character to factor
SCFA$Supplement_labels<- as_factor(SCFA$Supplement_labels)
## reorder factor 
SCFA$Supplement_labels <- factor(SCFA$Supplement_labels, levels=c("Pre-Inulin","Post-Inulin", "Pre-Calcium" , "Post-Calcium", "Pre-Combined","Post-Combined", "Post Final Washout"))

## rename 'supplement' values 
SCFA$Supplement <- recode(SCFA$Supplement, 
                                      "A1"="Inulin", "A2"="Inulin", 
                                      "B1"="Calcium", "B2"="Calcium", 
                                      "C1"="Combined", "C2"="Combined",
                                      "W"="Final Washout")
## convert 'supplement' from character to factor
SCFA$Supplement <- as_factor(SCFA$Supplement)

## Recode SCFA 
SCFA$SCFA <- recode(SCFA$SCFA, 
                    "propionic acid"="Propionic acid", "Acetic Acid"="Acetic acid")
## rename 'supplement' values 
SCFA$timepoint <- recode(SCFA$timepoint, "S1"="Before", "S3"="Before", "S5"="Before", 
                         "S2"="After", "S4"="After", "S6"="After", "S7"="Final Washout")
## convert 'supplement' from character to factor
SCFA$timepoint <- as_factor(SCFA$timepoint)


names(SCFA)

## Boxplots 
box.scfa <- ggplot(SCFA, aes(x=Supplement_labels, y=Concentration, fill=timepoint))+
  geom_boxplot() + 
  facet_grid(.~SCFA) + 
  theme_bw()+
  scale_colour_npg() + 
  scale_fill_npg() +
  xlab("Sample Collection Timepoint") + 
  ylab("Fecal Concentration of SCFA (ng/ng feces)") + 
  ggtitle("Fecal Concentrations of Short Chain Fatty Acids (SCFAs), Before and After each Intervention") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none")
box.scfa

# Adding Sample size to the bottom for each 

library(dplyr)
SCFAcount <- SCFA %>%
  group_by( SCFA, Supplement_labels) %>% 
  tally 
view(SCFAcount)
SCFAcount <- t(SCFAcount)




### Anova & Repeated Measures ANOVA for each SCFA

##Formic
formic <- subset(SCFA, SCFA %in% c("Formic acid"))
aovformic<- aov(Concentration~Supplement_labels, data=formic)
summary(aovformic)
TukeyHSD(aovformic)
##repeated measures anova
formicaov <- summary(aov(Concentration ~ Supplement_labels + Error(subjectID/Concentration), data=formic)) ##no significance
formicaov

## Acetic
acetic <- subset(SCFA, SCFA %in% c("Acetic acid"))
aovacetic<- aov(Concentration~Supplement_labels, data=acetic)
summary(aovacetic)
TukeyHSD(aovacetic)
##repeated measures anova
aceticaov <- summary(aov(Concentration ~ Supplement_labels + Error(subjectID/Concentration), data=acetic)) ##no significance
aceticaov

## Butyric
butyric <- subset(SCFA, SCFA %in% c("Butyric acid"))
aovbutyric<- aov(Concentration~Supplement_labels, data=butyric)
summary(aovbutyric)
TukeyHSD(aovbutyric)
##repeated measures anova
butyricaov <- summary(aov(Concentration ~ Supplement_labels + Error(subjectID/Concentration), data=butyric)) ##no significance
butyricaov

## Isovaleric
isovaleric <- subset(SCFA, SCFA %in% c("Isovaleric acid"))
aovisovaleric<- aov(Concentration~Supplement_labels, data=isovaleric)
summary(aovisovaleric)
TukeyHSD(aovisovaleric)
##repeated measures anova
isovalericaov <- summary(aov(Concentration ~ Supplement_labels + Error(subjectID/Concentration), data=isovaleric)) ##no significance
isovalericaov

## Propionic
propionic <- subset(SCFA, SCFA %in% c("Propionic acid"))
aovpropionic<- aov(Concentration~Supplement_labels, data=propionic)
summary(aovpropionic)
TukeyHSD(aovpropionic)
##repeated measures anova
propionicaov <- summary(aov(Concentration ~ Supplement_labels + Error(subjectID/Concentration), data=propionic)) ##no significance
propionicaov



##tried log transforming, shows lots of variation in the data.
SCFA_log <- log(SCFA$Concentration)
as.data.frame(SCFA_log)
SCFA <- cbind(SCFA_log, SCFA)
p<- ggplot(SCFA, aes(x=Supplement, y=SCFA_log, fill=SCFA))+geom_bar(stat="identity", position=position_dodge())
p


################## EXTRA CODE ###############
##  (Leah)
#################
# ## Barplots 
# bar.scfa<- ggplot(SCFA, aes(x=Supplement, y=Concentration, fill=SCFA)) + 
#   geom_bar(stat="identity", position='stack') ##stacks bars
# bar.scfa
# bar.scfa2<- ggplot(SCFA, aes(x=Supplement, y=Concentration, fill=SCFA)) + 
#   geom_bar(stat="identity", position=position_dodge()) ##side-by-side bars, can do position_dodge2() to have bars for all subjects (after doing this, it seems like there are some subjects with high concentrations of all SCFAs)
# bar.scfa2
# 
# ##function to create error bars
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# ##dataframe of mean and sd
# SCFA_sd <- data_summary(SCFA, varname="Concentration", 
#                                               groupnames=c("Supplement_labels", "SCFA"))
# SCFA_sd
# ##plot with error bars 
# bar.scfa.e<- ggplot(SCFA_sd, aes(x=Supplement_labels, y=Concentration, fill=SCFA)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=as.numeric(Concentration)-as.numeric(sd), ymax=Concentration+sd), width=.2,position=position_dodge(.9))
# bar.scfa.e
# ##pretty significant error bars (sd) for butyric acid and acetic acid. Suggests significant variation in the SCFA concentrations between individuals. 
# 
# ##add title, etc. 
# bar.scfa.e + 
#   labs(title="Short chain fatty acid concentration per intervention timepoint", x="Intervention Time Point ", y = "Concentration (ng/ng feces)") +
#   theme_classic()
# 
# ##line plot
# line.scfa <- ggplot(SCFA_sd, aes(x=Supplement_labels, y=Concentration, group=SCFA, color=SCFA)) + geom_line() + geom_point() +
#   geom_errorbar(aes(ymin=Concentration-sd, ymax=
#                       Concentration+sd), width=.2,
#                 position=position_dodge(0.05))
# line.scfa






#################
##  LBP 
################
RA_OTUs <- read_delim("C:/Users/laray/Box/Projects/2019-CRC/03_Raw data/RA_OTUs.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
## rename 'key' values 
RA_OTUs$key <- recode(RA_OTUs$key, "Inulin - 1"="Pre-Inulin", "Inulin - 2"="Post-Inulin", 
                                 "Calcium - 1"="Pre-Calcium", "Calcium - 2"="Post-Calcium", 
                                 "Inulin & Calcium - 1"="Pre-Combined", "Inulin & Calcium - 2"="Post-Combined", 
                                 "Washout"="Post Final Washout")
## convert 'key' from character to factor
RA_OTUs$key<- as_factor(RA_OTUs$key)
## reorder factor 
RA_OTUs$key <- factor(RA_OTUs$key, levels=c("Pre-Inulin","Post-Inulin", "Pre-Calcium" , "Post-Calcium", "Pre-Combined","Post-Combined", "Post Final Washout"))

## create 'supplement' values 
RA_OTUs$supp <- ifelse(RA_OTUs$key %in% c("Pre-Inulin","Post-Inulin"), "Inulin", 
                       ifelse(RA_OTUs$key %in% c("Pre-Calcium","Post-Calcium"), "Calcium", 
                              ifelse(RA_OTUs$key %in% c("Pre-Combined","Post-Combined"), "Combined", "Final Washout")))
                       
## convert 'supplement' from character to factor
RA_OTUs$supp <- as_factor(RA_OTUs$supp)


## convert LBP ng to ug 
RA_OTUs$LBPug <- RA_OTUs$LBP * 0.001
RA_OTUs <- relocate(RA_OTUs, "LBPug", .after=LBP)

## convert LBP to log
RA_OTUs$LBPlog <- log10(RA_OTUs$LBPug)
RA_OTUs <- relocate(RA_OTUs, "LBPlog", .after=LBPug)

## Weeks variable 
RA_OTUs$Weeks <- factor(RA_OTUs$Weeks )
RA_OTUs$Weeks 


### Boxplot
comparisonslbp <- list(c("Pre-Inulin","Post-Inulin"), 
                       c("Pre-Calcium","Post-Calcium" ), 
                       c("Pre-Combined","Post-Combined"))
comparisonslbpweek <- list(c("1","3"), 
                       c("1","5" ), 
                       c("1","7"))
  
  
box.lbp <- ggplot(RA_OTUs, aes(x=Weeks, y=LBPug, fill=BA))+
  geom_boxplot() + 
  theme_bw()+
  scale_colour_npg() + 
  scale_fill_npg() +
  xlab("Sample Collection Timepoint") + 
  ylab("Serum Concentration of LBP (ug/mL)") + 
  labs(title="Serum Concentrations of Lipopolysaccharide (LPS)-binding protein (LBP)", subtitle="Before and After each Intervention") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none") + 
  stat_compare_means(method="anova", label.y=-1, label.x.npc="left") +
  stat_compare_means(comparisons = comparisonslbpweek, label="p.format", tip.length = .01)
box.lbp

box.lbplog <- ggplot(RA_OTUs, aes(x=key, y=LBPlog, fill=BA))+
  geom_boxplot() + 
  theme_bw()+
  scale_colour_npg() + 
  scale_fill_npg() +
  xlab("Sample Collection Timepoint") + 
  ylab("Log10 Serum Concentration of LBP") + 
  labs(title="Serum Concentrations of Lipopolysaccharide (LPS)-binding protein (LBP)", subtitle="Before and After each Intervention") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none")
box.lbplog


### Summary 
table1(~LBPug | key, data=RA_OTUs, overall="Overall")
table1(~LBPug | ID, data=RA_OTUs, overall="Overall")

### Anova &  ANOVA for LBP

aovlbp.key <- aov(LBPug~key, data=RA_OTUs)
summary(aovlbp.key)
TukeyHSD(aovlbp.key)

aovlbp.key<- aov(LBPug~Weeks, data=RA_OTUs)
summary(aovlbp.key)
TukeyHSD(aovlbp.key) # no significance 
DunnettTest(x=RA_OTUs$LBPug, g=RA_OTUs$Weeks) # no significance

##After
after <- subset(RA_OTUs, BA %in% c("After"))
aovlbp<- aov(LBPug~key, data=after)
summary(aovlbp)
TukeyHSD(aovlbp)

##Inulin
inulin <- subset(RA_OTUs, supp %in% c("Inulin"))
aovlbp.in<- aov(LBPug~key, data=inulin)
summary(aovlbp.in)
TukeyHSD(aovlbp.in)

##Calcium
calcium <- subset(RA_OTUs, supp %in% c("Calcium"))
aovlbp.ca<- aov(LBPug~key, data=calcium)
summary(aovlbp.ca)
TukeyHSD(aovlbp.ca)

##Combined
combined <- subset(RA_OTUs, supp %in% c("Combined"))
aovlbp.co<- aov(LBPug~key, data=combined)
summary(aovlbp.co)
TukeyHSD(aovlbp.co)







#################
##  CARRYOVER EFFECT TESTING  
################

crc.co <- subset_samples(CRC, Weeks %in% (c("1", "3", "5", "7")))
crc.co

sample_data(crc.co)$cotime <- recode(sample_data(crc.co)$Weeks, "1"="Pre-Study (Baseline)", "3"="Pre-Intervention 2", 
                                     "5"="Pre-Intervention 3", "7"="Post- Final Washout")
view(sample_data(crc.co))

## ALPHA DIVERSITY 
comparisons <- list(c("Pre-Study (Baseline)", "Pre-Intervention 2"), 
                    c("Pre-Study (Baseline)", "Pre-Intervention 3"), 
                    c("Pre-Study (Baseline)", "Post- Final Washout"))



adiv.co <- plot_richness(crc.co, x="cotime", measures=c("Shannon")) + 
  geom_boxplot() +
  geom_jitter(width=0.15) +
  theme_bw() +
  scale_colour_npg() + 
  #scale_fill_npg() + 
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x = element_text(hjust = .5)) + 
  theme(legend.position="none", strip.text.x = element_blank()) +
  xlab("Sample Collection Timepoint") + ylab("Shannon Diversity Index") + 
  stat_compare_means(comparisons=comparisons, label="p.format", tip.length = .01)
adiv.co



## SCFA 
comparisons <- list(c("Pre-Study (Baseline)", "Pre-Intervention 2"), 
                    c("Pre-Study (Baseline)", "Pre-Intervention 3"), 
                    c("Pre-Study (Baseline)", "Post- Final Washout"))
view(SCFA)
names(SCFA)
SCFAco <- SCFA %>%  
  filter(str_detect(New_IDS, 'S1|S3|S5|S7')) %>%  
  mutate(cotime = if_else(str_detect(New_IDS, 'S1'), "Pre-Study (Baseline)", 
                          if_else(str_detect(New_IDS, 'S3'), "Pre-Intervention 2", 
                                  if_else(str_detect(New_IDS, 'S5'), "Pre-Intervention 3", "Post- Final Washout"))))

SCFAco$cotime <- as.factor(SCFAco$cotime)
SCFAco$cotimeord <- ordered(SCFAco$cotime, levels = c("Pre-Study (Baseline)", "Pre-Intervention 2", "Pre-Intervention 3", "Post- Final Washout"))
                                                  

scfa.co <- ggplot(SCFAco, aes(x=cotimeord, y=Concentration))+
  geom_boxplot() + 
  facet_grid(.~SCFA) + 
  theme_bw()+
  #scale_colour_npg() + 
  #scale_fill_npg() +
  xlab("Sample Collection Timepoint") + 
  ylab("Fecal Concentration of SCFA (ng/ng feces)") + 
  #ggtitle("Fecal Concentrations of Short Chain Fatty Acids (SCFAs), Before and After each Intervention") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none") + 
  #stat_compare_means(method="anova", label.y=-1, label.x.npc="left") + 
  stat_compare_means(comparisons=comparisons, label="p.format", tip.length = .01)
scfa.co


## LBP 
lbpco <- SCFAco %>%  
  select(c("New_IDS", "LBP", "cotimeord")) %>%  
  distinct()
lbpco$LBPug <- lbpco$LBP * 0.001

lbp.co <- ggplot(lbpco, aes(x=cotimeord, y=LBPug))+
  geom_boxplot() + 
  #facet_grid(.~SCFA) + 
  theme_bw()+
  #scale_colour_npg() + 
  #scale_fill_npg() +
  xlab("Sample Collection Timepoint") + 
  ylab("Serum Concentration of LBP (ug/mL)") + 
  #ggtitle("Fecal Concentrations of Short Chain Fatty Acids (SCFAs), Before and After each Intervention") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="none") + 
 # stat_compare_means(method="anova", label.y=-1, label.x.npc="left") + 
  stat_compare_means(comparisons=comparisons, label="p.format", tip.length = .01)
lbp.co


