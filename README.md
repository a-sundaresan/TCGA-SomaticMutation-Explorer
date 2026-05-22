# TCGA Somatic Mutation Explorer

![Language](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Data](https://img.shields.io/badge/Data-TCGA%20GDC-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)

Exploratory analysis of TCGA pan-cancer somatic mutation data using R.
Focuses on KRAS mutation frequency, amino acid changes, and clinical correlations across cancer types.

## Table of Contents
## Table of Contents
- [Data Sources](#data-sources)
- [KRAS Mutation Frequency Across Cancer Types](#kras-mutation-frequency-across-cancer-types)
- [KRAS Silent Mutation Frequency](#kras-silent-mutation-frequency)
- [KRAS Amino Acid Change Distribution](#kras-amino-acid-change-distribution)
- [G12D Mutation Analysis](#g12d-mutation-analysis)
- [Clinical Correlation in Pancreatic Cancer PAAD](#clinical-correlation-in-pancreatic-cancer-paad)
---

## Key Findings

- **COAD dominates KRAS mutations** — Colon Adenocarcinoma (COAD) has the highest frequency of KRAS mutations across all TCGA cancer types, followed by Lung Adenocarcinoma (LUAD)
- **UCEC leads silent mutations** — Uterine Corpus Endometrial Carcinoma (UCEC) shows the highest frequency of KRAS silent mutations, suggesting a distinct mutational mechanism
- **G12D is the most prevalent KRAS variant** — Among all amino acid changes, G12D is the most common KRAS mutation in the TCGA pan-cancer cohort
- **G12D is enriched in COAD and PAAD** — The G12D hotspot mutation is most prevalent in Colon Adenocarcinoma and Pancreatic Adenocarcinoma, consistent with its known role as a driver mutation in these cancer types
- **No significant age difference in PAAD** — Age at diagnosis is similar between KRAS mutated (mean 65.2y) and wild-type (mean 64.4y) patients (Wilcoxon p = 0.4522), suggesting KRAS mutation status does not influence age of onset in pancreatic cancer
---

**Data:** TCGA MAF file + clinical data (publicly available via GDC)  
**Tools:** R, ggplot2, dplyr, tidyr

## Data Sources
#### Here we are going to do some basic filtering on the TCGA somatic mutation data, combine with the clinical data and use simple visualizations to interrogate the data
TCGA somatic mutation data was downloaded from [here](https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc).\
This file is in Mutation Annotation Format. To know more about this use the info [here](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format). \
TCGA clinical data was downloaded from [here](https://api.gdc.cancer.gov/data/0fc78496-818b-4896-bd83-52db1f533c5c)


```r 
#Loading the necessary libaries 
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
```
```r
#Read the TCGA somatic mutation data
tcga_mutation_data<-read.delim("mc3.v0.2.8.PUBLIC.maf",
                               sep="\t",stringsAsFactors = F)
#Read TCGA  clinical data
tcga_clinical_data<-read.delim("clinical_PANCAN_patient_with_followup.tsv",
                               sep="\t",stringsAsFactors = F)

#Lets create a ID column in the tcga_mutation_data with first 12 characters of “Tumor_Sample_Barcode”
tcga_mutation_data$id<-substring(tcga_mutation_data$Tumor_Sample_Barcode,0,12)
```
In this guide we will focus on the KRAS gene. \
We can use the same data to investigate any gene/s of interest

```r
#Filtering the mutation data ONLY for KRAS mutations
kras_tcga_mutation_data<-tcga_mutation_data %>% filter(Hugo_Symbol=="KRAS")

#Filtering the clinical data ONLY for patients with KRAS mutations
KRAS_clinical<-filter(tcga_clinical_data,bcr_patient_barcode %in% kras_tcga_mutation_data$id)
```
## KRAS Mutation Frequency Across Cancer Types

Let us compute the frequency of KRAS mutations.\
Here I am grouping by the acronym column and computing the number of mutations and frequency of mutations for each specific cancer type.\
The acronym column is the short code for the different cancer types

```r
kras_mutation_freq_plot_df <- KRAS_clinical %>% group_by(acronym) %>% 
                                                dplyr::summarise(number_of_mutations = n()) %>% 
                                                mutate(freq=number_of_mutations/sum(number_of_mutations)) %>% 
                                                arrange(desc(freq))

#Lets us now create a basic bar plot with x axis as the cancer type and y axis for frequency of KRAS mutations
kras_freq_plot <- ggplot(kras_mutation_freq_plot_df, 
                         aes(x=reorder(acronym,(-freq)), fill=acronym,group=acronym))+
                  geom_bar(aes(y=freq),stat = "identity",position = "dodge")

# adding title,x and y labels as well as modifying the theme aesthetics of the plot
kras_freq_plot <- kras_freq_plot + 
                  theme_bw() +
                  theme(axis.text.x = element_text(size = 14,face = "bold"), 
                        axis.title.x = element_text(size = 16,face = "bold"),
                        axis.text.y = element_text(size = 14,face = "bold"), 
                        axis.title.y = element_text(size = 16,face = "bold"),
                        plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
                        plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
                        legend.text=element_text(size = 14,face = "bold"),
                        panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        legend.position="none",
                        legend.title = element_blank())+
                  labs(x = "Cancer Type",
                       y = "Frequency of KRAS mutations", 
                       subtitle = "KRAS",
                       title = "Mutation frequency plot for cancer types")
kras_freq_plot

rm(kras_mutation_freq_plot_df,kras_freq_plot)
```
![KRAS mutation frequency plot](/plots/KRAS_mutation_frequency_plot.png)

Based on the plot, we can see that "COAD" (Colon Adenocarcinoma) has the highest frequency of KRAS mutations followed by "LUAD" (Lung Adenocarcinoma).\
Let's now filter the mutation data ONLY for KRAS silent mutations

## KRAS Silent Mutation Frequency

```r 
kras_tcga_silent_mutation_data<-kras_tcga_mutation_data %>% filter(Variant_Classification=="Silent")

#Filtering the clinical data ONLY for patients with KRAS silent mutations
KRAS_clinical_silent_mutation<-filter(KRAS_clinical,bcr_patient_barcode %in% kras_tcga_silent_mutation_data$id)
```

Let us compute the frequency of KRAS silent mutations.\
As we did in the previous plot, I am grouping by the acronym column and 
computing the number of mutations and frequency of silent mutations for each specific cancer type.

```r
#frequency of KRAS silent mutations
kras_silent_mutation_freq_plot_df <- KRAS_clinical_silent_mutation %>% group_by(acronym) %>% 
                                                                       dplyr::summarise(number_of_mutations = n()) %>% 
                                                                       mutate(freq=number_of_mutations/sum(number_of_mutations)) %>% 
                                                                       arrange(desc(freq))

#Lets us now create a basic bar plot with x axis as the cancer type and y axis for frequency of KRAS silent mutations
kras_silent_mutation_freq_plot <- ggplot(kras_silent_mutation_freq_plot_df, 
                                         aes(x=reorder(acronym,(-freq)), fill=acronym,group=acronym))+
                                    geom_bar(aes(y=freq),stat = "identity",position = "dodge")
# adding title,x and y labels as well as modifying the theme aesthetics of the plot
kras_silent_mutation_freq_plot <- kras_silent_mutation_freq_plot + 
                                  theme_bw() +
                                  theme(axis.text.x = element_text(size = 14,face = "bold"), 
                                        axis.title.x = element_text(size = 16,face = "bold"),
                                        axis.text.y = element_text(size = 14,face = "bold"), 
                                        axis.title.y = element_text(size = 16,face = "bold"),
                                        plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
                                        plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
                                        legend.text=element_text(size = 14,face = "bold"),
                                        panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                        legend.position="none",
                                        legend.title = element_blank())+
                                  labs(x = "Cancer Type",
                                       y = "Frequency Of KRAS silent mutations", 
                                       subtitle = "KRAS",
                                       title = "Silent mutation frequency plot for cancer types")

kras_silent_mutation_freq_plot
rm(kras_tcga_silent_mutation_data,kras_silent_mutation_freq_plot,kras_silent_mutation_freq_plot_df)
```

![KRAS silent mutation frequency plot](/plots/KRAS_silent_mutation_frequency_plot.png)

Based on the plot, we can see that "UCEC" (Uterine Corpus Endometrial Carcinoma) has the highest frequency of KRAS silent mutations.

Let's now look at the amino acid change distribution within KRAS 

## KRAS Amino Acid Change Distribution

```r
#KRAS mutations amino acid change summary
kras_AA_change_freq_plot_df <- kras_tcga_mutation_data %>% group_by(HGVSp_Short) %>% 
                                                           dplyr::summarise(number_of_mutations = n()) %>% 
                                                           arrange(desc(number_of_mutations))
#removing ".p" from the amino acid short name
kras_AA_change_freq_plot_df$HGVSp_Short <- str_remove(kras_AA_change_freq_plot_df$HGVSp_Short, "p.")
```

Let us look at the distribution of KRAS amino acid changes.

```r
kras_AA_change_freq_plot_df <- kras_AA_change_freq_plot_df %>% filter(HGVSp_Short!=".")

#Lets us now create a basic bar plot with x axis for number of amino acid  changes and y axis as the the amio acid change within KRAS  
kras_AA_change_freq_plot <- ggplot(kras_AA_change_freq_plot_df, 
                                   aes(x=reorder(HGVSp_Short,(number_of_mutations)), fill=HGVSp_Short,group=HGVSp_Short))+
                            geom_bar(aes(y=number_of_mutations),stat = "identity",position = "dodge")+
                            coord_flip()
                            
# adding title,x and y labels as well as modifying the theme aesthetics of the plot
kras_AA_change_freq_plot <- kras_AA_change_freq_plot + 
                            theme(plot.title = element_text(hjust = 0.5, color = "black"),
                                                            legend.position = "none") +
                            theme_bw() +
                            theme(axis.text.x = element_text(size = 10,face = "bold"), 
                                  axis.title.x = element_text(size = 16,face = "bold"),
                                  axis.text.y = element_text(size = 10,face = "bold"), 
                                  axis.title.y = element_text(size = 16,face = "bold"),
                                  plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
                                  plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
                                  legend.text=element_text(size = 14,face = "bold"),
                                  panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                  legend.position="none",
                                  legend.title = element_blank())+
                            labs(x = "Amino Acid Change",
                                 y = "Number Of Mutations", 
                                 subtitle = "KRAS",
                                 title = "Amino acid change")

kras_AA_change_freq_plot
rm(kras_AA_change_freq_plot)
```
![KRAS amino acid change plot](/plots/KRAS_AA_change_plot.png)

Based on the plot, we can see that "G12D" is the most prominent among KRAS mutations.

Now lets see which cancer has the most frequent G12D KRAS mutations.

## G12D Mutation Analysis

```r 
#Filtering the mutation data ONLY for KRAS G12D mutations
kras_AA_G12D_mutation_data<-kras_tcga_mutation_data %>% filter(HGVSp_Short=="p.G12D")

#removing ".p" from the amino acid short name
kras_AA_G12D_mutation_data$HGVSp_Short <- str_remove(kras_AA_G12D_mutation_data$HGVSp_Short, "p.")

#Filtering the clinical data ONLY for patients with KRAS G12D mutations
KRAS_clinical_AA_G12D_mutation<-filter(KRAS_clinical,bcr_patient_barcode %in% kras_AA_G12D_mutation_data$id)
```

Let us look at the distribution of G12D KRAS mutations across cancer types.\
We will look at it in two ways, 
-	Raw numbers
- Percent 

```r
#frequency of KRAS silent mutations
KRAS_clinical_AA_G12D_mutation_freq_plot_df <- KRAS_clinical_AA_G12D_mutation %>% group_by(acronym) %>% 
                                                                                  dplyr::summarise(number_of_mutations = n())%>% 
                                                                                  arrange(desc(number_of_mutations))
                                                                                   
#Create a basic bar plot x axis as the cancer type and y axis for number of KRAS G12D mutations
#The number above each bar denotes the number of KRAS G12D mutations
kras_AA_G12D_mutation_count_plot <- ggplot(KRAS_clinical_AA_G12D_mutation_freq_plot_df, 
                                           aes(x=reorder(acronym,-number_of_mutations),y=number_of_mutations, fill=acronym))+
                                      geom_bar(stat = "identity",position = "dodge")+
                                      geom_text(aes(label = number_of_mutations), fontface = "bold",vjust = -0.25)

# adding title,x and y labels as well as modifying the theme aesthetics of the plot
kras_AA_G12D_mutation_count_plot <- kras_AA_G12D_mutation_count_plot + 
                                    theme_bw() +
                                    theme(axis.text.x = element_text(size = 10,face = "bold"), axis.title.x = element_text(size = 16,face = "bold"),
                                          axis.text.y = element_text(size = 10,face = "bold"), axis.title.y = element_text(size = 16,face = "bold"),
                                          plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
                                          plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
                                          legend.text=element_text(size = 14,face = "bold"),
                                          panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                          legend.position="none",
                                          legend.title = element_blank())+
                                    labs(x = "Cancer Type",
                                         y = "Number Of Mutations", 
                                         subtitle = "KRAS",
                                         title = "Amino acid change (G12D) count")
kras_AA_G12D_mutation_count_plot
```
![KRAS amino acid G12D mutation count plot](/plots/KRAS_AA_G12D_mutation_count.png)

From the plot we can see that the G12D mutation is more prevalent in COAD(Colon Adenocarcinoma) and PAAD (Pancreatic Adenocarcinoma).

Now lets look at the percent plot.
```r

#Create a basic bar plot x axis as the cancer type and y axis for number of KRAS G12D mutations
#The number above each bar denotes the percent of KRAS G12D mutations
#Percent plot
kras_AA_G12D_mutation_percent_plot <- ggplot(KRAS_clinical_AA_G12D_mutation,aes(x = factor(acronym),fill=acronym)) +
                                        geom_bar(aes(y = after_stat(count/sum(count))))+
                                        geom_text(aes(y = after_stat(count/sum(count)), 
                                                      label = scales::percent(after_stat(count)/sum(count))),
                                                      fontface = "bold",
                                                      stat = "count", vjust = -0.25) +
                                        scale_y_continuous(labels = scales::percent_format())

# adding title,x and y labels as well as modifying the theme aesthetics of the plot
kras_AA_G12D_mutation_percent_plot <- kras_AA_G12D_mutation_percent_plot + 
                                      theme_bw() +
                                      theme(axis.text.x = element_text(size = 10,face = "bold"), axis.title.x = element_text(size = 16,face = "bold"),
                                            axis.text.y = element_text(size = 10,face = "bold"), axis.title.y = element_text(size = 16,face = "bold"),
                                            plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
                                            plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
                                            legend.text=element_text(size = 14,face = "bold"),
                                            panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position="none",
                                            legend.title = element_blank())+
                                      labs(x = "Cancer Type",
                                            y = "", 
                                           subtitle = "KRAS",
                                           title = "Amino acid change(G12D) percent")
                                      
kras_AA_G12D_mutation_percent_plot

rm(kras_AA_G12D_mutation_percent_plot,kras_AA_G12D_mutation_count_plot,KRAS_clinical_AA_G12D_mutation,KRAS_clinical_AA_G12D_mutation_freq_plot_df)
```
![KRAS amino acid G12D mutation percent plot](/plots/KRAS_AA_G12D_mutation_percent_plot.png)

## Clinical Correlation in Pancreatic Cancer PAAD

Now lets focus on the pancreatic adenocarcinoma (PAAD) KRAS mutated samples.

Lets look at the distribution of Age at initial pathologic diagnosis and compare it between the KRAS mutated and KRAS non mutated samples in PAAD.

Filtering the clinical data ONLY for PANCREATIC patients with KRAS mutations (KRAS_mutated)
```r
#Here I am filtering the KRAS mutated PAAD samples, adding the group as "KRAS_mutated" and setting the status to 1
KRAS_mutated_clinical_pancreatic<-KRAS_clinical %>% filter(acronym=="PAAD")
KRAS_mutated_clinical_pancreatic$group<-"KRAS_mutated"
KRAS_mutated_clinical_pancreatic$status<-1

Filtering the clinical data ONLY for PANCREATIC patients without KRAS mutations(KRAS_WT)
#Here I am filtering the KRAS non mutated PAAD samples, adding the group as "KRAS_WT" and setting the status to 0
'%!in%' <- function(x,y)!('%in%'(x,y))
KRAS_WT_clinical<-filter(tcga_clinical_data,bcr_patient_barcode %!in% kras_tcga_mutation_data$id)
KRAS_WT_clinical_pancreatic <-KRAS_WT_clinical %>% filter(acronym=="PAAD")
KRAS_WT_clinical_pancreatic$group<-"KRAS_WT"
KRAS_WT_clinical_pancreatic$status<-0

#gathering the KRAS_WT dataframe with age_at_initial_pathologic_diagnosis,group and status columns
KRAS_WT_age_df<-data.frame(KRAS_WT_clinical_pancreatic$group,
                           as.numeric(KRAS_WT_clinical_pancreatic$age_at_initial_pathologic_diagnosis),
                           KRAS_WT_clinical_pancreatic$status)
colnames(KRAS_WT_age_df)<-c("Group","Age_at_initial_pathologic_diagnosis","Status")

#Summary of KRAS_WT age at initial pathologic diagnosis
summary(KRAS_WT_age_df$Age_at_initial_pathologic_diagnosis)

#Min.   1st Qu.  Median   Mean   3rd Qu.  Max. 
#39.00   58.00   65.00   64.36   71.00   88.00

#gathering the KRAS_mutated dataframe with age_at_initial_pathologic_diagnosis group and status columns
KRAS_mutant_age_df<-data.frame(KRAS_mutated_clinical_pancreatic$group,
                               as.numeric(KRAS_mutated_clinical_pancreatic$age_at_initial_pathologic_diagnosis),
                               KRAS_mutated_clinical_pancreatic$status)
colnames(KRAS_mutant_age_df)<-c("Group","Age_at_initial_pathologic_diagnosis","Status")

#Summary of KRAS_mutated age at initial pathologic diagnosis
summary(KRAS_mutant_age_df$Age_at_initial_pathologic_diagnosis)

#Min.    1st Qu.  Median  Mean   3rd Qu.  Max. 
#35.00   56.25    66.00   65.15   74.00   85.00 
```
Lets now combine the KRAS_WT and KRAS_mutated in a dataframe and calculate the mean of age at initial pathologic diagnosis

```r
#Combining the KRAS_WT and KRAS_mutated age at initial pathologic diagnosis
KRAS_WT_mutant_age_df=rbind(KRAS_WT_age_df, KRAS_mutant_age_df)
age_mean <- ddply(KRAS_WT_mutant_age_df, "Group", summarise, 
                  grp.mean=mean(Age_at_initial_pathologic_diagnosis))
print(age_mean)

|     Group     | grp.mean    |
|     -----     |    ---      |
| KRAS_mutated  |   65.15254  |
| KRAS_WT       |   64.35821  |


#Density Plot of age at initial pathologic diagnosis for KRAS_WT and KRAS_mutated
ggplot(KRAS_WT_mutant_age_df, aes(Age_at_initial_pathologic_diagnosis, fill=Group)) + 
       geom_density(color = "black", alpha = 0.3)+
       geom_vline(data=age_mean, aes(xintercept=grp.mean, color=Group),
                  linetype="dashed")+
       theme_bw() +
       theme(axis.text.x = element_text(size = 10,face = "bold"), axis.title.x = element_text(size = 16,face = "bold"),
             axis.text.y = element_text(size = 10,face = "bold"), axis.title.y = element_text(size = 16,face = "bold"),
             plot.title = element_text(size = 20, face = "bold", color = "darkblue",hjust = 0.5),
             plot.subtitle = element_text(size = 16, face = "bold", color = "maroon",hjust = 0.5),
             legend.text=element_text(size = 14,face = "bold"),
             panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank())
```
![KRAS WT mutant age density plot](/plots/KRAS_WT_mutant_age_df_density_plot.png)

```r
#Violin Plot of age at initial pathologic diagnosis for KRAS_WT and KRAS_mutated
xlabs=  paste(levels(KRAS_WT_mutant_age_df$Group),(paste("N=",table(KRAS_WT_mutant_age_df$Group),sep="")),sep="\n")
KRAS_WT_mutant_age_df_violin<-ggplot(KRAS_WT_mutant_age_df, aes(x=Group, y=Age_at_initial_pathologic_diagnosis, fill=Group)) + 
                              geom_violin(alpha=0.4,trim=F)+
                              ylab("Age at initial pathologic diagnosis")+
                              scale_x_discrete(labels=xlabs)+
                              geom_boxplot(width=0.1)
#Customizing the violin
KRAS_WT_mutant_age_df_violin+
                             theme_bw() +
                             theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(size=15,,face = "bold"),
                                   axis.text.y = element_text(size=15,,face = "bold"),
                                   axis.text.x = element_text(size=15,,face = "bold"),
                                   legend.text=element_text(size = 14,face = "bold"),
                                   panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.title = element_blank())

```
![KRAS WT mutant age violin plot](/plots/KRAS_WT_mutant_age_df_violin_plot.png)

#### Wilcoxon rank-sum test to compare age at diagnosis between KRAS mutated and WT
wilcox_test_result <- wilcox.test(Age_at_initial_pathologic_diagnosis ~ Group,
                                  data = KRAS_WT_mutant_age_df,
                                  exact = FALSE)

print(wilcox_test_result)

##### Wilcoxon rank-sum test with continuity correction
##### data:  Age_at_initial_pathologic_diagnosis by Group
##### W = 4216.5, p-value = 0.4522
##### alternative hypothesis: true location shift is not equal to 0

The Wilcoxon rank-sum test was used to compare age at initial pathologic 
diagnosis between KRAS mutated and KRAS wild-type PAAD patients. 

- KRAS mutated: mean age 65.2 years (N=~150)
- KRAS wild-type: mean age 64.4 years (N=~30)
- Wilcoxon p-value: update after running

The result suggests that KRAS mutation status does not significantly influence age of onset in pancreatic adenocarcinoma.

## Advanced Somatic Mutation Analysis with maftools

Having explored KRAS-specific mutations and clinical correlations, we now use the `maftools` package for a broader view of the somatic mutation landscape across the TCGA pan-cancer cohort. 
`maftools` provides a suite of functions specifically designed for MAF file analysis and produces publication-ready visualizations.

