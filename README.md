# TCGA Somatic Mutation Explorer

![Language](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Data](https://img.shields.io/badge/Data-TCGA%20GDC-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)

Exploratory analysis of TCGA pan-cancer somatic mutation data using R.
Focuses on KRAS mutation frequency, amino acid changes, clinical correlations and advanced somatic mutation analysis using maftools across TCGA cancer types.

## Table of Contents
- [Data Sources](#data-sources)
- [KRAS Mutation Frequency Across Cancer Types](#kras-mutation-frequency-across-cancer-types)
- [KRAS Silent Mutation Frequency](#kras-silent-mutation-frequency)
- [KRAS Amino Acid Change Distribution](#kras-amino-acid-change-distribution)
- [G12D Mutation Analysis](#g12d-mutation-analysis)
- [Clinical Correlation in Pancreatic Cancer PAAD](#clinical-correlation-in-pancreatic-cancer-paad)
- [Advanced Somatic Mutation Analysis with maftools](#advanced-somatic-mutation-analysis-with-maftools)
---

## Key Findings

- **COAD dominates KRAS mutations** — Colon Adenocarcinoma (COAD) has the highest frequency of KRAS mutations across all TCGA cancer types, followed by Lung Adenocarcinoma (LUAD)
- **UCEC leads silent mutations** — Uterine Corpus Endometrial Carcinoma (UCEC) shows the highest frequency of KRAS silent mutations, suggesting a distinct mutational mechanism
- **G12D is the most prevalent KRAS variant** — Among all amino acid changes, G12D is the most common KRAS mutation in the TCGA pan-cancer cohort
- **G12D is enriched in COAD and PAAD** — The G12D hotspot mutation is most prevalent in Colon Adenocarcinoma and Pancreatic Adenocarcinoma, consistent with its known role as a driver mutation in these cancer types
- **KRAS dominates PAAD** — KRAS is mutated in over 90% of pancreatic adenocarcinoma samples, confirming its role as the primary oncogenic driver in PAAD
- **No significant age difference in PAAD** — Age at diagnosis is similar between KRAS mutated (mean 65.2y) and wild-type (mean 64.4y) patients (Wilcoxon p = 0.4522), suggesting KRAS mutation status does not influence age of onset
---

**Data:** TCGA MAF file + clinical data (publicly available via GDC)  
**Tools:** R, ggplot2, dplyr, tidyr, maftools, stringr, scales

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
                           as.numeric(KRAS_WT_clinical_pancreatic$days_to_death),
                           as.numeric(KRAS_WT_clinical_pancreatic$days_to_last_followup),
                           KRAS_WT_clinical_pancreatic$status,
                           KRAS_WT_clinical_pancreatic$vital_status)
colnames(KRAS_WT_age_df)<-c("Group","Age_at_initial_pathologic_diagnosis","Days_to_death",
                             "Days_to_last_followup","Status","Vital_status")

#Summary of KRAS_WT age at initial pathologic diagnosis
summary(KRAS_WT_age_df$Age_at_initial_pathologic_diagnosis)

#Min.   1st Qu.  Median   Mean   3rd Qu.  Max. 
#39.00   58.00   65.00   64.36   71.00   88.00

#gathering the KRAS_mutated dataframe with age_at_initial_pathologic_diagnosis group and status columns
KRAS_mutant_age_df<-data.frame(KRAS_mutated_clinical_pancreatic$group,
                               as.numeric(KRAS_mutated_clinical_pancreatic$age_at_initial_pathologic_diagnosis),
                               as.numeric(KRAS_mutated_clinical_pancreatic$days_to_death),
                               as.numeric(KRAS_mutated_clinical_pancreatic$days_to_last_followup),
                               KRAS_mutated_clinical_pancreatic$status,
                               KRAS_mutated_clinical_pancreatic$vital_status)
colnames(KRAS_mutant_age_df)<-c("Group","Age_at_initial_pathologic_diagnosis","Days_to_death",
                                "Days_to_last_followup","Status","Vital_status")

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
- Wilcoxon p-value: 0.4522 (not significant)

The result suggests that KRAS mutation status does not significantly influence age of onset in pancreatic adenocarcinoma.

### Survival Analysis : KRAS Mutation Status in PAAD

Having established that KRAS mutation status does not influence age at diagnosis, we now 
ask a more clinically relevant question: does KRAS mutation status affect overall survival 
in pancreatic adenocarcinoma? We use Kaplan-Meier survival curves and the log-rank test 
to compare survival between KRAS mutated and wild-type PAAD patients.

```r
library(survival)
library(survminer)

# Prepare survival data
# Days_to_death for deceased patients, Days_to_last_followup for censored
KRAS_WT_mutant_age_df$OS_days <- ifelse(!is.na(KRAS_WT_mutant_age_df$Days_to_death),
  as.numeric(KRAS_WT_mutant_age_df$Days_to_death),
  as.numeric(KRAS_WT_mutant_age_df$Days_to_last_followup)
)

# Convert to months for readability
KRAS_WT_mutant_age_df$OS_months <- KRAS_WT_mutant_age_df$OS_days / 30.44

KRAS_WT_mutant_age_df$OS_event <- ifelse(
  KRAS_WT_mutant_age_df$Vital_status == "Dead", 1, 0
)
# Create survival object
surv_object <- Surv(
  time = KRAS_WT_mutant_age_df$OS_months,
  event = KRAS_WT_mutant_age_df$OS_event
)

# Fit Kaplan-Meier curves by KRAS group
km_fit <- survfit(surv_object ~ Group, 
                  data = KRAS_WT_mutant_age_df)

# Summary
summary(km_fit)$table
```

| Group | N | Events | Median Survival (months) | 95% CI Lower | 95% CI Upper |
|-------|---|--------|--------------------------|--------------|--------------|
| KRAS Mutated | 118 | 75 | 17.0 | 15.4 | 20.8 |
| KRAS Wild-type | 67 | 25 | 43.8 | 20.6 | NA* |

*Upper confidence limit not reached — insufficient events in KRAS wild-type group

#### Kaplan-Meier Survival Curves

```r
# Plot Kaplan-Meier curves
km_plot <- ggsurvplot(
  km_fit,
  data = KRAS_WT_mutant_age_df,
  pval = TRUE,                    # show log-rank p-value
  pval.method = TRUE,             # show test name
  conf.int = TRUE,                # show confidence intervals
  risk.table = TRUE,              # show number at risk table
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",        # show median survival line
  ggtheme = theme_bw(),
  palette = c("#E7B800", "#2E9FDF"),
  title = "Overall Survival by KRAS Mutation Status in PAAD",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  legend.labs = c("KRAS Mutated", "KRAS Wild-type"),
  legend.title = "Group",
  font.main = c(14, "bold", "darkblue"),
  font.x = c(12, "bold"),
  font.y = c(12, "bold"),
  font.tickslab = c(10, "bold"),
  risk.table.fontsize = 3.5
)

print(km_plot)
```

![Kaplan-Meier Survival Plot](/plots/KRAS_survival_KM.png)

#### Median Survival

```r
# Extract median survival for each group
surv_table <- summary(km_fit)$table
print(surv_table[, c("median", "0.95LCL", "0.95UCL")])
```
| Group | Median Survival (months) | 95% CI Lower | 95% CI Upper |
|-------|--------------------------|--------------|--------------|
| KRAS Mutated | 17.0 | 15.4 | 20.8 |
| KRAS Wild-type | 43.8 | 20.6 | NA* |

*Upper confidence limit not reached — insufficient events in KRAS wild-type group

#### Log-rank Test

```r
# Log-rank test for difference in survival between groups
logrank_test <- survdiff(surv_object ~ Group, 
                         data = KRAS_WT_mutant_age_df)
print(logrank_test)

Call:
survdiff(formula = surv_object ~ Group, data = KRAS_WT_mutant_age_df)

                     N Observed Expected (O-E)^2/E (O-E)^2/V
Group=KRAS_mutated 118       75     57.1      5.59      13.3
Group=KRAS_WT       67       25     42.9      7.45      13.3

 Chisq= 13.3  on 1 degrees of freedom, p= 3e-04 

# Extract p-value
p_value <- 1 - pchisq(logrank_test$chisq, 
                       length(logrank_test$n) - 1)
cat("Log-rank p-value:", round(p_value, 4), "\n")

Log-rank p-value: 3e-04 
```

The Kaplan-Meier analysis reveals a striking and statistically significant difference 
in overall survival between KRAS mutated and wild-type PAAD patients.

- **KRAS mutated:** median survival 17.0 months (95% CI: 15.4 – 20.8)
- **KRAS wild-type:** median survival 43.8 months (95% CI: 20.6 – NA)
- **Log-rank p-value: 0.0003**

KRAS wild-type patients survive more than **2.5 times longer** than KRAS mutated patients. 
The upper confidence limit for the wild-type group could not be estimated, indicating that 
a substantial proportion of wild-type patients remained alive at last follow-up.

These results suggest that KRAS mutation status is a significant prognostic factor in 
pancreatic adenocarcinoma, with KRAS mutations associated with markedly worse overall 
survival. This is consistent with the known role of oncogenic KRAS in driving aggressive 
disease biology through constitutive activation of proliferative and anti-apoptotic signaling 
pathways.

*Note: the KRAS wild-type group is relatively small (N=67) compared to the mutated group 
(N=118), reflecting the near-universal prevalence of KRAS mutations in PAAD. Results should 
be interpreted with this imbalance in mind.*


## Advanced Somatic Mutation Analysis with maftools

Having explored KRAS-specific mutations and clinical correlations, we now use the `maftools` package for a broader view of the somatic mutation landscape across the TCGA pan-cancer cohort. 
`maftools` provides a suite of functions specifically designed for MAF file analysis and produces publication-ready visualizations.

```r
library(maftools)
#Lets match the IDs in TCGA data to the Clinical data first

tcga_mutation_data$Tumor_Sample_Barcode <- substring(tcga_mutation_data$Tumor_Sample_Barcode, 1, 12)

# Verify it looks right
head(tcga_mutation_data$Tumor_Sample_Barcode)
# Should now show: "TCGA-02-0003"

# Rename clinical ID column to match
colnames(tcga_clinical_data)[colnames(tcga_clinical_data) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"


# Read the MAF file into maftools
# Clinical data is passed directly to annotate the MAF object
tcga_maf <- read.maf(maf = tcga_mutation_data,
                    clinicalData = tcga_clinical_data,verbose = FALSE)

# Quick summary of the MAF object
tcga_maf
```

### 1. Pan-cancer Mutation Summary

The summary plot provides an at-a-glance overview of the entire dataset:  
variant classifications, variant types, SNV classes, variants per sample, and the top 10 most frequently mutated genes.

```r

plotmafSummary(maf = tcga_maf,rmOutlier = TRUE,
               addStat = "median",dashboard = TRUE,
               titvRaw = FALSE)

```

![MAF Summary Plot](/plots/maf_summary_plot.png)

The summary reveals the overall mutational landscape across TCGA cancer types. 
Missense mutations dominate the variant classification, consistent with known cancer mutation patterns.  
TP53 and KRAS appear among the top frequently mutated genes pan-cancer. 

---

### 2. KRAS Lollipop Plot

The lollipop plot maps somatic mutations onto the KRAS protein structure, providing a visual 
summary of where mutations cluster along the protein domain. 
This is particularly informative for identifying hotspot residues.

```r

lollipopPlot(maf = tcga_maf,gene = "KRAS",
             AACol = "HGVSp_Short",showMutationRate = TRUE,
             showDomainLabel = FALSE,axisTextSize = c(1, 1))
             title(sub = "KRAS protein domains: H_N_Ras-like GTPase domain", 
                   cex.sub = 0.9)

```

![KRAS Lollipop Plot](/plots/KRAS_lollipop_plot.png)

The lollipop plot confirms that G12 (particularly G12D, G12V, G12C) is the dominant mutation hotspot in KRAS, concentrated in the GTPase domain.  
This is consistent with the known role of codon 12 mutations in disrupting GTP hydrolysis and locking KRAS in a constitutively active state.

---

### 3. Oncoplot — Top Mutated Genes in PAAD
The oncoplot reveals the somatic mutation landscape across the most frequently altered genes in pancreatic adenocarcinoma (PAAD). 
Each column represents a sample and each row a gene, with colors indicating mutation type.

```r
paad_maf <- subsetMaf(maf = tcga_maf, 
                      clinQuery = "acronym == 'PAAD'")

oncoplot(maf = paad_maf, top = 20, fontSize = 0.5,
         titleFontSize = 1,
         legendFontSize = 0.8,
         annotationFontSize = 0.8,
         SampleNamefontSize = 0.4 )

```

![Oncoplot](/plots/oncoplot_PAAD_top20.png)

The oncoplot reveals the somatic mutation landscape across 175 pancreatic adenocarcinoma (PAAD) samples. KRAS dominates as the most frequently mutated 
gene, altered in over 90% of samples, consistent with its established role as the primary oncogenic driver in pancreatic cancer. 
TP53 and SMAD4 appear as the next most frequently altered genes, reflecting the classical KRAS → TP53 → SMAD4 progression model of pancreatic tumorigenesis. 
The high frequency of co-mutation between KRAS and TP53 underscores the cooperative role of these alterations in driving aggressive disease. 

---

### 4. Somatic Interactions — Co-mutation Analysis

Somatic interaction analysis identifies gene pairs that are significantly co-mutated or mutually exclusive across samples. 
This is biologically important as it can reveal synthetic lethality relationships and pathway redundancies. 

```r

somaticInteractions(maf = tcga_maf,top = 25,
                    pvalue = c(0.05, 0.1))

```

![Somatic Interactions](/plots/somatic_interactions.png)

Significant co-mutations and mutual exclusivities are shown with p-values. 
Pairs shown in red are significantly co-mutated while blue pairs are mutually exclusive. 
These patterns can inform combination therapy strategies and help prioritise therapeutic targets.

---

### 5. Transition and Transversion Analysis

The transition/transversion (TiTv) plot summarises the ratio of transitions (A↔G, C↔T) to 
transversions (A↔C, A↔T, G↔C, G↔T) across samples. This ratio serves as a mutational 
fingerprint and can indicate exposure to specific mutagens such as UV radiation or tobacco smoke.

```r

titv_result <- titv(maf = tcga_maf,plot = TRUE,useSyn = TRUE)

```

![TiTv Plot](/plots/titv_plot.png)

The TiTv distribution across cancer types reflects known mutational processes. 
C>T transitions predominate across most cancer types, consistent with spontaneous deamination of methylated cytosines, a hallmark of aging-related mutagenesis. 
Cancer types with high C>A transversions (such as lung cancers) reflect tobacco smoke exposure.

---

### 6. Tumor Mutational Burden — TCGA Comparison

`tcgaCompare` compares the mutation load of our cohort against TCGA benchmark datasets, 
providing context for how the mutational burden in each cancer type compares to published results.

```r

tcga_maf_compare <- tcgaCompare(maf = tcga_maf,cohortName = "TCGA Pan-cancer",
                                logscale = TRUE,capture_size = 35.8,cohortFontSize=1,axisFontSize=1)

```

![TCGA Comparison](/plots/tcga_compare.png)

Tumor mutational burden (TMB) varies substantially across cancer types. 
Melanoma and lung cancers show the highest TMB, consistent with UV and tobacco-related mutagenesis respectively, 
while haematological malignancies and thyroid cancers show comparatively lower mutation loads.

