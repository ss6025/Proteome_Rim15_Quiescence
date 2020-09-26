---
title: "PCA and DE analysis for proteins"
output: 
  html_document: 
    keep_md: true
---

## Chrunks

* [Load libraries]
* [Read in tables]
* [Proteomic data preparation]
  * [QC.1 - Intensity and normalized intensity]
    * [Data preparation]
    * [Box plot]
    * [Scatter plot]
  * [QC.2 - Ratio and normalized ratio]
    * [Data preparation]
    * [Box plot]
    * [Scatter plot]
  * [Combine ratios and intensity into one table]
* [PCA analysis for each data type]
* [T-test for DE at each timepoint]

## Load libraries


## Read in tables 


## Proteomic data preparation
* This code extract the information from proteinGroup file generated in MaxQuant and put them in a tidy format
  * intensity
    * Run Anders & Huber, 2010 method to get normalized intensity for each protein
  * ratio
  * normalized_ratio
  
### QC.1 - Intensity and normalized intensity
#### Data preparation

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers as id variables
```

```
## `summarise()` regrouping output by 'Protein IDs', 'Gene names', 'Protein names' (override with `.groups` argument)
```

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
```

```
## `summarise()` regrouping output by 'nutrient', 'samplingTime', 'replicate' (override with `.groups` argument)
```

#### Box plot
![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-4-2.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-4-3.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-4-4.png)<!-- -->

#### Scatter plot
![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-5-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-5-2.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

### QC.2 - Ratio and normalized ratio
#### Data preparation

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers as id variables
```

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers", "channel", "nutrient", "samplingTime", "replicate", "ratio")
```

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers")
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers as id variables
```

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers", "channel", "nutrient", "samplingTime", "replicate", "normed_ratio")
```

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers", "channel", "nutrient", "samplingTime", "replicate", "genotype")
```

#### Box plot
![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-7-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-7-2.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-7-3.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-7-4.png)<!-- -->

#### Scatter plot
![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-8-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-8-2.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

### Combine ratios and intensity into one table 

```
## Joining, by = c("Protein IDs", "Gene names", "Protein names", "Fasta headers", "nutrient", "samplingTime", "replicate", "genotype")
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers, nutrient, samplingTime, replicate, genotype as id variables
```

```
## [1] 551664     10
```

```
##   Protein IDs Gene names                                 Protein names
## 1      Q9P305       IGO2                   mRNA stability protein IGO2
## 2      Q99385       VCX1              Vacuolar calcium ion transporter
## 3      Q99383       HRP1  Nuclear polyadenylated RNA-binding protein 4
## 4      Q99297       ODC2    Mitochondrial 2-oxodicarboxylate carrier 2
## 5      Q99258       RIB3 3,4-dihydroxy-2-butanone 4-phosphate synthase
## 6      Q99257      MEX67                      mRNA export factor MEX67
##                                                                                                                                             Fasta headers
## 1                   sp|Q9P305|IGO2_YEAST mRNA stability protein IGO2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=IGO2 PE=1 SV=3
## 2              sp|Q99385|VCX1_YEAST Vacuolar calcium ion transporter OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=VCX1 PE=1 SV=1
## 3  sp|Q99383|HRP1_YEAST Nuclear polyadenylated RNA-binding protein 4 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=HRP1 PE=1 SV=1
## 4    sp|Q99297|ODC2_YEAST Mitochondrial 2-oxodicarboxylate carrier 2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=ODC2 PE=1 SV=1
## 5 sp|Q99258|RIB3_YEAST 3,4-dihydroxy-2-butanone 4-phosphate synthase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=RIB3 PE=1 SV=1
## 6                    sp|Q99257|MEX67_YEAST mRNA export factor MEX67 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=MEX67 PE=1 SV=1
##   nutrient samplingTime replicate genotype      data     value
## 1        C          T00        R1 spike-in intensity  35919000
## 2        C          T00        R1 spike-in intensity  27087000
## 3        C          T00        R1 spike-in intensity 306210000
## 4        C          T00        R1 spike-in intensity  82491000
## 5        C          T00        R1 spike-in intensity 157850000
## 6        C          T00        R1 spike-in intensity  85041000
```

## PCA analysis for each data type
![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-2.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-3.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-4.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-5.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-6.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-7.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-8.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-9.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-10.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-11.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-12.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-13.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-14.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-15.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-10-16.png)<!-- -->

## T-test for DE at each timepoint 

```
## Using Protein IDs, Gene names, Protein names, Fasta headers, nutrient, samplingTime, data as id variables
## Using Protein IDs, Gene names, Protein names, Fasta headers, nutrient, samplingTime, data as id variables
```

```
## Joining, by = c("Protein IDs", "Protein names", "Gene names", "Fasta headers", "nutrient", "samplingTime", "data")
```

![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-11-1.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```
## Using Protein IDs, Gene names, Protein names, Fasta headers, samplingTime, genotype, data as id variables
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers, samplingTime, genotype, data as id variables
```

```
## Joining, by = c("Protein IDs", "Protein names", "Gene names", "Fasta headers", "genotype", "samplingTime", "data")
```

![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-11-3.png)<!-- -->![](Data-preperation-and-QC_files/figure-html/unnamed-chunk-11-4.png)<!-- -->
