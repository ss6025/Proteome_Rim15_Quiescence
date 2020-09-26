---
title: "proteomics analysis (based on normalized ratio)"
author: "Siyu Sun"
date: "5/4/2020"
output: 
  html_document:
    keep_md: true
---

## Chrunks

  * [Load libraries]
  * [Read in pre-processed protein table]
* [Two-way ANCOVA]
  * [Genotype wise]
    * [Modeling]
    * [Summary of model fitting]
    * [DE visualization]
    * [Kegg analysis visualization]
  * [Nutrient wise]
    * [Modeling]
    * [Summary of model fitting]
    * [DE visualization]
    * [Kegg analysis visualization]
* [Three way ANCOVA]
  * [Modeling]
  * [Summary of model fitting]
  * [DE visualization]
  * [Kegg analysis visualization]
* [Single protein check]

## Load libraries

 
## Read in pre-processed protein table


# Two-way ANCOVA
## Genotype wise
### Modeling


### Summary of model fitting

```
## Parsed with column specification:
## cols(
##   nutrient = col_character(),
##   `Gene names` = col_character(),
##   `Protein IDs` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   estimate = col_double(),
##   std.error = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double()
## )
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```
## [1] 429  11
```

### DE visualization

```r
summary <- sig_genotype %>% 
  filter(data %in% c("normed_intensity","normed_ratio")) %>%
  group_by(nutrient, data, type) %>%
  summarize(DE = n())
```

```
## `summarise()` regrouping output by 'nutrient', 'data' (override with `.groups` argument)
```

```r
print(summary)
```

```
## # A tibble: 8 x 4
## # Groups:   nutrient, data [4]
##   nutrient data             type                    DE
##   <chr>    <chr>            <chr>                <int>
## 1 C        normed_intensity upregulated in rim15     1
## 2 C        normed_intensity upregulated in wt       15
## 3 C        normed_ratio     upregulated in rim15   161
## 4 C        normed_ratio     upregulated in wt       25
## 5 P        normed_intensity upregulated in rim15    14
## 6 P        normed_intensity upregulated in wt       37
## 7 P        normed_ratio     upregulated in rim15    25
## 8 P        normed_ratio     upregulated in wt       19
```

```r
listInput <- list(
  `upregulated in wt (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in wt (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

library("SuperExactTest")
```

```
## Loading required package: grid
```

```
## 
## Attaching package: 'SuperExactTest'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, union
```

```r
set.seed(1234)
n=1074

sig_geno = supertest(listInput, n=n)

summary <- summary(sig_geno)$Table 

ol_genes <- summary %>% filter(P.value < 0.05) %>% dplyr::select(Elements)
  
plot(sig_geno, Layout="landscape", sort.by="size", keep=FALSE,
	bar.split=c(70,180), elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig_geno)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
upset(fromList(listInput), order.by = "freq")
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
###normalized intensity
listInput <- list(
  `upregulated in wt (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in wt (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

sig_geno = supertest(listInput, n=n)

summary <- summary(sig_geno)$Table 

ol_genes <- summary %>% filter(P.value < 0.05) %>% dplyr::select(Elements)
  
plot(sig_geno, Layout="landscape", sort.by="size", keep=FALSE,
	bar.split=c(70,180), elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig_geno)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
upset(fromList(listInput), order.by = "freq")
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

### Kegg analysis visualization

```r
sig_genotype <- read_csv("./tables/anova_genotype_sig.csv")
```

```
## Parsed with column specification:
## cols(
##   nutrient = col_character(),
##   `Gene names` = col_character(),
##   `Protein IDs` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   estimate = col_double(),
##   std.error = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double(),
##   type = col_character()
## )
```

```r
#convert protein IDs column into ENSEMBL, which is readable for kegg annotation
tmp<-bitr(sig_genotype$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", drop = TRUE)
```

```
## Loading required package: org.Sc.sgd.db
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:SuperExactTest':
## 
##     intersect, union
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     rename
```

```
## The following object is masked from 'package:plotly':
## 
##     rename
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     slice
```

```
## The following object is masked from 'package:plotly':
## 
##     slice
```

```
## The following object is masked from 'package:purrr':
## 
##     reduce
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## 
## Attaching package: 'AnnotationDbi'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     select
```

```
## The following object is masked from 'package:plotly':
## 
##     select
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## 
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(sig_genotype$`Protein IDs`, "UNIPROT", "ENSEMBL",
## "org.Sc.sgd.db", : 4.27% of input gene IDs are fail to map...
```

```r
sig <- left_join (tmp, sig_genotype %>% dplyr::rename(UNIPROT = `Protein IDs`), by = "UNIPROT")

sig_C_upWT <- sig %>% 
  filter(data == "normed_ratio" & estimate > 0) %>% 
  dplyr::select(nutrient, ENSEMBL, UNIPROT)

upWT <- compareCluster(ENSEMBL ~ nutrient, 
                           data=sig_C_upWT, 
                           organism = "sce",
                           fun="enrichKEGG",
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

dotplot(upWT) + 
  ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom") 
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
#ggsave("./plots/go_up_C_kegg.pdf", width=10,height=8)

sig_upRim15 <- sig %>% 
  filter(data == "normed_ratio" & estimate < 0) %>% 
  dplyr::select(nutrient, ENSEMBL, UNIPROT)

upRim15 <- compareCluster(ENSEMBL ~ nutrient, 
                           data=sig_upRim15,
                           organism = "sce",
                           #OrgDb=org.Sc.sgd.db::org.Sc.sgd.db,
                           #ont="BP", 
                           fun="enrichKEGG", 
                           #keyType = "UNIPROT", 
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

dotplot(upRim15) + 
  ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom") 
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
#ggsave("./plots/go_up_P.pdf", width=10,height=8)
```

## Nutrient wise
### Modeling


### Summary of model fitting

```
## Parsed with column specification:
## cols(
##   genotype = col_character(),
##   `Gene names` = col_character(),
##   `Protein IDs` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   estimate = col_double(),
##   std.error = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double()
## )
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```
## [1] 1797   11
```

```
## # A tibble: 6 x 11
##   genotype `Gene names` `Protein IDs` data  term  estimate std.error statistic
##   <chr>    <chr>        <chr>         <chr> <chr>    <dbl>     <dbl>     <dbl>
## 1 rim15KO  AAP1         P37898        norm… samp… -5.96e+6   2.19e+6     -2.72
## 2 rim15KO  AAP1         P37898        norm… samp… -1.15e-2   3.95e-3     -2.91
## 3 rim15KO  AAT2         P23542        norm… samp… -4.37e-2   1.35e-2     -3.24
## 4 rim15KO  AAT2         P23542        ratio samp… -1.98e-2   8.31e-3     -2.38
## 5 rim15KO  ABF1         P14164        norm… samp…  1.96e+6   7.94e+5      2.47
## 6 rim15KO  ABF2         Q02486        norm… samp… -3.27e-2   1.35e-2     -2.41
## # … with 3 more variables: p.value <dbl>, qvalue <dbl>, type <chr>
```

### DE visualization
#### normalized ratio

```r
summary <- sig_nutrient %>% 
  filter(data %in% c("normed_intensity","normed_ratio")) %>%
  group_by(genotype, data, type) %>%
  summarize(DE = n())
```

```
## `summarise()` regrouping output by 'genotype', 'data' (override with `.groups` argument)
```

```r
print(summary)
```

```
## # A tibble: 8 x 4
## # Groups:   genotype, data [4]
##   genotype data             type                DE
##   <chr>    <chr>            <chr>            <int>
## 1 rim15KO  normed_intensity upregulated in C   179
## 2 rim15KO  normed_intensity upregulated in P   162
## 3 rim15KO  normed_ratio     upregulated in C     8
## 4 rim15KO  normed_ratio     upregulated in P   507
## 5 WT       normed_intensity upregulated in C   180
## 6 WT       normed_intensity upregulated in P   142
## 7 WT       normed_ratio     upregulated in C     8
## 8 WT       normed_ratio     upregulated in P   419
```

```r
listInput <- list(
  `upregulated in C (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in P (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in C (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in P (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

library("SuperExactTest")
set.seed(1234)
n=1074

sig = supertest(listInput, n=n)

summary <- summary(sig)$Table 

ol_genes <- summary %>% filter(P.value < 0.05)
#%>% select(Elements)
  
plot(sig, Layout="landscape", sort.by="size", keep=FALSE,	#bar.split=c(70,180), 
     elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
upset(fromList(listInput), order.by = "freq")
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

#### normalized intensity

```r
listInput <- list(
  `upregulated in C (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in P (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in C (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in P (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

sig = supertest(listInput, n=n)

summary <- summary(sig)$Table 

ol_genes <- summary %>% filter(P.value < 0.05)
#%>% select(Elements)
  
plot(sig, Layout="landscape", sort.by="size", keep=FALSE,	#bar.split=c(70,180), 
     elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
upset(fromList(listInput), order.by = "freq")
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

### Kegg analysis & visualization

```r
sig_nutrient <- read_csv("./tables/ancova_nutrient_sig.csv")
```

```
## Parsed with column specification:
## cols(
##   genotype = col_character(),
##   `Gene names` = col_character(),
##   `Protein IDs` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   estimate = col_double(),
##   std.error = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double(),
##   type = col_character()
## )
```

```r
#convert protein IDs column into ENSEMBL, which is readable for kegg annotation
tmp<-bitr(sig_nutrient$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", drop = TRUE)
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(sig_nutrient$`Protein IDs`, "UNIPROT", "ENSEMBL",
## "org.Sc.sgd.db", : 5.61% of input gene IDs are fail to map...
```

```r
sig <- left_join (tmp, sig_nutrient %>% dplyr::rename(UNIPROT = `Protein IDs`), by = "UNIPROT")

sig_nutrient_genelist_upC <- sig %>% 
  filter(data == "normed_ratio" & estimate > 0) %>% 
  dplyr::select(genotype, ENSEMBL, UNIPROT)

wt_rim15_upC <- compareCluster(ENSEMBL ~ genotype, 
                           data=sig_nutrient_genelist_upC, 
                           organism = "sce",
                           fun="enrichKEGG",
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

dotplot(wt_rim15_upC) + 
  ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom") 
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
#ggsave("./plots/go_up_C_kegg.pdf", width=10,height=8)

sig_nutrient_genelist_upP <- sig %>% 
  filter(data == "normed_ratio" & estimate < 0) %>% 
  dplyr::select(genotype, ENSEMBL, UNIPROT)

wt_rim15_upP <- compareCluster(ENSEMBL ~ genotype, 
                           data=sig_nutrient_genelist_upP,
                           organism = "sce",
                           #OrgDb=org.Sc.sgd.db::org.Sc.sgd.db,
                           #ont="BP", 
                           fun="enrichKEGG", 
                           #keyType = "UNIPROT", 
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

dotplot(wt_rim15_upP) + 
  ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom") 
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```r
#ggsave("./plots/go_up_P.pdf", width=10,height=8)
```


## Three way ANCOVA
### Modeling


### Summary of model fitting

```
## Parsed with column specification:
## cols(
##   `Protein IDs` = col_character(),
##   `Gene names` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   df = col_double(),
##   sumsq = col_double(),
##   meansq = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double()
## )
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-13-3.png)<!-- -->

```
## [1] 119  10
```

```
## # A tibble: 6 x 10
##   `Protein IDs` `Gene names` data  term     df    sumsq   meansq statistic
##   <chr>         <chr>        <chr> <chr> <dbl>    <dbl>    <dbl>     <dbl>
## 1 P00360        TDH1         norm… samp…     1 1.45e+20 1.45e+20      7.57
## 2 P00445        SOD1         norm… samp…     1 2.26e- 1 2.26e- 1     14.6 
## 3 P00447        SOD2         norm… samp…     1 8.42e+17 8.42e+17      7.87
## 4 P00447        SOD2         norm… samp…     1 2.26e+ 0 2.26e+ 0     12.9 
## 5 P00830        ATP2         norm… samp…     1 4.31e+ 0 4.31e+ 0      9.17
## 6 P01120        RAS2         norm… samp…     1 2.41e- 1 2.41e- 1     10.3 
## # … with 2 more variables: p.value <dbl>, qvalue <dbl>
```

### DE expression visualization

```r
sig_3way <- read_csv("./tables/anova_nutrient_sig.csv")
```

```
## Parsed with column specification:
## cols(
##   `Protein IDs` = col_character(),
##   `Gene names` = col_character(),
##   data = col_character(),
##   term = col_character(),
##   df = col_double(),
##   sumsq = col_double(),
##   meansq = col_double(),
##   statistic = col_double(),
##   p.value = col_double(),
##   qvalue = col_double()
## )
```

```r
sig_3way_proteins <- sig_3way %>% filter(data == "normed_ratio")
gene_list <- sig_3way_proteins$`Gene names`

nr_df <- master %>% 
    #remove spike-in in the data frame
    filter(genotype != "spike-in") %>%
    #remove samples from nitrogen starved condition
    filter(nutrient != "N") %>%
    #filter(replicate != "R2") %>%
    filter(data == "normed_ratio") %>%
    dplyr::select(`Protein IDs`, `Gene names`, genotype, nutrient, replicate, samplingTime, value) %>%
    filter(!is.na(value)) %>%
    filter(is.finite(value)) %>%
  filter(`Gene names` %in% sig_3way_proteins$`Gene names`) %>%
    #filter(value > 0 ) %>%
    unite(cond, nutrient, genotype, samplingTime, replicate) %>%
    spread(key = cond, value = value) 

annotation_col <- names(nr_df)[-c(1:2)] %>% as.tibble() %>%
    separate(value, c("nutrient","genotype","samplingTime","rep"), sep = "_", remove = FALSE) %>%
    column_to_rownames(var = "value")
```

```
## Warning: `as.tibble()` is deprecated as of tibble 2.0.0.
## Please use `as_tibble()` instead.
## The signature and semantics have changed, see `?as_tibble`.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_warnings()` to see where this warning was generated.
```

```r
  breaksList = seq(-0.07, 0.07, by = 0.0001)
#make heatmaps
  pheatmap(nr_df[,-c(1:2)] %>% drop_na(), 
           colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(length(breaksList)),
           cluster_cols = T, scale = "row", 
           annotation_col = annotation_col)
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

### Kegg analysis visualization

```r
#GO term analysis
ego_3way <- clusterProfiler::enrichGO(gene = sig_3way_proteins$`Protein IDs`, 
             ont ="ALL", 
             universe = unique(master$`Protein IDs`),
             keyType = "UNIPROT", 
             OrgDb = "org.Sc.sgd.db", 
             pAdjustMethod = "BH", 
             qvalueCutoff = 0.05)
 
dotplot(ego_3way, ~ONTOLOGY)
```

![](Proteom-analysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
write_csv(summary(ego_3way),"./tables/GOterm/ego_3way_ANCOVA.csv")
```

```
## Warning in summary(ego_3way): summary method to convert the object to data.frame
## is deprecated, please use as.data.frame instead.
```
# Single protein check 
## For DE proteins from three-way ancova


