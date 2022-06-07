## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# download and install fdrci from github or CRAN
# https://uscbiostats.github.io/fdrci
# library(devtools)
# install_github("USCbiostats/fdrci")
library(fdrci)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("GEOquery")
library(GEOquery)
library(dplyr)
#library(stringr)
library(foreach)
library(ggplot2)

set.seed(052422)

## ---- message='hide', warning=FALSE-------------------------------------------
gset <- getGEO("GSE132342", GSEMatrix =TRUE, getGPL=TRUE)
gset = gset$GSE132342_series_matrix.txt.gz

pdat = phenoData(gset)
pdat.p = pData(pdat)
pdat.m = varMetadata(pdat)

fdat = featureData(gset)
fdat.n = featureNames(fdat)
fdat.m = varMetadata(fdat)
fdat.p = pData(fdat) # includes gene id's, target sequence, and 
#                      gene names (mostly HUGO gene symbols)

edat = exprs(gset) # normalized gene expression, genes in rows 
#                    samples in columns

rownames(edat) = fdat.p$Customer.Identifier

# Transform and combine phenotype data with expression data
pemat <- cbind(pdat.p, t(edat)) %>%
  mutate(time.os = as.numeric(`time.last.fu:ch1`),
         event.os = as.numeric(`status:ch1`),
         site = factor(`site:ch1`),
         stage = `Stage:ch1`,
         age = factor(`age:ch1`)) %>% 
  filter(stage %in% c(1, 2)) %>% 
  select(time.os, event.os, site, stage, age, fdat.p$Customer.Identifier[1]:fdat.p$Customer.Identifier[513])

## ---- message='hide', warning=FALSE-------------------------------------------
# create interaction term manually
pemat[, "ageStage"] = factor(paste(as.character(pemat[,"age"]), as.character(pemat[,"stage"]), sep="_"), levels=c("q4_2", "q1_1", "q1_2", "q2_1", "q2_2", "q3_1", "q3_2", "q4_1"))
# Run interaction scan
rslts <- foreach(j = 1:length(fdat.p$Customer.Identifier), .combine = rbind) %do% {
  fmla = paste(fdat.p$Customer.Identifier[j], "~ site + age + stage + ageStage")
  fit = lm(fmla, data=pemat)
  c(fdat.p$Customer.Identifier[j], anova(fit)["ageStage", "Pr(>F)"])
}
colnames(rslts) = c("gene", "p_value")
rslts = as.data.frame(rslts)
rslts[,"p_value"] = as.numeric(rslts[,"p_value"])
rslts.obs = rslts

## ---- message='hide', warning=FALSE-------------------------------------------
rslts.obs[,"q_value"] = p.adjust(rslts.obs[,"p_value"], method="BH")
knitr::kable(rslts.obs[rslts.obs$q_value < 0.1,], caption = "Parametic approach: BH FDR < 0.1", digits=4, row.names=FALSE)

## ---- message='hide', warning=FALSE-------------------------------------------
n.perm = 100
perml = vector('list', n.perm)
tmpvar = factor(paste(as.character(pemat[,"age"]), as.character(pemat[,"stage"]), sep="_"),
                levels=c("q4_2", "q1_1", "q1_2", "q2_1", "q2_2", "q3_1", "q3_2", "q4_1"))
for(perm in 1:n.perm){
  rslts <- foreach(j = 1:length(fdat.p$Customer.Identifier), .combine = rbind) %do% {
    pemat[, "ageStage"] = sample(tmpvar)
    fmla = paste(fdat.p$Customer.Identifier[j], "~ site + age + stage + ageStage")
    fit = lm(fmla, data=pemat)
    c(fdat.p$Customer.Identifier[j], anova(fit)["ageStage", "Pr(>F)"])
  }
  colnames(rslts) = c("gene", "p_value")
  rslts = as.data.frame(rslts)
  rslts[,"p_value"] = as.numeric(rslts[,"p_value"])
  perml[[perm]] = rslts
}

## ---- message='hide', warning=FALSE-------------------------------------------
mytbl = fdrTbl(rslts.obs$p_value,perml,"p_value",nrow(rslts.obs),2,5)
mytbl = mytbl[!is.na(mytbl$fdr), ]
knitr::kable(mytbl, caption = "FDR: no multiplicity adjustment",
             digits=3, row.names=FALSE)

## ---- message='hide', warning=FALSE-------------------------------------------
mytbl1 = fdrTbl(rslts.obs$p_value,perml,"p_value",nrow(rslts.obs),2,5,correct="BH")
mytbl1 = mytbl1[!is.na(mytbl1$fdr), ]
knitr::kable(mytbl1, caption = "FDR: with multiplicity adjustment",
             digits=3, row.names=FALSE)

## ---- message='hide', warning=FALSE-------------------------------------------
knitr::kable(rslts.obs[rslts.obs$p_value < 10^(-2.6),], caption = "Parametic approach: BH FDR < 0.2", digits=4, row.names=FALSE)

## ----fig.height=5, fig.width=7.5, message=FALSE, warning=FALSE----------------
p1 = ggplot( pemat, aes( x=stage, y=CAV1, fill=age) ) + 
  xlab("stage.age") + ylab("CAV1 expression") + theme_bw() +
  geom_violin(position=position_dodge(width = 0.8)) + 
  geom_jitter(aes(group=age), size=.02, position=position_jitterdodge(dodge.width=0.8, jitter.width=0.05)) + 
  labs(title="") + 
  geom_boxplot(width=0.1, notch=TRUE, outlier.shape = NA, position=position_dodge(width = 0.8)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
p1

## -----------------------------------------------------------------------------
sessionInfo()

