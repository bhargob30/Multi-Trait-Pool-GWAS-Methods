# Two-dimensional Scheme

#(1) install and load the packages.

# Note that, it is necessary to have 'snpStats' 
# to install 'PhenotypeSimulator'.
# So, if you don't have snpStats, install it using
# the following command.


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats", force = T)
install.packages('PhenotypeSimulator')

# now load the packages. If you don't have any of the below mentioned packages, install it accordingly.
library(PhenotypeSimulator)
library(dplyr)
library(ggplot2)
library(readr)
library(stabs)
getshared=function(char){
  return( substring(char, 28, nchar(char))       )
}
getindependent=function(char){
  return( substring(char, 38, nchar(char))       )
}

K=4
L=K
probs_1=seq(0,1, by=1/K)
probs_2=seq(0,1, by=1/L)


getclusters_trait2=function(df){
  temp_list=list()
  temp_quantile=quantile(df$Trait_2, probs=probs_2)
  for(i in 1:(length(probs_2)-1)){
    if(i< (length(probs_2)-1)){
      temp_list[[i]]=df %>%  filter( (Trait_2>=temp_quantile[i]) & (Trait_2<temp_quantile[i+1]) )
    }
    else{
      temp_list[[i]]=df %>%  filter( (Trait_2>=temp_quantile[i]) & (Trait_2<=temp_quantile[i+1]) )
    }
  }
  return(temp_list)
}

cls_genotype=function(clusters, genotype){
  k=length(clusters)
  l=list()
  for(i in 1:k){
    l[[i]]= genotype %>% filter(row.names(genotype) %in% clusters[[i]]        ) 
    
  } 
  
  return(l)  
}
pool_allele_freq1=function(clusters_genotype, snp){
  k=length(clusters_genotype)
  a=rep(0,k)
  for(i in 1:k){
    a[i]=sum(clusters_genotype[[i]][,snp])/(2*(nrow(clusters_genotype[[i]]))) 
  }
  return(a) 
}


Nsim=100
A=matrix(NA, nrow=Nsim, ncol=12) # columns=c(#snps affecting only trait_1, #snps affecting only trait_2, TotalDis_Trait1, CausalDis_Trait1, PleiDis_Trait1, IndDisTrait_1, TotalDis_Trait2,CausalDis_Trait2, PleiDis_Trait2, IndDisTrait_2, TotalDisTrait1+2, CausalDisTrait1+2 )

time1=Sys.time()
for(iter in 1:Nsim){
  phenotype=runSimulation(N=2000, P=2, genVar=0.6,                                                                                                                                                                                                                                                                                                                                                                                                                                  
                          h2s=0.5/0.6, rho=0, delta=0, tNrSNP=10000,                                                                                                                                     
                          cNrSNP=10, SNPfrequencies=c(0.05, 0.1,0.2),
                          pIndependentGenetic = 0.8,
                          pTraitIndependentGenetic = 0.5, seed=iter+100 )  
  
  
  fixedGen=phenotype[["phenoComponentsIntermediate"]][["genFixed"]]
  snpscombined=as.data.frame(t(fixedGen$cov_effect))
  pleitropicsnps=rownames(snpscombined %>% filter(  (!(Trait_1 ==0)) & (!(Trait_2==0)) ))         
  snps_causal_trait1=rownames(snpscombined %>% filter(  (!(Trait_1 ==0)) & (Trait_2==0)  ))         
  snps_causal_trait2=rownames(snpscombined %>% filter(  (Trait_1==0)    & (!(Trait_2 ==0))  ))            
  A[iter, c(1,2)]=c(length(snps_causal_trait1), length(snps_causal_trait2))
  
  pleitropicsnps=sapply(pleitropicsnps, getshared)
  
  snps_causal_trait1=unname(sapply(snps_causal_trait1, getindependent))
  snps_causal_trait2=unname(sapply(snps_causal_trait2, getindependent))
  causal_snps=colnames(fixedGen[["cov"]])
  Y=as.data.frame(phenotype[["phenoComponentsFinal"]][["Y"]])
  quantile_1=quantile(Y$Trait_1, probs=probs_1)
  
  clusters_trait1=list()
  for(i in 1:(length(probs_1)-1)){
    if(i < (length(probs_1)-1)){
      clusters_trait1[[i]]=Y %>%  filter( (Trait_1>=quantile_1[i]) & (Trait_1<quantile_1[i+1]) )
    }
    else{
      clusters_trait1[[i]]=Y %>%  filter( (Trait_1>=quantile_1[i]) & (Trait_1<=quantile_1[i+1]) )  
    }
  }
  
  clusters_traits=list()
  for(i in 1:length(clusters_trait1)){
    temp=getclusters_trait2(clusters_trait1[[i]])
    clusters_traits=append(clusters_traits, temp)
  }
  
  clusters=lapply(clusters_traits, rownames)
  
  genotypes=as.data.frame(phenotype[["rawComponents"]][["genotypes"]][["genotypes"]])
  
  clusters_genotype=cls_genotype(clusters, genotypes)
  SNPs=phenotype[["setup"]][["id_snps"]]
  allele_freq1=matrix(NA, nrow=length(clusters), ncol=length(SNPs))
  for(i in 1:length(SNPs)){
    allele_freq1[,i]=pool_allele_freq1(clusters_genotype = clusters_genotype, snp=SNPs[i])  
  }
  colnames(allele_freq1)=SNPs
  ind_percluster=sapply(clusters_genotype, nrow)
  
  clusters_genotype_re=list()
  set.seed(300+iter)
  for(i in 1:length(clusters)){
    a=ind_percluster[i]
    y=rep(i,a)
    
    d=matrix(NA, nrow=a, ncol=length(SNPs))
    for(j in 1:length(SNPs)){
      d[,j]= rbinom(a,2, allele_freq1[i,j])
      
    }
    
    clusters_genotype_re[[i]]=cbind(y,d)
    
  }
  
  data_regenerate=clusters_genotype_re[[1]]
  for(i in 2:length(clusters_genotype_re)){
    data_regenerate=rbind(data_regenerate,clusters_genotype_re[[i]] )
  }
  
  data_regenerate=as.data.frame(data_regenerate)
  colnames(data_regenerate)[2:ncol(data_regenerate)]=SNPs
  
  
  clus_pheno=list()
  for(i in 1:length(clusters)){
    clus_pheno[[i]]=Y %>% filter(rownames(Y) %in% clusters[[i]])
  }
  
  clusters_genphen_re=list()
  for(i in 1:length(clusters_genotype_re)){
    clusters_genphen_re[[i]]=cbind( clus_pheno[[i]]$Trait_1, clus_pheno[[i]]$Trait_2,  clusters_genotype_re[[i]][, 2:ncol(data_regenerate)]         )
    
  }
  
  data_combined_re=clusters_genphen_re[[1]]
  for(i in 2:length(clusters_genphen_re)){
    data_combined_re=rbind(data_combined_re,clusters_genphen_re[[i]] )
  }
  data_combined_re=as.data.frame(data_combined_re)
  colnames(data_combined_re)=c('y_1','y_2', SNPs)
  
  n0=nrow(data_combined_re)
  data_combined_re=data_combined_re[sample(1:n0, size=n0),]
  x=data_combined_re[,3:ncol(data_combined_re)]
  y1=data_combined_re$y_1
  y2=data_combined_re$y_2
  stab.lasso=stabsel(x=x, y=y1, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")
  
  temp1=stab.lasso[["selected"]]    
  A[iter, c(3,4,5,6)]=c(length(temp1),sum(parse_number(causal_snps)%in%temp1),sum(parse_number(pleitropicsnps) %in% temp1), sum(parse_number(snps_causal_trait1) %in% temp1) )
  stab.lasso=stabsel(x=x, y=y2, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")   
  
  temp2=stab.lasso[["selected"]]    
  A[iter, c(7,8,9,10)]=c(length(temp2),sum(parse_number(causal_snps)%in%temp2),sum(parse_number(pleitropicsnps) %in% temp2), sum(parse_number(snps_causal_trait2) %in% temp2) )  
  
  temp=unique(c(unname(temp1), unname(temp2)) )
  A[iter, c(11,12)]   = c(length(temp), sum(parse_number(causal_snps)%in%temp) )
  
}

colnames(A)=c('NoSNPs_Trait1', 'NoSNPs_Trait2', 'TotalDis_Trait1', 'CausalDis_Trait1', 'PleiDis_Trait1', 'IndDisTrait_1', 'TotalDis_Trait2', 'CausalDis_Trait2', 'PleiDis_Trait2', 'IndDisTrait_2', 'TotalDis_Trait1+2', 'CausalDis_Trait1+2' )

time2=Sys.time()
write.csv(A, 'c16_MB_0.8_0.5.csv')
print(time2-time1)