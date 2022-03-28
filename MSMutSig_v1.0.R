############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:

####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.

####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.

####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.

#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

######################################################################################################
####### MSMutSig - Finding significantly mutated MS loci
######################################################################################################
####### For details on the implementation 
####### see YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors 
####### Nat. Biotechnology (2017) DOI: doi:10.1038/nbt.3966
#
# Written by Yosef E. Maruvka and Gad Getz
#
#############################################


#########
#
# Input
# A table of loci by sample.
# Each entry should be 0,1,2.
# 0 -> not covered 
# 1 -> covered.
# 2 -> has a mutation.
# One for the A motif and for the C motif




library(stargazer)
library(survival)
library(ggplot2)
library(GGally)
library(gtable)
library(survival)
library(ggplot2)

rm(list=ls())

file_A=""
file_C=""
Type="STAD";



###
# Params 
r_ratio_C=0.5;#Variablity variable for the negative binomial model for the C motif
r_ratio_A=0.5;#Variablity variable for the negative binomial model for the A motif

diversity=0.9
UTR=0;
motif=1
mutation_rate_range=0;


loci_inf=read.table("loci_number_of_alleles_plus_mutation_count_plus_mutation_type.stat.splt",row.names=1)


#### Uploading the data

motif="A";
data_t=read.table(file_A,header=TRUE,row.names=1)

fi=which(loci_inf[row.names(data_t),"V7"]>diversity)
data=data_t[fi,]


#### Mutation rate per locus
la<-data>0
cov_sample=colSums(la)
cov_loci=rowSums(la)
la<-data>1
mut_sample=colSums(la)
mut_loci=rowSums(la)


Rate<-data.frame(rate=sort(mut_loci/cov_loci),rank=c(length(sort(mut_loci/cov_loci)):1))
Rate[,"cov"]=cov_loci[row.names(Rate)]
Rate[,"mut"]=mut_loci[row.names(Rate)]
Rate[,"length"]=loci_inf[row.names(Rate),13]
Rate[,"type"]=loci_inf[row.names(Rate),7]
Rate[,"gene"]=loci_inf[row.names(Rate),8]
Rate[,"diversity"]=loci_inf[row.names(Rate),6]




#############
#  Rate estimation based on all loci
#
fi=which(Rate[, "rate"]>0)
Rate3=Rate[fi,]
fm1=aggregate(Rate[, "rate"], list(ceiling(Rate$length)), mean)
#p<-ggplot(fm1,aes(x=fm1[,1],y=fm1[,2]));
#p+geom_point()+scale_y_log10() + xlab("# of repeats") +ylab("mutation rate") + ggtitle(Type)

#########
# P-value based on all loci
#
Rate3=Rate[fi,];
Rate3["p_vale"]=1-pbinom(Rate3[,"mut"],Rate3[,"cov"],fm1[ceiling(Rate3[,"length"])-4,2])

so=sort(Rate3[,"p_vale"])

#############
#  Rate estimation based on variable loci only
#

fm1=aggregate(Rate[fi, "rate"], list(ceiling(Rate[fi,"length"])), mean)
fm2=aggregate(Rate[fi, "rate"], list(ceiling(Rate[fi,"length"])), sd)



#########
# P-value based on variable loci only
# Negative binomial model 
#
r_ratio=r_ratio_A;
for (i in seq(1,dim(Rate3)[1]))
{
  lambda=fm1[round(Rate3[i,"length"])-4,2]*Rate3[i,"cov"]
  r=r_ratio*lambda;
  p=lambda/(r+lambda)
  
  q=1-p;
  Rate3[i,"p_vale"]=1-pnbinom(Rate3[i,"mut"],r,q)
  Rate3[i,"mean"]=lambda;
  
}

############### The C motif
Rate3_copy=Rate3;


motif=2
loci_inf=read.table("loci_number_of_alleles_plus_mutation_count_plus_mutation_type.stat.splt",row.names=1)

#### Uploading the data
motif="C";
data_t=read.table(file_C,header=TRUE,row.names=1)

fi=which(loci_inf[row.names(data_t),"V7"]>diversity)
data=data_t[fi,]


#### Mutation rate per locus
la<-data>0
cov_sample=colSums(la)
cov_loci=rowSums(la)
la<-data>1
mut_sample=colSums(la)
mut_loci=rowSums(la)


Rate<-data.frame(rate=sort(mut_loci/cov_loci),rank=c(length(sort(mut_loci/cov_loci)):1))
Rate[,"cov"]=cov_loci[row.names(Rate)]
Rate[,"mut"]=mut_loci[row.names(Rate)]
Rate[,"length"]=loci_inf[row.names(Rate),13]
Rate[,"type"]=loci_inf[row.names(Rate),7]
Rate[,"gene"]=loci_inf[row.names(Rate),8]
Rate[,"diversity"]=loci_inf[row.names(Rate),6]



#############
#  Rate estimation based on all loci
#
fi=which(Rate[, "rate"]>0)
Rate3=Rate[fi,]
fm1=aggregate(Rate[, "rate"], list(ceiling(Rate$length)), mean)

#########
# P-value based on all loci
#
Rate3=Rate[fi,];
Rate3["p_vale"]=1-pbinom(Rate3[,"mut"],Rate3[,"cov"],fm1[ceiling(Rate3[,"length"])-4,2])

so=sort(Rate3[,"p_vale"])

#############
#  Rate estimation based on variable loci only
#

fm1=aggregate(Rate[fi, "rate"], list(ceiling(Rate[fi,"length"])), mean)
fm2=aggregate(Rate[fi, "rate"], list(ceiling(Rate[fi,"length"])), sd)



#########
# P-value based on variable loci only
#
# Negative binomial model 
#

r_ratio=r_ratio_C;
for (i in seq(1,dim(Rate3)[1]))
{
  lambda=fm1[round(Rate3[i,"length"])-4,2]*Rate3[i,"cov"]
  r=r_ratio*lambda;
  p=lambda/(r+lambda)
  
  q=1-p;
  Rate3[i,"p_vale"]=1-pnbinom(Rate3[i,"mut"],r,q)
  Rate3[i,"mean"]=lambda;
  
}


###############################################
#Combining the A's and C's motifs

Rate_all=rbind(Rate3,Rate3_copy)
Rate3=Rate_all;

so=sort(Rate3[,"p_vale"])


####Significant loci

#knonw_cancer_genes=read.table("/Users/maruvka/Google Drive/MSMuTect/Cancer_gene_list_Serana.csv",header = T,sep=",")
#ms_cancer_gene=intersect(knonw_cancer_genes[,"Gene"],Rate3[,"gene"])
#for(i in seq(1,dim(Rate3)[1]))
#{Rate3[i,"cancer_gene"]=is.element(Rate3[i,"gene"],knonw_cancer_genes[,"Gene"])}



Rate3["q_vale"]=Rate3[,"p_vale"]*di[1];
Rate3["q_vale"]=p.adjust(Rate3[,"p_vale"], method = "BH")
di=dim(Rate3)


rate=Rate3[fi,]
or=order(rate["p_vale"])
di=dim(Rate3[fi,])
print (di)

temp=Rate3[fi,];
temp2=Rate3;
temp[Type]=row.names(temp)
temp2[Type]=row.names(temp2)
write.table(temp2,paste0("Results_MS_all_variable_loci_",Type,"_",motif,"_.csv"),sep=",",row.names=FALSE)
if (dim(temp)[1]>0){
  write.table(temp,paste0("Results_MS_hotspots_",Type,"_",motif,"_.csv"),sep=",",row.names=FALSE)
}

#Exon loci
fi_exon=which(loci_inf[,7]=="Frame_Shift_Del"|loci_inf[,7]=="In_Frame_Del")


Rate4=Rate3[order(Rate3["p_vale"]),]

par(mar=c(1,1,1,1))

Rate4=Rate3[order(Rate3["p_vale"]),]

Rate4[,"P_expected"]=seq(1:dim(Rate4)[1])/(dim(Rate4)[1]+1) ;



fi_exon=which(Rate4[,6]=="Frame_Shift_Del"|Rate4[,6]=="In_Frame_Del")
Rate5=Rate4[fi_exon,]
Rate4[fi_exon,"q_vale"]=p.adjust(Rate4[fi_exon,"p_vale"], method = "BH")


df1 <- data.frame(x = seq(0,log10(dim(Rate4)[1])), y=seq(0,log10(dim(Rate4)[1])))
ggplot(Rate4,aes(x=-log10(P_expected),y=-log10(p_vale)))+
  geom_point() +
  geom_line(data=df1,aes(x=x,y=y),colour="blue")+
  geom_text(data=Rate4[fi_exon,],aes(label=ifelse(q_vale<0.1,as.character(paste0(gene,"  ")),'')),hjust=1,vjust=0,size = 2, label.padding = unit(1.25, "lines"))+
  geom_point(data=Rate4[which( (Rate4[,6]=="Frame_Shift_Del"|Rate4[,6]=="In_Frame_Del" )&Rate4$q_vale<0.1),],aes(x=-log10(P_expected),y=-log10(p_vale)),colour="red")+
  theme_bw()  + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=bquote('Expected - '  ~  log[10]~'(p−value)')   ,y=bquote('Observed - '  ~  log[10]~'(p−value)'),parse=T)+  
  ggsave(paste0("QQplot_",Type,"_",r_ratio_A,"_",r_ratio_C,".pdf"))  

