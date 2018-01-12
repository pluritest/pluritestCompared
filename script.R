require(lumi)
load("./working/NMETH-BC11016B-pluritestsub.rdata")
featureNames(H9targetArray)<-fData(H9targetArray)[,1] # required for R3.x compatibility

unzip(zipfile = "./working/forLumi.txt.zip",exdir = "./working/")
# definitions for plots and computations
WA09<-c("WA09_r1_HS","WA09_r2_HS","WA09_r1_OO","WA09_r1_TEL","WA09_r2_OO","WA09_r2_TEL","WA09_r1_AE","WA09_r2_AE")
UKSCB<-c( "Shef3_r1_OO","Shef3_r2_OO","Oxford2_r1_OO","Oxford2_r2_OO","I9_r1_OO","I9_r2_OO")
Murdoch<- c("MEL1_INS_GFP_MG3_r1_AE","MEL1_INS_GFP_MG4_r1_AE","hES3_Mixl_GFP_r1_AE","hES3_Mixl_GFP_r2_AE","RM3","RM31" ) 
Kyoto<- c( "Tig108_4f3_r1_HS","Tig108_4f3_r2_HS","iPS201B7_r1_HS", "iPS201B7_r2_HS","KhES.1_r1_HS","KhES.1_r2_HS" )
Wicell<-  c("WA14_r1_JB","WA14_r2_JB","DF19.9.11T_4_r1_JB","DF19.9.11T_4_r2_JB","IMR90C4_r1_TEL","IMR90C4_r2_TEL")

# load the modified pluritest script 
# functionality has been removed to allow R2.11 - R 3.x compatibility 
# different output handling. Original is available in original workspace
source('./pluritest_mod_file.R')

# load the data using patched version of original pluritest script.
# NormOnly provides only a normalized lumi object. 
working.lumi1<-pluritest_mod("./working/forLumi.txt",NormOnly = TRUE)

# Add a simple Correction Factor to the Pluritest Modelmatrices
# This is very simple version of factor based batch effect removal as in RUV or SVA

M<-exprs(working.lumi1)
T<-exprs(H9targetArray)
T<-T[rownames(W15),]
M<-M[match(rownames(W15),rownames(M)),WA09]
shifttemp<-apply(M,1,mean)-T
shifttemp[is.na(shifttemp)] <- 0

W15<-cbind(W15[],shifttemp)
W12<-cbind(W12,shifttemp)

# run pluritest with correction
# Normalization is performed with a Target, so no cross sample effects other than through correction factor
res<-pluritest_mod("./working/forLumi.txt")

