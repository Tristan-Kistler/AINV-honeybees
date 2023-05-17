################################################################################
# set up A-matrix, M-matrix and D-matrix. Compute AINV.
#
# Input is output of pedigree.r, see its log-file. 
# There are two possible assumptions for the relatedness of a sire's drone-
# producing queens in the base population. Either unrelated (equi=0) or that
# their additive genetic relationship is in equilibrium (equi=1). The desired
# value has to be specified within the programme.
################################################################################

write(paste((date()),'start AMD-AINV'),file="log-AMD-AINV.txt")

memory.limit(size=50000)

numbers<-read.table("numbers.txt")
NTOT<-numbers[1,6]
NS<-numbers[1,2]
ND<-numbers[1,3]

################################################################################
# Assumptions for the base population: specify equi (0 or 1)
# When 0, drone-producing queens in the base population are taken unrelated.
# When 1, the additive genetic relationship between these queens is taken in equilibrium
################################################################################
equi<-0

################################################################################
# If pinv=0 AINV is computed using (co)variances re Mendelian sampling terms.
# If pinv=1 AINV is computed directly from A as solve(A)
################################################################################
pinv<-0

################################################################################
# If priF=0 inbreeding coefficients of possible matings are not written
# if priF=1 they are written for the last yip years of birth of colonies
# These colonies produce virgin queens and drone-producing queens.
################################################################################
priF<-0
yip<-3
if(priF!=0 & priF!=1){
  write(paste((date()),'priF must be 0 or 1'),file="log-AMD-AINV.txt", append=TRUE)}
if(priF==0 | priF==1){
################################################################################


if(equi!=0 & equi!=1){
  write(paste((date()),'equi must be 0 or 1'),file="log-AMD-AINV.txt", append=TRUE)}
if(equi==0 | equi==1){
 
assb<-0
varcolbas<-0.25+1/(2*ND)+1/(4*NS)
##### if equi equals 1 (additive genetic relatioship between DPQ in equilibrium)
if(equi==1){
assb<-(ND*NS+2*NS+ND)/(3*ND*NS+NS+ND)
varcolbas<-assb
varsb<-(1+(NS-1)*assb)/NS}

if(equi==0){
write(paste('In base population sires are unrelated. NS = ',NS,'ND =',ND,'and the variance Acolony in the base population equals',varcolbas),file="log-AMD-AINV.txt", append=TRUE)}
if(equi==1){
write(paste('In base population the relation between sires is in equilibrium. NS = ',NS,'ND =',ND,'and the variance Acolony in the base population equals',varcolbas),file="log-AMD-AINV.txt", append=TRUE)}
################################################################################


if(pinv!=0 & pinv!=1){
  write(paste((date()),'pinv must be 0 or 1'),file="log-AMD-AINV.txt", append=TRUE)}
if(pinv==0 | pinv==1){

ID<-read.table("pedigree.txt",sep="",header=FALSE)
blocks<-read.table("blocks.txt",sep="",header=FALSE)
ident<-read.table("ident.txt",sep="", header=FALSE, stringsAsFactors = FALSE)
if(numbers[1,8]==1){
ss<-read.table("singlesires.txt",sep="",header=FALSE)}

nb<-numbers[1,7]

A<-matrix(ncol=NTOT, nrow=NTOT)
A[,]<-0
D<-matrix(nrow=NTOT,ncol=NTOT)
D[,]<-0
I<-diag(1,NTOT)
M<-matrix(nrow=NTOT, ncol=NTOT)
M[,]<-0                   
R<-rep(0,NTOT)

################################################################################
# Check whether always both parents are known or both unknown  
################################################################################
bothparents<-1
for(i in 1:NTOT){
if((ID[i,6]!=0 & ID[i,7]==0) | (ID[i,6]==0 & ID[i,7]!=0)){
write('Not in all cases both parents known or both unknown',file="log-AMD-AINV.txt", append=TRUE)
bothparents<-0
}}
# This if has its bracket at the end of the programme
if(bothparents==1){
################################################################################



################################################################################
# Set up A-matrix
################################################################################

for(i in 1:NTOT){
gc()
dam<-ID[i,6]
sire<-ID[i,7]
sex<-ID[i,2]
  for(j in 1:i){
  damj<-ID[j,6]
  sirej<-ID[j,7]
    if(i==j){
################################################################################
# DIAGONAL ELEMENTS
# both parents known
################################################################################
        if(dam>0 & sire>0){

# sires with multiple drone producing queens
if(sex==2){
######## ns and nd relate to the number of drone producing queens in this sire
######## and the number of drones produced by this sire
ns<-ID[i,19]
nd<-ID[i,20]
######## in R[i], however, nd and ns relate to the numbers in the sire of this sire
nss<-ID[dam,19]
nds<-ID[dam,20]              
R[i]=ID[i,8]<-(1+ID[dam,9])/4+1/(2*nds)+(1+ID[sire,9])/(4*nss)+(nds*nss-nds-nss)*R[sire]/(4*nds*nss)+A[dam,sire]/2

######## If the sire is open mating R[i] is taken to be zero. Open sires are identified by ns>50
if(ns>50){R[i]<-0}


######## check whether a sire is a single sires with NS=1
    if(numbers[1,8]==1){
          seqs<-which(ss[,1]==i)
          single<-length(seqs)
          if(single>0){
            ns<-1
            R[i]=ID[i,8]<-(1+ID[dam,9])/4+1/(2*nds)+(1+ID[sire,9])/(4*nss)+(nds*nss-nds-nss)*R[sire]/(4*nds*nss)+A[dam,sire]/2
          }
    }
        
################################################################################
A[i,i]<-(1+0.5*A[dam,sire]+(ns-1)*R[i])/ns                   
ID[i,9]<-A[dam,sire]/2
vardi<-(1-ID[dam,9])/4+(1-ID[sire,9])/4+(nss-1)*(1+ID[sire,9]-R[sire])/(4*nss)
covdfs<-(2-R[sire])/(4*nds)
ID[i,10]<-vardi/ns+(ns-1)*covdfs/ns
}

else{
# dams
  if(sex==1){ 
  A[i,i]<-1+0.5*A[dam,sire]
  ID[i,9]<-A[dam,sire]/2
  nss<-ID[dam,19]
  nds<-ID[dam,20]
  ID[i,10]<-(1-ID[dam,9])/4+(1-ID[sire,9])/4+(nss-1)*(1+ID[sire,9]-R[sire])/(4*nss)
  }
# colonies
  else{
# sex==3
########  ns and nd relate to the number of drone producing queens in the sire
########  of this colony and the number of drones produced by these
  ns<-ID[dam,19]
  nd<-ID[dam,20]
  A[i,i]<-(1+ID[dam,9])/4+1/(2*nd)+(1+ID[sire,9])/(4*ns)+(nd*ns-nd-ns)*R[sire]/(4*nd*ns)+A[dam,sire]/2
  ID[i,9]<-A[dam,sire]/2            
  ID[i,10]<-(2-R[sire])/(4*nd)
  }
}
############ bracket of dam>0 & sire>0: dam and sire known
        } 
        else{
################################################################################
# dam and sire unknown, so for sexes 1 and 2 only.
################################################################################ 
          if(equi==0){
          A[i,i]<-1
          ID[i,9]<-0
          ID[i,10]<-1
            if(sex==2){
            R[i]=ID[i,8]<-0
#######   ns relates to the number of drone producing queens in this sire
            ns<-ID[i,19]
#######   probably not practically, but theoretically this may be a single sire
            if(numbers[1,8]==1){
              seqs<-which(ss[,1]==i)
              if(length(seqs)>0){
                ns<-1
              }
            }
            A[i,i]<-1/ns
            ID[i,9]<-0
            ID[i,10]<-1/ns
            }
          }
              else{
              A[i,i]<-1
              ID[i,9]<-0
              ID[i,10]<-1
                if(sex==2){
                R[i]=ID[i,8]<-assb
                A[i,i]<-varsb
                ID[i,9]<-0
                ID[i,10]<-varsb
                }
              }
        }

# bracket of if(i==j): diagonal elements
    } 
      else{

################################################################################
# OFF-DIAGONALS
################################################################################
          if(sire>0 & dam>0){
# sire and dam known
# check if parents are full sibs
            if(damj==dam & sirej==sire){
######      this concerns the covariance between fullsibs, and nd
######      relates to the number of drones produced by the sire of
######      the full sibs
            ns<-ID[dam,19]
            nd<-ID[dam,20]

            A[i,j]<-(1+ID[dam,9])/4+1/(2*nd)+(1+ID[sire,9])/(4*ns)+(nd*ns-nd-ns)*R[sire]/(4*nd*ns)+A[dam,sire]/2
            A[j,i]<-A[i,j]
            } 
              else{
              A[i,j]<-0.5*(A[dam,j]+A[sire,j])
              A[j,i]<-A[i,j]
              }  
          }
            else{
# Both parents unknown
            A[i,j]<-0
            A[j,i]<-0
            }
      }

  }
}

  
write(paste((date()),'A ready'),file="log-AMD-AINV.txt",append=TRUE)







if(pinv==0){
write(paste((date()),'Compute AINV using mendelian sampling'),file="log-AMD-AINV.txt",append=TRUE)
  
################################################################################
# Produce D and M
################################################################################

for(i in 1:nb){
start<-blocks[i,1]
end<-blocks[i,2]
ma<-ID[start,6]
pa<-ID[start,7]
nd<-ID[ma,20]
ass<-ID[pa,8]
vars<-(2-ass)/(4*nd)
D[start:end,start:end]<-vars
}
                          
for(i in 1:NTOT){
D[i,i]<-ID[i,10]
ma<-ID[i,6]
pa<-ID[i,7]
if(ma>0){M[i,ma]<-0.5}
if(pa>0){M[i,pa]<-0.5}
}

################################################################################
# test eigenvalues
################################################################################
sneg<-0
for(i in 1:nb){
start<-blocks[i,1]
end<-blocks[i,2]
nbb<-end-start+1
e<-eigen(t(D[start:end,start:end]))$values[]
neg<-0
  for(ii in 1:nbb){
    if(e[ii]<=0){
    neg<-neg+1
    sneg<-sneg+1}}
      
      if(neg>0){
#     print(paste('non-positive eigenvalues in block ',i))
      write(paste('non-positive eigenvalues in block ',i),file="log-AMD-AINV.txt",append=TRUE)
        for(jj in start:end){
        write(D[jj,start:end],file="log-AMD-AINV.txt", append=TRUE)}
#          for(iii in start:end){
#          print(D[iii,jjj=start:end])
#          }
      }
}
if(sneg==0){
write('There are no non-positive eigenvalues in the blocks',file="log-AMD-AINV.txt", append=TRUE)}

write(paste((date()),'D and M ready'),file="log-AMD-AINV.txt",append=TRUE)

################################################################################
# AINV<-t(I-M)%*%solve(D)%*%(I-M)
# Derive D-inverse by inverting blocks
# Diagonal elements, correct outside blocks
################################################################################
# Size of blocks is set to 100 max. Change if necessary. 
################################################################################

# print(paste('Start setting up D-matrix',(date())))

AINV<-matrix(nrow=NTOT,ncol=NTOT)
AINV[,]<-0

DINV<-matrix(nrow=NTOT,ncol=NTOT)
DINV[,]<-0
block<-diag(1,100)
blocki<-matrix(nrow=NTOT,ncol=NTOT)
blocki[,]<-0

################################################################################
# diagonals for DINV in case of non-full-sibs
################################################################################
for(i in 1:NTOT){DINV[i,i]<-1/D[i,i]}
#######################################

for(i in 1:nb){
start<-blocks[i,1]
end<-blocks[i,2]
ifs<-end-start+1
block[1:ifs,1:ifs]<-D[start:end,start:end]
blocki<-solve(block)
DINV[start:end,start:end]<-blocki[1:ifs,1:ifs]
block<-diag(1,100)}

write(paste((date()),'DINV ready'),file="log-AMD-AINV.txt",append=TRUE)
################################################################################

AINV<-t(I-M)%*%DINV%*%(I-M)

write(paste((date()),'AINV ready'),file="log-AMD-AINV.txt",append=TRUE)



################################################################################
# Test if inverse AINV is Amin1 indeed
# Comparison of the product of AINV and A with I is the logical way to
# make this comparison because direct inversion of A is a problem in large
# datasets. If the product is not identity the method doesn't give insight
# where the problem lies. To identify the problem it helps to know whigh 
# elements are incorrect by comparison of the direct inverse of A with AINV
# in a smaller dataset 
I<-AINV%*%A
dd<-sum(abs(diag(I)-1)>10^-8)
tel<-sum(abs(I)>10^-8)-NTOT
if(dd>0){
for(i in 1:NTOT){
if(abs(I[i,i]-1)>10^-8){print(i)}
#  for(j in i:NTOT){
#    if(i!=j){
#      if(I[i,j]>10^-8){print(paste(i,j))}}
#  }
}}

write(paste('A x A-inverse differs for ',dd,' diagonals and for ',tel,' off-diagonals from I'),file="log-AMD-AINV.txt", append=TRUE)

# This bracket concerns AINV using (co)variances on Mendelian sampling terms
}

if(pinv==1){
  write(paste((date()),'Compute AINV by direct inverse'),file="log-AMD-AINV.txt",append=TRUE)
  AINV<-matrix(nrow=NTOT,ncol=NTOT) 
  AINV<-solve(A)
}

################################################################################
# write AINV
################################################################################
dim<-(NTOT+1)*NTOT/2
AINVP<-matrix(nrow=dim,ncol=3)
nprint<-0
for(i in 1:NTOT){
  for(j in 1:i){
    if(abs(AINV[i,j])>10^-20){
      nprint<-nprint+1
      AINVP[nprint,1]<-i
      AINVP[nprint,2]<-j
      AINVP[nprint,3]<-AINV[i,j]}
  }
}

AINVPP<-matrix(nrow=nprint,ncol=3)
for(i in 1:nprint){
  for(j in 1:3){
    AINVPP[i,j]<-AINVP[i,j]
  }
}

write.table(AINVPP,file="AINV.giv",eol="\r\n",row.names=FALSE,col.names=FALSE)

write(paste((date()),'AINV written'),file="log-AMD-AINV.txt",append=TRUE)

################################################################################
# Write diagonal elements of A and identifications to ID
################################################################################
for(i in 1:NTOT){
ID[i,12]<-A[i,i]
id1<-ID[i,16]
ID[i,1]<-ident[id1,1]}
################################################################################

write.table(ID,file="pedigree_complete.txt",eol="\r\n",row.names=FALSE,col.names=FALSE)

if(numbers[1,8]==1){
write.table(ss,file="singlesires_ass.txt",row.names=FALSE,col.names=FALSE)}

write(paste((date()),'pedigree_complete written'),file="log-AMD-AINV.txt",append=TRUE)


################################################################################
# Write inbreeding coefficients of combinations of last yip years.
# Written is column 5 the inbreeding coefficient of workers when mating queens reared from
# queen A with drones of queens reared from queen B and vv. Column 6 the inbreeding
# coefficient of workers when mating queens reared from queen A with drones from
# queens B. Column 7 is the reverse: Drones from queen B mated to queens reared
# from queen A.
# ################################################################################
if(priF==1){

ncol<-numbers[1,5]
fxy<-matrix(nrow=ncol*ncol,ncol=7)
fxy[,]<-0
naxy<-0
year<-ID[NTOT,18]-yip
for(i in 1:NTOT){
  if(ID[i,2]==3 & ID[i,18]>year){
    for(j in 1:i){
      if(ID[j,2]==3 & ID[j,18]>year){
        naxy<-naxy+1
        fxy[naxy,1]<-i
        fxy[naxy,2]<-j
        ii<-ID[i,6]
        jj<-ID[j,6]
        fxy[naxy,3]<-ID[ii,1]
        fxy[naxy,4]<-ID[jj,1]
        fxy[naxy,5]<-50*A[i,j]
			   if(i==j){
			   fxy[naxy,5]<-50*(1+ID[i,9])}
        fxy[naxy,6]<-50*A[i,jj]
        fxy[naxy,7]<-50*A[ii,j]
		  }
    }
  }
}

ffxy<-matrix(nrow=naxy,ncol=7)
for(i in 1:naxy){
for(j in 1:7){
ffxy[i,j]<-fxy[i,j]}}

write.table(ffxy,file="F-future.txt",eol="\r\n",row.names=FALSE,col.names=FALSE)
write(paste((date()),'F of potential matings written'),file="log-AMD-AINV.txt",append=TRUE)
}
#####################################################################################


write(paste((date()),'ready'),file="log-AMD-AINV.txt",append=TRUE)

################################################################################
# Bracket for both parents known or both unknown
}

# Bracket for pinv is not 0 or 1
}

# bracket for equi is not 0 or 1
}

# bracket for priF is not 0 or 1
}