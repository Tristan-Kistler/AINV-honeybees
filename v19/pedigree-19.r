################################################################################
# 
# Produces pedigree file from input data of colonies with identifications
# of queens in colonies (=1a), their dam (=2a) the dam of the drone-producing
# queens mated to 1a (=4a) or the drone-producing queen(s) mated to 1a (=1b),
# and their years of birth, plus info about numbers of drone-producing queens
# and drones. Default values for NS and ND need to be provided in the programme.
# Also a value needs to be given for codetest. 1 = thorough test of the 
# sequence of the pedigree, less time consuming otherwise.
# 
# There is a description of input and output and uderpinning literature. 
# 
################################################################################

################################################################################         
# Set default values for NS and ND, and steering parameter codetest           # 
################################################################################
# Default values for NS and ND
NS<-100
ND<-10

# If codetest=1 a thorough sequence test is carried out, otherwise a less time
# consuming short test.
codetest<-0

################################################################################
write(paste('Logfile pedigree.r',(date())),file="log_pedigree.txt",append=FALSE)
write(paste('Default values in the programme for NS=',NS,'& ND=',ND),file="log_pedigree.txt",append=TRUE)
################################################################################

################################################################################

aa<-read.table("input-pedigree.txt",stringsAsFactors = FALSE,sep="",header=FALSE)
nrec<-length(aa[,1])

################################################################################
# Check whether for all records ns and nd are known
################################################################################
ccheck<-0
seqns<-which(aa[,14]==0)
seqnd<-which(aa[,15]==0)
if(length(seqns)>0 | length(seqnd>0)){ccheck<-1
  write('not all records in input-pedigree.txt have non-zero ns and or nd',file="log_pedigree.txt",append=TRUE)}
if(ccheck==0){

################################################################################
# Check whether queens with colonies only appear once
################################################################################
check<-matrix(nrow=nrec,ncol=18)
jcheck<-0
df<-as.data.frame(aa)
check<-df[order(df$V5),]
ch<-check[1,5]
for(i in 2:nrec){
  if(check[i,5]==ch){jcheck<-1
  write(paste(check[i,5],' occurs more than once in input-pedigree.txt'),file="log_pedigree.txt",append=TRUE)}
  ch<-check[i,5]}

if(jcheck==0){
################################################################################

numbers<-matrix(nrow=1,ncol=8)
numbers[,]<-0
numbers[1,1]<-nrec
numbers[1,2]<-NS
numbers[1,3]<-ND

ncol<-0

# singlesires<-matrix(nrow=100,ncol=3)
# singlesires[,]<-0
singlesires <- matrix(data=0, nrow=nrow(unique(aa[, c(17,18)])), ncol=3)
nsinglesires<-0
dsiredams<-matrix(nrow=nrec,ncol=3)
dsiredams[,]<-0
ndsiredams<-0
opensires<-c(1:nrec)
opensires[]<-0
nop<-0

################################################################################
# Check whether there are single sires with NS=1. Then set msinglesires=1, 
# otherwise it is 0.
################################################################################
msinglesires<-0
seqs<-which(aa[,14]==1)
if(length(seqs)>0){msinglesires<-1
write('There are sires in the dataset with only one drone-producing queen (NS=1)',file="log_pedigree.txt",append=TRUE)}
numbers[1,8]<-msinglesires

################################################################################
# Check whether all dams of single sires with NS=1 have a record
################################################################################
nodamss<-0
if(msinglesires==1){
  seqs<-which(aa[,14]==1)
  nseqs<-length(seqs)
  for(i in 1:nseqs){
  seqdam<-which(aa[,5]==aa[seqs[i],9])
    if(length(seqdam)==0){
    nodamss<-nodamss+1
    write(paste('The dam of single sire',aa[seqs[i],17],'should have a record.'),file="log_pedigree.txt",append=TRUE)}
    }
}

####### this if has its bracket at the end of the programme ####################
if(nodamss==0){

################################################################################
# Check whether all individuals have just one birth year
################################################################################
identyear<-matrix(nrow=3*nrec,ncol=3)
id<-0
identyear[,]<-0
for(i in 1:nrec){
  if(aa[i,5]%in%identyear[,1]==TRUE){
  seq<-which(identyear[,1]==aa[i,5])
    if(aa[i,4]!=identyear[seq,2]){identyear[seq,3]<-1}}
      else{
      id<-id+1
      identyear[id,1]<-aa[i,5]
      identyear[id,2]<-aa[i,4]}

  if(aa[i,7]!=0){
    if(aa[i,7]%in%identyear[,1]==TRUE){
    seq<-which(identyear[,1]==aa[i,7])
      if(aa[i,6]!=identyear[seq,2]){identyear[seq,3]<-1}}
        else{
        id<-id+1
        identyear[id,1]<-aa[i,7]
        identyear[id,2]<-aa[i,6]}
    }

  if(aa[i,9]!=0){
    if(aa[i,9]%in%identyear[,1]==TRUE){
    seq<-which(identyear[,1]==aa[i,9])
      if(aa[i,8]!=identyear[seq,2]){identyear[seq,3]<-1}}
        else{
        id<-id+1
        identyear[id,1]<-aa[i,9]
        identyear[id,2]<-aa[i,8]}
    }
}

if(id!=0){
  double<-which(identyear[,3]==1)
  jd<-length(double)
  if(jd>0){
    jdyear<-c(1:jd)
    for(i in 1:jd){
      jdyear[i]<-identyear[double[i],1]
    }
    
    write('There are individuals with more than one year of birth. See identyear.txt',file="log_pedigree.txt",append=TRUE)
    write.table(jdyear,file="identyear.txt",row.names=FALSE,col.names=FALSE)}
}

################################################################################
# file specific issues causing sequence problems:
# for example, suppose that the modified year of birth of the 2a with 
# identification ident needs to be one year earlier than in the input file 
# and changed to 2005
################################################################################
# seq<-which(aa[,5]=='4-9-174-2005')
# aa[seq,4]<-2005

################################################################################


################################################################################
# This brackets ends at the end of the programme(before the other one)
if(jd==0){
################################################################################
  
  
################################################################################
# Check whether each 4a has only one combination of NS and ND
# For that purpose, 4a is defined as the combi of column 9 and 18 that
# replaces column 18
# This does not concern records with known 1b
################################################################################
code<-c(1:3)

df<-as.data.frame(aa)
a<-df[order(df$V9,df$V18),]
for(i in 1:nrec){
code[1]<-a[i,9]
code[2]<-"-"
code[3]<-a[i,18]
if(a[i,9]!=0){a[i,18]<-paste(code, collapse="")}}

# find first record with known 4a and no 1b: istart
istart<-0
for(i in 1:nrec){
if(istart==0 & a[i,17]==0 & a[i,9]!=0){
istart<-i}}

# make list of 4a's, excluding parents of 1b that have an identification
# only if istart>0
if(istart>0){
  a4a<-matrix(nrow=nrec,ncol=3)
  a4a[,]<-0
  n4a<-1
  a4a[n4a,1]<-a[istart,18]
  a4a[n4a,2]<-a[istart,14]
  a4a[n4a,3]<-a[istart,15]
  for(i in (istart+1):nrec){
  if(a[i,17]==0){
    if(a[i,18]!=a4a[n4a,1]){
    n4a<-n4a+1
    a4a[n4a,1]<-a[i,18]
    a4a[n4a,2]<-a[i,14]
    a4a[n4a,3]<-a[i,15]}
      else{
        if((a4a[n4a,2]!=a[i,14]) | (a4a[n4a,3]!=a[i,15])){
        print(paste('unequal NS or ND for 4a',a[i,9],a4a[n4a,2],a[i,14],a4a[n4a,3],a[i,15]))}}
  }
}
b4a<-matrix(nrow=n4a,ncol=3)
for(i in 1:n4a){
for(j in 1:3){
b4a[i,j]<-a4a[i,j]}}
# write.table(b4a,file="4a.txt",row.names=FALSE,col.names=FALSE)
}
for(i in 1:nrec){
for(j in 1:18){
aa[i,j]<-a[i,j]}}
################################################################################

################################################################################
# sort data on 1a within year of birth 1a
################################################################################
df<-as.data.frame(aa)
a<-df[order(df$V4,df$V5),]
# write.table(a,file="input-sorted.txt",row.names=FALSE,col.names=FALSE)

# write(paste('number of colonies=',nrec),file="log_pedigree.txt")

################################################################################
# create file with all known individuals + colonies (if required) and dummy 
# sires. Numbers of colonies start with 300001 and of dummy sires 200001
################################################################################
col<-300000
sire<-200000
b<-matrix(nrow=4*nrec,ncol=20)
b[,]<-0
nped<-0
for(i in 1:nrec){
# this for ends after creation of a record for 4a
################################################################################
# A sire present? Create record for the sire. 
################################################################################
if(a[i,17]!=0){
mate<-a[i,17]
################################################################################
# if a[i,9]==0 it concerns open mating, where NS and ND are large
################################################################################
if(a[i,9]!=0){

################################################################################
# Singlesires (with NS=1) contains ident (=seq2) and ND. ND relates to the 
# number of drone-producing queens of that sire. Single sires can appear more 
# than once if they are used with different numbers of drones.
# For NS>1 and known 4a is assumed to represent the case where a 1b-ident is 
# fabricated and ND>1. This is a single-drone insemination from a mix of
# full sisters. A bit odd, but it happens.
################################################################################

  if(a[i,14]==1){
    if(a[i,17]%in%singlesires[,1]==FALSE){
    nsinglesires<-nsinglesires+1
    singlesires[nsinglesires,1]<-mate
    singlesires[nsinglesires,3]<-a[i,15]
    }
      else{
      seqq<-which(singlesires[,1]==a[i,17])
      msires<-length(seqq)
        if(msires==1){
        if(a[i,15]!=singlesires[seqq,3]){
        write(paste('single sire',a[i,17],'has more than one set of NS and ND. Check whether this is ok'),file="log_pedigree.txt",append=TRUE)}}
      nmsires<-1
        for(j in 1:msires){
        seqs<-seqq[j]
        if(singlesires[seqs,3]==a[i,15]){
        nmsires<-0}}
          if(nmsires>0){
          nsinglesires<-nsinglesires+1
          singlesires[nsinglesires,1]<-mate
          singlesires[nsinglesires,3]<-a[i,15]
          }
      }
  }

# this is a record for a single sire with known 4a, irrespective whether NS=1 or NS>1
if(a[i,17]%in%b[,1]==FALSE){
nped<-nped+1
b[nped,1]<-a[i,17]
b[nped,5]<-a[i,16]
b[nped,6]<-a[i,9]
b[nped,19]<-a[i,14]
b[nped,20]<-a[i,15]
b[nped,2]<-2
}
  else{
# this is the case where for this single sire already a queen-record exists. 
# If NS=1 the NS and ND of the mate of the queen should be retained, 
# otherwise replaced by those of the sire  
  seq<-which(b[,1]==a[i,17])
  b[seq,2]<-2
  if(a[i,14]>1){
  b[seq,19]<-a[i,14]
  b[seq,20]<-a[i,15]
  }
  }
}
# this is a record for a single sire that is open mating
  else{
  open<-which(opensires==a[i,17])
    if(length(open)==0){
    nop<-nop+1
    opensires[nop]<-a[i,17]
    nped<-nped+1
    b[nped,1]<-a[i,17]
    b[nped,2]<-2
    b[nped,5]<-a[i,16]
    b[nped,6]<-0
    b[nped,7]<-0
    b[nped,19]<-a[i,14]
    b[nped,20]<-a[i,15]}
  }
}
################################################################################
# if sire is not present, create record for the dummy sire, unless already done
################################################################################
if(a[i,17]==0 & a[i,18]!=0){
  if(a[i,18]%in%dsiredams[,1]==FALSE){
#
# 
  
  ndsiredams<-ndsiredams+1
  sire<-sire+1
  dsiredams[ndsiredams,1]<-a[i,18]
  dsiredams[ndsiredams,2]<-sire
  nped<-nped+1
  b[nped,1]<-sire
  mate<-b[nped,1]
  b[nped,2]<-2
# Because of sequence problems changed b[nped,5]<-a[i,8]+2 in +1 and back to +2
#  b[nped,5]<-a[i,8]+1
  b[nped,5]<-a[i,8]+2
  dsiredams[ndsiredams,3]<-b[nped,5]
  b[nped,6]<-a[i,9]
  b[nped,19]<-a[i,14]
  b[nped,20]<-a[i,15]}
     else{
     seq<-which(dsiredams[,1]==a[i,18])
     mate<-dsiredams[seq,2]}     
}
################################################################################
# create record for 1a (queen in the colony)
################################################################################
if(a[i,5]%in%b[,1]==FALSE){
nped<-nped+1
b[nped,1]<-a[i,5]
b[nped,2]<-1
if(a[i,9]!=0){b[nped,3]<-mate}
if(a[i,17]!=0 & a[i,9]==0){b[nped,3]<-mate}
b[nped,5]<-a[i,4]
b[nped,6]<-a[i,7]
b[nped,19]<-a[i,14]
b[nped,20]<-a[i,15]}

################################################################################
# create record for colony, if required
################################################################################
if(a[i,1]!=0){
nped<-nped+1
col<-col+1
ncol<-ncol+1
b[nped,1]<-col
b[nped,2]<-3
### Note: year of birth is year of birth of 1a
b[nped,5]<-a[i,4]
b[nped,6]<-a[i,5]
seq<-which(b[,1]==a[i,5])
b[nped,7]<-b[seq,3]

# this concerns the test location
# b[nped,11]<-a[i,3]
}
################################################################################
# create record for 2a (dam of queen in the colony)
################################################################################
if(a[i,7]!=0){
if(a[i,7]%in%b[,1]==FALSE){
nped<-nped+1
b[nped,1]<-a[i,7]
b[nped,2]<-1
b[nped,5]<-a[i,6]
# 2a is not in dataset as 1a and take default values for NS and ND
b[nped,19]<-NS
b[nped,20]<-ND
}}

################################################################################
# create record for 4a
################################################################################
if(a[i,9]!=0){
if(a[i,9]%in%b[,1]==FALSE){
nped<-nped+1
b[nped,1]<-a[i,9]
b[nped,2]<-1
b[nped,5]<-a[i,8]
# 4a is not in dataset as 1a and take default values for NS and ND
b[nped,19]<-NS
b[nped,20]<-ND
}}
}
################################################################################
# Above is the bracket of the loop going through all records in input-pedigree.txt
################################################################################
numbers[1,5]<-ncol

if(msinglesires==1){
ssinglesires<-matrix(nrow=nsinglesires,ncol=3)
if(nsinglesires>0){
for(i in 1:nsinglesires){
for(j in 1:3){
ssinglesires[i,j]<-singlesires[i,j]}}}
}

if(ndsiredams>0){
  siredams<-matrix(nrow=ndsiredams,ncol=3)
  for(i in 1:ndsiredams){
  for(j in 1:3){
  siredams[i,j]<-dsiredams[i,j]}}
}

################################################################################
# Check whether all animals have a unique number
################################################################################
check<-matrix(nrow=nped,ncol=20)
icheck<-0
bbb<-matrix(nrow=nped,ncol=20)
for(i in 1:nped){
for(j in 1:20){
bbb[i,j]<-b[i,j]}}
  df<-as.data.frame(bbb)
  check<-df[order(df$V1),]
  ch<-check[1,1]
  for(i in 2:nped){
  if(check[i,1]==ch){icheck<-1
  print(paste('ident',check[i,1],'occurs more than once, perhaps look at check-double-idents.txt'))}
  ch<-check[i,1]}

if(icheck==1){
write.table(b,file="check-double-idents.txt",row.names=FALSE,col.names=FALSE)}
  
################################################################################
# Replace idents by sequence numbers (sequence number 1)
################################################################################
c<-matrix(nrow=nped,ncol=20)
c[,]<-0

ident<-c(1:nped)
for(i in 1:nped){
ident[i]<-b[i,1]}

for(i in 1:nped){
c[i,1]<-i
  for(j in 2:5){
  c[i,j]<-b[i,j]}
  for(j in 8:20){
  c[i,j]<-b[i,j]}
  c[i,16]<-i
  c[i,18]<-b[i,5]
    if(b[i,3]!=0){
    seq<-which(ident[]==b[i,3])
    c[i,3]<-seq}
    if(b[i,6]!=0){
    seq<-which(ident[]==b[i,6])
    c[i,6]<-seq}
    if(b[i,7]!=0){
    seq<-which(ident[]==b[i,7])
    c[i,7]<-seq}
}

mode(c)<-'numeric'

################################################################################
# add sire by finding the mate of the dam and increase modified birthyear of
# colonies with 2
################################################################################

for(i in 1:nped){
if(c[i,2]==3){c[i,5]<-c[i,5]+2}
if(c[i,7]==0){
dam<-c[i,6]
if(dam>0){
c[i,7]<-c[dam,3]}}
}

# write.table(c,file="tussen-463.txt",row.names=FALSE,col.names=FALSE)

################################################################################
# To reach a situation where each individual has either two or zero parents
# base parents are created; additional dams 110001 and additional sires 210001
################################################################################

# make list of known 4a's to avoid that dams get different mates
    dams4a<-matrix(nrow=nrec,ncol=2)
    dams4a[,]<-0
    n4a<-0

nbdam<-0
nbsire<-0
for(i in 1:nped){
  if(c[i,6]==0 & c[i,7]>0){
  nbdam<-nbdam+1}

  if(c[i,6]>0 & c[i,7]==0){
## These three lines (and also 2 + 1 + 2 later) are needed because a dam can
## be mated to only one sire.
    if(c[i,6]%in%dams4a[,1]==FALSE){
    n4a<-n4a+1
    dams4a[n4a,1]<-c[i,6]
    nbsire<-nbsire+1
    }
  }
}

nbase<-nbdam+nbsire

base<-matrix(nrow=nbase,ncol=20)
identbase<-c(1:nbase)
base[,]<-0
identbase[]<-0

idam<-110000
isire<-210000
j<-0
seq1<-nped

for(i in 1:nped){
  if(c[i,6]==0 & c[i,7]>0){
  j<-j+1
  seq1<-seq1+1
  c[i,6]<-seq1
  idam<-idam+1
  identbase[j]<-idam
  base[j,1]<-seq1
  base[j,2]<-1
  base[j,5]<-c[i,5]-2
  base[j,16]<-seq1
  base[j,18]<-base[j,5]
  }

## These two lines are needed because a dam can be mated to only one sire.
  if(c[i,6]>0 & c[i,7]==0){
  seq<-which(c[i,6]==dams4a[,1])      

    if(dams4a[seq,2]==0){
    j<-j+1
    seq1<-seq1+1
    c[i,7]<-seq1

## These line is needed because a dam can be mated to only one sire.
dams4a[seq,2]<-seq1

    isire<-isire+1
    identbase[j]<-isire
    base[j,1]<-seq1
    base[j,2]<-2
# This is the modified birth year of an additional sire; make four years older 
# than progeny, instead of 2
# base[j,5]<-c[i,5]-4
    base[j,5]<-c[i,5]-2
    base[j,16]<-seq1
    base[j,18]<-base[j,5]
# check whether the dam has a record
    base[j,19]<-NS
    base[j,20]<-ND
    seqdam<-which(c[i,6]==c[,16])
    if(length(seqdam)>0){
    base[j,19]<-c[seqdam,19]
    base[j,20]<-c[seqdam,20]}
    }

## These two lines are needed because a dam can be mated to only one sire.
else{c[i,7]<-dams4a[seq,2]}
  }

}
################################################################################
# add base animals to pedigreefile
###############################################################################
nped<-nped+nbdam+nbsire
numbers[1,6]<-nped
d<-rbind(c,base)
identt<-c(ident,identbase)

################################################################################
#  Decrease modified birth year with one year for dams that have the same
#  birth year as offspring
#################################################################################
for(i in 1:nped){
dam<-d[i,6]
if(dam>0){
if(d[i,5]<d[dam,5]){
d[dam,5]<-d[i,5]-1}}}

################################################################################
# modify birthyear such that all animals in a full sib group have the same
# lowest birth year
################################################################################

df<-as.data.frame(d)
e<-df[order(df$V7,df$V6,df$V5),]

sire<-e[i,7]
dam<-e[i,6]
birth<-e[i,5]
nfs<-1
iblock<-0
for(i in 2:nped){
    if(e[i,7]>0){
      if(e[i,7]==sire & e[i,6]==dam){
      nfs<-nfs+1}
        else{
          if(nfs>1){
          iblock<-iblock+1
            for(j in (i-nfs):(i-1)){
            e[j,5]<-birth}
          }
        sire<-e[i,7]
        dam<-e[i,6]
        birth<-e[i,5]
        nfs<-1
        }
    }
}

if(nfs>1){
  iblock<-iblock+1
  for(j in (i-nfs+1):i){
  e[j,5]<-birth}
}

# write.table(e,file="tussen-610.txt",row.names=FALSE,col.names=FALSE)

numbers[1,7]<-iblock
################################################################################
# Sort modified birth year, sire, dam and create file blocks with full sib 
# groups
################################################################################
blocks<-matrix(nrow=iblock,ncol=2)

df<-as.data.frame(e)
f<-df[order(df$V5,df$V7,df$V6),]
sire<-f[1,7]
dam<-f[1,6]
nfs<-1
iblock<-0
for(i in 2:nped){
  if(f[i,7]==sire & f[i,6]==dam & f[i,7]!=0){
  nfs<-nfs+1}
    else{
      if(nfs>1){
      iblock<-iblock+1
      blocks[iblock,1]<-i-nfs
      blocks[iblock,2]<-i-1}
    sire<-f[i,7]
    dam<-f[i,6]
    nfs<-1
    }
}

if(nfs>1){
iblock<-iblock+1
blocks[iblock,1]<-i-nfs
blocks[iblock,2]<-i-1}

################################################################################
# New sequence numbers in final order of the pedigree.txt
################################################################################
for(i in 1:nped){
f[i,17]<-i}

df<-as.data.frame(f)
g<-df[order(df$V16),]
h<-matrix(nrow=nped,ncol=20)
for(i in 1:nped){
for(j in 1:20){
h[i,j]<-g[i,j]
  if(g[i,3]>0){
  mate<-g[i,3]
  h[i,3]<-g[mate,17]}
    if(g[i,6]>0){
    dam<-g[i,6]
    h[i,6]<-g[dam,17]}
  if(g[i,7]>0){
  sire<-g[i,7]
  h[i,7]<-g[sire,17]}
}}
df<-as.data.frame(h)
hh<-df[order(df$V17),]

################################################################################
# Write pedigree file with identifications in stead of sequence numbers.
# Identifications for individuals in the input file are those from the 
# input files, but for colonies, sires identified by their dam, and
# base animals the programme fabricates identifications
################################################################################

ide<-matrix(nrow=nped,ncol=20)
for(i in 1:nped){
  for(j in 1:20){
  ide[i,j]<-hh[i,j]}
seq1<-hh[i,16]
ide[i,1]<-identt[seq1]
dam<-hh[i,6]
  if(dam!=0){
  seq1<-hh[dam,16]
  ide[i,6]<-identt[seq1]}
    if(hh[i,2]==1 & hh[i,3]!=0){
    mate<-hh[i,3]
    seq1<-hh[mate,16]
    ide[i,3]<-identt[seq1]
    matedam<-hh[mate,6]
      if(matedam!=0){
      seq1<-hh[matedam,16]
      ide[i,4]<-identt[seq1]}
    }
sire<-hh[i,7]
  if(sire!=0){
  seq1<-hh[sire,16]
  ide[i,7]<-identt[seq1]
  siredam<-hh[sire,6]
    if(siredam!=0){
    seq1<-hh[siredam,16]
    ide[i,8]<-identt[seq1]}
  }
}
write.table(ide,file="pedigree_ident.txt",row.names=FALSE,col.names=FALSE)

################################################################################
# check whether all parents are earlier in the data than offspring             #
################################################################################
# print('check whether all parents are earlier in the data than offspring')
itel<-0
write('        Check sequence issues in the pedigree file',file="log_pedigree.txt",append=TRUE)
write('        Listed are sequence number 2, identification, sex and year of birth of individuals and dams, sires or aunts',file="log_pedigree.txt",append=TRUE)
write('        Check whether all parents are earlier in the data than offspring',file="log_pedigree.txt",append=TRUE)

for(i in 1:nped){
  if(hh[i,17]<hh[i,6]){
  itel<-itel+1
  dam<-hh[i,6]
  ident1<-identt[hh[i,1]]
  ident2<-identt[hh[dam,1]]
  write(paste('ident= ',hh[i,17],ident1,hh[i,2],hh[i,5],'dam= ',hh[dam,17],ident2,hh[dam,2],hh[dam,5]),file="log_pedigree.txt",append=TRUE)
  }
    if(hh[i,17]<hh[i,7]){
    itel<-itel+1
    dad<-hh[i,7]
    ident1<-identt[hh[i,1]]
    ident2<-identt[hh[dad,1]]
    write(paste('ident= ',hh[i,17],ident1,hh[i,2],hh[i,5],'dad= ',hh[dad,17],ident2,hh[dad,2],hh[dad,5]),file="log_pedigree.txt",append=TRUE)
    }
}

if(codetest!=1){write("codetest is not 1 and therefore thorough testing of the sequence in the pedigree is suppressed",file="log_pedigree.txt",append=TRUE )}
if(codetest==1){

################################################################################
#
# check whether full aunts (uncles) are earlier in the data than offspring
################################################################################
# print('Check whether full aunts are earlier in the data than offspring')
  write('        Check whether full aunts are earlier in the data than offspring',file="log_pedigree.txt",append=TRUE)
  for(ind in 1:(nped-2)){
  ma<-hh[ind,6]
  pa<-hh[ind,7]
    if(ma>0 & pa>0){
    mama<-hh[ma,6]
    pama<-hh[ma,7]
      if(mama>0 & pama>0){
        for(iaunt in (ind+1):nped){
        maaunt<-hh[iaunt,6]
        paaunt<-hh[iaunt,7]
          if(mama==maaunt & pama==paaunt){
          itel<-itel+1
          ident1<-identt[hh[ind,1]]
          ident2<-identt[hh[iaunt,1]]
          write(paste('ident= ',hh[ind,17],ident1,hh[ind,2],hh[ind,5],'dam= ',hh[ind,6],'aunt= ',hh[iaunt,17],ident2,hh[iaunt,2],hh[iaunt,5]),file="log_pedigree.txt",append=TRUE)
  }}}}}

################################################################################
# Extra: Check whether full aunts are earlier in the data than offspring, through sire
################################################################################
  write('        Extra: Check whether full aunts are earlier in the data than offspring, through sire',file="log_pedigree.txt",append=TRUE)

  for(ind in 1:(nped-2)){
  ma<-hh[ind,6]
  pa<-hh[ind,7]
    if(ma>0 & pa>0){
    mapa<-hh[pa,6]
    papa<-hh[pa,7]
      if(mapa>0 & papa>0){
        for(iaunt in (ind+1):nped){
        maaunt<-hh[iaunt,6]
        paaunt<-hh[iaunt,7]
          if(mapa==maaunt & papa==paaunt){
          itel<-itel+1
          write(paste('ident= ',hh[ind,17],hh[ind,1],hh[ind,2],hh[ind,5],'sire= ',hh[ind,7],'aunt= ',hh[iaunt,17],hh[iaunt,1],hh[iaunt,2],hh[iaunt,5]),file="log_pedigree.txt",append=TRUE)
  }}}}}
}

  if(itel==0){
  write('The sequence of the pedigree is ok',file="log_pedigree.txt",append=TRUE)}
  if(itel>0){
  write('The sequence of the pedigree is NOT ok, pedigree_ident,txt may help to identify problems',file="log_pedigree.txt",append=TRUE)
  write.table(ide,file="pedigree_ident.txt",row.names=FALSE,col.names=FALSE)
  }

write('   ',file="log_pedigree.txt",append=TRUE)
write('blocks.txt, ident.txt, numbers.txt and pedigree.txt are input for AMD-AINV.r',file="log_pedigree.txt",append=TRUE)
if(msinglesires==1){
write('because there are single sires, singlesires.txt is also input for AMD-AINV.r',file="log_pedigree.txt",append=TRUE)} 

write('',file="log_pedigree.txt",append=TRUE) 
write('In the script of AMD-AINV.r four parameters need to be defined',file="log_pedigree.txt",append=TRUE)
write('equi: about the relationshsips between DPQ in the base population;',file="log_pedigree.txt",append=TRUE)
write('pinv: whether you like to compute AINV utilizing mendelian sampling (co)variances, or directly;',file="log_pedigree.txt",append=TRUE)
write('priF: whether or not you like to write inbreeding coefficients of possible future matings;',file="log_pedigree.txt",append=TRUE)
write('yip:  and if you like to compute inbreeding coefficients, for how many years back',file="log_pedigree.txt",append=TRUE)
write('See manual for details',file="log_pedigree.txt",append=TRUE)


################################################################################
# Add sequence number 2 to ssinglesires
################################################################################
if(msinglesires==1){
  for(i in 1:nsinglesires){
  ssinglesires[i,2]<-ssinglesires[i,1]
  seq<-which(ide[,1]==ssinglesires[i,1])
  ssinglesires[i,1]<-seq
  ssinglesires[i,3]<-singlesires[i,3]
  }
write.table(ssinglesires,file="singlesires.txt",row.names=FALSE,col.names=FALSE)
}
################################################################################
# write files. Some are suppressed but may be useful in special cases
################################################################################
# write.table(siredams,file="siredams.txt",row.names=FALSE,col.names=FALSE)
# write.table(f,file="pedigree-check.txt",row.names=FALSE,col.names=FALSE)
write.table(identt,file="ident.txt",row.names=FALSE,col.names=FALSE)
write.table(hh,file="pedigree.txt",row.names=FALSE,col.names=FALSE)
write.table(blocks,file="blocks.txt",row.names=FALSE,col.names=FALSE)
write.table(numbers,file="numbers.txt",row.names=FALSE,col.names=FALSE)

###### This bracket is for not executing most of the programme because #########
###### There are queens with more than one year of birth
}

###### This bracket is for not executing most of the programme because #########
###### there are single sires of which the dam doesn't have a record   #########
}

###### This bracket is for not executing the programme because double queens
###### in input-pedigree.txt
}


}








