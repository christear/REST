library(data.table)
#library(parallel)
###
# 1, read the output from the prediction
# 2, read the motif position file based on predction output
# 3, read the sequence file baed on predction output
# 4, relabel the true:1 and false:0 as the last coloum
###
#pasmotifs = c('AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA', 'AATACA', 'CATAAA', 'GATAAA','AATGAA', 'TTTAAA','ACTAAA','AATAGA')
dcut = 100
###

argv=commandArgs(TRUE)
tsv = as.character(argv[1])
pasf = as.character(argv[2])
pasfda = read.table(pasf)
pasmotifs = unlist(pasfda[,1])

# read files
cat("#reading files ...\n")
predout = fread(text = tsv,sep = "\t")
motifposf = paste(tsv,".pas.motif.pos",sep = "")
motifpos = fread(input = motifposf,sep = "\t")
motifpos = motifpos[order(motifpos[,1]),]
faf = paste(tsv,".fa",sep = "")
fa = fread(input = faf,sep = "\t",head = F)
faseq = fa[seq(2,nrow(fa),2),1]
#faseq = fa[-grep(">",fa[,1]),1]

##
cat("#processing PAS motif and internal primming...\n")
# PAS have any of the first three PAS motif
pasmkl = apply(motifpos[,2:4],1,sum) > 0
# PAS have any of the PAS motif
pasmkl2 = apply(motifpos[,-1],1,sum) > 0
### determin internal primming
determinip = function(seql){
    subs = substr(seql,101,110)
    ip = ifelse(length(grep("A{7}",subs,perl = T)) > 0 | sum(unlist(strsplit(subs,"")) == "A") > 7,1,0)
    ip
}
#
ipmkl = apply(faseq,1,function(x) determinip(x[1]))
# relabel to get the updated true/false
cat("#relabeling ...\n")
newlab = predout[,2]
if(FALSE){
    ### relabel PAS in ground truth without PAS motif, but with internal primming as negative
    newlab[predout[,2] > 0 & pasmkl2 == FALSE & ipmkl == 1] = 0
}
## relabel PAS in FN without PAS motif and with internal primming as negative
newlab[predout[,2] == 1 & predout[,5] == 0 & pasmkl2 == FALSE & ipmkl == 1] = 0
## relabel predicted PAS not from ground truth (FP) as negative
newlab[predout[,2] == 0 & predout[,5] == 1] = 0
###  relabel FP with PAS motif as positive
#newlab[predout[,2] == 0 & predout[,5] == 1 & pasmkl == TRUE] = 1
###  relabel FP with PAS motif as positive and without internal primming as positive
newlab[predout[,2] == 0 & predout[,5] == 1 & pasmkl == TRUE & ipmkl == 0] = 1
newmat = cbind(predout,newlab)
colnames(newmat)[ncol(newmat)] = "relabel"
write.table(newmat,file = gsub("\\.tsv",".relabeled.tsv",tsv),sep = "\t",quote = F,row.names = F)

# summary of the prediction
cat("#summarizing ...\n")
tp = sum(predout[,2] == 1 & predout[,5] == 1)
fp = sum(predout[,2] == 0 & predout[,5] == 1)
tn = sum(predout[,2] == 0 & predout[,5] == 0)
fn = sum(predout[,2] == 1 & predout[,5] == 0)
cat("#TP:FP:TN:FN\t",tp,fp,tn,fn,"\n")
gtinf = rbind(c(sum(predout[,2] == 1),sum(predout[,2] == 1 & pasmkl2 == TRUE),sum(predout[,2] == 1 & ipmkl == 1)),c(sum(predout[,2] == 0),sum(predout[,2] == 0 & pasmkl2 == TRUE),sum(predout[,2] == 0 & ipmkl == 1)))
predinf = rbind(c(sum(predout[,5] == 1),sum(predout[,5] == 1 & pasmkl2 == TRUE),sum(predout[,5] == 1 & ipmkl == 1)),c(sum(predout[,5] == 0),sum(predout[,5] == 0 & pasmkl2 == TRUE),sum(predout[,5] == 0 & ipmkl == 1)))
newinf = rbind(c(sum(newlab == 1),sum(newlab == 1 & pasmkl2 == TRUE),sum(newlab == 1 & ipmkl == 1)),c(sum(newlab == 0),sum(newlab == 0 & pasmkl2 == TRUE),sum(newlab == 0 & ipmkl == 1)))
allinf = rbind(gtinf,predinf,newinf)
rownames(allinf) = c("ground_truth_P","ground_truth_N","predicted_P","predicted_N","newlabel_P","newlabel_N")
colnames(allinf) = c("total_number","PAS","internal_primming")
write.table(allinf,file = paste(tsv,".pas.internal.primming.inf.txt",sep = ""),sep = "\t",quote = F)

calfreq = function(mat,dcut = 100){
    #mat2 = apply(mat,2,function(x) replacev(x,100-dcut,100))
    mat2 = as.data.frame(mat)
    #ttfreq = apply(mat,2,function(x) sum(x > 100-dcut & x < 100))
    f1 = sum(mat2[,1] > 100-dcut & mat2[,1] < 100)
    f2 = sum(mat2[,1] == 0 & mat2[,2] > 100-dcut & mat2[,2] < 100)
    f313 = sapply(3:ncol(mat2),function(i) sum(apply(mat2[,1:(i-1)],1,max) == 0 & mat2[,i] > 100-dcut & mat2[,i] < 100))
    freq = c(f1,f2,f313)/nrow(mat)
    names(freq) = pasmotifs
    freq
}
gtids = (1:nrow(motifpos))[predout[,2] == 1]
gtmotif = calfreq(motifpos[gtids,-1])
#
predids = (1:nrow(motifpos))[predout[,5] == 1]
predmotif = calfreq(motifpos[predids,-1])
#
newids = (1:nrow(motifpos))[newlab == 1]
newmotif = calfreq(motifpos[newids,-1])
motifmat = rbind(gtmotif,predmotif,newmotif)
rownames(motifmat) = c("ground_truth","predicted_truth","new_truth")
write.table(motifmat,file = paste(tsv,".pas.summary.txt",sep = ""),sep = "\t",quote = F)




