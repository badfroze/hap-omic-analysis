argv <- commandArgs(TRUE)
start <- as.numeric(argv[1])
# end <- as.numeric(argv[2])
#####
# load('D:/工作/三体/figure/paper_final/chr21_snv_pama.Rdata')
# load('chr21_snv_pama.Rdata')
#
z = start
# pav = read.table(paste0('/data5/QYYang/yyzou/down/snp/case_pa_chr',z,'.vcf'))
# mav = read.table(paste0('/data5/QYYang/yyzou/down/snp/case_ma_chr',z,'.vcf'))
# csv = read.table(paste0('/data5/QYYang/yyzou/down/snp/case_chr',z,'.vcf'))
# upa = setdiff(pav[,2],mav[,2])
# uma = setdiff(mav[,2],pav[,2])
# up = intersect(upa,csv[,2])
# um = intersect(uma,csv[,2])

# tt = paste0('/data5/QYYang/yyzou/down/haplo-rs/chr',z,'.case_hu.mpileup.txt')
#pa
####
# pos = paste0(up,collapse = '|')
# name = paste0('chr',z,'_pa_unique_snp.txt')
# sl = paste('egrep -w','\"',pos,'\"',tt,'>',name)
# re = system(sl)

# up = read.table(paste0('/data5/QYYang/yyzou/down/haplo-rs/uniq.pa.chr',z,'.rs.mpileup.txt'),fill = T,col.names = c(1:6))
# mapa = as.data.frame(array(,dim = c(nrow(up),5)))
# for (i in 1:nrow(up)){
	# hs = c('A','T','C','G')
	# ng = toupper(as.character(up[i,5]))
	# nt = strsplit(ng,'')
	# nt = nt[[1]]
	# ref = toupper(as.character(up[i,3]))
	# ale = hs[-which(hs == ref)]
	# count1 = length(which(nt == ale[1])) + length(which(nt == ale[2])) + length(which(nt == ale[3]))
	# count2 = length(which(nt == '.'))
	# count3 = length(which(nt == ','))
	# mapa[i,1] = as.numeric(up[i,2])
	# mapa[i,2] = count1
	# mapa[i,3] = count2
	# mapa[i,4] = count3
	# mapa[i,5] = as.numeric(up[i,4])
# }

# um = read.table(paste0('/data5/QYYang/yyzou/down/haplo-rs/uniq.ma.chr',z,'.rs.mpileup.txt'),fill = T,col.names = c(1:6))
# mama = as.data.frame(array(,dim = c(nrow(um),5)))
# for (i in 1:nrow(um)){
	# hs = c('A','T','C','G')
	# ng = toupper(as.character(um[i,5]))
	# nt = strsplit(ng,'')
	# nt = nt[[1]]
	# ref = toupper(as.character(um[i,3]))
	# ale = hs[-which(hs == ref)]
	# count1 = length(which(nt == ale[1])) + length(which(nt == ale[2])) + length(which(nt == ale[3]))
	# count2 = length(which(nt == '.'))
	# count3 = length(which(nt == ','))
	# mama[i,1] = as.numeric(um[i,2])
	# mama[i,2] = count1
	# mama[i,3] = count2
	# mama[i,4] = count3
	# mama[i,5] = as.numeric(um[i,4])
# }
# save(mapa,mama,file = paste0('hap_rs_chr',z,'.Rdata'))
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
gtf = read.table('gencode.v19.annotation.gtf',fill = T,col.names = c(1:50))
tf = gtf[which(gtf[,3] == 'gene'),c(1:5,7,10,16,19,22)]
colnames(tf) = paste0('V',1:ncol(tf))
for (i in 1:nrow(tf)){
	temp = as.character(tf[i,7])
	tem = substr(temp,1,15)
	tf[i,11] = tem
}
ge = tf[,c(1:6,8:11)]
colnames(ge) = paste0('V',1:ncol(ge))
#############################################################################
# setwd("D:/工作/三体/figure/haplo-rnaseq/rep_hu")
load('/data5/QYYang/yyzou/down/haplo-rs/ge.Rdata')
ch = c(1:22,'X','Y')
# x = 1
for (x in 1:24){
	z = ch[x]
	load(paste0('hu-2.hap_rs_chr',z,'.Rdata'))
	#######################################
	mapa = mapa[,c(1:4)]
	mama = mama[,c(1:4)]
	a = merge(mapa,mama,by = 'V1',all = T)
	a[is.na(a)] = 0
	a[,8] = a[,2] + a[,6] + a[,7]
	a[,9] = a[,5] + a[,3] + a[,4]
	#######################################
	#mapa
	te = a
	# te = a[,c(1,8)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,c(1,8)]
	colnames(mat) = c('V1','V2')
	for (i in 1:(length(num4) - 1)){
		n1 = num4[i] + 1
		n2 = num4[(i + 1)]
		va = te[n1:n2,8]
		matd = as.data.frame(array(,dim = c(1,2)))
		matd[1,1] = num[n1]
		matd[1,2] = max(va)
		mat = rbind(mat,matd)
	}
	mp = mat
	#######################################
	#mama
	te = a
	# te = a[,c(1,9)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,c(1,9)]
	colnames(mat) = c('V1','V2')
	for (i in 1:(length(num4) - 1)){
		n1 = num4[i] + 1
		n2 = num4[(i + 1)]
		va = te[n1:n2,9]
		matd = as.data.frame(array(,dim = c(1,2)))
		matd[1,1] = num[n1]
		matd[1,2] = max(va)
		mat = rbind(mat,matd)
	}
	mm = mat
	######################################
	# sum(mp[which(mp[,1]>12000000),2])
	# sum(mm[which(mm[,1]>12000000),2])
	######################################
	chr = paste0('chr',z)
	gg = ge[which(ge[,1] == chr),]
	for (i in 1:nrow(gg)){
		num1 = as.numeric(as.character(gg[i,4]))
		num2 = as.numeric(as.character(gg[i,5]))
		region = c(num1:num2)
		ln = abs(num1 - num2)
		ctp = mp[mp[,1]%in%intersect(region,mp[,1]),2]
		ctp = sum(ctp)
		ctm = mm[mm[,1]%in%intersect(region,mm[,1]),2]
		ctm = sum(ctm)
		#
		gg[i,11] = ctp
		gg[i,12] = ctm
		gg[i,13] = ln
		gg[i,14] = ctp/ln*1000000
		gg[i,15] = ctm/ln*1000000
	}
	mat = gg
	#
	assign(paste0('mat',x),mat)
	# temp = as.data.frame(array(0,dim = c(1,7)))
	# mlgb = rbind(mlgb,mat)
}
mlgb = as.data.frame(array(,dim = c(0,15)))
# ch = c(1:22,'X','Y')
for (x in 1:24){
	temp = get(paste0('mat',x))
	mlgb = rbind(mlgb,temp)
}
#
save(mlgb,file = 'rep_hu_1_hap_rs.Rdata')
#
load('rep_hu_1_hap_rs.Rdata')
m1 = mlgb
load('rep_hu_2_hap_rs.Rdata')
m2 = mlgb
#
mt = m1
mt[,14] = mt[,11]/mt[,13]
mt[,15] = mt[,12]/mt[,13]
mt[,16] = mt[,14] + mt[,15]
mt[,17] = mt[,14]*1000000/sum(mt[,16])
mt[,18] = mt[,15]*1000000/sum(mt[,16])
num1 = which(mt[,17] > 0)
num2 = which(mt[,18] > 0)
num3 = union(num1,num2)
mt[,19] = 0
mt[num3,19] = 1
mb = mt[which(mt[,19] == 1),]
mb[,20] = apply(mb[,17:18], 1, function(x) (x[1]-x[2]))/apply(mb[,17:18],1,max)
b1 = mb
colnames(b1) = c('chr','datasets','gene','start','end','strand','type','known','name','ensg_name','count_pa','count_ma','gene_length','cvlpa','cvlma','cvl','tpm_pa','tpm_ma','expressed_or_not','change_ratio')
#
mt = m2
mt[,14] = mt[,11]/mt[,13]
mt[,15] = mt[,12]/mt[,13]
mt[,16] = mt[,14] + mt[,15]
mt[,17] = mt[,14]*1000000/sum(mt[,16])
mt[,18] = mt[,15]*1000000/sum(mt[,16])
num1 = which(mt[,17] > 0)
num2 = which(mt[,18] > 0)
num3 = union(num1,num2)
mt[,19] = 0
mt[num3,19] = 1
mb = mt[which(mt[,19] == 1),]
mb[,20] = apply(mb[,17:18], 1, function(x) (x[1]-x[2]))/apply(mb[,17:18],1,max)
b2 = mb
colnames(b2) = c('chr','datasets','gene','start','end','strand','type','known','name','ensg_name','count_pa','count_ma','gene_length','cvlpa','cvlma','cvl','tpm_pa','tpm_ma','expressed_or_not','change_ratio')
save(m1,m2,file = 'count_hap.Rdata')
#
t1 = b1[,c(1:10,17,18,20)]
t2 = b2[,c(1:10,17,18,20)]
t3 = merge(t1,t2,by = c('chr','datasets','gene','start','end','strand','type','known','name','ensg_name'))
type = unique(t3[,7])
###################################deseq2
load('count_hap.Rdata')
library(stringi)
library(DESeq2)
#
de1 = m1[,c(1,7,9,10,11,12)]
de2 = m2[,c(1,7,9,10,11,12)]
mat = merge(de1,de2,by = c('V1','V7','V9','V10'))
mat[,9] = mat[,7]
mat[,10] = mat[,6]
mat[,6:7] = mat[,9:10]
mat = mat[,1:8]
colnames(mat)[5:8] = c('pa1','pa2','ma1','ma2')
#
source('diff.express.seq2.R')
matd = mat[,5:8]
abet = 'hap'
diff.express.seq2(matd,abet)
#
de3 = read.csv('hap.csv')
de4 = cbind(mat,de3)
de5 = na.omit(de4)
#
de5[,21] = -log(de5[,15],10)
mat = de5
for (i in 1:nrow(mat)){
	#
	t1 = mat[i,12]
	t2 = mat[i,21]
	if (t1 >= log(1.5,2) & t2 >= -log(0.05,10)){
		mat[i,22] = 'a'
	}
	else if (t1 <= -log(1.5,2) & t2 >= -log(0.05,10)){
		mat[i,22] = 'b'
	}
	else{
		mat[i,22] = 'c'
	}
	
}
# colnames(mat) = c('chr','gene','esg','logfc','pv','logpv','type')
p = ggplot(mat, aes(x=log2FoldChange, y=V21,color = V22)) + geom_point(size = 2) + xlim(-15,15) + ylim(0,300)
p = p + scale_discrete_manual(values=c('#CC0000',"#2f5688","#BBBBBB"),aesthetics = 'colour',labels = c('a','b','c'))
p = p + geom_hline(yintercept = -log(0.05,10),linetype = 'dashed')
p = p + geom_vline(xintercept = c(-log(1.5,2),log(1.5,2)),linetype = 'dashed')
#
p = p + theme_bw() +
	theme(panel.grid.major=element_line(colour=NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank())
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
ma21 = mat[which(mat[,1] == 'chr21'),]
p = ggplot(ma21, aes(x=log2FoldChange, y=V21,color = V22)) + geom_point(size = 2) + xlim(-15,15) + ylim(0,100)
p = p + scale_discrete_manual(values=c('#CC0000',"#2f5688","#BBBBBB"),aesthetics = 'colour',labels = c('a','b','c'))
p = p + geom_hline(yintercept = -log(0.05,10),linetype = 'dashed')
p = p + geom_vline(xintercept = c(-log(1.5,2),log(1.5,2)),linetype = 'dashed')
#
p = p + theme_bw() +
	theme(panel.grid.major=element_line(colour=NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank())
##











