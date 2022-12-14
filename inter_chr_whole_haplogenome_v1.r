argv <- commandArgs(TRUE)
av1 <- as.numeric(argv[1])
av2 <- as.numeric(argv[2])
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
#part1
####################################################################################################################
####################################################################################################################
###################################################################################################################
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case/case.var_raw.vcf  > case.vcf
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case_maternal_2/case_ma.var_raw.vcf  > case_ma.vcf
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case_paternal_2/case_pa.var_raw.vcf  > case_pa.vcf
###################################################################################################################
csv = read.table('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/case.vcf')
pav = read.table('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/case_pa.vcf')
mav = read.table('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/case_ma.vcf')
###################################################################################################################
###################################################################################################################
ch = c(1:22,'X','Y')
i = av1
ii = ch[i]
j = av2
jj = ch[j]
###################################################################################################################
# bp = read.table(paste0('/home/ccwu/haplotype/haplotype/whole_genome/chr/chr',ii,'.uniq.bedpe'))
#awk '$1=="chr1" && $4=="chr2" {print}' ../inter_chr.uniq.bedpe > temp.uniq.bedpe
line0 = paste0('temp',ii,'_',jj,'.uniq.bedpe')
line1 = paste0('\'$1==\"chr',ii,'\" && $4==\"chr',jj,'\" {print}\'')
line2 = paste0('\'$1==\"chr',jj,'\" && $4==\"chr',ii,'\" {print}\'')
line3 = paste('awk',line1,'/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/chr/inter_chr.uniq.bedpe >',line0)
line4 = paste('awk',line2,'/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/chr/inter_chr.uniq.bedpe >>',line0)
system(line3)
system(line4)
####
bp = read.table(line0)
####
line5 = paste('rm',line0)
system(line5)
#####################################################################
temp1 = as.character(bp[1,1])
temp2 = as.character(bp[1,4])
iii = substr(temp1,4,nchar(temp1))
jjj = substr(temp2,4,nchar(temp2))
###################################################################################################################
chr = paste0('chr',iii)
temp.cs = csv[which(csv[,1] == chr),]
temp.pa = pav[which(pav[,1] == chr),]
temp.ma = mav[which(mav[,1] == chr),]
temp.pos = unique(c(temp.cs[,2],temp.pa[,2],temp.ma[,2]))
n = floor(max(temp.pos)/1000000)
for (x in 0:n){
	gap1 = x*1000000 + 1
	gap2 = (x + 1)*1000000
	gapn = intersect(c(gap1:gap2),temp.pos)
	assign(paste0('region',x),gapn)
}
n1 = apply(bp,1,function(x) length(intersect(get(paste0('region',floor(as.numeric(x[2])/1000000))),c((as.numeric(x[2]) + 1):(as.numeric(x[3]))))))
###################################################################################################################
chr = paste0('chr',jjj)
temp.cs = csv[which(csv[,1] == chr),]
temp.pa = pav[which(pav[,1] == chr),]
temp.ma = mav[which(mav[,1] == chr),]
temp.pos = unique(c(temp.cs[,2],temp.pa[,2],temp.ma[,2]))
n = floor(max(temp.pos)/1000000)
for (x in 0:n){
	gap1 = x*1000000 + 1
	gap2 = (x + 1)*1000000
	gapn = intersect(c(gap1:gap2),temp.pos)
	assign(paste0('regionsec',x),gapn)
}
n2 = apply(bp,1,function(x) length(intersect(get(paste0('regionsec',floor(as.numeric(x[5])/1000000))),c((as.numeric(x[5]) + 1):(as.numeric(x[6]))))))
#################################################################################################################################################
n1 = as.numeric(n1)
n2 = as.numeric(n2)
#######################
num1 = which(n1 > 0)
num2 = which(n2 > 0)
num3 = union(num1,num2)
bt = bp[num3,7:8]
btseq = as.data.frame(array(,dim = c(nrow(bt),2)))
######################################################################################################
for (y in 1:24){
	ch = c(1:22,'X','Y')
	yy = ch[y]
	name1 = paste0('/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/pet/chr',yy,'.pet1.uniq.part.sam')
	name2 = paste0('/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/pet/chr',yy,'.pet2.uniq.part.sam')
	pet1_1 = read.table(name1)
	pet1_2 = read.table(name2)
	pet = rbind(pet1_1,pet1_2)
	save(pet,file = paste0('/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/pet/pet_chr',yy,'.Rdata'))
}
#/home/ccwu/haplotype/haplotype/whole_genome##########################################################
load(paste0('/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/pet/pet_chr',iii,'.Rdata'))
pet1 = pet
#####
load(paste0('/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/pet/pet_chr',jjj,'.Rdata'))
pet2 = pet
#####
pet = merge(pet1,pet2,by='V1',all = F,suffixes = c('pet1','pet2'))
colnames(pet) = c('title','pos1','seq1','pos2','seq2')
colnames(bt) = c('title','none')
pet_bt = merge(pet,bt,by='title',all = F,suffixes = c('pet','bt'))
#
# save(pet_bt,file = paste0('part2_pair_chr',ii,'.Rdata'))
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
#part3
####################################################################################################################
####################################################################################################################
# csv = read.table('case.vcf')
# pav = read.table('case_pa.vcf')
# mav = read.table('case_ma.vcf')
###################################################################################################################
# ch = c(1:22,'X','Y')
# i = start
# ii = ch[i]
# chr = paste0('chr',ii)
# temp.cs = csv[which(csv[,1] == chr),]
# temp.pa = pav[which(pav[,1] == chr),]
# temp.ma = mav[which(mav[,1] == chr),]
# temp.pos = unique(c(temp.cs[,2],temp.pa[,2],temp.ma[,2]))
##################################################################################################
##################################################################################################
#############################################################
# chr = read.table('/home/ccwu/haplotype/haplotype/whole_genome/homo_sapiens/chr21.fa')
chr = read.table(paste0('/home/ccwu/data1/homo_sapiens/chr',iii,'.fa'))
chrn = chr[-1,]
assign(paste0('chr',iii),chrn)
chr = read.table(paste0('/home/ccwu/data1/homo_sapiens/chr',jjj,'.fa'))
chrn = chr[-1,]
assign(paste0('chr',jjj),chrn)
source('pickup.sequence.R')
#############################################################
#############pet_bt
source('list_string_diff.R')
##############################################################################
paternal = c()
maternal = c()
intra = c()
bt = pet_bt
mat = as.data.frame(array(,dim = c(nrow(bt),2)))
##
for (j in 1:nrow(bt)){
	##############################################################################################################
	tem1 = as.character(bt[j,3])
	n1 = bt[j,2]
	len = nchar(tem1)
	n2 = n1 + len - 1
	#left# 0 both 1 pa 2 ma
	seq1 = pickup.sequence(n1,n2,iii)
	nn = 0
	if (seq1 != tem1){
		if (nchar(seq1) == nchar(tem1)){
			temp = list_string_diff(seq1, tem1, only.position = FALSE)
			pos = n1 + temp[1,1] - 1
			n1 = length(which(temp.cs[,2] == pos))
			n2 = length(which(temp.pa[,2] == pos))
			n3 = length(which(temp.ma[,2] == pos))
			n4 = c(n1,n2,n3)
			n = n1 + n2 + n3
			if (n == 3){
				nn = 0
			}
			if (n == 0){
				nn = 0
			}
			if (n == 1){
				num = which(n4 == 1)
				if (num == 1){
					nn = 0
				}
				if (num == 2){
					nn = 1
				}
				if (num == 3){
					nn = 2
				}
			}
			if (n == 2){
				num = which(n4 == 0)
				if (num == 1){
					nn = 0
				}
				if (num == 2){
					nn = 2
				}
				if (num == 3){
					nn = 1
				}
			}
		}
	}
	te1 = nn
	##############################################################################################################
	tem1 = as.character(bt[j,5])
	n1 = bt[j,4]
	len = nchar(tem1)
	n2 = n1 + len - 1
	#left# 0 both 1 pa 2 ma
	seq1 = pickup.sequence(n1,n2,jjj)
	nn = 0
	if (seq1 != tem1){
		if (nchar(seq1) == nchar(tem1)){
			temp = list_string_diff(seq1, tem1, only.position = FALSE)
			pos = n1 + temp[1,1] - 1
			n1 = length(which(temp.cs[,2] == pos))
			n2 = length(which(temp.pa[,2] == pos))
			n3 = length(which(temp.ma[,2] == pos))
			n4 = c(n1,n2,n3)
			n = n1 + n2 + n3
			if (n == 3){
				nn = 0
			}
			if (n == 0){
				nn = 0
			}
			if (n == 1){
				num = which(n4 == 1)
				if (num == 1){
					nn = 0
				}
				if (num == 2){
					nn = 1
				}
				if (num == 3){
					nn = 2
				}
			}
			if (n == 2){
				num = which(n4 == 0)
				if (num == 1){
					nn = 0
				}
				if (num == 2){
					nn = 2
				}
				if (num == 3){
					nn = 1
				}
			}
		}
	}
	te2 = nn
	##############################################################
	mat[j,1] = te1
	mat[j,2] = te2
	te = paste(te1,te2)
	if (te != '0 0'){
		#
		if (te == '1 0'){
			paternal = c(paternal,j)
		}
		if (te == '0 1'){
			paternal = c(paternal,j)
		}
		if (te == '1 1'){
			paternal = c(paternal,j)
		}
		#
		if (te == '2 0'){
			maternal = c(maternal,j)
		}
		if (te == '0 2'){
			maternal = c(maternal,j)
		} 
		if (te == '2 2'){
			maternal = c(maternal,j)
		}
		#
		if (te == '1 2'){
			intra = c(intra,j)
		} 
		if (te == '2 1'){
			intra = c(intra,j)
		}
	}
}
bd = bt[,c(2,4)]
save(iii,jjj,paternal,maternal,intra,bd,mat,file = paste0('inter_results_chr',ii,'_',jj,'.Rdata'))
##############################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
#part4
library(RColorBrewer)
library(pheatmap)
cl = brewer.pal(11,'Spectral')
cl2 = rev(cl)
###################################
ch = c(1:22,'X','Y')
n = c()
for (k in 1:24){
	kk = ch[k]
	load(paste0('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/part3_results_chr',kk,'.Rdata'))
	####
	x = ceiling(max(bd[,1])/10000000)
	pa = matrix(0,x,x)
	for (i in 1:length(paternal)){
		temp = paternal[i]
		pos1 = bd[temp,1]
		p1 = ceiling(pos1/10000000)
		pos2 = bd[temp,2]
		p2 = ceiling(pos2/10000000)
		pa[p1,p2] = pa[p1,p2] + 1
		pa[p2,p1] = pa[p2,p1] + 1
	}
	ma = matrix(0,x,x)
	for (i in 1:length(maternal)){
		temp = maternal[i]
		pos1 = bd[temp,1]
		p1 = ceiling(pos1/10000000)
		pos2 = bd[temp,2]
		p2 = ceiling(pos2/10000000)
		ma[p1,p2] = ma[p1,p2] + 1
		ma[p2,p1] = ma[p2,p1] + 1
	}
	n = c(n,x)
	assign(paste0('pa',k),pa)
	assign(paste0('ma',k),ma)
}
m = sum(n)
wpa = matrix(0,m,m)
for (i in 1:24){
	#
	num1 = sum(n[1:i]) - n[i] + 1
	num2 = sum(n[1:i])
	xj = get(paste0('pa',i))
	wpa[num1:num2,num1:num2] = xj
}
wma = matrix(0,m,m)
for (i in 1:24){
	#
	num1 = sum(n[1:i]) - n[i] + 1
	num2 = sum(n[1:i])
	xj = get(paste0('ma',i))
	wma[num1:num2,num1:num2] = xj
}
####################################################
pa = wpa
ma = wma
for (i in 1:23){
	ii = ch[i]
	for (j in (i + 1):24){
		jj = ch[j]
		#iii,jjj,paternal,maternal,intra,bd,mat
		load(paste0('inter_results_chr',ii,'_',jj,'.Rdata'))
		###
		i4 = which(ch == iii)
		j4 = which(ch == jjj)
		###
		for (k in 1:length(paternal)){
			temp = paternal[k]
			pos1 = bd[temp,1]
			p1 = ceiling(pos1/10000000) + sum(n[1:i4]) - n[i4]
			pos2 = bd[temp,2]
			p2 = ceiling(pos2/10000000) + sum(n[1:j4]) - n[j4]
			pa[p1,p2] = pa[p1,p2] + 1
			pa[p2,p1] = pa[p2,p1] + 1
		}
		for (k in 1:length(maternal)){
			temp = maternal[k]
			pos1 = bd[temp,1]
			p1 = ceiling(pos1/10000000) + sum(n[1:i4]) - n[i4]
			pos2 = bd[temp,2]
			p2 = ceiling(pos2/10000000) + sum(n[1:j4]) - n[j4]
			ma[p1,p2] = ma[p1,p2] + 1
			ma[p2,p1] = ma[p2,p1] + 1
		}
	}
}
pm = ma
for (i in 1:nrow(pm)){
	for (j in i:nrow(pm)){
		pm[i,j] = pa[i,j]
	}
}
wpm = pm
###################################################
cl = brewer.pal(11,'Spectral')
cl2 = rev(cl)
# pm[which(pm > 10,arr.ind = T)] = 10
test = wpm
test[75,] = 0
test[,75] = 0
test[which(test > 100,arr.ind = T)] = 100
pheatmap(test,cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = cl2[c(1:4,7)])(100))
##############################################################################################
##############################################################################################
kk = 1
load(paste0('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/part3_results_chr',kk,'.Rdata'))
x1 = ceiling(max(bd[,1])/100000)
kk = 21
load(paste0('/home/ccwu/gpu/2021_downsy/haplotype/whole_genome_2/part3_results_chr',kk,'.Rdata'))
x2 = ceiling(max(bd[,1])/100000)
#
pa = matrix(0,x1,x2)
ma = matrix(0,x2,x1)
ch = c(1:22,'X','Y')
i = 1
j = 21
# for (i in 1:23){
ii = ch[i]
# for (j in (i + 1):24){
jj = ch[j]
#iii,jjj,paternal,maternal,intra,bd,mat
load(paste0('inter_results_chr',ii,'_',jj,'.Rdata'))
###
i4 = which(ch == iii)
j4 = which(ch == jjj)
###
for (k in 1:length(paternal)){
	temp = paternal[k]
	pos1 = bd[temp,1]
	p1 = ceiling(pos1/100000) #+ sum(n[1:i4]) - n[i4]
	pos2 = bd[temp,2]
	p2 = ceiling(pos2/100000) #+ sum(n[1:j4]) - n[j4]
	pa[p1,p2] = pa[p1,p2] + 1
	# pa[p2,p1] = pa[p2,p1] + 1
}
for (k in 1:length(maternal)){
	temp = maternal[k]
	pos1 = bd[temp,1]
	p1 = ceiling(pos1/100000) #+ sum(n[1:i4]) - n[i4]
	pos2 = bd[temp,2]
	p2 = ceiling(pos2/100000) #+ sum(n[1:j4]) - n[j4]
	# ma[p1,p2] = ma[p1,p2] + 1
	ma[p2,p1] = ma[p2,p1] + 1
}
# }
# }
pm = ma
for (i in 1:nrow(pm)){
	for (j in i:nrow(pm)){
		pm[i,j] = pa[i,j]
	}
}
wpm = pm
pheatmap(pa,cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = c("white","red"))(100))












