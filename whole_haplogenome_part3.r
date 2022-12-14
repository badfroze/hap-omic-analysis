argv <- commandArgs(TRUE)
start <- as.numeric(argv[1])
####haplotype
###################################################################################################################
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case/case.var_raw.vcf  > case.vcf
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case_maternal_2/case_ma.var_raw.vcf  > case_ma.vcf
# awk -F" " '{print $1,$2}' /home/ccwu/gpu/2021_downsy/haplotype/case_paternal_2/case_pa.var_raw.vcf  > case_pa.vcf
###################################################################################################################
csv = read.table('case.vcf')
pav = read.table('case_pa.vcf')
mav = read.table('case_ma.vcf')
###################################################################################################################
ch = c(1:22,'X','Y')
i = start
ii = ch[i]
chr = paste0('chr',ii)
temp.cs = csv[which(csv[,1] == chr),]
temp.pa = pav[which(pav[,1] == chr),]
temp.ma = mav[which(mav[,1] == chr),]
temp.pos = unique(c(temp.cs[,2],temp.pa[,2],temp.ma[,2]))
##############################################################################################################################################################################
###############################################################################################################################################################################
# temp.pos = unique(c(temp.cs[,5],temp.pa[,5],temp.ma[,5]))
# bp = read.table(paste0('/home/ccwu/haplotype/haplotype/whole_genome/chr/chr',ci,'.uniq.bedpe'))
#############################################################
chr = read.table(paste0('/home/ccwu/haplotype/haplotype/whole_genome/homo_sapiens/chr',ii,'.fa'))
chrn = chr[-1,]
assign(paste0('chr',ii),chrn)
source('pickup.sequence.R')
#############################################################
load(paste0('hap_pair_chr',ii,'.Rdata'))
#############bt btseq
source('list_string_diff.R')
#####################
paternal = c()
maternal = c()
intra = c()
mat = as.data.frame(array(,dim = c(nrow(bt),2)))
##
for (j in 1:nrow(bt)){
	##############################################################################################################
	n1 = bt[j,2] + 1
	n2 = bt[j,3]
	#left# 0 both 1 pa 2 ma
	seq1 = pickup.sequence(n1,n2,ii)
	tem1 = btseq[j,1]
	if (seq1 != tem1){
		if (nchar(seq1) == nchar(tem1)){
			temp = list_string_diff(seq1, tem1, only.position = FALSE)
			pos = bt[j,2] + temp[1,1] - 1
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
		else{
			nn = 0
		}
	}
	else{
		nn = 0
	}
	te1 = nn
	##############################################################################################################
	n3 = bt[j,5] + 1
	n4 = bt[j,6]
	#right#
	seq2 = pickup.sequence(n3,n4,ii)
	tem2 = btseq[j,2]
	if (seq2 != tem2){
		if (nchar(seq2) == nchar(tem2)){
			temp = list_string_diff(seq2, tem2, only.position = FALSE)
			pos = bt[j,5] + temp[1,1] - 1
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
		else{
			nn = 0
		}
	}
	else{
		nn = 0
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
save(paternal, maternal, intra, mat, bt, file = paste0('hap_v3_chr',ii,'.Rdata'))

