argv <- commandArgs(TRUE)
start <- as.numeric(argv[1])
####haplotype
for (i in c(2:20,22:24)){
	ch = c(1:22,'X','Y')
	ii = ch[i]
	#
	load(paste0('D:/工作/三体/figure/figure4_haplohic/part3_results_chr',ii,'.Rdata'))
	###################################
	library(RColorBrewer)
	library(pheatmap)
	cl = brewer.pal(11,'Spectral')
	cl2 = rev(cl)
	###################################
	fb = 10000
	dx = ceiling(max(bd)/fb)
	###################################
	pa = matrix(0,dx,dx)
	for (i in 1:length(paternal)){
		temp = paternal[i]
		pos1 = bd[temp,1]
		p1 = ceiling(pos1/fb)
		pos2 = bd[temp,2]
		p2 = ceiling(pos2/fb)
		pa[p1,p2] = pa[p1,p2] + 1
		pa[p2,p1] = pa[p2,p1] + 1
	}
	####
	ma = matrix(0,dx,dx)
	for (i in 1:length(maternal)){
		temp = maternal[i]
		pos1 = bd[temp,1]
		p1 = ceiling(pos1/fb)
		pos2 = bd[temp,2]
		p2 = ceiling(pos2/fb)
		ma[p1,p2] = ma[p1,p2] + 1
		ma[p2,p1] = ma[p2,p1] + 1
	}
	save(pa,ma,file = paste0('hap_results_chr',ii,'.Rdata'))
}
# ma[which(ma > 10,arr.ind = T)] = 10
# pheatmap(ma,cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = cl2[c(1:4,7)])(100))
####
# load('D:/工作/三体/figure/paper_final/haplo_mat_fb50K.Rdata')
# load('D:/工作/三体/figure/paper_final/haplo_mat_imputed2_fb10K.Rdata')
# load('D:/工作/三体/figure/figure3_haplotype/whole_hic_data_chr21/haplo_pama_matrix_20k.Rdata')
# fb = 10000
# pm = ma
# for (i in 1:nrow(pm)){
	# for (j in i:nrow(pm)){
		# pm[i,j] = pa[i,j]
	# }
# }
# library(RColorBrewer)
# library(pheatmap)
# cl = brewer.pal(11,'Spectral')
# cl2 = rev(cl)
yc = ma
yc[which(yc > 5,arr.ind = T)] = 5
####
dz = 1000000/fb
num1 = 280*dz + 1
num2 = 282*dz
pheatmap(yc[num1:num2,num1:num2],cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = c("white","red"))(100))
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
s1 = read.table('D:/工作/三体/figure/figure3_haplotype/whole_hic_data_chr21/case1_chr21_10000.ginteractions.tsv')
############################
i1 = 1
############################
# region1 = 1
# region2 = 48000000
fb = 10000
#########################################################
# a = read.table(paste0('s',i1,'_chr21'))
a = get(paste0('s',i1))
#24,300,000-24,800,000
# mat1 = as.data.frame(array(0,dim=c(ncol(a),3)))
# ma = c()
# for (i in 1:nrow(a)){
	# num1 = a[i,2]
	# num2 = a[i,5]
	# if (num1 >= region1 & num1 <= region2){
		# if (num2 >= region1 & num2 <= region2){
			# ma = c(ma,i)
		# }
	# }
# }
mat = a
#########
ut = nrow(pm)
matd = as.data.frame(array(0,dim=c(ut,ut)))
for (i in 1:nrow(mat)){
	num1 = (mat[i,2] - 0)/fb + 1
	num2 = (mat[i,5] - 0)/fb + 1
	matd[num1,num2] = mat[i,7]
	matd[num2,num1] = mat[i,7]
}
########################################################
# mato = matd
# library(RColorBrewer)
# mato[which(mato > 10,arr.ind = T)] = 10
# pheatmap(mato,cluster_rows = F,cluster_cols = F,border = F,labels_row = ' ',labels_col = ' ',color = colorRampPalette(colors = c("white","red"))(100))
#
##################################
# cl = brewer.pal(11,'Spectral')
# cl2 = rev(cl)
# pheatmap(mato,cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = cl2[c(1:4,7)])(100))
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
mt1 = as.matrix(matd)
mt2 = pm
dm1 = diag(mt1)
dm2 = diag(mt2)
sz = (sum(mt1) - sum(dm1))/sum(mt2)
temp = mt1/sz
mt3 = mt2 + temp
# mt1[which(mt1 > 10,arr.ind = T)] = 10
##############
cl = brewer.pal(11,'Spectral')
cl2 = rev(cl)
yc = mt3
yc[which(yc > 5,arr.ind = T)] = 5
# yc = yc + mt1
# pheatmap(yc[521:660,521:660],cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = cl2[c(1:4,7)])(100))
dz = 1000000/fb
num1 = 35.5*dz + 1
num2 = 36.5*dz
pheatmap(yc[num1:num2,num1:num2],cluster_rows = F,cluster_cols = F,border = F,labels_row = 0,labels_col = 0,color = colorRampPalette(colors = cl2[c(1:4,7)])(100))
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
####pa ma
#########
pm = ma
region1 = 26000000
region2 = 33000000
height = 3000000
range1 = (region1 - height)/fb
range2 = (region2 + height)/fb
matm = as.matrix(pm[range1:range2,range1:range2])
ut = nrow(matm)
ut2 = (region2 - region1)/fb + 1
ut3 = height/fb
matsq = matrix(nrow=ut3,ncol=ut2)
for (i in 1:nrow(matsq)){
	sql = ut2
	p1 = ut - ut2 - (i - 1)
	p2 = p1 + ut2 - 1
	p3 = i
	p4 = p3 + ut2 - 1
	matt = matm[p3:p4,p1:p2]
	temp = diag(matt)
	matsq[i,] = temp
}
#
mato = matsq
library(RColorBrewer)
mato[which(mato > 10,arr.ind = T)] = 10
##############
cl = brewer.pal(11,'Spectral')
cl2 = rev(cl)
pheatmap(mato,cluster_rows = F,cluster_cols = F,border = F,labels_row = ' ',labels_col = ' ',color = colorRampPalette(colors = cl2[c(1:4,7)])(100),cellwidth=1,cellheight=1)
















