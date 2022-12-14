#
setwd("D:/工作/三体/figure/haplo-cnv")
#
# z = 1
ch = c(1:22,'X','Y')
for (x in 1:22){
	z = ch[x]
	load(paste0('hap_cnv_chr',z,'.Rdata'))
	#######################################
	mapa = mapa[,c(1:4)]
	mama = mama[,c(1:4)]
	a = merge(mapa,mama,by = 'V1',all = T)
	a[is.na(a)] = 0
	a[,8] = a[,2] + a[,6] + a[,7]
	a[,9] = a[,5] + a[,3] + a[,4]
	#######################################
	#mapa
	te = mapa
	# te = a[,c(1,8)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,1:2]
	colnames(mat) = c('V1','V2')
	for (i in 1:(length(num4) - 1)){
		n1 = num4[i] + 1
		n2 = num4[(i + 1)]
		va = te[n1:n2,2]
		matd = as.data.frame(array(,dim = c(1,2)))
		matd[1,1] = num[n1]
		matd[1,2] = max(va)
		mat = rbind(mat,matd)
	}
	mp = mat
	#######################################
	#mama
	te = mama
	# te = a[,c(1,9)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,1:2]
	colnames(mat) = c('V1','V2')
	for (i in 1:(length(num4) - 1)){
		n1 = num4[i] + 1
		n2 = num4[(i + 1)]
		va = te[n1:n2,2]
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
	cl = c(250000000,243000000,198000000,192000000,182000000,171000000,159000000,146000000,141000000,136000000,135000000,134000000,115000000,107000000,103000000,90000000,83000000,80000000,59000000,64000000,48000000,52000000,156000000,57000000)
	ge = cl[x]
	bin = 1000000
	num = ge/bin
	cap = c()
	cam = c()
	cop = c()
	com = c()
	for (i in 1:num){
		num1 = (i - 1)*bin + 1
		num2 = i*bin
		region = c(num1:num2)
		ctp = mp[mp[,1]%in%intersect(region,mp[,1]),2]
		c1 = length(ctp)
		ctp = sum(ctp)
		ctm = mm[mm[,1]%in%intersect(region,mm[,1]),2]
		c2 = length(ctm)
		ctm = sum(ctm)
		cap = c(cap,ctp)
		cam = c(cam,ctm)
		cop = c(cop,c1)
		com = c(com,c2)
	}
	#
	mat = as.data.frame(array(,dim = c(length(cap),2)))
	mat[,1] = z
	mat[,2] = c(1:num)
	mat[,3] = cap/cop
	mat[,4] = cam/com
	mat[,5] = mat[,3]/mat[,2]
	#
	# is.nan(mat[,2])
	#
	assign(paste0('mat',x),mat)
	# temp = as.data.frame(array(0,dim = c(1,7)))
	# mlgb = rbind(mlgb,mat)
}
mlgb = as.data.frame(array(,dim = c(0,5)))
# ch = c(1:22,'X','Y')
for (x in 1:21){
	temp = get(paste0('mat',x))
	tem = as.data.frame(array(0,dim = c(7,5)))
	mlgb = rbind(mlgb,temp,tem)
	#
}
x = 22
temp = get(paste0('mat',x))
# tem = as.data.frame(array(0,dim = c(7,4)))
mlgb = rbind(mlgb,temp)
mlgb[,6] = c(1:nrow(mlgb))
################################################
mat = mlgb
mat[is.nan(mat[,5]),3:5] = 0
ln = length(which(mat[,5] != 0))
p1 = sum(mat[,3])/ln
p2 = sum(mat[,4])/ln
mat[,7] = mat[,3]/p1
mat[,8] = mat[,4]/p2
mat[,9] = 1
mat[which(mat[,3] == 0),7] = 0
################################################
s1 = sort(mat[,7],decreasing  = T)
ss1 = s1[152]
s2 = sort(mat[,8],decreasing  = T)
ss2 = s2[152]
n1 = c()
n2 = c()
for (i in 1:22){
	temp = mat[which(mat[,1] == i),]
	t1 = length(which(temp[,7] >= ss1))/nrow(temp)
	t2 = length(which(temp[,8] >= ss2))/nrow(temp)
	n1 = c(n1,t1)
	n2 = c(n2,t2)
}
################################################
mat[which(mat[,5] > 4),5] = 4
mat[which(mat[,6] > 4),6] = 4
p = ggplot(mat,aes(x=V1,y=V7,fill = V5)) + geom_bar(stat="identity",width = 1) + xlim(0,nrow(mlgb))
p = p + scale_fill_gradient2(low = '#2f5688',high = '#CC0000',midpoint = 1)
#######################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
##hap复刻多少CNV
cnv = read.table('D:/工作/三体/figure/figure1_CNV/20211023_CNV/raw_cnv_case_parent12_control/case_hu.cnvnator')
cv = cnv[which(cnv[,10] <= 0.05),]
for (i in 1:nrow(cv)){
	temp = as.character(cv[i,3])
	for (j in 1:nchar(temp)){
		tem = substr(temp,j,j)
		if (tem == '-'){
			break
		}
	}
	num1 = substr(temp,1,(j - 1))
	num2 = substr(temp,(j + 1),nchar(temp))
	cv[i,11] = num1
	cv[i,12] = num2
}
##
mt = as.data.frame(array(,dim = c(0,14)))
for (i in 1:22){
	ch = paste0('chr',i)
	cno = cv[which(cv[,2] == ch),]
	#
	load(paste0('hap_cnv_chr',i,'.Rdata'))
	#######################################
	mapa = mapa[,c(1:4)]
	mama = mama[,c(1:4)]
	a = merge(mapa,mama,by = 'V1',all = T)
	a[is.na(a)] = 0
	a[,8] = a[,2] + a[,6] + a[,7]
	a[,9] = a[,5] + a[,3] + a[,4]
	#######################################
	#mapa
	te = mapa
	# te = a[,c(1,8)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,1:2]
	colnames(mat) = c('V1','V2')
	for (j in 1:(length(num4) - 1)){
		n1 = num4[j] + 1
		n2 = num4[(j + 1)]
		va = te[n1:n2,2]
		matd = as.data.frame(array(,dim = c(1,2)))
		matd[1,1] = num[n1]
		matd[1,2] = max(va)
		mat = rbind(mat,matd)
	}
	mp = mat
	#######################################
	#mama
	te = mama
	# te = a[,c(1,9)]
	# te = te[which(te[,2] > 0),]
	num = te[,1]
	num1 = num[1:(length(num) - 1)]
	num2 = num[2:length(num)]
	num3 = num2 - num1
	num4 = which(num3 >= 150)
	mat = te[1,1:2]
	colnames(mat) = c('V1','V2')
	for (j in 1:(length(num4) - 1)){
		n1 = num4[j] + 1
		n2 = num4[(j + 1)]
		va = te[n1:n2,2]
		matd = as.data.frame(array(,dim = c(1,2)))
		matd[1,1] = num[n1]
		matd[1,2] = max(va)
		mat = rbind(mat,matd)
	}
	mm = mat
	######################################################
	for (j in 1:nrow(cno)){
		num1 = as.numeric(cno[j,11])
		num2 = as.numeric(cno[j,12])
		#
		region = c(num1:num2)
		ctp = mp[mp[,1]%in%intersect(region,mp[,1]),2]
		c1 = length(ctp)
		ctp = sum(ctp)
		ctm = mm[mm[,1]%in%intersect(region,mm[,1]),2]
		c2 = length(ctm)
		ctm = sum(ctm)
		cap = c(cap,ctp)
		cam = c(cam,ctm)
		cop = c(cop,c1)
		com = c(com,c2)
		#
		cno[j,13] = ctp/c1
		cno[j,14] = ctm/c2
	}
	mt = rbind(mt,cno)
}
#
# ln = length(which(mat[,5] != 0))
# p1 = sum(mat[,3])/ln
# p2 = sum(mat[,4])/ln
#
mt[is.nan(mt[,13]),13]=0
mt[is.nan(mt[,14]),14]=0
#
ln = nrow(mt)
p1 = sum(mt[,13])/ln
p2 = sum(mt[,14])/ln
#
mt[,15] = mt[,13]/p1
mt[,16] = mt[,14]/p2
mt[,17] = apply(mt[,15:16],1,mean)
#
plot(mt[,5],mt[,16],ylim = c(0,5))
mr<-lm(V16~V5, data=mt)
summary(mr)
abline(mr)




