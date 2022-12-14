argv <- commandArgs(TRUE)
start <- as.numeric(argv[1])
end <- as.numeric(argv[2])
# tt = read.table('test_fit_1.txt')
# for (i in 2:33){
	# nb = (i - 1)*50000 + 1
	# name = paste0('test_fit_',nb,'.txt')
	# te = read.table(name)
	# tt = rbind(tt,te)
# }
# fit = unique(tt)
# save(fit,file = 'fit.Rdata')
###################################################
br = read.table('test/bk_region.txt')
###################################################
load('fit.Rdata')
b = fit
#
for (i in start:end){
	#
	mat = as.data.frame(array(,dim = c(1,3)))
	#
	mat[1,1] = i
	#
	pos1 = b[i,7]
	pos2 = b[i,8]
	line1 = which(br[,2] <= pos1)
	line2 = which(br[,3] >= pos1)
	lin = intersect(line1,line2)
	if (length(lin) == 0){
		po1 = 0
		posf1 = line2[1]
	}
	else{
		po1 = lin
	}
	line1 = which(br[,2] <= pos2)
	line2 = which(br[,3] >= pos2)
	lin = intersect(line1,line2)
	if (length(lin) == 0){
		po2 = 0
		posf2 = line2[1]
	}
	else{
		po2 = lin
	}
	po = max(c(po1,po2))
	re1 = po
	if (posf1 != posf2){
		re1 = posf1
	}
	mat[1,2] = re1
	#############################################
	pos1 = b[i,11]
	pos2 = b[i,12]
	line1 = which(br[,2] <= pos1)
	line2 = which(br[,3] >= pos1)
	lin = intersect(line1,line2)
	if (length(lin) == 0){
		po1 = 0
		posf1 = line2[1]
	}
	else{
		po1 = lin
	}
	line1 = which(br[,2] <= pos2)
	line2 = which(br[,3] >= pos2)
	lin = intersect(line1,line2)
	if (length(lin) == 0){
		po2 = 0
		posf2 = line2[1]
	}
	else{
		po2 = lin
	}
	po = max(c(po1,po2))
	re2 = po
	if (posf1 != posf2){
		re2 = posf1
	}
	mat[1,3] = re2
	write.table(mat,'hit_pos_2.txt',append = T,quote =F,row.names =F,col.names =F)
}
###########################################################################################
a = read.table('hit_pos_2.txt')
hp = a[order(a[,1]),]
rownames(hp) = hp[,1]
num1 = which(hp[,2] > 0)
num2 = which(hp[,3] > 0)
num = intersect(num1,num2)
temp = hp[num,]
##########################################################################################
#RCAN1:CHR21:35888740-35987412
#
# pos1 = 35888740
# pos2 = 35987412
# line1 = which(br[,2] <= pos1)
# line2 = which(br[,3] >= pos1)
# lin = intersect(line1,line2)
# if (length(lin) == 0){
	# po1 = 0
# }
# else{
	# po1 = lin
# }
# line1 = which(br[,2] <= pos2)
# line2 = which(br[,3] >= pos2)
# lin = intersect(line1,line2)
# if (length(lin) == 0){
	# po2 = 0
# }
# else{
	# po2 = lin
# }
#

























