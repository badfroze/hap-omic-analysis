#
br = read.table('test/bk_region.txt')
source('pickup.sequence.R')
chr = read.table('/home/ccwu/data1/homo_sapiens/chr21.fa')
chr21 = chr[-1,]
for (i in 1:nrow(br)){
	pos1 = br[i,2]
	num1 = pos1 - 22
	num2 = pos1 + 22
	seq1 = pickup.sequence(num1,num2)
	n = 0
	for (j in 1:(nchar(seq1) - 3)){
		tem = substr(seq1,j,(j+3))
		if (tem == 'TTAA'){
			n = n + 1
		}
		if (tem == 'AATT'){
			n = n + 1
		}
	}
	##################################
	pos2 = br[i,3]
	num3 = pos2 - 22
	num4 = pos2 + 22
	seq2 = pickup.sequence(num3,num4)
	m = 0
	for (j in 1:(nchar(seq2) - 3)){
		tem = substr(seq2,j,(j+3))
		if (tem == 'TTAA'){
			m = m + 1
		}
		if (tem == 'AATT'){
			m = m + 1
		}
	}
	##################################
	seq3 = pickup.sequence(pos1,pos2)
	k = 0
	for (j in 1:(nchar(seq3) - 3)){
		tem = substr(seq2,j,(j+3))
		if (tem == 'TTAA'){
			k = k + 1
		}
		if (tem == 'AATT'){
			k = k + 1
		}
	}
	
	
	##################################
	fin = paste(br[i,1],pos1,pos2,seq1,seq2,n,m,k)
	write.table(fin,'MseI_pos.txt',append =T,quote =F,row.names= F,col.names =F)
}



