pickup.sequence <- function(num1,num2,ch){
	###num1 refers to start
	###num2 refers to end
	###wrote by chengchao wu###
	###for design in situ primer###
	###############################
	###for test, just using 21 chromosome
	chr.num = ch
	#################################################################################
	#pick up the sequence
	ge = get(paste0('chr',chr.num))
	initial = num1
	final = num2
	centre = floor((initial + final)/2)
	#######################
	width = final - initial
	num = width + 1
	initial.temp = initial%%50
	if (initial.temp == 0){
		initial.line = floor(initial/50) - 1
		initial.point = 50
	}
	else{
		initial.line = floor(initial/50)
		initial.point = initial%%50
	}
	final.temp = final%%50
	if (final.temp == 0){
		final.line = floor(final/50) - 1
		final.point = 50
	}
	else{
		final.line = floor(final/50)
		final.point = final%%50
	}
	if (initial.line == final.line){
		frag = ge[(initial.line + 1)]
		frag = toupper(substr(frag,initial.point,final.point))
	}
	if (initial.line < final.line){
		ge.temp = ge[(initial.line + 1):(final.line + 1)]
		chr.temp = as.character(ge.temp[1])
		for (n in 2:(final.line - initial.line + 1)){
			chr.temp = paste(chr.temp,as.character(ge.temp[n]),sep="")
		}
		ini = initial.point
		fin = ini + width
		frag = toupper(substr(chr.temp,ini,fin))
	}
	return(frag)
}
