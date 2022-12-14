list_string_diff = function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE, only.position = TRUE){
    if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ")
    if(ignore.case){
        a = toupper(a)
        b = toupper(b)
    }
    
    split_seqs = strsplit(c(a, b), split = "")
    only.diff = split_seqs[[1]] != split_seqs[[2]]
    only.diff[
        (split_seqs[[1]] %in% exclude) |
        (split_seqs[[2]] %in% exclude)
    ] = NA
    
    diff.info = data.frame(which(is.na(only.diff)|only.diff), 
                                 split_seqs[[1]][only.diff], split_seqs[[2]][only.diff])
    names(diff.info) = c("position", "seq.a", "seq.b")
    
    if(!show.excluded) diff.info = na.omit(diff.info)
    if(only.position){
        diff.info$position
    }else diff.info
}
