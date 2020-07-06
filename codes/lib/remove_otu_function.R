# remove_otu_function -
# Input:
#   1) mothur-generated 16S matrix
#   2) mothur-generated taxonomy 
#   3) list of genera to be removed
# Returns an 16S matrix without unwantd OTUs

remove_otu = function(otu_mat, taxonomy, remove_genus){
  otu_to_remove = taxonomy[taxonomy$genus %in% remove_genus, "OTU"]
  remove_otu_mat = otu_mat[,!colnames(otu_mat) %in% otu_to_remove]
  return(remove_otu_mat)
} 
# End remove_otu
