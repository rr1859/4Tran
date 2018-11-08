valid_primer = function(chr, start, end, strand){
	if(as.character(strand) == "+"){
		primer_valid = NULL
		primer2_valid = which(primer2[,1] == chr & primer2[,2] > (start-1000) & primer2[,2] < start & primer2[,6] == "-")
		if(length(primer2_valid) > 0){
			primer2_start = primer2[primer2_valid[1],2]
			secondary_enz_site = tail(secondary_enz[secondary_enz[,1] == chr & secondary_enz[,2] < primer2_start,],1)
			no_secondary_enz = NULL
			no_secondary_enz = which(secondary_enz[,1] == chr & secondary_enz[,2] > primer2_start & secondary_enz[,2] < start)
			no_primary_enz = NULL
			no_primary_enz = which(primary_enz[,1] == chr & primary_enz[,2] < start & primary_enz[,2] > secondary_enz_site[,2])
			if(length(no_primary_enz) == 0 & length(no_secondary_enz) == 0) 
				return (1)
			else
				return(0)
		}
		else
			return(0)	
	}
	if(as.character(strand) == "-"){
		primer_valid = NULL
		primer2_valid = which(primer2[,1] == chr & primer2[,2] < (start+1000) & primer2[,2] > start & primer2[,6] == "+")
		if(length(primer2_valid) >0){
			primer2_start = primer2[primer2_valid[1],2]
			secondary_enz_site = head(secondary_enz[secondary_enz[,1] == chr & secondary_enz[,2] > primer2_start,],1)
			no_secondary_enz = NULL			
			no_secondary_enz = which(secondary_enz[,1] == chr & secondary_enz[,2] < primer2_start & secondary_enz[,2] > start)
			no_primary_enz = NULL
			no_primary_enz = which(primary_enz[,1] == chr & primary_enz[,2] > start & primary_enz[,2] < secondary_enz_site[,2])
			if(length(no_primary_enz) == 0 & length(no_secondary_enz) == 0)
				return (1)
			else
				return(0)
		}
		else
			return(0)
	}
}
primary_enz = read.table("/scratch/rr1859/Genome_files/mm10/mm10_hind3_restriction_sites_oligomatch.bed",stringsAsFactors=FALSE)
secondary_enz = read.table("/scratch/rr1859/Genome_files/mm10/mm10_csp6_restriction_sites_oligomatch.bed",stringsAsFactors = FALSE)
primer_names = read.table("new_primer_names_unique.txt")
for(primer in primer_names[,1]){
	print(primer)
	if(file.exists(paste(primer,"_2_mm10.bed",sep = "")))
		primer1 = read.table(paste(primer,"_2_mm10.bed",sep = ""), stringsAsFactors = FALSE)
	else
		next
	if(file.exists(paste(primer,"_1_mm10.bed",sep = "")))
		primer2 = read.table(paste(primer,"_1_mm10.bed",sep = ""),stringsAsFactors = FALSE)
	else
		next
	output = NULL
	output = unlist(lapply(1:nrow(primer1), function(i) valid_primer(primer1[i,1],primer1[i,2],primer1[i,3],primer1[i,6])))
	if(length(which(output == 1)) > 0)
		write.table(primer1[which(output == 1),], paste(primer,"_valid_mm10.bed", sep = ""), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
	else
		print("no valid primers")
}
