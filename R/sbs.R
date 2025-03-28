sbs_palette = c("C>A" = "#23c9f2",
								"C>G" = "#000000",
								"C>T" = "#eb4030",
								"T>A" = "#d4d3d3",
								"T>C" = "#afd476",
								"T>G" = "#f2d3d1")   

trinucs_function = function(){
	complement = c("A" = "T", "C" = "G", "T" = "A", "G" = "C")
	mutated_base = c("C", "T")
	trinucs = rep("", times = 96)
	it = 1
	for(i in 1:4){
		base1 = complement[i]
		for(j in 1:2){
			base2 = mutated_base[j]
			poss = complement[complement!=base2]
			for(k in 1:4){
				base3 = complement[k]
				for(l in 1:3){
					change = poss[l]
					trinucs[it] = paste0(base1, base2, base3, ">", change)
					it = it+1
				}
			}
		}
	}
	return(trinucs)
}

#' @export
annotate_contexts = function(dat, genome = BSgenome.Hsapiens.UCSC.hg19){
	trinucs = trinucs_function()
	trinucs_blank = data.table(trinuc = tstrsplit(trinucs, ">")[[1]], change = tstrsplit(trinucs, ">")[[2]])
	dat = dat[nchar(REF) == 1 & nchar(ALT) == 1]
	if(substr(dat$CHROM[1], 1, 3) != "chr"){
		dat[,CHROM := paste0("chr", CHROM)]
	}
	dat[,up:=POS-1]
	dat[,down:=POS+1]
	gr = GenomicRanges::makeGRangesFromDataFrame(df = dat, seqnames.field = "CHROM", start.field = "POS", end.field = "POS")
	mutpos = BSgenome::getSeq(genome, gr)
	flip = which(data.table::data.table(data.frame(mutpos))$mutpos %in% c("G", "A"))
	dat[,strand:="+"]
	dat[flip,strand:="-"]
	gr = GenomicRanges::makeGRangesFromDataFrame(df = dat, seqnames.field = "CHROM", start.field = "up", end.field = "down")
	trinucs = BSgenome::getSeq(genome, gr)
	dat$trinuc = as.character(trinucs)
	dat[,comp_alt:=ALT]
	dat[strand == "-",comp_alt:=complement[ALT]]
	return(dat[])
}

#' @export
plot_variants = function(dat, genome = BSgenome.Hsapiens.UCSC.hg19){
	trinucs = trinucs_function()
	trinucs_blank = data.table(trinuc = tstrsplit(trinucs, ">")[[1]], change = tstrsplit(trinucs, ">")[[2]])
	dat = dat[nchar(REF) == 1 & nchar(ALT) == 1]
	dat = annotate_contexts(dat, genome = genome)
	catalog = dat[,c("trinuc", "comp_alt")]
	catalog = catalog[,.N, by = c("trinuc", "comp_alt")]
	catalog[order(comp_alt, trinuc)]
	
	setkey(catalog, trinuc, comp_alt)
	## Fill in to 96
	catalog = catalog[trinucs_blank]
	catalog[,ref:=substr(trinuc, 2, 2)]
	catalog = catalog[order(ref, comp_alt, trinuc)]
	catalog[is.na(N)]$N = 0
	catalog[,xticks:=paste0(trinuc, comp_alt)]
	catalog[,xticks:=factor(xticks, levels = xticks)]
	catalog[,change:=paste0(ref, ">", comp_alt)]
	p = ggplot(catalog, aes(x = xticks, y = N, fill = change)) + geom_bar(stat = "identity", width = 0.6) + scale_fill_manual(values = sbs_palette) + 
		scale_x_discrete(labels = catalog$trinuc) + 
		theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
		xlab("Trinucleotide Context")
	return(p)
}

#' @export
plot_frequencies_sbs = function(sbs_catalog, sample_name, ylab = "N"){
	sbs_catalog = sbs_catalog[,.SD, .SDcols = c("MutationType", sample_name)]
	names(sbs_catalog)[2] = ylab
	sbs_catalog[,change:=gsub("[A-Z]\\[(.*)\\][A-Z]", "\\1", MutationType)]
	sbs_catalog[,trinuc:=gsub("([A-Z])\\[([A-Z])>[A-Z]\\]([A-Z])", "\\1\\2\\3", MutationType)]
	sbs_catalog[,ref:=gsub("([A-Z])\\[([A-Z])>[A-Z]\\]([A-Z])", "\\2", MutationType)]
	sbs_catalog[,comp_alt:=gsub("([A-Z])\\[([A-Z])>([A-Z])\\]([A-Z])", "\\3", MutationType)]
	sbs_catalog = sbs_catalog[order(ref, comp_alt, trinuc)]
	sbs_catalog[,MutationType:=factor(MutationType, levels = MutationType)]
	p = ggplot(sbs_catalog, aes(x = MutationType, y = N, fill = change)) + geom_bar(stat = "identity", width = 0.6) + scale_fill_manual(values = sbs_palette) + 
		#scale_x_discrete(labels = catalog$trinuc) + 
		theme_bw() + theme(axis.text.x = element_blank(), legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
		xlab("Trinucleotide Context")
	return(p)
}
