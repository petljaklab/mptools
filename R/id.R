complement = c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

#' @export
indel_contexts = function(var, genome = BSgenome.Hsapiens.UCSC.hg19){
	vars = var[,c("CHROM", "POS", "REF", "ALT")]
	vars[,length:=nchar(ALT) - nchar(REF)]
	vars = vars[abs(length) == 1]
	vars[,affected_base:=""]
	vars[length>0,affected_base:=substr(ALT, 2, 2)]
	vars[length<0,affected_base:=substr(REF, 2 ,2)]
	vars[,type:="ins"]
	vars[length<0,type:="del"]
	vars[,comp_base:=affected_base]
	vars[affected_base %in% c("G", "A"),comp_base:=complement[affected_base]]
	vars[,strand:="+"]
	vars[affected_base %in% c("G", "A"),strand:="-"]
	vars[,start:=POS-1]
	vars[,end:=POS+3]
	vars[strand == "-",start:=POS-1]
	vars[strand == "-",end:=POS+3]
	vars[,CHROM:=paste0("chr", CHROM)]
	gr = makeGRangesFromDataFrame(vars[,c("CHROM", "start", "end")], seqnames.field = "CHROM")
	vars[,seq := as.character(BSgenome::getSeq(genome, gr))]
	vars[,comp_seq:=seq]
	vars[strand == "-",comp_seq:=comp_seq(seq)]
	vars[,mut_type:=""]
	vars[type == "ins" & comp_base == "T", mut_type:="insT"]
	vars[type == "ins" & comp_base == "C", mut_type:="insC"]
	vars[type == "del" & comp_base == "T", mut_type:="delT"]
	vars[type == "del" & comp_base == "C", mut_type:="delC"]
	vars[,pentanuc:=comp_seq]
	vars[,comp_seq:=substr(pentanuc, 2, 4)]
	vars[,p5:=substr(comp_seq, 1, 1)]
	vars[,p3:=substr(comp_seq, 3, 3)]
	return(vars[])
}

comp_seq = function(seq){
	split = strsplit(seq, split = "")
	unlist(lapply(split, function(x){
		paste0(complement[rev(x)], collapse = "")
	}))
}

#' @export
plot_indels = function(indat, matrix = F, genome = BSgenome.Hsapiens.UCSC.hg19){
	indat = indat[nchar(REF) != nchar(ALT) & (nchar(REF) == 1 | nchar(ALT) == 1)]
	indat[,typ:="insertion"]
	indat[nchar(REF) > nchar(ALT),typ:="deletion"]
	indat[,size:=abs(nchar(REF)-nchar(ALT))]
	## We have a number of categories to deal with
	### 1bp del/ins of a C/T, homopolymer length 1/0-6/5+
	### 2-5+bp del/ins of 2/1-6/5+ repeats
	### 2-5+bp del of 1-5+ microhomologies
	## Firstly we do the 1bp indels
	ones = indat[size == 1]
	if(nrow(ones) > 0){
		ones[,modded:=ifelse(nchar(REF) > nchar(ALT), yes = substr(REF, 2, 2), no = substr(ALT, 2, 2))]
		ones[,up:=POS-6]
		ones[,down:=POS+6]
		ones[,chrchrom:=paste0("chr", CHROM)]
		gr = makeGRangesFromDataFrame(df = ones, seqnames.field = "chrchrom", start.field = "up", end.field = "down")
		seq_context  = as.character(getSeq(genome, gr))
		ones[,seqs:=seq_context]
		homopolymers = function(search_seq, base){
			## First get a vector of the sequence
			search_seq = unlist(strsplit(search_seq, ""))
			pos = 8
			direction = 1
			hom_length = 0
			done = F
			while(!done){
				query_seq = search_seq[pos]
				if(query_seq == base){
					hom_length = hom_length + 1
				}else if(query_seq != base){
					## We didn't find a continuation of the homopolymer
					if(direction == 1){
						direction = -1
						pos = 8
					}else{
						done = T
					}
				}
				if(pos == length(search_seq)){
					direction = -1
					pos = 8
				}
				pos = pos + direction
			}
			return(hom_length)
		}
		ones[,homp_len:=homopolymers(seqs, modded), by = 1:nrow(ones)]
	}
	## Now we look at >1bp indels
	## Strategy: Assume everything is a repeat indel
	## Then, we take the 1-rep indels and check them for microhomologies
	large = indat[size > 1]
	nlarge = nrow(large)
	large[,modded_seq:=ifelse(nchar(REF) > nchar(ALT), yes = substr(REF, 2, nchar(REF)), no = substr(ALT, 2, nchar(ALT)))]
	large[,chrchrom:=paste0("chr", CHROM)]
	large[,up:=POS - size + 1]
	large[,down:=POS+size]
	large[,pos2:=POS]
	large[,pos3:=POS+1]
	large[,upseq:="N"]
	large[,downseq:="N"]
	large[,reps:=0]
	large[,c("left","right"):=0]
	n_repeats = 0
	# Do the initial iteration out here, because we need to look both up and down. Afterward we set a direction and iteratively look in that direction+
	## Get upstream and downstream sequences
	gr = makeGRangesFromDataFrame(df = large, seqnames.field = "chrchrom", start.field = "up", end.field = "pos2")
	large[,upseq:=as.character(getSeq(genome, gr))]
	gr = makeGRangesFromDataFrame(df = large, seqnames.field = "chrchrom", start.field = "pos3", end.field = "down")
	large[,downseq:=as.character(getSeq(genome, gr))]
	## I think this should never trigger, but put this here to be sure
	if(nrow(large[modded_seq == upseq & modded_seq == downseq]) > 0){
		stop("Edge case detected - improve the program")
	}
	## Decide if we're looking upstream or downstream for more repeats
	large[,direction:=ifelse(modded_seq == upseq, yes = -1, no = 1)]
	## Append the non-repetitive seqs to this table, filter down table
	large_d = large[0]
	large_d = rbind(large_d, large[!(modded_seq == upseq | modded_seq == downseq)])
	large = large[(modded_seq == upseq | modded_seq == downseq)]
	## Define left and right bounds of new intervals
	large[,left:=ifelse(direction > 0, yes = pos3 + size, no = up - size)]
	large[,right:=ifelse(direction > 0, yes = down+size, no = pos2 - size)]
	n_repeats = 1
	while(nrow(large) > 1 & n_repeats < 6){
		large[,reps:=n_repeats]
		gr = makeGRangesFromDataFrame(df = large, seqnames.field = "chrchrom", start.field = "left", end.field = "right")
		large[,downseq:=as.character(getSeq(genome, gr))]
		large_d = rbind(large_d, large[!(modded_seq == upseq | modded_seq == downseq)])
		large = large[(modded_seq == upseq | modded_seq == downseq)]
		large[,left:=left+(size*direction)]
		large[,right:=right+(size*direction)]
		n_repeats = n_repeats + 1
	}
	large[,reps:=n_repeats]
	large_d = rbind(large_d, large)
	## Now we process the microhomology deletions
	microhom_candidates = large_d[typ == "deletion" & reps == 1]
	n_microhom_candidates = nrow(microhom_candidates)
	## If l = length(del), we begin with length= l-1 microhomologies
	microhom_d = microhom_candidates[0]
	micro_len = 1
	while(nrow(microhom_candidates) > 0){
		microhom_candidates[,micro_search:=micro_len]
		microhom_candidates[,check_del:=substr(modded_seq, 1, nchar(modded_seq) - micro_len)]
		microhom_candidates[,check_upseq:=substr(upseq, 1 + micro_len, nchar(upseq))]
		microhom_candidates[,check_downseq:=substr(downseq, 1, nchar(downseq) - micro_len )]
		microhom_d = rbind(microhom_d, microhom_candidates[check_del == check_upseq | check_del == check_downseq], fill = TRUE)
		microhom_candidates = microhom_candidates[!(check_del == check_upseq | check_del == check_downseq)]
		micro_len = micro_len + 1
	}
	if(nrow(microhom_d) > 0){
		microhom_d[,micro_size:=nchar(check_del)]
		microhom_failed = microhom_d[nchar(check_del) == 0]
		microhom_d = microhom_d[nchar(check_del) > 0]
		microhom_d[,discrete_len:=size]
		microhom_d[discrete_len >=5,discrete_len:=5]
		microhom_d[,discrete_micro:=micro_size]
		microhom_d[discrete_micro >= 5,discrete_micro:=5]
		microhom_sum = microhom_d[order(discrete_len, discrete_micro)][discrete_micro > 0,.(n = .N), by = c("discrete_len", "discrete_micro")]
		microhom_sum[,name:=paste0("micro", discrete_len, discrete_micro)]
		setkey(microhom_sum, name)
		ref_microhom = data.table("Var1" = "micro", "Var2" = c(rep(2, 1), rep(3, 2), rep(4, 3), rep(5,5)), "Var3" = c(1, 1:2, 1:3, 1:5))
		ref_microhom[,name:=paste0(Var1, Var2, Var3)]
		setkey(ref_microhom, name)
		microhom_sum = microhom_sum[,c("name", "n")][ref_microhom[,"name"]]
	}else{
		ref_microhom = data.table("Var1" = "micro", "Var2" = c(rep(2, 1), rep(3, 2), rep(4, 3), rep(5,5)), "Var3" = c(1, 1:2, 1:3, 1:5))
		ref_microhom[,name:=paste0(Var1, Var2, Var3)]
		ref_microhom[,n:=NA]
		microhom_sum = ref_microhom[,c("name", "n")]
		microhom_failed = data.table()
	}
	## Now summarize it all into a vector
	### 1bpdel[C,T], 1bpins[C,T]
	
	ref_ones = rbind(data.table(expand.grid(c("del"), c("C", "T"), 1:6)), data.table(expand.grid(c("ins"), c("C", "T"), 0:5)))
	ref_ones[,name:=paste0(Var1, Var2, Var3)]
	ref_ones = ref_ones[,"name"]
	setkey(ref_ones, name)
	if(nrow(ones) > 0){
		ones[,pmodded:=modded]
		ones[modded == "G",pmodded:="C"]
		ones[modded == "A",pmodded:="T"]
		ones[typ == "insertion" & homp_len == 6,homp_len:=5]
		ones_sum = ones[order(typ, pmodded, homp_len)][,.(n = .N), by = c("typ", "pmodded", "homp_len")]
		ones_sum[,name:=paste0(substr(typ, 1, 3), pmodded, homp_len)]
		setkey(ones_sum, name)
		ones = ones_sum[,c("name", "n")][ref_ones]
	}else{
		ones = ref_ones
		ones$n = NA
	}
	
	
	### 2bpdel[rep1-6+], 3bpdel[rep1-6+], 4bpdel[rep1-6+],5+bpdel[rep1-6+], then insertions, same thing
	large_d = rbind(large_d[!(typ == "deletion" & reps == 1)], microhom_failed, fill = T)
	large_d[,discrete_len:=size]
	large_d[discrete_len >=5,discrete_len:=5]
	large_d[typ == "insertion" & reps == 6,reps:=5]
	d_sum = large_d[order(typ, discrete_len, reps)][,.(n = .N), by = c("typ", "discrete_len", "reps")]
	d_sum[,name:=paste0("rep", substr(typ, 1, 3), discrete_len, reps)]
	setkey(d_sum, name)
	ref_large = rbind(data.table(expand.grid("repdel", 2:5, 1:6)), data.table(expand.grid("repins", 2:5, 0:5)))
	ref_large[,name:=paste0(Var1, Var2, Var3)]
	ref_large = ref_large[,"name"][order(name)]
	setkey(ref_large, name)
	d_sum = d_sum[,c("name", "n")][ref_large]
	
	## MicrohomD[2]_len1, MicrohomD[3]_len1-2, MicrohomD[4]_len1-3, Microhom[5+]_len1-5+, 
	indel_catalog = rbind(ones, d_sum, microhom_sum)
	indel_catalog[is.na(n),n:=0]
	indel_catalog[,name:=factor(name, levels = name)]
	if(sum(indel_catalog$n) != nrow(indat)){
		stop("You lost or gained indels during processing")
	}
	colors = c(rep("#FEC984", 6), 
						 rep("#FF9100", 6), 
						 rep("#BCDF9C", 6), 
						 rep("#43AC3D", 6), 
						 rep("#FED3C1", 6), 
						 rep("#FF9A7F", 6),
						 rep("#F55D43", 6),
						 rep("#CF2C1E", 6),
						 rep("#D9E6F3", 6),
						 rep("#A4CEE7", 6),
						 rep("#5CA6D2", 6),
						 rep("#1678BA", 6),
						 rep("#E8E6F3", 1),
						 rep("#C4C3E1", 2),
						 rep("#9799C7", 3),
						 rep("#7557AA", 5))
	if(matrix){
		return(indel_catalog)
	}
	base_plot = ggplot(indel_catalog, aes(x = name, y = n)) + geom_bar(stat = "identity", fill = colors) + theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")) + scale_y_continuous(expand = c(0,0))
	blocks = seq(from = 1, to = 73, by = 6)
	blocks = blocks - 0.5
	blocks2 = blocks + 6
	blocks = c(blocks, 73.5, 75.5, 78.5)
	blocks2 = c(blocks2[-13], 73.5, 75.5, 78.5, 83.5)
	blockdf = data.frame(x1 = blocks, x2 = blocks2, colors = unique(colors), labs = c("C","T", "C", "T", "2", "3", "4", "5+", "2", "3", "4", "5+", "2", "3", "4", "5+"))
	boxes = ggplot(blockdf, aes(xmin = x1, xmax = x2, ymin=0, ymax = max(indel_catalog$n))) + geom_rect(fill = unique(colors)) + theme_cowplot() + ylab("n") + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	#boxes + geom_fit_text(grow = F)
	boxes_top = boxes + geom_fit_text(aes(label = labs), grow = F)
	
	annotations_df = data.frame(x1 = c(0, 12.5, 24.5, 48.5, 72.5), x2 = c(12.5, 24.5, 48.5, 72.5, 83), label = c("\n1bp Insertion", "\n1bp Deletion", "2+bp Deletion\n(Deletion Length)", "2+bp Insertion\n(Insertion Length)", "Microhomology Del\n(Deletion Length)"))
	annos = ggplot(annotations_df, aes(xmin = x1, xmax = x2, ymin = 0, ymax = max(indel_catalog$n), label = label)) + geom_rect(fill = "#FF000000") + theme_cowplot() + geom_fit_text(grow = F) + ylab("n") + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	bottext = data.frame(x = 1:83, lab = c(rep(c(1:5, "6+"), times = 2), rep(c(0:4, "5+"), times = 2), rep(c(1:5, "6+"), times = 4), rep(c(0:4, "6+"), times = 4), 1, 1:2, 1:3, 1:4, "5+"))
	bot_box = ggplot(bottext, aes(xmin = x-0.5, xmax = x+0.5, ymin = 0, ymax = max(indel_catalog$n), label = lab)) + geom_text() + theme_cowplot()+ ylab("n") + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	#anno2 = ggplot(annotations_df, aes(x = (x1+x2)/2, y = 0, label = label)) + geom_text() + theme_cowplot() + ylab("n") + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	p = plot_grid(annos, boxes_top, base_plot, boxes, ncol = 1, rel_heights = c(0.2, 0.15, 1, 0.15))
	
	#p = ggplot(indel_catalog, aes(x = name, y = n)) + geom_bar(stat = "identity", fill = colors) + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_y_continuous(expand = c(0,0))
	return(p)
}

#' @export
plot_frequencies_indel = function(indel_catalog, sample_name, ylab = "N"){
	colors = c(rep("#FEC984", 6), 
						 rep("#FF9100", 6), 
						 rep("#BCDF9C", 6), 
						 rep("#43AC3D", 6), 
						 rep("#FED3C1", 6), 
						 rep("#FF9A7F", 6),
						 rep("#F55D43", 6),
						 rep("#CF2C1E", 6),
						 rep("#D9E6F3", 6),
						 rep("#A4CEE7", 6),
						 rep("#5CA6D2", 6),
						 rep("#1678BA", 6),
						 rep("#E8E6F3", 1),
						 rep("#C4C3E1", 2),
						 rep("#9799C7", 3),
						 rep("#7557AA", 5))
	indel_catalog = indel_catalog[,.SD, .SDcols = c("MutationType", sample_name)]
	indel_catalog[,MutationType:=factor(MutationType, levels = unique(MutationType))]
	setnames(indel_catalog, old = c("MutationType", sample_name), new = c("MutationType", "y"))
	base_plot = ggplot(indel_catalog, aes(x = MutationType, y = y)) + geom_bar(stat = "identity", fill = colors) + theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")) + ylab(ylab) +  scale_y_continuous(expand = c(0,0))
	blocks = seq(from = 1, to = 73, by = 6)
	blocks = blocks - 0.5
	blocks2 = blocks + 6
	blocks = c(blocks, 73.5, 75.5, 78.5)
	blocks2 = c(blocks2[-13], 73.5, 75.5, 78.5, 83.5)
	blockdf = data.frame(x1 = blocks, x2 = blocks2, colors = unique(colors), labs = c("C","T", "C", "T", "2", "3", "4", "5+", "2", "3", "4", "5+", "2", "3", "4", "5+"))
	boxes = ggplot(blockdf, aes(xmin = x1, xmax = x2, ymin=0, ymax = max(indel_catalog$y))) + geom_rect(fill = unique(colors)) + theme_cowplot() + ylab(ylab) + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	#boxes + geom_fit_text(grow = F)
	boxes_top = boxes + geom_fit_text(aes(label = labs), grow = F)
	
	annotations_df = data.frame(x1 = c(0, 12.5, 24.5, 48.5, 72.5), x2 = c(12.5, 24.5, 48.5, 72.5, 83), label = c("\n1bp Insertion", "\n1bp Deletion", "2+bp Deletion\n(Deletion Length)", "2+bp Insertion\n(Insertion Length)", "Microhomology Del\n(Deletion Length)"))
	annos = ggplot(annotations_df, aes(xmin = x1, xmax = x2, ymin = 0, ymax = max(indel_catalog$y), label = label)) + geom_rect(fill = "#FF000000") + theme_cowplot() + geom_fit_text(grow = F) + ylab(ylab) + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	bottext = data.frame(x = 1:83, lab = c(rep(c(1:5, "6+"), times = 2), rep(c(0:4, "5+"), times = 2), rep(c(1:5, "6+"), times = 4), rep(c(0:4, "6+"), times = 4), 1, 1:2, 1:3, 1:4, "5+"))
	bot_box = ggplot(bottext, aes(xmin = x-0.5, xmax = x+0.5, ymin = 0, ymax = max(indel_catalog$y), label = lab)) + geom_text() + theme_cowplot()+ ylab(ylab) + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	#anno2 = ggplot(annotations_df, aes(x = (x1+x2)/2, y = 0, label = label)) + geom_text() + theme_cowplot() + ylab("n") + theme(axis.line = element_line(color = "#FF000000"), axis.text = element_text(color = "#FF000000"), axis.ticks = element_line(color = "#FF000000"), axis.title = element_text(color = "#FF000000"), plot.margin = unit(c(0,0,0,0), "cm")) + scale_x_continuous(expand = c(0,0))
	
	p = plot_grid(annos, boxes_top, base_plot, boxes, ncol = 1, rel_heights = c(0.2, 0.15, 1, 0.15))
	return(p)
}
