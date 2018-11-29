suppressPackageStartupMessages({
    library(plyr)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(Cairo)
    library(gridExtra)
})

# Implementation of HRD scores based on PMID: 26015868

# Copy-number based stats -----------------------------------------------------------------------------------------
calc_cn_stats = function(facets_rdata,
                         genome = 'hg19',
                         algorithm = 'em',
                         sample_name = NULL) {

    # Load data   
    load(facets_rdata)
    if (is.null(sample_name)) sample_name = out$IGV$ID[1]

    # Centromere locations
    chrom_info = switch(genome,
                        hg19 = get_chrom_info('hg19'),
                        hg18 = get_chrom_info('hg18'))

    seg = get_seg(fit, algorithm)
    
    # Create chrom_info for sample
    sample_chrom_info = get_sample_genome(seg, chrom_info)

    # Calculated length of interrogated genome
    interrogated_genome = sum(sample_chrom_info$size)
    autosomal_genome = sum(sample_chrom_info$size[sample_chrom_info$chr %in% 1:22])

    # Check for whole-genome duplication // PMID 30013179
    wgd_treshold = 0.5 # treshold 
    frac_elevated_mcn = sum(seg$length[which(seg$mcn >= 2 & seg$chrom %in% 1:22)])/autosomal_genome
    wgd = frac_elevated_mcn > wgd_treshold
    
    # Calculate fraction of genome altered
    sample_ploidy = ifelse(wgd, round(fit$ploidy), 2)
    if (wgd) {
        diploid_length = sum(seg$length[which(seg$tcn == sample_ploidy & seg$lcn == 1)])
    } else if (!wgd) {
        diploid_length = sum(seg$length[which(seg$mcn == sample_ploidy & seg$lcn >= 1)]) 
    }
    frac_altered = (interrogated_genome-diploid_length)/interrogated_genome
    
    # Altered arms, slightly adjusted to PMID 29622463
    seg = left_join(seg, sample_chrom_info[, c('chr', 'centromere')], by = c('chrom' = 'chr'))
    altered_arms = sapply(unique(seg$chrom), function(x) {
        seg_p = seg[which(seg$chrom == x & seg$start < seg$centromere),]
        seg_p$end = ifelse(seg_p$end > seg_p$centromere, seg_p$centromere, seg_p$end)
        seg_p$length = seg_p$end - seg_p$start
        seg_q = seg[which(seg$chrom == x & seg$end > seg$centromere),]
        seg_q$start = ifelse(seg_q$start < seg_q$centromere, seg_q$centromere, seg_q$start)
        seg_q$length = seg_q$end - seg_q$start
        
        if (wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn == 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn == 1)])
        } else if (!wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$mcn == sample_ploidy & seg_p$lcn >= 1)]) 
            seg_q_unaltered = sum(seg_q$length[which(seg_q$mcn == sample_ploidy & seg_q$lcn >= 1)]) 
        }
        paste0(paste0(x, c('p', 'q'))[c((sample_chrom_info$plength[sample_chrom_info$chr == x]-seg_p_unaltered)/sample_chrom_info$plength[sample_chrom_info$chr == x]>.8,
                                        (sample_chrom_info$qlength[sample_chrom_info$chr == x]-seg_q_unaltered)/sample_chrom_info$qlength[sample_chrom_info$chr == x]>.8) %>% 
                                          map_lgl(~ifelse(is.na(.), FALSE, .))], collapse = ',')
    }) %>% strsplit(., ',') %>% 
        unlist %>% 
        discard(. %in% c('', '13p', '14p', '15p', '21p', '22p'))

    # Weigthed fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>% 
        gather(arm, length, -chr) %>% 
        filter(!paste0(chr, arm) %in% c('13p', '14p', '15p', '21p', '22p')) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms])/sum(length)) %>%
        as.numeric
    
    # # Calculate number of CN events and their length
    # if (wgd) {
    #     cn_segments = rbind(seg[which(seg$tcn != 2),], # any non-diploid
    #         seg[which(seg$tcn == 2 & seg$lcn == 0),]) # unless CN-LOH
    #     cn_segments_no = nrow(cn_segments)
    #     cn_segments_mean_length = mean(cn_segments$end - cn_segments$start)
    #     cn_segments_median_length = median(cn_segments$end - cn_segments$start)
    # } else if (wgd) {
    #     cn_segments = rbind(seg[which(seg$tcn != 4),], # any non-diploid
    #         seg[which(seg$tcn == 4 & seg$lcn != 2),]) # unless CN-LOH
    #     cn_segments_no = nrow(cn_segments)
    #     cn_segments_mean_length = mean(cn_segments$end - cn_segments$start)
    #     cn_segments_median_length = median(cn_segments$end - cn_segments$start)
    # }

    data.frame(
        sample = sample_name,
        wgd = wgd,
        fcna = frac_altered,
        weigthed_fcna = frac_altered_w,
        aneuploidy_score = length(altered_arms),
        altered_arms = paste0(altered_arms, collapse = ',') 
    )
}

# Telomeric allelic imbalance - NtAI ------------------------------------------------------------------------------
# PMID: 22576213

calc_ntai = function(facets_rdata,
                     genome = 'hg19',
                     min_size = 0,
                     min_probes = 250,
                     return_loc = FALSE,
                     return_plot = FALSE,
                     algorithm = 'em',
                     sample_name = NULL) {
    
    # Load data   
    load(facets_rdata)
    if (is.null(sample_name)) sample_name = out$IGV$ID[1]
    
    # Centromere locations
    chrom_info = switch(genome, hg19 = get_chrom_info('hg19'), hg18 = get_chrom_info('hg18'))
    
    # Choose CNCF or EM algorithm
    if (algorithm == 'em') {
        fit$cncf$tcn = fit$cncf$tcn.em
        fit$cncf$lcn = fit$cncf$lcn.em
    }
    seg = fit$cncf
    if (!is.null(fit$start)) {
        seg$start = fit$start
        seg$end = fit$end
    } else {
        seg$start = out$IGV$loc.start
        seg$end = out$IGV$loc.end
    }
    seg$mcn = seg$tcn - seg$lcn
    seg[which(seg$tcn == 1),'mcn'] = 1 # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
    seg[which(seg$tcn == 1),'lcn'] = 0
    
    # Filter segments that do not pass filters
    seg = seg[which(seg$num.mark >= min_probes),]
    seg = seg[which(seg$end - seg$start >= min_size),]

    # Should be TRUE
    
    # Shrink segments with identical copy number
    new_seg = data.frame()
    for(i in unique(seg$chrom)){
        chrom_seg = seg[which(seg$chrom == i),]
        if (nrow(chrom_seg) > 1) { chrom_seg = shrink_seg(chrom_seg) }
        new_seg = rbind(new_seg, chrom_seg)
    }
    seg = new_seg
    
    # Add a column to call AI
    # Codes for telomeric/interstitial/whole chromosome:
    # 0 = no AI
    # 1 = telomeric
    # 2 = interstitial
    # 3 = whole chromosome
    seg$AI = NA
    seg$chrom_ploidy = NA
    sample_ploidy = fit$ploidy # sample ploidy from FACETS
    
    # Loop through chromosomes
    for(i in unique(seg$chrom)){
        
        # Subset on chromosome, proceed to next if no segment
        chrom_seg = seg[which(seg$chrom == i),]
        if(nrow(chrom_seg) == 0){ next }
        
        # Determine major ploidy for chromosome
        ploidy = vector()
        for(k in unique(seg$tcn)){
            tmp = chrom_seg[which(chrom_seg$tcn == k),]
            ploidy = c(ploidy, setNames(sum(tmp$end-tmp$start), k))
            ploidy = ploidy[!names(ploidy) %in% 0] #Remove any ploidy calls of zero
        }
        ploidy = as.numeric(names(ploidy)[which.max(ploidy)])
        chrom_seg$chrom_ploidy = ploidy # update "ploidy" column, so the new calculated value can be returned
        
        if (ploidy%%2 == 0) { # if even
            chrom_seg$AI = c(0,2)[match(chrom_seg$mcn == chrom_seg$lcn, c('TRUE', 'FALSE'))]
            chrom_seg$AI[which(is.na(chrom_seg$AI))] = 0 # adjust NAs
        } else if (ploidy%%2 != 0) { # if odd
            chrom_seg$AI = c(0,2)[match(chrom_seg$mcn + chrom_seg$lcn == ploidy & chrom_seg$lcn != 0, c('TRUE', 'FALSE'))]
            chrom_seg$AI[which(is.na(chrom_seg$AI))] = 0 # adjust NAs
        }
        
        # Put back into original seg
        seg$chrom_ploidy[which(seg$chrom == i)] = ploidy
        seg$AI[which(seg$chrom == i)] = chrom_seg$AI

        # Check relative position to centromere
        if(chrom_seg$AI[1] == 2 & nrow(chrom_seg) != 1 & chrom_seg$end[1] < (chrom_info$centromere[i])){
            seg$AI[which(seg$chrom == i)][1] = 1 # if the first segment of chromosome is AI and does not extend to centromere --> telomeric AI
        }
        if(chrom_seg$AI[nrow(chrom_seg)] == 2 & nrow(chrom_seg) != 1 & chrom_seg$start[nrow(chrom_seg)] > (chrom_info$centromere[i])){
            seg$AI[which(seg$chrom == i)[nrow(chrom_seg)]] = 1 # if the last segment of chromosome is AI and starts beyond the centromere --> telomeric AI
        }
        if(nrow(seg[which(seg$chrom == i),]) == 1 & seg$AI[which(seg$chrom == i)][1] != 0){
            seg$AI[which(seg$chrom == i)[1]] = 3 # if only one segment on chromosome and AI --> chromosomal AI
        }
    }
    
    # Prepare return 
    seg_loh = seg[which(seg$lcn == 0),]
    ntai_out = data.frame(Sample = sample_name,
        NtAI = nrow(seg[which(seg$AI == 1),]), # telomeric AI
        NtAI_Mean_Size = mean(seg$end[which(seg$AI == 1)] - seg$start[which(seg$AI == 1)]),
        NiAI = nrow(seg[which(seg$AI == 2),]), # interstitial AI
        NiAI_Mean_Size = mean(seg$end[which(seg$AI == 2)] - seg$start[which(seg$AI == 2)]),
        NcAI = nrow(seg[which(seg$AI == 3),]), # chromosomal AI
        NtLOH = nrow(seg_loh[which(seg_loh$AI == 1),]), # telomeric LOH
        NtLOH_Mean_Size = mean(seg_loh$end[which(seg_loh$AI == 1)] - seg_loh$start[which(seg_loh$AI == 1)]),
        NiLOH = nrow(seg_loh[which(seg_loh$AI == 2),]), # interstitial LOH
        NiLOH_Mean_Size = mean(seg_loh$end[which(seg_loh$AI == 2)] - seg_loh$start[which(seg_loh$AI == 2)]),
        NcLOH = nrow(seg_loh[which(seg_loh$AI == 3),])) # chromosomal LOH
    
    if (!return_loc & !return_plot) {
        return(ntai_out)
    } else if (return_loc & !return_plot) {
        return(list(ntai_out = ntai_out, ntai_seg = seg))
    } else if (return_loc & return_plot)  {
        print(paste0('Plotting function not implemented'))
        return(list(ntai_out = ntai_out, ntai_seg = seg))
    } else if (!return_loc & return_plot)  {
        print(paste0('Plotting function not implemented'))
        return(ntai_out)
    }
}

# Calculate HRD-LOH score -----------------------------------------------------------------------------------------
# PMID: 23047548
calc_hrdloh = function(facets_rdata,
                       algorithm = 'em',
                       return_loc = FALSE,
                       return_plot = FALSE,
                       sample_name = NULL) {

    # Load data   
    load(facets_rdata)
    if (is.null(sample_name)) sample_name = out$IGV$ID[1]
      
    # Choose CNCF or EM algorithm
    if (algorithm == 'em') {
        fit$cncf$tcn = fit$cncf$tcn.em
        fit$cncf$lcn = fit$cncf$lcn.em
    }
    seg = fit$cncf
    if (!is.null(fit$start)) {
        seg$start = fit$start
        seg$end = fit$end
    } else {
        seg$start = out$IGV$loc.start
        seg$end = out$IGV$loc.end
    }
    seg$mcn = seg$tcn - seg$lcn
    seg[which(seg$tcn == 1),'mcn'] = 1 # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
    seg[which(seg$tcn == 1),'lcn'] = 0

    chr_del =vector()
    for (j in unique(seg$chrom)){
        if (all(!is.na(seg$lcn[which(seg$chrom == j)]) & all(seg$lcn[which(seg$chrom == j)] == 0))) {
            chr_del = c(chr_del, j) # one parental copy of chromosome completely lost
        }
    }
    
    # Check for LOH
    seg_loh = seg[!is.na(seg$lcn),] # remove segments with lcn = NA
    seg_loh = seg_loh[which(seg_loh$lcn == 0 & seg_loh$mcn != 0),]
    seg_loh = seg_loh[which(seg_loh$end - seg_loh$start > 15e6),] # check if segment is long enough
    seg_loh = seg_loh[!which(seg_loh$chrom %in% chr_del),] # remove if whole chromosome lost

    # Mark HRD-LOH segments in original seg
    seg$HRD_LOH = NA
    for (k in 1:nrow(seg_loh)) {
        ii = which(seg$chrom == seg_loh$chrom[k] & seg$start == seg_loh$start[k] & seg$end == seg_loh$end[k])
        seg$HRD_LOH[ii] = TRUE
    }    

    loh_out = data.frame(Sample = sample_name, HRD_LOH = nrow(seg_loh))
    if (!return_loc & !return_plot) {
        return(loh_out)
    } else if (return_loc & !return_plot) {
        return(list(loh_out = loh_out, loh_seg = seg))
    } else if (return_loc & return_plot)  {
        print(paste0('Plotting function not implemented'))
        return(list(ntai_out = ntai_out, ntai_seg = seg))
    } else if (!return_loc & return_plot)  {
        print(paste0('Plotting function not implemented'))
        return(loh_out)
    }
}


# Large-scale state transitions - LST -----------------------------------------------------------------------------
# PMID: 22933060
calc_lst = function(
    facets_rdata,
    genome = 'hg19',
    min_probes = 50,
    lst_length = 10e6,
    algorithm = 'em',
    return_loc = FALSE,
    return_plot = FALSE,
    pdf = FALSE,
    sample_name = NULL)
{
	# Load data	
	load(facets_rdata)
	if (is.null(sample_name)) sample_name = out$IGV$ID[1]

	# Centromere locations
	chrom_info = switch(genome, hg19 = get_chrom_info('hg19'), hg18 = get_chrom_info('hg18'))

	# Choose CNCF or EM algorithm
    if (algorithm == 'em') {
        fit$cncf$tcn = fit$cncf$tcn.em
        fit$cncf$lcn = fit$cncf$lcn.em
    }
    seg = fit$cncf
    if (!is.null(fit$start)) {
        seg$start = fit$start
        seg$end = fit$end
    } else {
        seg$start = out$IGV$loc.start
        seg$end = out$IGV$loc.end
    }
	seg$mcn = seg$tcn - seg$lcn
    seg[which(seg$tcn == 1),'mcn'] = 1 # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
    seg[which(seg$tcn == 1),'lcn'] = 0

	# Return loci of LSTs?
	if (return_loc | return_plot) {
		lst_seg = matrix(0,0,20)
        colnames(lst_seg) = c(colnames(seg), 'LST_breakpoint', 'chr_arm')
	}

	# Count LSTs
	lst = c()
	for (chr in unique(seg$chrom)) {
		chrom_seg = seg[seg$chrom == chr,]
		if (nrow(chrom_seg) < 2) next
		
		# Split into chromosome arms
		p_arm = chrom_seg[chrom_seg$start <= chrom_info$centstart[chrom_info$chr == chr],] # segments starting on p arm, these might overlap centromere
        q_arm = chrom_seg[chrom_seg$end >= chrom_info$centend[chrom_info$chr == chr],] # segments ending on q arm
        q_arm = shrink_seg(q_arm)  # shrink segments with same CN
        p_arm = shrink_seg(p_arm)
        if (nrow(p_arm) > 0) p_arm[nrow(p_arm),'end'] = chrom_info$centstart[chrom_info$chr == chr]  # cut p-arm segment spanning centromere at centromere start
        if (nrow(q_arm) > 0) q_arm[1,'start'] = chrom_info$centend[chrom_info$chr == chr]  # set first q arm segment to start at centromere end

        # P arm
        # Smoothen 3-Mb segments
        n_3mb = which((p_arm$end-p_arm$start) < 3e6)
        while(length(n_3mb) > 0) {
            p_arm = p_arm[-(n_3mb[1]),] # non-juxtaposed segments will be removed, which is fine
            p_arm = shrink_seg(p_arm)
        	n_3mb = which((p_arm$end - p_arm$start) < 3e6)
        }

        # Now check for LST
        if (nrow(p_arm) >= 2) { # if more than one segment
			p_arm = cbind(p_arm, c(0,1)[match((p_arm$end-p_arm$start) >= lst_length, c('FALSE','TRUE'))]) # mark segments that pass length test
            for (k in 2:nrow(p_arm)) {
                if (p_arm[k,ncol(p_arm)] == 1 & p_arm[(k-1),ncol(p_arm)]==1 & (p_arm[k,'start']-p_arm[(k-1),'end']) < 3e6){ # if two juxtaposed segments are 10 Mb and the space between them is less than 3 Mb...
                    lst = c(lst, 1) # ...then add to LST
                    if (return_loc | return_plot) {
                        ## Number indicates if start (1) or end (2) defines the breakpoint
                        a = cbind(p_arm[(k-1),1:18], 2, 'p-arm')
                        b = cbind(p_arm[k,1:18], 1, 'p-arm')
                    	colnames(a)[19:20] = colnames(b)[19:20] = c('LST_breakpoint', 'chr_arm')
                        lst_seg = rbind(lst_seg, a, b)
                    }
                }
            }
        }

        # Q arm
        # Smoothen 3-Mb segments
        n_3mb = which((q_arm$end-q_arm$start) < 3e6)
        while(length(n_3mb) > 0) {
            q_arm = q_arm[-(n_3mb[1]),] # non-juxtaposed segments will be removed, which is fine
            q_arm = shrink_seg(q_arm)
        	n_3mb = which((q_arm$end - q_arm$start) < 3e6)
        }

        # Now check for LST
        if (nrow(q_arm) >= 2) { # if more than one segment
			q_arm = cbind(q_arm, c(0,1)[match((q_arm$end-q_arm$start) >= lst_length, c('FALSE','TRUE'))]) # mark segments that pass length test
            for (k in 2:nrow(q_arm)) {
                if (q_arm[k,ncol(q_arm)] == 1 & q_arm[(k-1),ncol(q_arm)]==1 & (q_arm[k,'start']-q_arm[(k-1),'end']) < 3e6){ # if two juxtaposed segments are 10 Mb and the space between them is less than 3 Mb...
                lst = c(lst, 1) # ...then add to LST
                    if (return_loc | return_plot) {
                        ## Number indicates if start (1) or end (2) defines the breakpoint
                        a = cbind(q_arm[(k-1),1:18], 2, 'q-arm')
                        b = cbind(q_arm[k,1:18], 1, 'q-arm')
                    	colnames(a)[19:20] = colnames(b)[19:20] = c('LST_breakpoint', 'chr_arm')
                        lst_seg = rbind(lst_seg, a, b)
                    }
                }
            }
        }
	}

    # Return values
    lst_out = data.frame(Sample = sample_name, LST = sum(lst))
   
    if (!return_loc & !return_plot) { return(lst_out) }
    if (return_loc) { return(list(lst_out = lst_out, lst_seg = lst_seg)) }
    if (return_plot) {
        if (pdf == TRUE) {
            out_name = paste0(sample_name, '-lst-plot-', algorithm, '.pdf')
            CairoPDF(out_name, w = 12, h = 3)
            plot_lst(out, fit, lst_seg, sum(lst), chrom_info, sample_name)
            dev.off()
            print(paste0('Plot saved to ', out_name))
        } else if (pdf == FALSE) {
            out_name = paste0(sample_name, '-lst-plot-', algorithm, '.png')
            CairoPNG(out_name, w = 1200, h = 300)
            plot_lst(out, fit, lst_seg, sum(lst), chrom_info, sample_name)
            dev.off()
            print(paste0('Plot saved to ', out_name))
        }
    }
}


# Genome-wide LOH  ------------------------------------------------------------------------------------------------
calc_gloh = function(
    facets_rdata,
    genome = 'hg19',
    algorithm = 'em',
    return_loc = FALSE,
    return_plot = FALSE,
    pdf = FALSE,
    sample_name = NULL)
{
   
    # Load data   
    load(facets_rdata)
    if (is.null(sample_name)) sample_name = out$IGV$ID[1]

    # Centromere locations
    chrom_info = switch(genome, hg19 = get_chrom_info('hg19'), hg18 = get_chrom_info('hg18'))

    # Choose CNCF or EM algorithm
    if (algorithm == 'em') {
        fit$cncf$tcn = fit$cncf$tcn.em
        fit$cncf$lcn = fit$cncf$lcn.em
    }
    seg = fit$cncf
    if (!is.null(fit$start)) {
        seg$start = fit$start
        seg$end = fit$end
    } else {
        seg$start = out$IGV$loc.start
        seg$end = out$IGV$loc.end
    }
    seg$mcn = seg$tcn - seg$lcn
    seg[which(seg$tcn == 1),'mcn'] = 1 # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
    seg[which(seg$tcn == 1),'lcn'] = 0
    seg$length = seg$end - seg$start # calculate length of copy-number segments
    seg = subset(seg, chrom < 23) # gLOH only calculated for autosomes

    # Add arm length to chromosome info
    chrom_info$plength = chrom_info$centromere
    chrom_info$qlength = chrom_info$size - chrom_info$centromere

    # Track cumulative lengths
    loh_length = 0 # regions of loh
    exclude_length = 0 # regions of loh that either span 90% of whole chromosome or chromosome arm
    sample_loh_centromere = 0 # loh regions spanning centromere, but less than 90% of whole chromosome and either arm
    seg$gloh = FALSE # mark regions for plotting

    # Create chrom_info for sample
    sample_chrom_info = ldply(unique(seg$chrom), function(x) {
        chr = x
        centromere = chrom_info[which(chrom_info$chr == x), 'centromere']
        chrstart = min(seg[which(seg$chrom == x), 'start'])
        chrend = max(seg[which(seg$chrom == x), 'end'])
        size = chrend
        plength = centromere - chrstart
        qlength = chrend - centromere
        c(chr = chr, centromere = centromere, chrstart = chrstart, chrend = chrend, size = size, plength = plength, qlength = qlength)
    })
    
    # Calculated length of interrogated genome
    interrogated_genome = sum(sample_chrom_info$size)
    autosomal_genome = sum(sample_chrom_info$size[sample_chrom_info$chr %in% 1:22])
    
    # Check for whole-genome duplication // PMID: 30013179
    wgd_treshold = 0.5 # treshold 
    frac_elevated_mcn = sum(seg$length[which(seg$mcn >= 2 & seg$chrom %in% 1:22)])/autosomal_genome
    wgd = frac_elevated_mcn > wgd_treshold

    # Regions of LOH
    # If duplicated genome, lcn = 1 also counts toward LOH
    if (!wgd) {
        seg$loh = with(seg, ifelse(lcn == 0 & mcn > 0, TRUE, FALSE))
    } else if (wgd) {
        seg$loh = with(seg, ifelse(lcn < 2 & mcn > 1, TRUE, FALSE))
    }
    seg$loh[is.na(seg$loh)] = FALSE
    
    # Loop through segments, check exclusion criteria
    if (nrow(seg) > 0) {
        for (i in 1:nrow(seg)) {
            if (seg$loh[i] == FALSE) { next }
            chrom = seg[i, 'chrom'] # segment is on this chromosome
            seg_start = seg[i, 'start'] # starts here
            seg_end = seg[i, 'end'] # ends here
            seg_length = seg[i, 'length'] # is this long
            chrom_size = sample_chrom_info[which(sample_chrom_info$chr == chrom), 'size'] # chrom size
            chrom_centromere = sample_chrom_info[which(sample_chrom_info$chr == chrom), 'centromere'] # centromere position
            p_arm = sample_chrom_info[which(sample_chrom_info$chr == chrom), 'plength'] # length of p arm
            q_arm = sample_chrom_info[which(sample_chrom_info$chr == chrom), 'qlength'] # length of q arm

            # Calculate span of LOH region
            if (seg_length >= 0.9 * chrom_size) { # exclude if >= 90% of chromosome length
                exclude_length = exclude_length + seg_length
            } else if (seg_end <= chrom_centromere & seg_length < 0.9 * p_arm) { # on p arm, less than 90% of length
                loh_length = loh_length + seg_length
                seg$gloh[i] = TRUE
            } else if (seg_start >= chrom_centromere & seg_length < 0.9 * q_arm) { # ditto for q arm
                loh_length = loh_length + seg_length
                seg$gloh[i] = TRUE
            } else if (seg_start < chrom_centromere & seg_end > chrom_centromere) { # spans centromere
                p_length = chrom_centromere - seg_start
                q_length = seg_end - chrom_centromere
                if (p_length < 0.9 * p_arm) { # check if length of part on p arm less than 90% of p arm
                    loh_length = loh_length + p_length
                    seg$gloh[i] = TRUE
                } else { exclude_length = exclude_length + p_length }
                if (q_length < 0.9 * q_arm) { # ditto for q arm
                    loh_length = loh_length + q_length
                    seg$gloh[i] = TRUE
                } else { exclude_length = exclude_length + q_length }
            } 
        }
    }

    # Return values
    gloh_out = data.frame(Sample = sample_name,
        gLOH = loh_length/(interrogated_genome-exclude_length),
        wgd = wgd,
        gLOH_No = length(which(seg$loh == TRUE)),
        gLOH_Median_Length = median(seg$length[which(seg$loh == TRUE)]),
        gLOH_Mean_Length = mean(seg$length[which(seg$loh == TRUE)]))

    if (!return_loc & !return_plot) { return(gloh_out) }
    if (return_loc) { return(list(list(gloh_out = gloh_out, gloh_seg = seg))) }
    if (return_plot) {
        if (pdf == TRUE) {
            out_name = paste0(sample_name, '-gloh-plot-', algorithm, '.pdf')
            CairoPDF(out_name, w = 12, h = 6)
            plot_gloh(out, seg, gloh_out, chrom_info, sample_name)
            dev.off()
            print(paste0('Plot saved to ', out_name))
        } else if (pdf == FALSE) {
            out_name = paste0(sample_name, '-gloh-plot-', algorithm, '.png')
            CairoPNG(out_name, w = 1200, h = 600)
            plot_gloh(out, seg, gloh_out, chrom_info, sample_name)
            dev.off()
            print(paste0('Plot saved to ', out_name))
        }
    }
}


# Helper functions ------------------------------------------------------------------------------------------------

# Get segmentation from fit object
get_seg = function(fit, algorithm) {
    if (algorithm == 'em') {
        fit$cncf$tcn = fit$cncf$tcn.em
        fit$cncf$lcn = fit$cncf$lcn.em
    }
    
    seg = fit$cncf
    seg$length = seg$end-seg$start
    seg[which(seg$tcn <= 1),'lcn'] = 0 # correct mcn, lcn for cases of tcn = 1 // sometimes FACETS set lcn = NA when tcn = 1, when it clearly has to be 0
    seg$mcn = seg$tcn - seg$lcn
    seg
}

# Sample-specific genome info
get_sample_genome = function(seg, chrom_info) {
    ldply(unique(seg$chrom), function(x) {
        chr = x
        centromere = chrom_info[which(chrom_info$chr == x), 'centromere']
        chrstart = min(seg[which(seg$chrom == x), 'start'])
        chrend = max(seg[which(seg$chrom == x), 'end'])
        size = chrend-chrstart
        plength = centromere - chrstart
        if (plength < 0) plength = 0 # acrocentric chromosomes
        qlength = chrend - centromere
        c(chr = chr,
          centromere = centromere,
          chrstart = chrstart,
          chrend = chrend,
          size = size,
          plength = plength,
          qlength = qlength)
    })
}

# Shrink segments with identical copy number that are separated for some reason
# Logic: if different tcn ==> don't shrink
# if any NAs in mcn/lcn ==> don't shrink
# if identical mcn/lcn content ==> do shrink
shrink_seg = function(chrom_seg) {
    if (nrow(chrom_seg) < 2) {
        chrom_seg
    } else {
        new_chr = chrom_seg
        seg_class = c(1)
        for(j in 2:nrow(new_chr)){
        	# if adjacent segments have same allelic content, assign to same class
        	if (new_chr[(j-1),'tcn'] != new_chr[j,'tcn']) { seg_class = c(seg_class, seg_class[j-1]+1) }# if tcn differs, definitely don't condense
        	else if (any(is.na(c(new_chr[(j-1),'mcn'], new_chr[(j),'mcn'])))) { seg_class = c(seg_class, seg_class[j-1]+1) } # but if tcn the same, check for difference in allelic content
            else if (new_chr[(j-1),'mcn'] == new_chr[j,'mcn'] & new_chr[(j-1),'lcn'] == new_chr[j,'lcn']) { seg_class = c(seg_class, seg_class[j-1]) }
            else { seg_class = c(seg_class, seg_class[j-1]+1) }
        }
        for(j in unique(seg_class)){
        	# condense segments belonging to same class
            new_chr[seg_class %in% j,'end'] = max(new_chr[seg_class %in% j,'end'])
            new_chr[seg_class %in% j,'num.mark'] = sum(new_chr[seg_class %in% j,'num.mark'])
            }
        new_chr[!duplicated(seg_class),]
    }
}

### Read in gzip file from url
get_gz_from_url = function(url, ...) {
    # http://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
    con = gzcon(url(url))
    txt = readLines(con)
    dat = read.delim(textConnection(txt), ...)
    return(dat)
    close.connectio(con)
}

### Get chromosome info for given genome assembly
get_chrom_info = function(genome = 'hg19') {
	# gaps = genome
	# if (genome == 'hg19') { genome = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz'
	# } else if (genome == 'hg18') { genome = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz' } 
	
 #    # Get chromInfo table from UCSC
 #    chrom = get_gz_from_url(genome, header = FALSE)
 #    chrom = subset(chrom, grepl("^chr[0-9XY]{1,2}$", chrom[,1]))
 #    # Get gap table from UCSC
 #    if (gaps == 'hg19') {
 #    	gaps = get_gz_from_url('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz', header = FALSE)
 #    	centro = subset(gaps, gaps[,8] == "centromere")
	# } else if (gaps == 'hg18') {
	# 	gaps = ldply(c(seq(1:22), 'X', 'Y'), function(x) {
	# 			get_gz_from_url(paste0('http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr', x, '_gap.txt.gz'), header = F)
	# 		})
	# 	centro = subset(gaps, gaps[,8] == "centromere")
	# 	centro = centro[!duplicated(centro$V2),]
	# }
    
 #    # Merge the relevant info from the two tables
 #    chrominfo = merge(chrom[,1:2], centro[,2:4], by.x = 1, by.y = 1) # merge size and centromere location
 #    chrominfo$centromere = rowMeans(chrominfo[,3:4]) # convert centromere start and end into one location (the mean)
 #    chrominfo = chrominfo[,c(1,2,5,3,4)] # keep chromosome, size and centromere location
 #    colnames(chrominfo) = c("chr", "size", "centromere", "centstart", "centend")
 #    chrominfo[,1] = as.character(chrominfo[,1])
 #    chrominfo$chr = sub("chr", "", chrominfo$chr)
 #    chrominfo$chr = sub("X", "23", chrominfo$chr)
 #    chrominfo$chr = sub("Y", "24", chrominfo$chr)
 #    chrominfo[,1] = as.numeric(chrominfo[,1])
 #    chrominfo = chrominfo[order(chrominfo$chr), ]  # order by chromosome number
 #    rownames(chrominfo) = as.character(chrominfo[,1])
 #    chrominfo = as.matrix(chrominfo)
 #    return(chrominfo)
    
 	if (genome == 'hg19') {
    	return(data.frame(chr = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
    		size = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566),
    		centromere = c(123035434,93826171,92004854,51160117,47905641,60330166,59554331,45338887,48867679,40754935,53144205,36356694,17500000,17500000,18500000,36835801,23763006,16960898,26181782,27869569,12788129,14500000,60132012,11604553),
	    	centstart = c(121535434,92326171,90504854,49660117,46405641,58830166,58054331,43838887,47367679,39254935,51644205,34856694,16000000,16000000,17000000,35335801,22263006,15460898,24681782,26369569,11288129,13000000,58632012,10104553),
    		centend = c(124535434,95326171,93504854,52660117,49405641,61830166,61054331,46838887,50367679,42254935,54644205,37856694,1.9e+07,1.9e+07,2e+07,38335801,25263006,18460898,27681782,29369569,14288129,1.6e+07,61632012,13104553)))
		} else if (genome == 'hg18') {
			return(data.frame(chr = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
			size = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954),
			centromere = c(122356957,93189898,92037544,50854874,47941398,60438125,59558273,45458052,48607499,40434941,52950781,35445461,16934000,16570000,16760000,36043302,22237133,16082897,28423622,27150399.5,11760000,12830000,60098737,11453954),
			centstart = c(121236957,91689898,90587544,49354874,46441398,58938125,58058273,43958052,47107499,39244941,51450781,34747961,1.6e+07,15070000,15260000,35143302,22187133,15400898,26923622,26267569,10260000,11330000,58598737,11253954),
			centend = c(123476957,94689898,93487544,52354874,49441398,61938125,61058273,46958052,50107499,41624941,54450781,36142961,17868000,18070000,18260000,36943302,22287133,16764896,29923622,28033230,13260000,14330000,61598737,11653954)))
		}	
}

# Plot functions --------------------------------------------------------------------------------------------------
# Get chromosome position, cumulative over whole genome
get_cumulative_chr_maploc_mod = function(mat, chrom_info) {
    chrom_lengths = chrom_info$size
    cum_chrom_lengths = cumsum(as.numeric(chrom_lengths))
    mid = cum_chrom_lengths - (chrom_lengths/2)
    names(mid) = 1:22
    chr_maploc_to_gen_maploc = function(x) { mat[mat$chrom==x,]$maploc + cum_chrom_lengths[x-1] }
    chr_maploc = sapply(2:22, chr_maploc_to_gen_maploc)
    chr_maploc = unlist(chr_maploc)
    chr_maploc = c(mat[mat$chrom==1,]$maploc, chr_maploc)
    mat = cbind(mat[mat$chrom < 23,], chr_maploc)
    list(mat=mat, mid=mid)
}

# Plot HRD-specific CN plot
plot_lst = function(out,
                    fit,
                    lst_seg,
                    lst,
                    chrom_info,
                    sample_name) {
    
    col1 ="lightgrey"
    col2 = "darkgrey"
    
    mat = out$jointseg
    mat = subset(mat, chrom < 23)
    mat = get_cumulative_chr_maploc_mod(mat, chrom_info)
    mid = mat$mid
    mat = mat$mat
    cncf = fit$cncf
    cncf = subset(cncf, chrom < 23)
    dipLogR = out$dipLogR
    cnlr_median = rep(cncf$cnlr.median, cncf$num.mark)
    mat = cbind(mat, cnlr_median)

    starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
    ends = cumsum(c(cncf$num.mark))
    my_starts = mat[starts,c('chr_maploc','cnlr_median')]
    my_ends = mat[ends,c('chr_maploc','cnlr_median')]

    centr_pos = rbind(centromere_start = get_cumulative_chr_maploc_mod(data.frame(chrom = chrom_info$chr, maploc = chrom_info$centstart), chrom_info)$mat,
        centromere_end = get_cumulative_chr_maploc_mod(data.frame(chrom = chrom_info$chr, maploc = chrom_info$centend), chrom_info)$mat)

    col_rep = 1 + rep(mat$chrom - 2 * floor(mat$chrom/2))

    # Tease out LST breakpoints
    if (nrow(lst_seg) > 0) {
        break_pos = rbind(
            data.frame(chrom = lst_seg[lst_seg$LST_breakpoint==1,]$chrom, maploc = lst_seg[lst_seg$LST_breakpoint==1,]$start),
            data.frame(chrom = lst_seg[lst_seg$LST_breakpoint==2,]$chrom, maploc = lst_seg[lst_seg$LST_breakpoint==2,]$end)
        )
        break_pos = break_pos[order(break_pos[,1], break_pos[,2], decreasing = F),]
        break_pos = do.call('rbind', lapply(seq(1, nrow(break_pos), by = 2), function(x) {
            data.frame(chrom = break_pos$chrom[x], maploc = mean(break_pos[x:x+1, 'maploc']))
        }))
        break_pos = get_cumulative_chr_maploc_mod(break_pos, chrom_info)
    } else { break_pos = NULL }

    # Where is the allelic ratio undetermined?
    cncf_na = which(is.na(cncf$lcn))
    
    # Plot
    cnlr = ggplot(mat, environment = environment()) +
        geom_point(aes(y = cnlr, x = chr_maploc), colour = c(col1, col2)[col_rep], size = .8) +
        scale_x_continuous(breaks = mid[!is.na(names(mid))], labels = names(mid)[!is.na(names(mid))]) +
        xlab('') +
        ylim(-3,3) +
        ylab('Log-Ratio') +
        geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) +       
        theme(axis.text.x  = element_text(angle=90, vjust=0),
              axis.text.y = element_text(angle=90, vjust=0),
              text = element_text(size=8)) +
        ggtitle(paste0(sample_name, '\nLST ', lst)) +
        geom_segment(x = 0, xend = 10e6, y = 2.5, yend = 2.5, col = 'black') + # show scale for 10 Mb segment
        geom_point(data = centr_pos, aes(x = chr_maploc, y = dipLogR), col = 'darkgreen', size = 1) + # add centromere location
        annotate(geom = 'text', x = 10e6/2, y = 2.6, label = '10 Mb', hjust = .5, vjust = 0, size = 2) +
        theme_bwmin()
    if (length(cncf_na) == 0) {
        cnlr = cnlr + geom_segment(data=cncf,aes(x=my_starts$chr_maploc, xend=my_ends$chr_maploc, y=my_starts$cnlr_median, yend=my_ends$cnlr_median), col='darkorchid3', size=1)
    } else {
        cnlr = cnlr + 
        geom_segment(data=cncf[-cncf_na,],aes(x=my_starts$chr_maploc[-cncf_na], xend=my_ends$chr_maploc[-cncf_na], y=my_starts$cnlr_median[-cncf_na], yend=my_ends$cnlr_median[-cncf_na]), col='darkorchid3', size=1) +
        geom_segment(data=cncf[cncf_na,],aes(x=my_starts$chr_maploc[cncf_na], xend=my_ends$chr_maploc[cncf_na], y=my_starts$cnlr_median[cncf_na], yend=my_ends$cnlr_median[cncf_na]), col='cyan3', size=1)
    }
    if (!is.null(break_pos)) { cnlr = cnlr + geom_vline(xintercept = break_pos$mat$chr_maploc, color = 'red', lty = 'dashed') }
    print(cnlr)
}

# Plot HRD-specific CN plot
plot_gloh = function(out,
                     seg,
                     gloh_out,
                     chrom_info,
                     sample_name) {
    
    col1 ="lightgrey"
    col2 = "darkgrey"
    
    mat = out$jointseg
    mat = subset(mat, chrom < 23)
    mat = get_cumulative_chr_maploc_mod(mat, chrom_info)
    mid = mat$mid
    mat = mat$mat
    mat = subset(mat, chrom < 23)
    cncf = seg
    cncf = subset(cncf, chrom < 23)
    dipLogR = out$dipLogR
    cnlr_median = rep(cncf$cnlr.median, cncf$num.mark)
    mat = cbind(mat, cnlr_median)
    mafR = rep(sqrt(abs(cncf$mafR)), cncf$num.mark)
    mat = cbind(mat, mafR)

    starts = cumsum(c(1,cncf$num.mark))[1:length(cncf$num.mark)]
    ends = cumsum(c(cncf$num.mark))
    my_starts = mat[starts,c('chr_maploc','cnlr_median', 'mafR')]
    my_ends = mat[ends,c('chr_maploc','cnlr_median', 'mafR')]

    # Mark centromeres
    centr_pos = rbind(centromere_start = get_cumulative_chr_maploc_mod(data.frame(chrom = chrom_info$chr, maploc = chrom_info$centstart), chrom_info)$mat,
        centromere_end = get_cumulative_chr_maploc_mod(data.frame(chrom = chrom_info$chr, maploc = chrom_info$centend), chrom_info)$mat)

    col_rep = 1 + rep(mat$chrom - 2 * floor(mat$chrom/2))

    # Log-ratio plot
    cnlr = ggplot(mat, environment = environment()) +
        geom_point(aes(y = cnlr, x = chr_maploc), colour = c(col1, col2)[col_rep], size = .8) +
        scale_x_continuous(breaks = mid[!is.na(names(mid))], labels = names(mid)[!is.na(names(mid))]) +
        xlab('') +
        ylim(-3,3) +
        ylab('Log-Ratio') +
        geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) +
        theme(axis.text.x  = element_text(angle=90, vjust=0),
              axis.text.y = element_text(angle=90, vjust=0),
              text = element_text(size=8)) +
        geom_point(data = centr_pos, aes(x = chr_maploc, y = dipLogR), col = 'darkgreen', size = 1) + # add centromere location
        theme_bwmin()
    if (any(cncf$gloh) == T) {
        cnlr = cnlr + geom_segment(data=cncf[which(cncf$gloh==T),],aes(x=my_starts$chr_maploc[which(cncf$gloh==T)], xend=my_ends$chr_maploc[which(cncf$gloh==T)], y=my_starts$cnlr_median[which(cncf$gloh==T)], yend=my_ends$cnlr_median[which(cncf$gloh==T)]), col='red', size=1) +
        geom_segment(data=cncf[which(cncf$gloh==F),],aes(x=my_starts$chr_maploc[which(cncf$gloh==F)], xend=my_ends$chr_maploc[which(cncf$gloh==F)], y=my_starts$cnlr_median[which(cncf$gloh==F)], yend=my_ends$cnlr_median[which(cncf$gloh==F)]), col='darkorchid3', size=1)
    } else {
        cnlr = cnlr + geom_segment(data=cncf,aes(x=my_starts$chr_maploc, xend=my_ends$chr_maploc, y=my_starts$cnlr_median, yend=my_ends$cnlr_median), col='darkorchid3', size=1)
    }
    
    # Log-odds ratio plot
    valor = ggplot(mat, environment = environment()) +
        geom_point(aes(y=valor,x=chr_maploc), colour=c(col1,col2)[col_rep], size=.4) +
        scale_x_continuous(breaks = mid[!is.na(names(mid))], labels = names(mid)[!is.na(names(mid))]) +
        xlab('') +
        ylim(-4,4) +
        ylab('Variant allele log odds ratio') +
        theme(axis.text.x = element_text(angle=90, vjust=0, size=8),
              axis.text.y = element_text(angle=90, vjust=0, size=8),
              text = element_text(size=10),
              panel.grid.minor.x=element_line(colour='white', size=.5),
              panel.grid.major.x=element_line(colour='white', size=0)) + 
        theme_bwmin() + geom_vline(xintercept = my_starts$chr_maploc[which(!duplicated(cncf$chrom))], col = 'grey', lty = 'dashed')
    if (any(cncf$gloh) == T) {
        valor = valor + geom_segment(data=cncf[which(cncf$gloh==F),], aes(x=my_starts$chr_maploc[which(cncf$gloh==F)], xend=my_ends$chr_maploc[which(cncf$gloh==F)], yend=my_ends$mafR[which(cncf$gloh==F)], y=my_starts$mafR[which(cncf$gloh==F)]), col='darkorchid3', size=1, lineend='butt') +
        geom_segment(data=cncf[which(cncf$gloh==F),], aes(x=my_starts$chr_maploc[which(cncf$gloh==F)], xend=my_ends$chr_maploc[which(cncf$gloh==F)], yend=-my_ends$mafR[which(cncf$gloh==F)], y=-my_starts$mafR[which(cncf$gloh==F)]), col='darkorchid3', size=1, lineend='butt') +
        geom_segment(data=cncf[which(cncf$gloh==T),], aes(x=my_starts$chr_maploc[which(cncf$gloh==T)], xend=my_ends$chr_maploc[which(cncf$gloh==T)], yend=my_ends$mafR[which(cncf$gloh==T)], y=my_starts$mafR[which(cncf$gloh==T)]), col='red', size=1, lineend='butt') +
        geom_segment(data=cncf[which(cncf$gloh==T),], aes(x=my_starts$chr_maploc[which(cncf$gloh==T)], xend=my_ends$chr_maploc[which(cncf$gloh==T)], yend=-my_ends$mafR[which(cncf$gloh==T)], y=-my_starts$mafR[which(cncf$gloh==T)]), col='red', size=1, lineend='butt')
    } else {
        valor = valor + geom_segment(data=cncf, aes(x=my_starts$chr_maploc, xend=my_ends$chr_maploc, yend=my_ends$mafR, y=my_starts$mafR), col='darkorchid3', size=1, lineend='butt') +
        geom_segment(data=cncf, aes(x=my_starts$chr_maploc, xend=my_ends$chr_maploc, yend=-my_ends$mafR, y=-my_starts$mafR), col='darkorchid3', size=1, lineend='butt')
    }
    
    # Concatenate plots
    grid.arrange(cnlr, valor, ncol = 1, top = grid::textGrob(paste0(sample_name, '\ngLOH ', signif(gloh_out$gLOH, 2)*100, '% | WGD=', gloh_out$WGD), gp=grid::gpar(fontsize=18, fontfamily = 'Helvetica', font=1)))
}

# Minimalistic ggplot2 theme in black and white
theme_bwmin = function (base_size = 12, base_family = "ArialMT") 
{
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(panel.grid.major = element_blank(), panel.border = element_rect(fill = NA, 
            colour = "black", size = 1), panel.grid.minor = element_blank(), 
            strip.background = element_rect(fill = NA, colour = NA))
}
