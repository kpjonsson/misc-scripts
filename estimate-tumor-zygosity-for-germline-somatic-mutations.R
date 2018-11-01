library(binom)

### Calculate expected allele frequency of mutation given major/minor copy number and tumor purity estimate
expected_vaf_germline = function(purity, mcn, lcn) {

    tcn = mcn + lcn # total number of alleles in tumor
    tumor_normal_alleles = purity*tcn + (1-purity)*2 # total number of alleles in sample

    mcn_mut_alleles = (purity*mcn + (1-purity)) # number of mutated alleles in tumor if mcn, plus mut alleles in normal cells
    lcn_mut_alleles = (purity*lcn + (1-purity)) # number of mutated alleles in tumor if lcn, plus mut alleles in normal cells

    c(vaf_mcn = mcn_mut_alleles/tumor_normal_alleles, vaf_lcn = lcn_mut_alleles/tumor_normal_alleles)
}

expected_vaf_somatic = function(purity, mcn, lcn) {

    tcn = mcn + lcn # total number of alleles in tumor
    tumor_normal_alleles = purity*tcn + (1-purity)*2 # total number of alleles in sample

    mcn_mut_alleles = purity*mcn # number of mutated alleles in tumor if mcn, plus mut alleles in normal cells
    lcn_mut_alleles = purity*lcn # number of mutated alleles in tumor if lcn, plus mut alleles in normal cells

    c(vaf_mcn = mcn_mut_alleles/tumor_normal_alleles, vaf_lcn = lcn_mut_alleles/tumor_normal_alleles)
}

# Calculate binomial proportion CI
# calc_ci = function(t_alt_count, t_depth, end) {

#     p_hat = t_alt_count/t_depth
#     c = 1/(1+(1.96^2)/t_depth)
#     d = p_hat + (1.96^2)/(2*t_depth)
#     e = 1.96 * sqrt(p_hat*(1-p_hat)/t_depth + (1.96^2)/(4*t_depth^2))

#     if (end == 0) { c*(d-e) } else { c*(d+e) }
# }

### From Alex Penson
expected_mutant_copies_in_tumor = function(VAF, purity, tcn, germline = FALSE){
  mutant_tumor_copies <- VAF * (2*(1 - purity) + purity * tcn)
  if(germline == TRUE) mutant_tumor_copies <- mutant_tumor_copies - (1 - purity)
  mutant_tumor_copies <- mutant_tumor_copies / purity
}

estimate_loh_ci = function(purity, mcn, lcn, t_ref_count, t_alt_count, germline = T, n_vaf = NULL) {

    t_depth = t_ref_count+t_alt_count
    t_vaf = t_alt_count/t_depth

    if (germline & n_vaf %nin% c(NULL, NA)) {
        ref_bias = n_vaf / 0.5
        if (ref_bias > 1) ref_bias = 1
    } else { # when somatic or germline but n_vaf not given
        n_vaf = 0.5
        ref_bias = 1
    }

    if (t_depth < 10) { # Depth threshold
        expected_vaf=NA
        expected_vaf_lower_ci=NA
        expected_vaf_upper_ci=NA
        ref_cn=NA
        alt_cn=NA
        observed_vaf_lower_ci=NA
        observed_vaf_upper_ci=NA
        cn_vaf_concordance=NA
        allelic_imbalance=NA
        loss_of_heterozygosity=NA
    } else {

        # Get expected VAFs given mcn/lcn
        if (germline) {
            expected_vafs = calc_expected_vaf_germline(purity, mcn, lcn)
        } else {
            expected_vafs = calc_expected_vaf_somatic(purity, mcn, lcn)
        }
        vaf_mcn = expected_vafs[['vaf_mcn']]
        vaf_lcn = expected_vafs[['vaf_lcn']]

        # Distance between observed and expected VAFs
        delta_vaf_mcn = abs(t_vaf-vaf_mcn)
        delta_vaf_lcn = abs(t_vaf-vaf_lcn)

        # Calculate CI for observed VAF
        observed_vaf_ci = binom.confint(t_alt_count, t_depth, methods = 'wilson')
        observed_vaf_upper_ci = observed_vaf_ci[['upper']]
        observed_vaf_lower_ci = observed_vaf_ci[['lower']]

        # Determine which is the mutant allele is more likely to be on mcn/lcn
        if(delta_vaf_lcn < delta_vaf_mcn) {
            expected_vaf = vaf_lcn
            # expected_vaf_lower_ci = calc_ci(expected_vaf*t_depth, t_depth, 0)
            # expected_vaf_upper_ci = calc_ci(expected_vaf*t_depth, t_depth, 1)
            expected_vaf_ci = binom.confint(expected_vaf * t_depth * ref_bias, t_depth, methods = 'wilson')
            expected_vaf_upper_ci = expected_vaf_ci[['upper']]
            expected_vaf_lower_ci = expected_vaf_ci[['lower']]
            ref_cn = mcn
            alt_cn = lcn
        } else  {
            expected_vaf = vaf_mcn
            # expected_vaf_lower_ci = calc_ci(expected_vaf*t_depth, t_depth, 0)
            # expected_vaf_upper_ci = calc_ci(expected_vaf*t_depth, t_depth, 1)
            expected_vaf_ci = binom.confint(expected_vaf * t_depth * ref_bias, t_depth, methods = 'wilson')
            expected_vaf_upper_ci = expected_vaf_ci[['upper']]
            expected_vaf_lower_ci = expected_vaf_ci[['lower']]
            ref_cn = lcn
            alt_cn = mcn
        }

        # boundary_vaf = 0
        # boundary_vaf_lower_ci = 0
        # boundary_vaf_upper_ci = 0
        # if (t_vaf > 0.5) {
        #     if (germline) {
        #         boundary_vaf = calc_expected_vaf_germline(purity, 2, 1)[['vaf_mcn']]
        #     } else {
        #         boundary_vaf = calc_expected_vaf_somatic(purity, 2, 1)[['vaf_mcn']]
        #     }
        #     # boundary_vaf_lower_ci = calc_ci(t_depth*boundary_vaf, t_depth, 0)
        #     # boundary_vaf_upper_ci = calc_ci(t_depth*boundary_vaf, t_depth, 1)
        #     boundary_vaf_ci = binom.confint(t_depth*boundary_vaf, t_depth, methods = 'wilson')
        # } else {
        #     if (germline) {
        #         boundary_vaf = calc_expected_vaf_germline(purity, 2, 1)[['vaf_lcn']]
        #     } else {
        #         boundary_vaf = calc_expected_vaf_somatic(purity, 2, 1)[['vaf_lcn']]
        #     }
        #     # boundary_vaf_lower_ci = calc_ci(t_depth*boundary_vaf, t_depth, 0)
        #     # boundary_vaf_upper_ci = calc_ci(t_depth*boundary_vaf, t_depth, 1)
        #     boundary_vaf_ci = binom.confint(t_depth*boundary_vaf, t_depth, methods = 'wilson')
        # }

        ### Given purity and CN state, does the observed and expected VAFs agree?
        # cn_vaf_concordance = F
        # if ((t_vaf > 0.5 & t_vaf > boundary_vaf_ci$upper) | (t_vaf <= 0.5 & t_vaf < boundary_vaf_ci$lower))
        # if (observed_vaf_lower_ci <= expected_vaf_upper_ci && expected_vaf_lower_ci <= observed_vaf_upper_ci) {
        #     cn_vaf_concordance = TRUE
        # }
        cn_vaf_concordance = between(t_vaf, expected_vaf_lower_ci, expected_vaf_upper_ci)

        imbalance_treshold = ifelse(germline, n_vaf, purity/2)
        ### Is the tumor non-diploid?
        # if (t_vaf >= 0.5) {
        if (t_vaf > imbalance_treshold) {
            if (germline) {
                imbalanced_vaf_ci = binom.confint(calc_expected_vaf_germline(purity, 2, 1)[['vaf_mcn']] * t_depth * ref_bias,
                    t_depth,
                    methods = 'wilson')
            } else {
                imbalanced_vaf_ci = binom.confint(calc_expected_vaf_somatic(purity, 2, 1)[['vaf_mcn']] * t_depth * ref_bias,
                    t_depth,
                    methods = 'wilson')
                # print(imbalanced_vaf_ci)
            }
            if (imbalanced_vaf_ci$lower <= t_vaf) {
                allelic_imbalance = "ALT_GAIN"
            } else { allelic_imbalance = 'N/A' }
        } else {
            if (germline) {
                imbalanced_vaf_ci = binom.confint(calc_expected_vaf_germline(purity, 2, 1)[['vaf_lcn']] * t_depth * ref_bias,
                    t_depth,
                    methods = 'wilson')
            } else {
                imbalanced_vaf_ci = binom.confint(calc_expected_vaf_somatic(purity, 2, 1)[['vaf_lcn']] * t_depth * ref_bias,
                    t_depth,
                    methods = 'wilson')
            }
            if (imbalanced_vaf_ci$upper >= t_vaf){
                allelic_imbalance = "REF_GAIN"
            } else { allelic_imbalance = 'N/A' }
        }

        ### If non-diploid and ref CN == 0, then LOH w.r.t to EITHER allele
        if (cn_vaf_concordance &
            (allelic_imbalance == 'ALT_GAIN' && ref_cn == 0) |
            (allelic_imbalance=='REF_GAIN' && alt_cn == 0)) {
            loss_of_heterozygosity = TRUE
        } else { loss_of_heterozygosity = FALSE}
        # if ((t_vaf > 0.5 & observed_vaf_ci$upper > expected_vaf_ci$lower) | (t_vaf <= 0.5 & observed_vaf_ci$lower < expected_vaf_ci$upper)) {
        #     do_intervals_overlap = 1
        # }
    }

    data.frame(expected_vaf=expected_vaf, expected_vaf_lower_ci=expected_vaf_lower_ci, expected_vaf_upper_ci=expected_vaf_upper_ci,
      ref_cn=ref_cn, alt_cn=alt_cn,
      observed_vaf_lower_ci=observed_vaf_lower_ci, observed_vaf_upper_ci=observed_vaf_upper_ci,
      cn_vaf_concordance=cn_vaf_concordance,
      allelic_imbalance=allelic_imbalance,
      loss_of_heterozygosity=loss_of_heterozygosity
      )
}
