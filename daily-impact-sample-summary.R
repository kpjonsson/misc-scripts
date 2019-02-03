#!/opt/common/CentOS_6-dev/R/R-3.4.1/bin/Rscript
suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(data.table)
    library(janitor)
    library(kpjmisc)
})

fread_ = function(x) { suppressMessages(suppressWarnings(read_tsv(x, comment = '#', progress = F))) }
fwrite_ = function(fn, p) { fwrite(fn, p, quote = F, sep = '\t', col.names = T) }

args = commandArgs(TRUE)
if (args[1] %in% c('-h', '-help', '--h', '--help')) {
    message(paste0('Usage: daily_impact_samples_summary.R YYYY-MM-DD',
                   'If no date give, will compare today\'s pull to yesterday\'s',
                   'Otherwise compares the IMPACT repo from given date and the day before.')
            )
    q()
}

if (!is.na(args[1])) {
    message(args[1])
    input_date = ymd(str_replace_all(args[1], '-', '_'))
    yesterday = input_date-1
} else {
    input_date = today()
    yesterday = today()-1
}

# Load old and new list of samples, get the difference
old_samples = fread_(paste0('~/log/impact_logs/sample-list-', yesterday,'.tsv')) %>% 
    clean_names
clin = suppressMessages(read_impact_samples())

if (input_date == today()) {
    fwrite_(clin, paste0('~/log/impact_logs/sample-list-', today(),'.tsv'))
    new_samples = clean_names(clin) %>% 
        filter(!sample_id %in% old_samples$sample_id)
} else {
    fwrite_(clin, paste0('~/log/impact_logs/sample-list-', today(),'.tsv'))
    new_samples = fread_(paste0('~/log/impact_logs/sample-list-', input_date,'.tsv')) %>% 
        clean_names %>% 
        filter(!sample_id %in% old_samples$sample_id)
}

message(paste('New samples:', nrow(new_samples)))
if (nrow(new_samples) == 0) {
    q()
}

maf = suppressMessages(read_impact_maf(germline = T)) %>% 
    filter(Tumor_Sample_Barcode %in% new_samples$sample_id)
cnas = suppressMessages(read_impact_cna()) %>% 
    filter(sample_id %in% new_samples$sample_id)
rearr = fread('~/res/dmp/mskimpact/data_SV.txt') %>% 
    clean_names %>% 
    filter(sample_id %in% new_samples$sample_id)

# Check for new gliomas -------------------------------------------------------------------------------------------
gliomas = filter(new_samples, tolower(cancer_type) %like% 'Glioma') %>% 
    select(patient_id, sample_id, cancer_type, cancer_type_detailed)

if (nrow(gliomas) > 0) { 
    new_gliomas = T
    new_gliomas_url = portal_link(gliomas$patient_id)
} else {
    new_gliomas = F
}

# Check for new BRCA1/2-mutated samples ---------------------------------------------------------------------------
gml = filter(maf,
             Hugo_Symbol %in% c('BRCA1', 'BRCA2'),
             Mutation_Status == 'GERMLINE')

som = filter(maf,
             Hugo_Symbol %in% c('BRCA1', 'BRCA2'),
             Mutation_Status != 'GERMLINE')

cnas = filter(cnas,
              Hugo_Symbol %in% c('BRCA1', 'BRCA2'),
              cna < 0)

rearr = filter(rearr,
               site1_gene %in% c('BRCA1', 'BRCA2') |
               site2_gene %in% c('BRCA1', 'BRCA2'))

brca = filter(new_samples, sample_id %in% c(gml$Tumor_Sample_Barcode,
                                            som$Tumor_Sample_Barcode,
                                            cnas$sample_id,
                                            rearr$sample_id))

if (nrow(brca) > 0) { 
    new_brca = T
    new_brca_url = portal_link(brca$patient_id)
} else {
    new_brca = F
}

# Compose e-mail summary ------------------------------------------------------------------------------------------
# Compose mail
glioma_string = ifelse(
    new_gliomas==F,
    'There were no new glioma samples\n',
    paste0('There were ', nrow(gliomas), ' new glioma samples\n',
           'Link to cBioPortal:\n',
           new_gliomas_url, '\n')
    )
brca_string = ifelse(
    new_brca==F,
    'There were no new BRCA1/2-mutated samples\n',
    paste0('There were ', nrow(brca), ' BRCA1/2-mutated samples\n',
           'Link to cBioPortal:\n',
           new_brca_url))

mail_body = paste0(
    'Number of new IMPACT samples: ', nrow(new_samples), '\n\n',
    '***BRCA/Gray***\n',
    brca_string,
    '\n',
    '***Glioma***\n',
    glioma_string
)

subject_line = paste0('MSK-IMPACT update ', input_date)
system(paste0("echo '", mail_body ,"' | mail -s '", subject_line ,"' 'jonssonp@mskcc.org'"))
