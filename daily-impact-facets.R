#!/opt/common/CentOS_6-dev/R/R-3.4.1/bin/Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(lubridate)
    library(dplyr)
    library(stringr)
    library(glue)
    library(jsonlite)
})

args = commandArgs(TRUE)

if (args[1] %in% c('-h', '-help', '--h', '--help')) {
    message("Usage: daily_impact_facets.R YYYY-MM-DD\nIf no date give, will use latest.")
    q()
}

if (!is.na(args[1])) {
    message(args[1])
    input_date = ymd(str_replace_all(args[1], '-', '_'))
} else {
    input_date = today()-1
}
# Read bam file key
facets_dir = '/ifs/res/taylorlab/impact_facets/all'
bam_dir = '/ifs/dmpshare/share/irb12_245'
key_base = '/ifs/res/taylorlab/dmp_mirror_key/{date}/dmp_key_paired_{date}.txt'
message(paste('Finding new IMPACT samples from BAM file keys, comparing', input_date, 'to', input_date-1), '.')
new_key = fread(glue(key_base, date = str_replace_all(input_date, '\\-', '\\_')))
old_key = fread(glue(key_base, date = str_replace_all(input_date-1, '\\-', '\\_')))

# Find new pairs and unpaired tumor samples
new_pairs = filter(new_key,
                   !tumor_sample %in% old_key$tumor_sample,
                   !is.na(normal_sample)) %>% 
    mutate_at(vars(matches('bamname')), funs(str_c(bam_dir, '/', ., '.bam')))

unpaired = filter(new_key,
                  !tumor_sample %in% old_key$tumor_sample,
                  is.na(normal_sample))
message(paste('Found:\n', nrow(new_pairs), 'new tumor-normal pairs\n', nrow(unpaired), 'new unmatched tumors'))

if (nrow(unpaired) > 0) {
    write(paste(unpaired$tumor_sample, collapse = '\n'),
          file = paste0(facets_dir, '/log/', 'unpaired_tumor_samples.txt'), append = T)
}

# Run Facets on new pairs
if (nrow(new_pairs) == 0) {
    q()
}
    
facets_cmd = paste(
    'cmoflow_facets --R_lib 0.5.6 --normal-bam {normal_bam} --tumor-bam {tumor_bam}',
    '--normal-name {normal_name} --tumor-name {tumor_name} --output-dir {dir}',
    '--cval 50 --purity_cval 150 --min_nhet 15 --purity_min_nhet 15 --seed 100'
)

log = select(new_pairs, tumor_sample, normal_sample) %>% 
    mutate(run_status = '',
           run_command = '',
           counts_file = '')

for (i in 1:nrow(new_pairs)) {
    output_dir = str_c(facets_dir, '/', new_pairs$tumor_sample[i], '_', new_pairs$normal_sample[i])
    sample_facets = glue(facets_cmd,
                         normal_bam = new_pairs$normal_bamname[i],
                         tumor_bam = new_pairs$tumor_bamname[i],
                         normal_name = new_pairs$normal_sample[i],
                         tumor_name = new_pairs$tumor_sample[i],
                         dir = output_dir) 
    frun = system(sample_facets, intern = T)
    log$run_command[i] = sample_facets
    
    if (is.null(attr(frun, 'status'))) attr(frun, 'status') = 0
    
    if (attr(frun, 'status') == 1) {
        log$run_status[i] = 'error'
    } else {
        log$run_status[i] = 'ok'
        run_info = as.character(fromJSON(frun[[3]])$spec$`_tasks`[[5]]$script)
        log$counts_file[i] = str_extract(run_info, '(?<=\ )[A-Za-z0-9\\-\\_\\/\\.]+countsMerged[A-Za-z0-9\\-\\_\\/\\.]+(?=\ )')
    }
}

write.table(log,
            paste0(facets_dir, '/log/facets_runs_', str_replace_all(input_date, '\\-', '\\_'), '.txt'),
            quote = F, sep = '\t', col.names = T, row.names = F
            )
