#!/usr/bin/env python
################################################################################
descr = 'Wrapper for the cmo_facets doFacets wrapper'
################################################################################

import argparse
import subprocess
import os
import sys
import re

def rerun_facets():
    parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-l', '--lib-version', help = 'Version of FACETS R pacakge', required = False, default = '0.5.14')
    parser.add_argument('-f', '--counts-file', help = 'Merged SNP counts file (*.dat.gz)', required = True)
    parser.add_argument('-c', '--cval', help = 'cval parameter', required = False)
    parser.add_argument('-pc', '--purity-cval', help = 'Purity cval, to run 2-pass mode', required = False)
    parser.add_argument('-m', '--min-nhet', help = 'Heterozygous SNPs required for segmentation',
        required = False, default = 15)
    parser.add_argument('-pm', '--purity-min-nhet', help = 'Purity min. nhet, to run 2-pass mode',
        required = False, default = 15)
    parser.add_argument('-d', '--diplogr', help = 'Diploid log-ratio', required = False)
    parser.add_argument('-D', '--directory', help = 'Output directory', required = False)
    parser.add_argument('-t', '--tag', help = 'Output file tag', required = False)
    parser.add_argument('-s', '--seed', help = 'Set seed parameter', required = False, default = 100)

    args = parser.parse_args()
    facets_args = {
        'lib-version': args.lib_version,
        'counts_file': args.counts_file,
        'min_nhet': args.min_nhet,
        'seed': args.seed
    }

    if args.diplogr is not None:
        facets_args['dipLogR'] = args.diplogr

    if args.cval is not None:
        facets_args['cval'] = args.purity_cval

    if args.purity_cval is not None:
        facets_args['purity_cval'] = args.purity_cval

    if args.purity_min_nhet is not None:
        facets_args['purity_min_nhet'] = args.purity_min_nhet 

    if args.directory is None:
        base_dir = os.path.dirname(args.counts_file)
        directory = base_dir + '/facets_' + args.lib_version + 's' + str(args.seed) + 'm' + str(args.min_nhet)
        if 'purity_min_nhet' in facets_args: 
            directory += 'pm' + str(facets_args['purity_min_nhet'])
        if 'cval' in facets_args: 
            directory += 'c' + str(facets_args['cval'])
        if 'purity_cval' in facets_args: 
            directory += 'pc' + str(facets_args['purity_cval'])
        if 'dipLogR' in facets_args: 
            directory += '_diplogr' + str(facets_args['dipLogR'])
        facets_args['directory'] = directory
    else:
        facets_args['directory'] = args.directory

    if os.path.exists(facets_args['directory']):
        sys.exit('Output directory exists, speciffy a diffferent one.')
    else:
        print(facets_args['directory'])
        subprocess.call('mkdir -p ' + facets_args['directory'], shell = True)

    if args.tag is None:
        tumor = re.search(r'P-[0-9]{7}-T0[0-9]{1}-IM[0-9]{1}', args.counts_file)
        normal = re.search(r'P-[0-9]{7}-N0[0-9]{1}-IM[0-9]{1}', args.counts_file)
        if bool(tumor):
            d = tumor.group()
            if bool(normal):
                d += '_' + normal.group()
            facets_args['TAG'] = d
        else:
            facets_args['TAG'] = 'facets_rerun'
    else:
        facets_args['TAG'] = args.tag

    call = 'cmo_facets --lib-version ' + facets_args['lib-version'] + ' doFacets'
    for key, value in facets_args.items():
        if key != 'lib-version':
            call += ' --' + key + ' ' + str(value) 
    print('Running...' + '\n' + call)

    subprocess.call(call, shell = True)

if __name__ == '__main__':
    rerun_facets()
