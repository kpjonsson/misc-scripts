#!/usr/bin/env python
################################################################################
descr = 'Wrapper for the cmo_facets doFacets wrapper'
################################################################################

import argparse
import subprocess
import os
import sys
import re

bam_dir = '/ifs/dmpshare/share/irb12_245/'
key = '/home/jonssonp/keys/key.txt'

def run_impact_facets():
    parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-l', '--lib-version', help = 'Version of FACETS R pacakge', required = False, default = '0.5.6')
    parser.add_argument('-t', '--tumor-sample', help = 'MSK-IMPACT tumor sample name', required = True)
    parser.add_argument('-n', '--normal-sample', help = 'MSK-IMPACT normal sample name', required = False)
    parser.add_argument('-c', '--cval', help = 'cval parameter', required = False, default = 50)
    parser.add_argument('-pc', '--purity-cval', help = 'Purity cval, to run 2-pass mode', required = False, default = 150)
    parser.add_argument('-m', '--min-nhet', help = 'Heterozygous SNPs required for segmentation',
        required = False, default = 15)
    parser.add_argument('-pm', '--purity-min-nhet', help = 'Purity min. nhet, to run 2-pass mode',
        required = False, default = 15)
    parser.add_argument('-s', '--seed', help = 'Set seed parameter', required = False, default = 100)

    args = parser.parse_args()
    tumor_sample = args.tumor_sample
    normal_sample = args.normal_sample
    patient = tumor_sample[:9]
    platform = tumor_sample[16:]

    facets_args = {
        'v': args.lib_version,
        'c': str(args.cval),
        'pc': str(args.purity_cval),
        'm': str(args.min_nhet),
        'pm': str(args.purity_min_nhet),
        's': str(args.seed)
    }

    # If full sample name is provided
    if bool(re.match(r'P-[0-9]{7}-T[0-9]{2}-IM[0-9]{1}', tumor_sample)):
        query = query_key(patient)
        normal_sample, normal_bam, tumor_bam = pair_tumor_sample(tumor_sample, normal_sample, query)
        run_facets(normal_sample, normal_bam, tumor_sample, tumor_bam, facets_args)

    # If only partial sample name is provided, look for match
    elif bool(re.match(r'P-[0-9]{7}-T[0-9]{2}', tumor_sample)):
        query = query_key(patient)
        name_match = re.match(tumor_sample+'-IM[0-9]{1}', query)
        
        if name_match is not None: 
            print('Partial sample ID provided, found ' + name_match)
            tumor_sample = name_match.group()
            normal_sample, normal_bam, tumor_bam = pair_tumor_sample(tumor_sample, normal_sample, query)
            run_facets(normal_sample, normal_bam, tumor_sample, tumor_bam, facets_args)

        else:
            print('Partial sample ID provided, found no match')
        
    # If only patient ID is input
    elif bool(re.match(r'P-[0-9]{7}', tumor_sample)): 
        print('Patient ID provided, will look for all tumor samples from patient, ok [y/n]?')
        user_in = raw_input()
        if bool(re.match(r'(y|Y|yes)', user_in)):
            query = query_key(patient)
            all_tumors = re.findall(r'P-[0-9]{7}-T[0-9]{2}-IM[0-9]{1}', query)

            if len(all_tumors) > 0:
                print('Found: ' + ', '.join(all_tumors))
                for sample in all_tumors:
                    normal_sample = None
                    normal_sample, normal_bam, tumor_bam = pair_tumor_sample(sample, normal_sample, query)
                    run_facets(normal_sample, normal_bam, sample, tumor_bam, facets_args)

            else:
                sys.exit('Found no tumor samples from patient ' + patient)

        else:
            sys.exit()    
    else:
        sys.exit('Invalid tumor sample name provided, try again.')
 
def query_key(patient):
    """Grep bam file key for patient ID"""

    query = subprocess.Popen(['grep', patient, key], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    output, err = query.communicate()

    if err != '':
        sys.exit(err)
    elif output is '' or bool(re.match(patient, output)) is not True:
        sys.exit(patient + ' not found in BAM file key')
    else:
        return output

def pair_tumor_sample(tumor_sample, normal_sample, query_output):
    """Match tumor sample to normal based on platform, pick highest numbered sample"""
    
    sample_list = {}
    for k in [s for s in query_output.split('\n') if s != '']:
        v = k.split(',')
        if bool(re.match(r'P-[0-9]{7}-N[0-9]{2}-IM[0-9]{1}', v[0])):
            type = 'normal'
        else:
            type = 'tumor'
        platform = re.search(r'(?<=IM)[0-9]{1}', v[0]).group()
        number = re.search(r'(?<=[A-Z]{1}0)[0-9]{1}', v[0]).group()
        sample_list[v[0]] = {'name': v[1], 'type': type, 'platform': platform, 'number': number}

    tumor_bam = ''.join([bam_dir, sample_list[tumor_sample].get('name'), '.bam'])

    if normal_sample is not None:
        if normal_sample in sample_list:
            normal_bam = sample_list[normal_sample].get('name')
        else:
            sys.exit(normal_sample + 'not found in BAM file key')
    else:
        tumor_platform = sample_list[tumor_sample].get('platform')
        normals = [x for x in sample_list.keys() if sample_list[x].get('type') == 'normal' and sample_list[x].get('platform') == tumor_platform]
        if len(normals) == 0:
            sys.exit('Could not find an appropriate normal in BAM file key')
        normal_numbers = map(int, [(v.get('number')) for (k, v) in sample_list.items() if k in normals])
        best_normal = normals[normal_numbers.index(max(normal_numbers))]
        normal_bam = ''.join([bam_dir, sample_list[best_normal].get('name'), '.bam'])

    return best_normal, normal_bam, tumor_bam

def run_facets(normal_sample, normal_bam, tumor_sample, tumor_bam, facets_args):
    """Construct and run cmoflow_facets command"""

    facets_cmd = ' '.join(
        ['cmoflow_facets',
        '--R_lib', facets_args['v'],
        '--normal-bam', normal_bam,
        '--tumor-bam', tumor_bam,
        '--normal-name', normal_sample,
        '--tumor-name', tumor_sample,
        '--cval', facets_args['c'],
        '--purity_cval', facets_args['pc'],
        '--min_nhet', facets_args['m'],
        '--purity_min_nhet', facets_args['pm']])

    print 'Running Facets:\nTumor: ' + tumor_sample + '\nNormal: ', normal_sample
    subprocess.call(facets_cmd, shell = True)

if __name__ == '__main__':
    run_impact_facets()
