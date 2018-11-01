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
    facets_args = {
        'lib-version': args.lib_version,
        'tumor_sample': args.tumor_sample,
        'min_nhet': args.min_nhet,
        'seed': args.seed
    }

    patient = args.tumor_sample[:9]
    platform = args.tumor_sample[16:]

    # Look for all patient BAM files
    query = subprocess.Popen(['grep', patient, key], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    output, err = query.communicate()

    # Check, unpack query
    if err != '':
        sys.exit(err)

    if output == '' or args.tumor_sample not in output:
        sys.exit(args.tumor_sample + 'not found in BAM file key.')

    # output = [s for s in output.split('\n') if s != '']
    # output = {sample:name for sample, name in (x.split(',')[:2] for x in output)}
    sample_list = {}
    for k in [s for s in output.split('\n') if s != '']:
        v = k.split(',')
        if bool(re.match(r'P-[0-9]{7}-N[0-9]{2}-IM[0-9]{1}', v[0])):
            type = 'normal'
        else:
            type = 'tumor'
        platform = re.search(r'(?<=IM)[0-9]{1}', v[0]).group()
        number = re.search(r'(?<=[A-Z]{1}0)[0-9]{1}', v[0]).group()
        sample_list[v[0]] = {'name': v[1], 'type': type, 'platform': platform, 'number': number}

    tumor_bam = ''.join([bam_dir, sample_list[args.tumor_sample].get('name'), '.bam'])

    if args.normal_sample is not None:
        if args.normal_sample in sample_list:
            normal_bam = sample_list['P-0029357-N01-IM6'].get('name')
        else:
            sys.exit(args.normal_sample + 'not found in BAM file key.')
    else:
        tumor_platform = sample_list[args.tumor_sample].get('platform')
        normals = [x for x in sample_list.keys() if sample_list[x].get('type') == 'normal' and sample_list[x].get('platform') == tumor_platform]
        if len(normals) == 0:
            sys.exit('Could not find an appropriate normal in BAM file key.')
        normal_numbers = map(int, [(v.get('number')) for (k, v) in sample_list.items() if k in normals])
        best_normal = normals[normal_numbers.index(max(normal_numbers))]
        normal_bam = ''.join([bam_dir, sample_list[best_normal].get('name'), '.bam'])

    facets_cmd = ' '.join(['cmoflow_facets', '--lib-version', args.lib_version, '--normal-bam', normal_bam, '--tumor-bam', tumor_bam, '--normal-name', best_normal, '--tumor-sample', args.tumor_sample, '--cval', str(args.cval), '--purity_cval', str(args.purity_cval), '--min_nhet', str(args.min_nhet), '--purity_min_nhet', str(args.purity_min_nhet)])

    print 'Running Facets:\nTumor: ' + args.tumor_sample + '\nNormal: ', best_normal

    subprocess.call(facets_cmd, shell = True)

if __name__ == '__main__':
    run_impact_facets()
