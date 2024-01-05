#!/usr/bin/env python
"""Querying isolate MLST types from Pubmlst.org and sending results to REDcap"""

import argparse
import requests
from pyfaidx import Fasta
import json
import sys

#------------------------------------------------
def get_args():
    """get command line args"""

    parser = argparse.ArgumentParser(
            description='Three positional arguments',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
            metavar='file',
            type=str,
            help="Isolate fasta file")

    parser.add_argument('species',
            metavar='species',
            type=str,
            help="species, ex: calbicans")

    parser.add_argument('ref_genome',
            metavar='ref_genome',
            type=str,
            help="Abbreviation of ref genome, ex: SC5314 A21")

    return parser.parse_args()


#------------------------------------------------
def main():
    """main"""

    args = get_args()
    with open(args.file) as fasta:
        seqs = fasta.read()

    # mlst loci variables
    if args.species == "calbicans":
        species_loci=["AAT1a", "ACC1", "ADP1", "MPIb", "SYA1", "VPS13", "ZWF1b"]
    elif args.species == "cglabrata":
        species_loci = ["FKS", "LEU2", "NMT1", "TRP1", "UGP1", "URA3"]
    elif args.species == "ckrusei":
        species_loci=["HIS3", "LEU2", "NMT1", "TRP1", "ADE2", "LYS2D" ]
    elif args.species == "ctropicalis":
        species_loci=["ICL1", "MDR1", "SAPT2", "SAPT4", "XYR1", "ZWF1a"]
    else:
        species_loci = []

    # url variables
    base_url = "https://rest.pubmlst.org/db/"
    seq_db_url = f"pubmlst_{args.species}_seqdef/"
    scheme_url = "schemes/1/sequence"

    # params
    isolate_st_params = {"base64" : False, "sequence" : seqs, "details": True}
    isolate_locus_params = {"base64" : False, "sequence" : seqs, "details" : True}

    # query functions for MLST sequence type or loci
    def isolate_st_query(params):
        return (requests.post((base_url + seq_db_url + scheme_url),
            json = isolate_st_params)).json()

    def isolate_locus_query(params):
        return (requests.post((base_url + seq_db_url + locus_url),
            json = isolate_locus_params)).json()


    # Perform initial query for exact matches
    isolate_seqs = isolate_st_query(isolate_st_params)

    # start dict
    isolate_data = {'primary_id' : args.file.removesuffix("_mlst.fa")}
    isolate_data['typing_scheme'] = '0'
    isolate_data['st_ref_genome'] = args.ref_genome
    if 'fields' in isolate_seqs:
        isolate_data['ST'.lower()] = isolate_seqs.get('fields').get('ST')

    # add locus exact matches
    for locus in isolate_seqs.get('exact_matches'):
        isolate_data[locus.lower() + '_exact_match'] = isolate_seqs.get('exact_matches', {}).get(locus)[0].get('allele_id')
        locus_contig = isolate_seqs.get('exact_matches').get(locus)[0].get('contig')
        isolate_data[locus.lower() + '_sequence'] = isolate_seqs.get('exact_matches').get(locus)[0].get('contig') + str(Fasta(args.file)[locus_contig][:])

    # if ST is present, then all loci are too; send dict as json to redcap api
    # otherwise search for best matches to finish dict and then send to redcap

    if 'st' in isolate_data:
        isolate_data['sequence_typing_data_complete'] = '2'
    else:
        for locus in species_loci:
            if locus in isolate_seqs:
                pass
            else:
                locus_url = f"loci/{locus}/sequence"
                locus_seqs = isolate_locus_query(isolate_locus_params)
                if 'best_match' in locus_seqs:
                    isolate_data[locus.lower() + '_best_match'] =str(["allele ID",
                            locus_seqs.get('best_match').get('allele_id'),
                            "identity",
                            locus_seqs.get('best_match').get('identity'),
                            "mismatches",
                            locus_seqs.get('best_match').get('mismatches')])
                    locus_contig =locus_seqs.get('best_match').get('contig')
                    isolate_data[locus.lower() + '_sequence'] = locus_seqs.get('best_match').get('contig') + str(Fasta(args.file)[locus_contig][:])
        isolate_data['sequence_typing_data_complete'] = '1'

    data = json.dumps([isolate_data])
    fields = {
        'token': '',
        'content': 'record',
        'action': 'import',
        'format': 'json',
        'type': 'flat',
        'data': data,
        'overwriteBehavior': 'normal',
        'forceAutoNumber': 'false',
        'returnContent': 'count',
        'returnFormat': 'json'
    }
    r = requests.post('https://redcap.ahc.umn.edu/api/',data=fields)
    print(isolate_data['primary_id'] + ' HTTP Status: ' + " " + str(r.status_code), file=sys.stderr)
    print(r.json(), file=sys.stderr)


#------------------------------------------------
if __name__ == '__main__':
    main()
