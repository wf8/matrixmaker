#! /usr/bin/python
"""
matrix_maker.py

A simple utility that finds the NCBI taxid for a list of taxa, and
then searches GenBank for a set of genes using the taxids. The script
then downloads the sequences and produces a MAFFT alignment.

Copyright 2015 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import argparse
import csv
import os
import sys
import time
import urllib
import re
from Bio import Entrez
from Bio import SeqIO


class Gene(object):
    """
    Class responsible for managing the search terms of genes we are searching for in NCBI.
    """


    def __init__(self, name=""):
        self.name = name         # name of gene region
        self.gene_names = []     # list of search terms for this gene
        self.exclusions = []     # list of search terms to exclude for this gene



class Taxon(object):
    """
    Class responsible for managing the data for each taxon.
    """


    def __init__(self, binomial, taxid=""):
        """
        Optionally accept the NCBI taxid.
        """
        self.binomial = binomial    # genus_species
        self.taxid = taxid          # NCBI taxid
        self.synonyms = []          # list of synonym
        self.sequences = {}         # a dictionary holding lists of Bio.SeqRecords, each key is a gene name


    def get_taxid(self, email):
        """
        Gets the NCBI taxid from entrez taxonomy.
        """
        toolname = "matrix_maker"
        params = {
            'db': 'taxonomy',
            'tool': toolname,
            'email': email,
            'term': self.binomial.replace(' ', '_'),
            'rettype': 'xml',
        }
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
        url = url + urllib.urlencode(params)
        data = urllib.urlopen(url).read()
        if re.search('<Id>(\S+)</Id>', data):
            self.taxid = re.search('<Id>(\S+)</Id>', data).group(1)
        else:
            # taxid was not found using the binomial,
            # so now check for synonyms...
            for synonym in self.synonyms:
                params = {
                    'db': 'taxonomy',
                    'tool': toolname,
                    'email': email,
                    'term': synonym,
                    'rettype': 'xml',
                }
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
                url = url + urllib.urlencode(params)
                data = urllib.urlopen(url).read()
                if re.search('<Id>(\S+)</Id>', data):
                    self.taxid = re.search('<Id>(\S+)</Id>', data).group(1)
                    break
            # if taxid is still not found
            if self.taxid == '':
                self.taxid='not found'
        return self.taxid


    def get_sequences(self, email, gene):
        """
        Searches Entrez Nucleotide database for taxid and and a list of gene names and
        downloads results. Appends the resulting list of Bio.SeqRecords to self.sequences.
        Exclude any search term exclusions found in the records description.
        """
        Entrez.email = email
        # (txid202994[Organism] AND (rbcL[All Fields] OR internal transcribed spacer[All Fields])
        term = "txid" + self.taxid + "[Organism] AND ("
        for i, name in enumerate(gene.gene_names):
            if i == 0:
                term = term + name + "[All Fields]"
            else:
                term = term + " OR " + name + "[All Fields]"
        term = term + ")"
        #print("Using search term: " + term)
        try:
            handle = Entrez.esearch(db="nuccore", term=term)
            records = Entrez.read(handle)
            gi_list = records["IdList"]
            gi_str = ",".join(gi_list)
            #print("Found GenBank GIs: " + gi_str)
            handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "gb")
            final_records = []
            for record in records:
                exclude = False
                # check to make sure exclusions are not in the description
                for exclusion in gene.exclusions:
                    if exclusion in record.description:
                        exclude = True
                        break
                # check to make sure the search terms are actually in the description 
                # (they might have been in other parts of the record)
                include = False
                for name in gene.gene_names:
                    if name in record.description:
                        include = True
                        break
                if not exclude and include:
                    final_records.append(record)
        except Exception as e:
            print("\nCaught error while downloading from NCBI: " + str(e))
            print("Trying again...")
            time.sleep(0.5)
            return self.get_sequences(email, gene)
        self.sequences[gene.name] = final_records
        return final_records


    def get_longest_seq(self, gene_name, max_seq_length = -1):
        """
        Returns the longest SeqRecord shorter than max_seq_length for a given gene.
        """
        try:
            records = self.sequences[gene_name]
        except:
            records = []
        longest_len = 0
        longest_seq = None
        for record in records:
            if len(record) < max_seq_length or max_seq_length == -1:
                if len(record) > longest_len:
                    longest_len = len(record)
                    longest_seq = record
        return longest_seq


def main():

    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", "-e", help="Email address for NCBI database searches.")
    parser.add_argument("--genes", "-g", help="Text file that contains a list of all gene names.")
    parser.add_argument("--max_seq_length", "-m", help="Optional. Sets the maximum sequence length to include. Use this to exclude genomes.")
    parser.add_argument("--species", "-s", help="Text file that contains a list of all species binomials and their synonyms.")
    parser.add_argument("--taxids", "-t", help="Optional. Text file that contains a list of all taxids. Use this to avoid repeating the NCBI taxid lookups.")
    args = parser.parse_args()

    print("\n\nmatrix_maker.py\n\n")

    if not args.email:
        print("NCBI requires an email address for database searches. Please use the --email flag to specify an email address.\n")
        sys.exit(0)
    else:
        email = args.email

    if not args.species or not os.path.isfile(args.species):
        print("Please specify a valid list of taxa to search for.\n")
        sys.exit(0)

    if args.max_seq_length:
        max_seq_length = int(args.max_seq_length)
    else:
        max_seq_length = -1

    genes = []
    if not args.genes or not os.path.isfile(args.genes):
        print("Please specify a valid list of genes to search for.\n")
        sys.exit(0)
    else:
        # read in gene names....
        # format of file:
        # gene_name,include,rbcL,RBCL
        # gene_name,exclude,RRRBCL
        with open(args.genes, 'rb') as csvfile:
            genereader = csv.reader(csvfile, delimiter=",")
            for row in genereader:
                if row[1] == "include":
                    gene = Gene(row[0])
                    for i in range(2, len(row)):
                        if row[i] != "":
                            gene.gene_names.append(row[i])
                    genes.append(gene)
                if row[1] == "exclude":
                    for gene in genes:
                        if gene.name == row[0]:
                            for i in range(2, len(row)):
                                if row[i] != "":
                                    gene.exclusions.append(row[i])
                    
    # list of all taxon objects
    taxa = []

    # check for taxid
    print("Checking for taxids csv file...")
    if args.taxids and os.path.isfile(args.taxids):
        with open(args.taxids, 'rb') as csvfile:
            print("Found taxids csv file, reading taxids...\n")
            taxidsreader = csv.reader(csvfile, delimiter=",")
            for row in taxidsreader:
                taxa.append(Taxon(row[0], row[1]))
    else:
        print("No taxids csv file found.\n")

    # open species list file, get synonyms and any missing taxids
    with open(args.species, 'rb') as csvfile:
        print("Checking list of species, getting missing taxids from NCBI...")
        taxids_file = open("taxids.csv", "w")
        namesreader = csv.reader(csvfile, delimiter=",")
        i = 1
        num_lines = sum(1 for line in open(args.species))
        for row in namesreader:
            # update status
            percent = str(round(100 * i/float(num_lines), 2))
            sys.stdout.write('\r' + 'Completed: ' + str(i) + '/' + str(num_lines) + ' (' + percent + '%)')
            sys.stdout.flush()
            i += 1
            # check to see if we already have a taxid for this species
            found = False
            for taxon in taxa:
                if taxon.binomial == row[0]:
                    found = True
                    taxids_file.write(taxon.binomial + "," + taxon.taxid + "\n")
                    # add synonyms
                    for j in range(1, len(row)):
                        taxon.synonyms.append(row[j])
                    break
            if not found:
                taxon = Taxon(row[0])
                # add synonyms
                for j in range(1, len(row)):
                    taxon.synonyms.append(row[j])
                taxa.append(taxon)
                # get the taxid from NCBI
                taxon.get_taxid(email)
                # dont overload genbank
                time.sleep(0.1)
                taxids_file.write(taxon.binomial + "," + taxon.taxid + "\n")
        taxids_file.close()
        print("\nWriting all taxids to file taxids.csv...")

    print("\nDownloading sequences from NCBI...") 
    for gene in genes:
        print("\nSearching for gene: " + gene.name)
        i = 1
        for taxon in taxa:
            # update status
            percent = str(round(100 * i/float(len(taxa)), 2))
            sys.stdout.write('\r' + 'Completed: ' + str(i) + '/' + str(num_lines) + ' (' + percent + '%)')
            sys.stdout.flush()
            i += 1
            if taxon.taxid != "not found":
                taxon.get_sequences(email, gene)
                # dont overload genbank
                time.sleep(0.2)

        print("\nGenerating unaligned FASTA file...")
        unaligned_file = open(gene.name + ".fasta", "w")
        for taxon in taxa:
            record = taxon.get_longest_seq(gene.name, max_seq_length)
            if record != None:
                # output format: >binomial_accession_description
                description = taxon.binomial + "_" + record.id + "_" + record.description
                description = description.replace(" ", "_")
                unaligned_file.write(">" + description + "\n")
                unaligned_file.write(str(record.seq) + "\n\n")
        unaligned_file.close()

        print("Making alignment with MAFFT...")
        try:
            from Bio.Align.Applications import MafftCommandline
            mafft_cline = MafftCommandline(input=gene.name + ".fasta")
            mafft_cline.set_parameter("--auto", True)
            mafft_cline.set_parameter("--adjustdirection", True)
            print(str(mafft_cline))
            stdout, stderr = mafft_cline()
            print("Writing alignment to FASTA file...")
            with open("aligned_" + gene.name + ".fasta", "w") as handle:
                handle.write(stdout)
        except:
            print("Problem finding MAFFT, alignment skipped.")

    print("\nGenerating summary results spreadsheet...\n")
    summary = open("result.csv", "w")
    header = "taxon,"
    for gene in genes:
        header += gene.name + ","
    summary.write(header + "\n")
    for taxon in taxa:
        accessions = taxon.binomial + ","
        for gene in genes:
            # each column will be the longest sequences accession
            record = taxon.get_longest_seq(gene.name, max_seq_length)
            if record != None:    
                accessions += record.id + ","
            else:
                accessions += ","
        summary.write(accessions + "\n")
    summary.close()
    print("Done!\n")


        


if __name__ == "__main__":
    main()
