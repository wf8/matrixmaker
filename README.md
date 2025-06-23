
# Matrix Maker

Matrix Maker is a simple but handy utility to download DNA sequences from GenBank
and construct alignments when the taxa (with or without synonyms) and gene regions of interest are already known.

### Citation:

Freyman, W.A. and A.H. Thornhill. 2016. Matrix Maker [Computer software]. Retrieved from https://github.com/wf8/matrixmaker

### Other stuff:

Copyright 2016 Will Freyman - willfreyman@gmail.com
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

### Requirements:

    Python 3+
    Biopython
    MAFFT v6.9+

### Quick usage:

    python matrix_maker.py -e your_email@here.edu -g example/genes.csv -m 5000 -s example/synonyms.csv

### Detailed usage:

    usage: matrix_maker.py [-h] [--email EMAIL] [--genes GENES]
                           [--max_seq_length MAX_SEQ_LENGTH] [--species SPECIES]
                           [--taxids TAXIDS]

    optional arguments:
      -h, --help            show this help message and exit
      --email EMAIL, -e EMAIL
                            Email address for NCBI database searches.
      --genes GENES, -g GENES
                            Text file that defines the gene regions of interest
                            using both include and exclude terms.
      --max_seq_length MAX_SEQ_LENGTH, -m MAX_SEQ_LENGTH
                            Optional. Sets the maximum sequence length to include.
                            Use this to exclude genomes.
      --species SPECIES, -s SPECIES
                            Text file that contains a list of all species
                            binomials and their synonyms.
      --taxids TAXIDS, -t TAXIDS
                            Optional. Text file that contains a list of all
                            taxids. Use this to avoid repeating the NCBI taxid
                            lookups.

