# cLCP-mACS
Efficient lightweight strategy to solve the multi-string Average Common Substring (ACS) problem, that consists in the pairwise comparison of a single string against a collection of m strings simultaneously, in order to obtain m ACS induced distances.

This software is an implementation of the algorithm described in:

>F. Garofalo, G. Rosone, M. Sciortino and D. Verzotto
*The colored longest common prefix array computed via sequential scans.*
Proceedings of String Processing and Information Retrieval - 25th International Symposium, SPIRE 2018, Lecture Notes in Computer Science, Springer <sup id="a1">[1](#f1)</sup>


### Install

```sh
git clone https://github.com/giovannarosone/cLCP-mACS.git
cd cLCP-mACS
make all
```

### Run

```sh
./cLCP-mACS [-h] [-v] [-p] [-l] [-f input_format] [-Q amount] ref_seq target_seqs ref_color output
```

##### Input
<!-- cLCP-mACS needs LCP and ID array (array of colors) of a collection of sequences (_target collection_) and the LCP and ID (color) of one of the sequence in the collection (_reference sequence_). In current impementation these are provided by the output of eGSA tool. -->

cLCP-mACS takes as input the GESA (_Generalized Enhanced Suffix Arrays_) of a collection of sequences (_target collection_) and the GESA of one of the sequence in the collection (_reference sequence_). At the current state, these files are assumed as resulting by [eGSA tool][240fb5f5] computation:

  [240fb5f5]: https://github.com/felipelouza/egsa "eGSA: Generalized enhanced suffix array construction in external memory [CPM'13, ALMOB 2017]"

```sh
ref_seq       GESA file name of reference sequence (without .gesa extension)
target_seqs   GESA file name of target collection (without .gesa extension)
ref_color     ID of reference sequence in the target collection
```
##### Output
cLCP-mACS computes the multi-ACS measure between the reference sequence and the remaining sequences of the target collection. The computation collaterally produces in addition to the the file with the distance values (`.acs`) two working files containing the array D (`.d`) and a partial cLCP (`.xclcp`).
```sh
output       Output/Working files name
```
##### Options:

```sh
-h    help message
-v    verbose textual output
```
The option `-f    input_format` specifies the format of the input to cLCP-mACS. For now, only `-f    1` option is admitted (by default) corresponding to the GESA computed by eGSA.

The option `-p` tells if the preprocessing step can be skipped. You should use this option only if the target collection has been already processed in a previous execution (this is indicated by the presence of a `.info` file, and files `.bwt`,`.lcp`,`.id`).

Option `-l` also skips preprocessing step provided that are available a ` .lenSeqs.aux` file and `.bwt`,`.lcp`,`.id` files (as produced by [BCR tool](https://github.com/giovannarosone/BCR_LCP_GSA)).

The option `-Q amount` dictates the amount of RAM (in Bytes) available to accomodate partial cLCP page. `amount` has to be at least 2 x _m_, where _m_ is the number of sequences in target collection.

#### Contributors

Fabio Garofalo,  University of Palermo

Marinella Sciortino,  University of Palermo

Giovanna Rosone, University of Pisa (project lead)

Davide Verzotto, University of Pisa

#### Citation
If you use cLCP-mACS in an academic setting you could cite the following paper:

    @inproceedings{GRSV_Spire18,
    author    = {Garofalo, Fabio and Rosone, Giovanna and Sciortino, Marinella and Verzotto, Davide},
    title     = {The colored longest common prefix array computed via sequential scans},
    booktitle = {String Processing and Information Retrieval - 25th International Symposium,
               {SPIRE} 2018, Proceedings},
    pages     = {},
    year      = {2018},
    series    = {Lecture Notes in Computer Science},
    volume    = {},
    publisher = {Springer}
    }
---
1. <small id="f1"> Supported by the project Italian MIUR-SIR CMACBioSeq ("_Combinatorial methods for analysis and compression of biological sequences_")
grant n.~RBSI146R5L.</small> [↩](#a1)
