HIFIBarcode V1.0 (201707)

DESCRIPTION
	HIFIBarcode is used to produce full-length COI barcodes from pooled PCR
amplicons generated by individual specimens.

INSTALLATION
	tar -zxvf HIFIBarcode.v1.0.tar.gz

Source code
(1) Wrappper script:
	HIFIBarcode.v1.0.pl: the wrapper script to run other Perl scripts to do the work you choose.
	
(2) Main scripts:
	Seven Perl scripts in folder HIDIBarcode.V1.0/bin/, including:
	1_split_extract.pl
	2_uniqu_sort_cluster.Pro.pl
	3_sep_extract_overlap.pl
	4_cluster_fromend.pl
	5_forgap_filling.pl
	6_rename_kmer.pl
	7_final.pl
(3) Published softwares:
	VSEARCH v2.4.4. Rognes, Torbjørn, et al. "VSEARCH: a versatile open source tool for metagenomics." PeerJ 4 (2016): e258.
	COPE CMR v1.0.3. Liu, Binghang, et al. "COPE: an accurate k-mer-based pair-end reads connection tool to facilitate genome assembly." Bioinformatics 28.22 (2012): 2870-2874.
	SOAPBarcode. Liu, Shanlin, et al. "SOAPBarcode: revealing arthropod biodiversity through assembly of Illumina shotgun sequences of PCR amplicons." Methods in Ecology and Evolution 4.12 (2013): 1142-1150.

Pre-requisites
	PERL v5

EXAMPLES
1: run wrapper script to get shell text and then sh runHIFIBarcode.sh to run HIFIBarcode
	perl HIFIBarcode.v1.0.pl  --fq1 test_1.fq --fq2 test_2.fq --index index_primer.txt  --length 5 --cpunum 10 --outdir test  --outpre testout
	sh test/runHIFIBarcode.sh

NOTE
The proposed method assemble short-read Illumina sequences based on k-mer sequence matches and such misassembly was not observed in their real data, but it’s still possible. Thus, we provide a function to hightlight NON-HIGHEST-ABUNDANCE sequences within file which has a suffix of "note.txt".

CONTACT US

Email:
liushanlin at genomics dot cn
zhouchengran at genomics dot cn

CITATION
Liu, Yang, Zhou, Zhou.(2017). Filling reference gaps via assembling DNA barcodes using high-throughput sequencing - moving toward to barcoding the world. Unpublished.

LATEST RELEASE
Version 1.1 201709

Copyright
This package is released under version 3 of the GNU General Public License (GPLv3). Please refer to https://www.gnu.org/licenses/gpl-3.0.html.
