## *HIFIBarcode*
---

### DESCRIPTION
HIFIBarcode is used to produce full-length COI barcodes from pooled PCR amplicons generated by individual specimens.


### Change logs
- HIFIBarcode v2.0.0, 201901. <b>HUGE CHANGE</b>: 
	- We change programming language form PERL to Python3. We will not update old version(perl) anymore. 
	- This python version is about 15 times fast than perl version, highly recommending to use this version!
	- It will not neet cmr to connect paired reads.
- HIFIBarcode v1.3.0, 201805. fix a bug on step5 mid.lis
- HIFIBarcode v1.0.0, 201707, fist release.

### INSTALLATION
```
git clone https://github.com/comery/HIFI-barcode-hiseq.git
python3 HIFIBarcode.py -h
```
	
```text
usage: HIFI-hiseq [-h] [-v]
                  {all,filter,assign,buildend,chain,gapfill,mkout,polish,bold_identification}
                  ...

Description

    An automatic pipeline for HIFI-hiseq project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences).

Version

    2.0.0 2018-12-13 The first python version.

Author

    yangchentao at genomics.cn, BGI.
    zhouchengran at genomics.cn, BGI.
    liushanlin at genomics.cn, BGI.

positional arguments:
  {all,filter,assign,buildend,chain,gapfill,mkout,polish,bold_identification}
    all                 run filter, assign, buildend, chain, and gapfill
    filter              filter raw reads by quality or expected_err
    assign              assign clean reads to samples
    buildend            buildends for each sample, output DNA fragment with tag
    chain               connect middle meta-paried reads to longer contigs
    gapfill             gap filling to generate raw contigs
    mkout               rename raw contigs to final COI barcodes
    polish              polish COI barcode assemblies, output confident barcodes
    bold_identification
                        do taxa identification on BOLD system

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

#### (3) Published softwares:
>VSEARCH v2.4.4. Rognes, Torbjørn, et al. "VSEARCH: a versatile open source tool for metagenomics." PeerJ 4 (2016): e258.   

>SOAPBarcode. Liu, Shanlin, et al. "SOAPBarcode: revealing arthropod biodiversity through assembly of Illumina shotgun sequences of PCR amplicons." Methods in Ecology and Evolution 4.12 (2013): 1142-1150.

### Pre-requisites
Operating system:  
HIFIBarcode is designed to run on most platforms, including UNIX, Linux and MacOS/X. Microsoft Windows. We have tested on Linux and on MacOS/X, because these are the machines we develop on. HIFIBarcode is written in python language, and a version 3.5 or higher is required.

#### Dependencies:
- biopython version 1.5 or higher (required). Please check https://biopython.org/ and https://pypi.org/project/biopython/#description for more details on installation of biopython.
- Another python package - bold_identification is also required for getting complete function of HIFI-SE. See https://pypi.org/project/bold-identification/
- HIFIBarcode supposed you have installed the VSEARCH on your device, and its path in your $PATH. See https://github.com/torognes/vsearch
- HIFIBarcode needs software - barcode to achieve gap-filling procee, you can find it in bin/

### EXAMPLES
- step 1, run all

```shell
python3 HIFIBarcode.py all -outpre hifi -index 5 -q1 test1.fq.gz -q2 test2.fq.gz -primer indexed_primer.txt
```
- step 2, run gap-filling  

```text
# at first, go to shell dir
cd hifi_gapfill/shell
# because this step will cost more memory, so I isolate this step.  
# then, you can run all *.sh file in this fold by way you like.  
# e.g. you can use qsub if you are a cluster user, or you can just  
# run on your local computer.  
ls *.sh|while read a; do qsub -cwd -l vf=1.5g $a;done
or:
sh *.sh
```

- step 3, generate output contigs  

```shell
python3 HIFIBarcode.py mkout -outpre hifi -d hifi_gapfill
```

- [optional] step 4, polish assemblies

```shell
python3 HIFIBarcode.py -outpre hifi -i hifi_barcodes.fa -index 5 
```

### PERL version
if you still use old versions, you can see [PERL VERSION](./Perl_version/PERL_Version.md) for help.

### CONTACT US

>yangchentao at genomics dot cn  
liushanlin at genomics dot cn  
zhouchengran at genomics dot cn  

### CITATION
Liu, Yang, Zhou, Zhou.(2017). Filling reference gaps via assembling DNA barcodes using high-throughput sequencing - moving toward to barcoding the world. Gigascience. 2017 Oct 25. doi: 10.1093/gigascience/gix104.


### Copyright
This package is released under version 3 of the GNU General Public License (GPLv3). Please refer to https://www.gnu.org/licenses/gpl-3.0.html.
