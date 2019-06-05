# frequency-scanner
Reports frequency of all single nucleotide changes for a group of fastq files 

# Dependencies
This program requires that you have the following programs installed and working on your path

1. Python
2. fastq_quality_filter[http://hannonlab.cshl.edu/fastx_toolkit/] 
3. bwa-mem[http://bio-bwa.sourceforge.net/] (reccomend installing through a package manager)
4. samtools[http://samtools.sourceforge.net/] (reccomend installing through a package manager)
5. bam-readcount[https://github.com/genome/bam-readcount] You have to build this one yourself. I reccomend downloading directly into the folder you're installing this program into. Then make and build right in this folder. 


# Running
Once you've got the above dependencies working and on your path just run frequency_scanner.py with the only argument being your reference fasta file. It'll run on all the fastq files in your folder. You can get a list of configurable options by typing `python frequency_scanenr.py -h`

The output will be a big csv with all the frequencies in it. This can be easily imported into R and you can make pretty pictures. 
