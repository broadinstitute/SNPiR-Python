# SNPiR-Python
<p align="justify">
<strong> SNPiR-Python</strong> is a Python implementation of the original SNPiR code. SNPiR identifies Genomic Variants (SNPs) in RNA-seq data with high accuracy using several hard filtering methods. 
</p>

The rationale behind reimplementing SNPiR into python is due to speed and modernization. 

<p align="justify">
To help speed up the analysis, PBLAT is used in substitution of BLAT. PBLAT is aparalalized method of BLAT, allowing it to run on multiple cpu threads. 
</p>

<p align="justify">
The original SNPiR method utilizes GATK3’s UnifiedGenotyper to identify variants in RNA-seq data. UnifiedGenotyper is now deprecated in favor of the updated GATK4 HaplotypeCaller method for variant identification. Currently SNPiR-Python does not run HaplotypeCaller, but instead takes a VCF file as input, allowing the user freedom to run their own HaplotypeCaller with the GATK best practice or SNPiR recommended parameters. 
The original SNPiR was written in Perl while this reimplementation is written in Python. Perl is becoming less popular as each day passes, giving rise to the widely used Python language. It is subjective wether python is truly better than Perl, but it is more advantageous to utilize the more popular and widely used language Python. 
</p>

<p align="justify">
The following are precision-recall plots of SNPiR and SNPiR-Python variant identification on the gold standard NA12878 dataset from Genome in a Bottle consortium (Zook et al., 2014; Zook et al., 2019). These plots ensure that SNPiR-Python output does not differ from the original SNPiR. 
</p>

TL/DR: SNPiR-Python is a faster and more updated version of the original SNPiR code 

**Reference:**
<ul>Piskol R, Ramaswami G, Li JB. Reliable identification of genomic variants from RNA-seq data. Am J Hum Genet. 2013 Oct 3;93(4):641-51. doi: 10.1016/j.ajhg.2013.08.008. Epub 2013 Sep 26. PMID: 24075185; PMCID: PMC3791257.</ul>
