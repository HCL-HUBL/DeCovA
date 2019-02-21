### DeCovA_1.6.0



#### Requirements

DeCovA requires at least R and bedtools/GATK softwares to be installed;
additionnally, it can use picard-tools (for deduplication), samtools (for mapq filter), and GATK (alternatively to bedtools; required if a base-q filter is needed; GATK is also aware of pair reads overlap).
The script will first attempt to run programs installed as root with the following names: samtools, bedtools, picard-tools, GenomeAnalysisTK;
if not found, it will try to find them according to the paths provided in the command-line.

DeCovA also requires perl modules: IO::Compress::Gzip.

An annotation file needs to be provided (-r option), for all the options that use gene coordinates:
UCSC refgene.txt  or Ensembl .gtf/gff files are OK.

ex:

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
```

```bash
wget ftp://ftp.ensembl.org/pub/grch37/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
```

```bash
wget ftp://ftp.ensembl.org/pub/grch37/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
```



#### Running

DeCovA can be executed via command-line execution of the main perl script:

```bash
perl path/to/DeCovA/bin/DeCovA [options]
```

the script can be changed to executable:

```bash
chmod 755 path/to/DeCovA/bin/DeCovA
```

then the command is :

```bash
./path/to/DeCovA/bin/DeCovA [options]
```

The DeCovA script directory can be added to the $PATH :

```bash
echo 'export PATH=$PATH:/home/me/path/to/DeCovA/bin/DeCovA' >> /home/me/.bashrc
```

then the command is only :

```bash
DeCovA [options]
```



DeCovA can also be installed:

```bash
cd path/to/DeCovA/
perl Makefile.PL
make
```

then, as root:

```bash
sudo make install
```


then just enter 

```bash
DeCovA [options]
```



#### List of parameters

inputs:
​	-f / --file [file]: list of bam files (comma separated, or set several times)
​	-F / --fList [file]: file with such a list of bam files (one bam per line)
​	-d / --dir [dir]: directory(ies) where to find bam files (comma separated, or set several times)
​	-s / --suffix [str]: suffix to add before opening bam files
​	-r / --ref [file]: gene annotation file (can be .gz)
​	--fmt [gtf/gff3/ucsc] : gene annotation file format (ucsc <=> UCSC refGene) ; if not provided, determined from extension (txt => UCSC refGene)
​	-b / --bed [file]: bed file, used to analyse depth coverage
​	-m / --mut [file]: mut file, used to plot known mutations ; format: "chr<tab>pos(1-based)<tab>info" (vcf files are ok ; can be .gz)
​	-i / --id [str]: list of of genes/transcripts ids (comma separated, or set several times)
​	-I / --idList [file]: file with a list of of genes/transcripts ids (one id per line)
​	-g / --genome [file]: path to genome.fa file, if available (required if using GATK)
​	--sex_file [file]: format: patient<tab>sex
​	--raw_cov [file]: use this coverage tool output .cov file (to skip bam analysis)
​	--bed_cov [file]: use this DeCovA's output .cov.txt file (to skip cov bed analysis in CNV detect)

outputs:
​	-O / --outdir [dir]: out directory (default: folder named with date)
​	-S / --graphSum : will perform graphSums (sum of covered samples by position)
​	-A / --allSample : will perform graphAllSample (depthline by gene and by sample, all samples graph on same .png file)
​	-X / --bySample : will perform graphBySample (depthline by gene and by sample, one sample by .png file)
​	-M / --noDepthMut : does not print, foreach file, depth at known mutations provided by opt -m (default: yes if -opt m)
​	-P / --covPlot : will perform covPlots
​	-B / --covBed : will output cov of bed intervals
​	-C / --CNV : will output CNV foreach bed intervals
​	--Reseq : [float,0-1] : print bed interval if cov < value (def: do not print)
​	--geneReport : will print all uncovered genomic intervals (within gene region) in 1 txt file per sample (default: no)
​	--bedReport : will print all uncovered intervals (within bed intervals) in 1 txt file per sample (default: no)
​	--summary [Y/N] : to print summary txt file (default: yes if -S -A -X)
​	-k / --keepCov : do not erase coverage file at the end of the process
​	-K / --keepBed : do not erase bed file inferred from gene list, at the end of the process (and eventually rename)

parameters:
  *gene/transcript regions analysis param.:
​	-N / --nonCoding : analyse also Non coding transcripts (default: no)
​	-U / --noUTR : does not take into account UTR regions, for graphs (default: yes)
​	-u / --noUTRinTxt : does not take into account UTR regions, for summary txt file and plots (default: yes)
​	-t / --depthThreshold [int]: depth thresholds (comma separated, or set several times)
​	-T / --printThreshold [int]: depth threshold used for txt outputs (must be one of those in opt -t; default : the smallest one)
​	--noGraphThreshold : all graphs will be printed, whatever the coverage 
​		(default: only the genes not fully covered at threshold in -opt -T will be drawn)
​	--noAllTranscripts : does not print All transcripts on same file, in graphBySample (default: yes)
​	--maxDepth [int]: max depth value when printing graph (optional)
​	-l / --expand2val [int]: length to add at each ends of exons, on graphs (default: 0) ; or [int1,int2] : lengths to add in 5' and 3'
​	--UDstream [int]: length to add at each ends of genes, on graphs ; or [int1,int2] : lengths to add upstram and downstream
​	--splitBedFromId : if padding creates overlapping exons, take the mid between them (for report)
​	--mergeBedFromId : merge overlapping exons
​	-L / --expand2bed : expand length of gene analysed regions to bed coord, if -l < bed , on graphs (default: no)
​	--Ltxt [+/-int]: does take into account expanded length (from -l and -L) for txt outputs (default: no), or add a different length
​	--UDtxt [+/-int]: does take into account up/downStream length for txt outputs (default: no), or add a different length
​	-R / --noReverse : does not reverse regions if sens of transcript = (-) (default: yes)
​	--nGraph : max nber of graphs per sheet (default : all samples or all transcripts)

  *plot param:
​	--binPlot [int]: bin width for covPlot (default=10)
​	--maxPlot [int]: max depth for covPlot (default=100)
​	--genePlot : will perform plots for regions extracted from genes coord, not only for bed intervals (default: no)
​	--interPlot : will produce intersection covPlot (default: no)

  *bam filters
 	--dedup : do not take in account dup reads (default keep all reads; enter "do" to perform Picard deduplication) 
​	--mbq : minimum base quality (default 0; requires gatk)
​	--mmq : minimum mapping quality (default 0)

  *cov_bed param:
​	--cov_fields [min/max/tot/mean/median/cov]: fields foreach intervals in covBed (comma separated) (default: min,mean,cov)
​	--Lbed [int]: length added out of bed interval ends (default: 0)
​	--split_bed : splits overlapping bed intervals for Cov and CNV analyses
​	--no_overlap_bed : removes overlapping bed intervals for Cov and CNV analyses
​	--cut_bed [+/-cutL:x,minL:y,maxL:z,keepLast:s]: cut bed intervals in shorter fragments:
​				cutL : length of segmentation (def: 150)
​				minL : min length required to keep the last interval, after segmentation (def: --cutL/2)
​				maxL : length above which bed intervals will be segmented, in N segments of "cutL" length (def: as --cutL)
​				keepLast : if last interval shorter than minL :
​						enter m (merge) if want that last two ones are simply merged
​						enter h (half) if want that last two ones are output with length = half of their sum
​						enter n if want to through it out
​	--reAnnot_bed : removes and replaces 4th column of bed file with gene info
​		(optional args: g,t,e,i,o : indicates to annotate with gene/transcript/exon/intron/intergenic infos; default: all)

  *CNV_detect param:
​	--level2 : "avg"/"med" : use average/median as center of depths of a region (def: med)(if spread2 is set, level2 is unset, unless explicitedly)
​	--spread2 : "std"/"qtile" : use standard deviation/deviation from quartile as dispersion of depths of a region (def: none)(std forces avg, qtile forces med)
​	--level_del [float [0-1]] (def: 0.8)
​	--level_dup [float >1] (def: 1.2)
​	--spread_del [float <0] (def: none)
​	--spread_dup [float >0]  (def: none)
​	--range [float]: samples kept for avg-std calculation if within mediane+/-range*quartile (def: none, ie all samples used)
​	--highQual [li:float/ls:float/si:float/ss:float/c:int]: flag as high qual if one of following criteria, comma separated :
​		li=level inf, ls=level sup, si=spread inf, ss=spread sup, c=consecutive :
​		ex : li:0.25,ls:1.75,si:-5,ss:5,c:2
​	--ex_region [float [0-1]] : region excluded from analysis if CNVs/N_samples >value (def: 1)
​	--ex_sample [float [0-1]] : sample excluded from analysis if CNVs/N_regions >value (def: 1)
​	--ex_cov [float [0-1]] : region excluded from analysis if none of the samples have cov >=value (def: 0)
​	--ex_DP [int] : region excluded from analysis if avg depth <=value (def: 0)
​	--max_nonCNVcons [int]: max nber of non-CNV consecutive intervals tolerated within a CNV (def: 0)
​	--max_nonCNVrate [int]: max rate of non-CNV intervals tolerated within a CNV (def: 0)
​	--ratioByGender [a/g/no]:	enter "a" : foreach region from all chrom, depth ratio computed separately for F and M ; 
​				enter "g" : foreach region from gonosomes only, depth ratio computed separately for F and M
​				def: no (depth ratio for F and M together)
​	--normAllChr : total depth used to norm sample depths = sum on all chr, whatever the sex (def: double the depth for chrX if male, and skip chrY in the sum)
​	--normDepth [mean/tot] : total depth used to norm sample depths = sum of total depths of each region or sum of mean depths of each region (def)
​	--graph_byGene : to enable graph for gene affected by a CNV (def: no)
​	--graph_byChr : to enable graph by chromosome (def: no)
​	--graph_byCNV : to enable graph around each CNV (def: yes)
​	--CNV_fields [min/max/med/avg/std/Q1/Q3]: list of fields foreach region (comma separated) (default: none)

  *external tools path:
​	--bedtools [dir/file]: enter path to executable, if not installed as root or not in path
​	--samtools [dir/file]: enter  path to executable, if not installed as root
​	--picard [dir/file]: enter path to executable .jar, if not installed as root
​	--gatk [+/-dir/file]: cov analysis will be performed by gatk (default:bedtools; enter path to executable, if not installed as root)

  *general:
​	-x / --ram [int]: memory for gatk (in Go)
​	--cpu [int]: multi-thread for gatk (def: 1)
​	-v / --version : current version 
​	-h / --help : help



#### examples



```bash
./path/to/DeCovA -d path/to/bam_dir/ -r path/to/Refseq.txt -b path/to/targets.bed -M path/to/mut.list -t 20,50,100 -A -S -P -B -C
```



```bash
./path/to/DeCovA -f path/to/file1.bam -f path/to/file2.bam -r path/to/Refseq.txt -b path/to/targets.bed -M path/to/mut.list -t 20,50,100 -A -S -P -B -C
```



```bash
./path/to/DeCovA -f path/to/file.list -r path/to/Refseq.txt -b path/to/targets.bed -M path/to/mut.list -t 20,50,100 -A -S -P -B -C
```



```bash
./path/to/DeCovA -d path/to/bam_dir/ -r path/to/Refseq.txt -i GENE1,GENE2,NM_xxx1,NM_xxx2 -M path/to/mut.list -t 20,50,100 -A -S -P -B -C
```



```bash
./path/to/DeCovA -d path/to/bam_dir/ -r path/to/Refseq.txt -i genes.list -b path/to/targets.bed -t 20,50,100 -A -S -P -B -C
```

