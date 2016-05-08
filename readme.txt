DeCovA_1.2g.pl

DeCovA script requires at least R and bedtools softawres to be installed;
additionnally, it can use picard tools (for deduplication), samtools (for mapq filter), and GATK (alternatively to bedtools; required if a base-q filter is needed).


a refSeq file need to be provided (-r option), for all the options that use gene coordinates:
you can download it at UCSC:
$ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
	(may change hg19 to hg38
	and /refGene.txt.gz to /ensGene.txt.gz)
then unzip 
$ gunzip refGene.txt.gz

ftp://ftp.ensembl.org/pub/release-79/gtf/homo_sapiens/

DeCovA can be executed via command-line execution of the perl DeCovA-x.pl script:
$ perl path/to/script/DeCovA-x.pl [options]
if the script has been changed to executable:
$ chmod 755 path/to/script/DeCovA-x.pl
then:
$ path/to/script/DeCovA-x.pl [options]
(in the current directory: $ ./DeCovA-x.pl [options])
and if the dir has been added to the path:
$ echo 'export PATH=$PATH:/home/me/path/to/script' >> /home/me/.bashrc
then
$ DeCovA-x.pl [options]



list of parameters:

usage: perl /path/to/DeCovA-x [option]
	-f / --file : list of bam files, either command line (comma separated), or as a file (\".list\")
	-d / --dir : directory where to find bam files (optional)
	-s / --suffix : suffix to add before opening bam files (optional)
	-r / --refseq : RefSeq file (required)
	-b / --bed : bed file, used to analyse depth coverage
	-m / --mut : mut file, used to plot known mutation: format: chr	pos(1-based)	info
	-o / --outfile : text outfile name (optional, default: covBySample.txt)
	-O / --outdir : out directory (optional ; default: folder named with date)
	-i / --id : list of of refseq genes/transcripts id, either command line (comma separated), or as a file (\".list\")
	-N / --nonCoding : analyse also Non coding transcripts (default: no)
	-U / --noUTR : does not take into account UTR regions for all analyses (default: yes)
	-u /--noUTRinTxt : does not take into account UTR regions for summary txt file (default: yes)
	-t / --depthThreshold : depth thresholds, comma separated (at least 1 required)
	-T / --printThreshold : depth threshold below which graph will be printed (one of the thresholds in opt -t)
	--maxDepth : max depth value when printing graph (optional)
	-l / --expand2val : length to add out of exons, to inform on intron-exon junction (default: 20)
	-L / --expand2bed : expand length of analysed regions to bed coord, if -l < bed (default: no)
	--expandForTxt : does take into account expanded length (from -l and -L) for summary txt file (default: no)
	-R / --noReverse : does not reverse regions if sens of transcript = (-) (default: yes)
	-S / --noGraphSum : does not perform graphSum (default: yes)
	-A / --noAllSample : does not perform graphAllSample (default: yes)
	-X / --noBySample : does not perform graphBySample (default: yes)
	-E / --noEachTranscript : does not print All transcripts on same file, in graphBySample (default: yes)
	-M / --noDepthMut : does not print, foreach file, depth at known mutations provided by opt -M (default: yes)
	-P / --noCovPlot : does not perform covPlots (default: yes)
	-B / --noCovBed : does not output cov of bed intervals (default: yes)
	--printReports : print all uncovered intervals in 1 txt file per transcript (default: no)
	--noSummary : does not print summary txt file (default: yes)
	--binPlot : bin width for covPlot (default=10)
	--maxPlot : max depth for covPlot (default=100)
	--genePlot : performs plots for regions extracted from genes coord, not only for bed intervals
	--interPlot : intersection covPlot (default: no; enter yes to perform)
	--bedtools : enter path to executable, if not installed as root or not in path
	--samtools : enter path to executable, if not installed as root
	--picard : enter path to executable .jar, if not installed as root
	--gatk : cov analysis will be performed by gatk (default:bedtools; enter path to executable, if not installed as root)
	-g / --genome : path to genome.fa file, if available (required if using GATK)
 	--dedup : do not take in account dup reads (optional; default keep all reads; enter \"do\" to perform deduplication) 
	--mbq : minimum base quality (optional; default 0; requires gatk)
	--mmq : minimum mapping quality (optional; default 0)
	--nGraph : max nber of graphs per sheet (default : all samples or all transcripts)
	-z / --compress : archive output folder
	-v / --version : current version 
	-h / --help : help


examples:

$ perl DeCovA-x.pl -d path/to/dir/ -r path/to/Refseq.txt -bed path/to/targets.bed -M path/to/mut.list -t 20,50,100

$ perl DeCovA-x.pl -f path/to/file1.bam,path/to/file2.bam -r path/to/Refseq.txt -bed path/to/targets.bed -M path/to/mut.list -t 20,50,100

$ perl DeCovA-x.pl -f path/to/file.list -r path/to/Refseq.txt -bed path/to/targets.bed -M path/to/mut.list -t 20,50,100

$ perl DeCovA-x.pl -d path/to/dir/ -r path/to/Refseq.txt -i GENE1,GENE2,NM_xxx1,NM_xxx2 -M path/to/mut.list -t 20,50,100

$ perl DeCovA-x.pl -d path/to/dir/ -r path/to/Refseq.txt -i GENE1,GENE2,NM_xxx1,NM_xxx2 -bed path/to/targets.bed -M path/to/mut.list -t 20,50,100



