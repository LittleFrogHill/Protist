# Protist
Kenny Data

## 1.Trinity assembly
	LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	Trinity --seqType fq --max_memory 60G --CPU 20  \
	  	--left ./A006200108_131270_S158_L003_R1_001_val_1.fq.gz,./A006200108_131272_S159_L003_R1_001_val_1.fq.gz,./A006200108_131274_S160_L003_R1_001_val_1.fq.gz,./A006200108_131276_S161_L003_R1_001_val_1.fq.gz,./A006200108_131278_S162_L003_R1_001_val_1.fq.gz,./A006200108_131280_S163_L003_R1_001_val_1.fq.gz,./A006200108_131282_S164_L003_R1_001_val_1.fq.gz,./A006200108_131284_S165_L003_R1_001_val_1.fq.gz,./A006200108_131286_S166_L003_R1_001_val_1.fq.gz,./A006200108_131288_S167_L003_R1_001_val_1.fq.gz \
		--right ./A006200108_131270_S158_L003_R2_001_val_2.fq.gz,./A006200108_131272_S159_L003_R2_001_val_2.fq.gz,./A006200108_131274_S160_L003_R2_001_val_2.fq.gz,./A006200108_131276_S161_L003_R2_001_val_2.fq.gz,./A006200108_131278_S162_L003_R2_001_val_2.fq.gz,./A006200108_131280_S163_L003_R2_001_val_2.fq.gz,./A006200108_131282_S164_L003_R2_001_val_2.fq.gz,./A006200108_131284_S165_L003_R2_001_val_2.fq.gz,./A006200108_131286_S166_L003_R2_001_val_2.fq.gz,./A006200108_131288_S167_L003_R2_001_val_2.fq.gz \
	  --output ./trinity.out

## 2.Filter
### get longest transcription
	/home/shangao/.conda/pkgs/trinity-2.1.1-6/opt/trinity-2.1.1/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >unigene.fasta
	
### cd-hit cluster 0.9 identy
	/home/shangao/software/cd-hit-v4.8.1-2019-0228/cd-hit -i unigene.fasta -o unigene.cdhit.fasta -d 0 -M 100000 -T 48 -c 0.9
	
### check length
	/NVME/Software/Assembly/trinityrnaseq-v2.11.0/util/misc/fasta_seq_length.pl unigene.cdhit.fasta > unigene.cdhit.fasta.length
	R
	data <- read.table("unigene.cdhit.fasta.length",header=T)
	data[,2][which(as.numeric(data[,2])>=2000)]<-2000

	library(ggplot2)
	pdf("length_distribution.pdf",height=7,width=10)
	ggplot(as.data.frame(data), aes(x = as.numeric(data[,2])))+geom_histogram(binwidth =100)+
	xlab("Transcripts Length Interval")+
	ylab("Number ofTranscripts")+
	labs(title="Transcripts Length Distribution")+
	scale_x_continuous(breaks=seq(100,2000,by=100),
	labels=c("100","200","300","400","500","600","700","800","900","1000","1100","1200","1
	300","1400","1500","1600","1700","1800","1900",">=2000"))
	dev.off()

## 3.RSEM mapping
	LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	for i in 131270_S158 \
		131272_S159 \
		131274_S160 \
		131276_S161 \
		131278_S162 \
		131280_S163 \
		131282_S164 \
		131284_S165 \
		131286_S166 \
		131288_S167
	do
	echo ''perl /home/shangao/Software/Assembly/trinityrnaseq-v2.11.0/util/align_and_estimate_abundance.pl \
		--transcripts ../unigene.cdhit0.9.fasta --seqType fq \
		--left /home/shangao/Data/protist/clean_data/A006200108_${i}_L003_R1_001_val_1.fq.gz \
		--right /home/shangao/Data/protist/clean_data/A006200108_${i}_L003_R2_001_val_2.fq.gz \
		--est_method RSEM \
		--aln_method bowtie2 \
		--trinity_mode  \
		--prep_reference  \
	       	--output_dir RSEMResult/${i} \
		--thread_count 5'' > rsem_shell/$i.sh
	mkdir RSEMResult/$i
	nohup sh rsem_shell/$i.sh > rsem_shell/$i.log &
	done




## 1.Annotation with docker KOBAS

	docker run -it -v /home/shangao:/gpfs -v /home/shangao/Data/kobas_DB/seq_pep:/opt/kobas-3.0/seq_pep -v /home/shangao/Data/kobas_DB/sqlite3:/opt/kobas-3.0/sqlite3 -v /NVME2/Scratch/gaoshan/:/NVME2/Scratch/gaoshan/ kobas
	
	annotate.py -i ./DEGs.adj.gene.filter500.fa.transdecoder.pep -s hsa -t fasta:pro -o anno_seq.tsv -e 1e-5 -r 1 -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep
	
	identify.py -f anno_seq.tsv -o identify.out -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep
