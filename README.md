# Protist
Kenny Data

## 1.Trinity assembly
	LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	Trinity --seqType fq --max_memory 60G --CPU 20  \
	  	--left ./A006200108_131270_S158_L003_R1_001_val_1.fq.gz,./A006200108_131272_S159_L003_R1_001_val_1.fq.gz,./A006200108_131274_S160_L003_R1_001_val_1.fq.gz,./A006200108_131276_S161_L003_R1_001_val_1.fq.gz,./A006200108_131278_S162_L003_R1_001_val_1.fq.gz,./A006200108_131280_S163_L003_R1_001_val_1.fq.gz,./A006200108_131282_S164_L003_R1_001_val_1.fq.gz,./A006200108_131284_S165_L003_R1_001_val_1.fq.gz,./A006200108_131286_S166_L003_R1_001_val_1.fq.gz,./A006200108_131288_S167_L003_R1_001_val_1.fq.gz \
		--right ./A006200108_131270_S158_L003_R2_001_val_2.fq.gz,./A006200108_131272_S159_L003_R2_001_val_2.fq.gz,./A006200108_131274_S160_L003_R2_001_val_2.fq.gz,./A006200108_131276_S161_L003_R2_001_val_2.fq.gz,./A006200108_131278_S162_L003_R2_001_val_2.fq.gz,./A006200108_131280_S163_L003_R2_001_val_2.fq.gz,./A006200108_131282_S164_L003_R2_001_val_2.fq.gz,./A006200108_131284_S165_L003_R2_001_val_2.fq.gz,./A006200108_131286_S166_L003_R2_001_val_2.fq.gz,./A006200108_131288_S167_L003_R2_001_val_2.fq.gz \
	  --output ./trinity.out
	  
![image](https://user-images.githubusercontent.com/34407101/144247486-95676f5d-89c2-4b6b-852a-ef534322aeec.png)


## 2.Filter
### get longest transcription
	/home/shangao/.conda/pkgs/trinity-2.1.1-6/opt/trinity-2.1.1/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >unigene.fasta
	
### cd-hit cluster 0.9 identy
	/home/shangao/software/cd-hit-v4.8.1-2019-0228/cd-hit -i unigene.fasta -o unigene.cdhit.fasta -d 0 -M 100000 -T 48 -c 0.9
	
### remove the seq length less than 300bp
	seqkit seq -m 300

![image](https://user-images.githubusercontent.com/34407101/144271665-d36cc2c3-df98-4bab-9568-8c32e0a70d3c.png)

	
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
	
[length_distribution.pdf](https://github.com/LittleFrogHill/Protist/files/7635508/length_distribution.pdf)

	
### busco
 ![image](https://user-images.githubusercontent.com/34407101/144270793-111b81aa-ca99-4cfe-9b3a-eb1f8246f9cb.png)

### blobtools 
 ![protist cumulative](https://user-images.githubusercontent.com/34407101/144270893-a68d6898-4e72-4a78-9a6a-e809db1cde66.png)

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

## 4.merge
	find * -name '*.genes.results'> merge_result/merged.genes.quant.file
	find * -name '*.isoforms.results'> merge_result/merged.isoform.quant.file
	/home/shangao/Software/Assembly/trinityrnaseq-v2.11.0/util/misc/merge_RSEM_output_to_matrix.pl --rsem_files merge_result/merged.genes.quant.file --mode counts > merge_result/merged.genes.result.file

## 5.DEseq
	R
	library('DESeq2')
	colData <- read.table('sample.list', header=T,row.names=1)
	countData <-read.table('./merged.genes.result.file', row.names='id',header=T)
	dds <- DESeqDataSetFromMatrix(countData = round(countData),colData = colData, design = ~ group)
	dds <- dds[ rowSums(counts(dds)) > 50, ]
	dds2 <- DESeq(dds)
	res <- results(dds2)
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds2, normalized=TRUE)),by='row.names',sort=FALSE)
	write.table(file='DEseq.result.gene',resdata,col.names=T,row.names=F,sep=',',quote=F)
	
## 6.get q<0.05 DEGs
	head -n 1 DEseq.result.gene |sed 's/,/\t/g' > gene.title
	sed 's/,/\t/g' DEseq.result.gene |awk '$7<= 0.05 {print$0}' > gene
	cat gene.title gene > DEseq.adj.gene
	
	rld <- rlog(dds, blind=FALSE)
	vsd <- vst(dds, blind=FALSE)
	
	CCG Sample ID    Sample Name
	131270         Fisculla_asex_1
	131272         Fisculla_asex_2
	131274         Fisculla_asex_3
	131276         Fisculla_asex_4
	131278         Fisculla_asex_5
	131280         Fisculla_sex_6
	131282         Fisculla_sex_7
	131284         Fisculla_sex_8
	131286         Fisculla_sex_9
	131288         Fisculla_sex_10	
	
## 7.Plot

### MA_plot

#### In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

	pdf("Ma.pdf")
	plotMA(res, ylim=c(-10,10),alpha =0.005)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/144275542-4497977e-23de-4ada-96aa-c172c727a0c1.png)

### Plot counts
#### It can also be useful to examine the counts of reads for a single gene across the groups. A simple function for making this plot is plotCounts, which normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in intgroup, where more than one variable can be specified. Here we specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by rowname or by numeric index.

	pdf("count.pdf")
	plot( assay(rld)[, 5:6], col="#00000020", pch=20, cex=0.3 )
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/144279240-df58206b-88f5-47da-8aa4-6a57d1c8da3c.png)

### rlog_plot
	
	pdf("rlog.pdf")
	plot( assay(rld)[, 5:6], col="#00000020", pch=20, cex=0.3 )
![image](https://user-images.githubusercontent.com/34407101/144280051-ff83d978-5909-4b0a-b1a2-26a46299a61b.png)

### sample_distances_plot
	
	sampleDistMatrix <- as.matrix( sampleDists )
	rownames(sampleDistMatrix) <- colnames(rld)
	colnames(sampleDistMatrix) <- c(rep('asex', times=5), rep('sex', times=5))
	library( "gplots" )
	library( "RColorBrewer" )
	colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	
	pdf("heatmap_sample_distances.pdf")
	heatmap.2( sampleDistMatrix, trace="none", col=colours)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/144289021-d7d3694a-fe62-44e5-9fc9-08c83e1e6077.png)

### PCA_plot
	pdf("pca_sample_distances.pdf")
	plotPCA( rld, intgroup = ("group") )
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/144292107-ae651dd5-4082-48e4-827f-062389da59ec.png)

### SD_top500_plot
	
	pdf("heatmap_top500SDgenes.pdf")
	topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )
	heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[colData(rld)$group ] )
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/144292473-5675ec31-3e2f-440d-8110-b1384caf3264.png)

#### SD_top35_plot
![image](https://user-images.githubusercontent.com/34407101/144296237-7cbff44a-9a8b-4f0e-b3d9-9aed9dd0321f.png)
#### SD_top50_plot
![image](https://user-images.githubusercontent.com/34407101/144296431-341370eb-2278-4302-8f3c-638a6eee3674.png)
#### SD_top100_plot
![image](https://user-images.githubusercontent.com/34407101/144296381-fcd0ff76-3e81-48f2-bcd1-68dc455855dd.png)

https://jokergoo.github.io/ComplexHeatmap-reference/book/other-tricks.html

### Volcano
	pdf("volcano.pdf")	
	EnhancedVolcano(res,lab = rownames(res),x = 'log2FoldChange',y = 'pvalue',title = 'sex vs asexual', FCcutoff = 0.5, pointSize = 1.0,labSize = 1.0)
	dev.off()	
![image](https://user-images.githubusercontent.com/34407101/144298869-4b265913-7ddd-426b-9148-46fe01bfaf8c.png)


## 8.Annotation with docker KOBAS

	docker run -it -v /home/shangao:/gpfs -v /RAID/Data/databases/kobas_DB/seq_pep:/opt/kobas-3.0/seq_pep -v /RAID/Data/databases/kobas_DB/sqlite3:/opt/kobas-3.0/sqlite3 -v /NVME2/Scratch/gaoshan/:/NVME2/Scratch/gaoshan/ -v /RAID/Data/gaoshan/:/RAID/Data/gaoshan/  kobas
	
	annotate.py -i ./DEGs.adj.gene.filter500.fa.transdecoder.pep -s hsa -t fasta:pro -o anno_seq.tsv -e 1e-5 -r 1 -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep
	
	identify.py -f anno_seq.tsv -o identify.out -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep

### If you don’t know the code for your species it can be found here: 
https://www.kegg.jp/kegg/catalog/org_list.html
### If your species of interest is not available then you should choose the code for the closest-related species available
https://agbase-docs.readthedocs.io/en/latest/kobas/using_kobas_cmd.html
https://agbase-docs.readthedocs.io/en/latest/kobas/using_kobas_cmd.html
![image](https://user-images.githubusercontent.com/34407101/154854832-d13f301d-d1c3-4c9c-9727-091d0ad76597.png)


## 9.clusterpr
http://yulab-smu.top/clusterProfiler-book/chapter12.html#bar-plot

https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html

Human_database

	R
	d <- read.table('hsa_anno.out')
	geneList <- d[,2]
	names(geneList) <- as.character(d[,1])
	geneList <- sort(geneList, decreasing = TRUE)
	gene <- names(geneList)[abs(geneList) > 0.5]
	library(clusterProfiler)
	library(org.Hs.eg.db)
	kk <- enrichKEGG(gene         = gene,organism     = 'hsa',pvalueCutoff = 0.05)
	head(kk)
	
### barplot
![image](https://user-images.githubusercontent.com/34407101/123746111-7c3b3c80-d8b1-11eb-8b65-5d1c10a29a63.png)
### dotplot
![image](https://user-images.githubusercontent.com/34407101/123746159-8cebb280-d8b1-11eb-9e27-9d8bd082c41a.png)
### Network
	> pdf('Network1.pdf',pointsize=2,width=100,height=100)
	> cnetplot(edox, foldChange=geneList,showCategory=20)
	> dev.off()
	
	> pdf('Network2.pdf',pointsize=2,width=100,height=100)
	> cnetplot(kk, categorySize="pvalue", foldChange=geneList,)
	x=             showCategory=  foldChange=    layout=        ...=           
	> cnetplot(kk, categorySize="pvalue", foldChange=geneList,showCategory=20)
	> dev.off()
	
	> pdf('Network3.pdf',pointsize=2,width=100,height=100)
	> cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE,showCategory=20)
	Warning message:
	ggrepel: 1 unlabeled data points (too many overlaps). Consider increasing max.overlaps 
	> dev.off()
	
	> pdf('Network_Enriched1.pdf',pointsize=1,width=50,height=50)
	> cnetplot(edox, node_label="all",showCategory=20)
	> dev.off()

### Heatmap
	> pdf('Heatmap.pdf',pointsize=1,width=80,height=10)
	> heatplot(edox, foldChange=geneList)
	> dev.off()
	
### Enrich
	> pdf('Enrich.pdf',pointsize=5,width=50,height=50)
	> emapplot(kk, pie_scale=1.5,layout="kk")
	> dev.off()

	
Ptm_database
	
	d <- read.table('ptm.annoDEGs.out')
	geneList <- d[,2]
	names(geneList) <- as.character(d[,1])
	geneList <- sort(geneList, decreasing = TRUE)
	gene <- names(geneList)[abs(geneList) > 0.5]
	library(clusterProfiler)
	
	
	
	pdf('bar1.pdf')
	barplot(kk1, showCategory=20)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/123746304-bf95ab00-d8b1-11eb-97d8-0d40c73ab091.png)
	
	pdf('dot.pdf')
	dotplot(kk, showCategory=30) + ggtitle("dotplot for ptm_database")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/123746738-582c2b00-d8b2-11eb-90ab-7df079da3067.png)

## 10.	Phylogenetic
### blast and prepare seqs from different species

### aln
	mafft rad50.fa > rad50.aln
### IQ-tree
	iqtree-omp -s rad50.aln -st AA -nt 16 -quiet -bb 1000 -m TESTNEW -msub nuclear


[results.pdf](https://github.com/LittleFrogHill/Protist/files/7640959/results.pdf)
![image](https://user-images.githubusercontent.com/34407101/144400820-cf8f612e-c734-43e1-a975-c534ccfa1ba6.png)
![image](https://user-images.githubusercontent.com/34407101/144400849-d8594782-cc34-43cf-afe5-c979d7417e13.png)

	
	
	
	
	
	
