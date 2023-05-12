# Protist
## *De novo Fisculla terrestris* Transcription Data Analysis Pipeline

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
![image](https://user-images.githubusercontent.com/34407101/171874537-00e64280-cf09-4e9a-95a3-33abbbc214b4.png)
[length_distribution.pdf](https://github.com/LittleFrogHill/Protist/files/7635508/length_distribution.pdf)
	
### busco
	busco -i ../new_reference/unigene.cdhit0.9.over300.fasta \
		--auto-lineage-euk \
		-o busco \
		-m tran \
		--cpu 40

 ![image](https://user-images.githubusercontent.com/34407101/144270793-111b81aa-ca99-4cfe-9b3a-eb1f8246f9cb.png)

### blobtools 
	/home/shangao/Software/blobtoolkit/blobtools2/blobtools create \
	--fasta ../new_reference/unigene.cdhit0.9.fasta \
	/home/shangao/Data/clean_data_protsist/trinity/longest_transcription/blobtools/protist

	/home/shangao/Software/blobtoolkit/blobtools2/blobtools add \
	--taxrule bestsumorder \
	--taxdump /RAID/Data/databases/taxdump/ \
	--hits /home/shangao/Data/clean_data_protsist/trinity/longest_transcription/blobtools/longest.ncbi1.blastn.out \
	--busco busco/run_eukaryota_odb10/full_table.tsv \
	/home/shangao/Data/clean_data_protsist/trinity/longest_transcription/blobtools/protist
	
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

## 4.Merge transcription
	find * -name '*.genes.results'> merge_result/merged.genes.quant.file
	find * -name '*.isoforms.results'> merge_result/merged.isoform.quant.file
	/home/shangao/Software/Assembly/trinityrnaseq-v2.11.0/util/misc/merge_RSEM_output_to_matrix.pl --rsem_files merge_result/merged.genes.quant.file --mode counts > merge_result/merged.genes.result.file

## 5.DEseq calculate DEGs
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
	
## 6. Filter q<0.05 DEGs
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
	
## 7.Plots

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
	
	pdf("heatmap_sample_distances_withbar.pdf")
	heatmap.2( sampleDistMatrix, trace="none", col=colours,dendrogram = "row",ColSideColors =plot_color,srtCol = 45, offsetCol = -0.5,cexRow = 0.8, cexCol = 1.0)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171907670-1a1c0172-1961-43f9-8294-a5f62913d28c.png)


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


## 9.clusterProfiler
http://yulab-smu.top/clusterProfiler-book/chapter12.html#bar-plot
https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html

### Human_database annotated

	R
	d <- read.table('hsa_anno.out')
	geneList <- d[,2]
	names(geneList) <- as.character(d[,1])
	geneList <- sort(geneList, decreasing = TRUE)
	gene <- names(geneList)[abs(geneList) > 0.5]
	library(clusterProfiler)
	library(org.Hs.eg.db)
	
### KEGG enrichment
#### KEGG pathway over-representation analysis
	kk <- enrichKEGG(gene         = gene,organism     = 'hsa',pvalueCutoff = 0.05)
	head(kk)
	
	pdf("KEGG_enrichment.pdf",width=8)
	dotplot(kk, showCategory=30) + ggtitle("KEGG pathway over-representation analysis")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171643777-3fb91524-2ff7-41e2-a41f-9031c7497585.png)


#### KEGG pathway gene set enrichment analysis
	 kk2 <- gseKEGG(geneList = geneList,organism = 'hsa',minGSSize = 20,maxGSSize=50,pvalueCutoff = 0.05,verbose = FALSE)
	 kk2$Description
	 [1] "Epstein-Barr virus infection"                 
	 [2] "Proteasome"                                   
	 [3] "Toxoplasmosis"                                
	 [4] "DNA replication"                              
	 [5] "Citrate cycle (TCA cycle)"                    
	 [6] "Gap junction"                                 
	 [7] "Ribosome biogenesis in eukaryotes"            
	 [8] "Cell cycle"                                   
	 [9] "Lipid and atherosclerosis"                    
	[10] "Chemical carcinogenesis - receptor activation"
	[11] "Cysteine and methionine metabolism"           
	[12] "Glutathione metabolism"                       
	[13] "N-Glycan biosynthesis"                        
	[14] "Tight junction"                               
	[15] "cAMP signaling pathway"                       
	[16] "Pathogenic Escherichia coli infection"        
	[17] "Arginine and proline metabolism"              
	[18] "Measles"                                      
	[19] "Oocyte meiosis"                               
	[20] "Thyroid hormone signaling pathway"            
	[21] "Phagosome"                                    
	[22] "Adrenergic signaling in cardiomyocytes"

	pdf("KEGG_gse_enrichment.pdf",width=7.5)
	dotplot(kk2, showCategory=30) + ggtitle("KEGG pathway gene set enrichment analysis")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171643378-9bd8a86a-98f8-451e-b1c5-712ff5cd4baa.png)

	pdf("KEGG_GSEA.pdf")
	gseaplot2(kk2, geneSetID = c(4,19,8), pvalue_table = TRUE,color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171684961-66c7ad25-8159-41c9-9bb4-001827d90828.png)


#### KEGG module over-representation analysis
	mkk <- enrichMKEGG(gene = gene,organism = 'hsa',pvalueCutoff = 1,qvalueCutoff = 1,keyType= "kegg")
![image](https://user-images.githubusercontent.com/34407101/171635548-13a962b6-c849-421b-ad35-df0ec1ce8c9a.png)

#### KEGG module gene set enrichment analysis
	mkk2 <- gseMKEGG(geneList = geneList,organism = 'hsa',pvalueCutoff = 1)
![image](https://user-images.githubusercontent.com/34407101/171635836-418d2848-a1a2-407b-a7f7-aa01b89b8063.png)

#### KEGG Pathway
	library("pathview")
	hsa04114 <- pathview(gene.data  = geneList,pathway.id = "hsa04114",species    = "hsa", limit=list(gene = 2, cpd = 1))
![image](https://user-images.githubusercontent.com/34407101/171658956-03787493-36c4-428f-833f-f5183fb99eec.png)
	hsa03030 <- pathview(gene.data  = geneList,pathway.id = "hsa03030",species    = "hsa", limit=list(gene = 2, cpd = 1))
	hsa03430 <- pathview(gene.data  = geneList,pathway.id = "hsa03430",species    = "hsa", limit=list(gene = 2, cpd = 1))
	hsa03420 <- pathview(gene.data  = geneList,pathway.id = "hsa03420",species    = "hsa", limit=list(gene = 2, cpd = 1))
	hsa04110 <- pathview(gene.data  = geneList,pathway.id = "hsa04110",species    = "hsa", limit=list(gene = 2, cpd = 1))
	
### GO enrichment
#### GO classification
	ggo <- groupGO(gene = gene,OrgDb= org.Hs.eg.db,ont= "bp",level= 3, readable = TRUE)
	pdf("GO_group.pdf",height=20)
	barplot(ggo, showCategory=100,drop=TRUE) + ggtitle("GO classification level=3 100 items")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171648035-3f8a58c4-5000-4cf7-92bc-eca5bb623e67.png)
![image](https://user-images.githubusercontent.com/34407101/171654203-42944571-f30d-4374-bf6a-eb604575c5cd.png)
	
	ggo2 <- groupGO(gene= gene,OrgDb= org.Hs.eg.db,ont= "bp",level= 2, readable = TRUE)
	pdf("GO_group_level2.pdf")
	barplot(ggo2, showCategory=100,drop=TRUE) + ggtitle("GO classification level=2")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171659558-f8bac429-971c-43b3-b878-8b4886c8ed6a.png)

#### GO over-representation analysis
	ego <- enrichGO(gene= gene,universe= names(geneList),OrgDb= org.Hs.eg.db,ont= "bp",pAdjustMethod = "BH",qvalueCutoff = 0.05,readable= TRUE)
	pdf("GO_enrichment.pdf",width=30)
	goplot(ego)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171630414-c20bc4ec-939e-4708-9df1-acf4971d3344.png)
	pdf("GO_enrich_cnetplot.pdf")
	cnetplot(ego, categorySize="pvalue", foldChange=geneList,showCategory="negative regulation of cell cycle G2/M phase transition")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171679011-f100607a-cf5e-4e05-8517-18148573a8b5.png)
	pdf("GO_enrich_graph.pdf",width=50,height=50)
	plotGOgraph(ego)
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171681172-7f13546b-b614-4c2c-989f-d660b6f2a2ac.png)
	
#### GO Gene Set Enrichment Analysis (GSEA)
	ego3 <- gseGO(geneList=geneList,OrgDb=org.Hs.eg.db,ont= "bp",minGSSize=100,maxGSSize= 500,pvalueCutoff = 0.05,verbose = FALSE)
	pdf("GO_GSEA.pdf")
	gseaplot2(ego3, geneSetID = (1,13,118), pvalue_table = TRUE,color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
	dev.off()
![image](https://user-images.githubusercontent.com/34407101/171684448-b3bc8cc6-0222-4d05-9d2f-efc386baf1d2.png)

### Ptm_database

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

## 10. Meiosis-related genes Blast
### download the Meiosis-related genes and blastp to our database
	cat blastp2pro.sh
	for i in Acrispumella_msimbaziensis \
		Apoikiospumella_mondseeiensis \
		Cornospumella_fuschlensis \
		Dinobryon_sp_1 \
		Dinobryon_sp_2 \
		Epipyxis_sp \
		Ochromonas_or_Spumella_sp \
		Pedospumella_encystans \
		Pedospumella_sinomuralis \
		Poterioochromonas_malhamensis \
		Poteriospumella_lacustris_1 \
		Poteriospumella_lacustris_2 \
		Poteriospumella_lacustris_3 \
		Spumella_bureschii \
		Spumella_lacusvadosi \
		Spumella_vulgaris \
		Synura_sp \
		Uroglena_sp
	do
	cd $i
	cat *.fasta > $i.fasta
	sed -i 's/\ /-/g' $i.fasta
	#makeblastdb -in $i.fasta -dbtype prot -out $i
	blastp -db /home/shangao/Data/clean_data_protsist/trinity/longest_transcription/trandecoder/unigene.cdhit0.9.fasta.transdecoder.pep.cdhit -query $i.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-4 -num_threads 10 > ../$i.2pro.blastn.out
	cd ..
	done
![image](https://user-images.githubusercontent.com/34407101/172344425-d04b5658-1eda-4788-8243-cb331d41aa32.png)

## 11.	Phylogenetic
### blast and prepare seqs from different species

### aln
	mafft rad50.fa > rad50.aln
### IQ-tree
	iqtree-omp -s rad50.aln -st AA -nt 16 -quiet -bb 1000 -m TESTNEW -msub nuclear

[results.pdf](https://github.com/LittleFrogHill/Protist/files/7640959/results.pdf)
![image](https://user-images.githubusercontent.com/34407101/144400820-cf8f612e-c734-43e1-a975-c534ccfa1ba6.png)
![image](https://user-images.githubusercontent.com/34407101/144400849-d8594782-cc34-43cf-afe5-c979d7417e13.png)

## 12.	clustertree(https://mp.weixin.qq.com/s?__biz=MzI5NjUyNzkxMg==&mid=2247485611&idx=1&sn=69990a569c623730583b56e54d55b58b&scene=21#wechat_redirect)

	conda activate BRAKER
	R
	library(tximport)
	library(DESeq2)
	library(ggplot2)
	library(ggrepel)
	library(dendextend)
	library(grid)
	library(gridExtra)
	library(tibble)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(forcats)
	library(ochRe)
	library(ggpubr)
	library(goseq)
	library(GO.db)
	library(GOfuncR)
	library(qvalue)
	library(stringr)
	library(ComplexHeatmap)
	library(RColorBrewer)
	library(circlize)
	library(seqinr)
	library(UpSetR)
	library(purrr)
	# Calculate DEGs
	colData <- read.table('../../../sample.list', header=T,row.names=1)
	countData <-read.table('../../../merged.genes.result.file', row.names='id',header=T)
	dds <- DESeqDataSetFromMatrix(countData = round(countData),colData = colData, design = ~ group)
	dds <- dds[ rowSums(counts(dds)) > 50, ]
	dds2 <- DESeq(dds)
	res <- results(dds2)
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds2, normalized=TRUE)),by='row.names',sort=FALSE)
	rld <- rlog(dds, blind=FALSE)
	
	# Extracting transformed count values
	vsd <- vst(dds, blind=FALSE)
	vsd_matrix <- assay(vsd)
	
	# Cluster
	# Cutoff readcount>50, padj 0.05
	as.data.frame(res) %>%
	filter(padj < 0.05) -> res_diff
	as.data.frame(res) %>% filter(padj < 0.05) -> res_diff
	head(res)
	head(res_diff)
	dge_set <- c(rownames(res_diff))
	dge_set <- unique(dge_set)
	
	# Generate z-score matrix for heatmaps
	vsd_matrix_dge <- vsd_matrix[dge_set,]
	heat <- t(scale(t(vsd_matrix)))
	heat <- heat[dge_set,]
	
	# Clusters rows by Pearson correlation as distance method
	hc <- hclust(as.dist(1 - cor(t(as.matrix(vsd_matrix_dge)))))
	my_transcript_partition_assignments <- cutree(hc, h = 80/100*max(hc$height))
	
	# Visualise dendrogram with clusters
	clust.cutree <- dendextend:::cutree(as.dendrogram(hc), h = 80/100*max(hc$height), order_clusters_as_data = FALSE)
	idx <- order(names(clust.cutree))
	clust.cutree <- clust.cutree[idx]
	idx <- order(names(clust.cutree))
	clust.cutree <- clust.cutree[idx]
	df.merge <- merge(my_transcript_partition_assignments, clust.cutree, by = "row.names")
	df.merge.sorted <- df.merge[order(df.merge$y),]
	lbls <- unique(df.merge.sorted$x)
	dend1 <-color_branches(as.dendrogram(hc), h = 80/100*max(hc$height), groupLabels = lbls)
	pdf(file = paste0("./clustering_dendrogram.1.pdf"), width = 12, height = 5)
	lot(dend1, leaflab = "none")
	bline(h=80/100*max(hc$height), lty = 2, col="grey")
	dev.off()
	
![image](https://user-images.githubusercontent.com/34407101/216330344-aa3f8ee0-5775-4f98-8b05-321ebbfdf2c5.png)

	# Make factor_labeling.txt for goseq
	data.frame(my_transcript_partition_assignments) %>%
	rownames_to_column(var = "transcripts") -> factor_labeling
	colnames(factor_labeling) <- c("transcript", "cluster")
	#write.table(factor_labeling, 
	#            file = paste0(mydir, "/Module_5/GO_analysis/factor_labeling.txt"),
	#            row.names = FALSE,
	#            col.names = FALSE,
	#            quote = FALSE,
	#            sep = "\t")
	
	# Make list of clusters
	clusterlist <- list()
	for (i in c(1:max(my_transcript_partition_assignments))) {
	  cluster <- heat[(my_transcript_partition_assignments == i),]
	  clusterlist[[i]] <- cluster
	}

	# Plot clusters as boxplots
	p <- list()
	splan <- 3 - 1L
	for (i in 1:max(my_transcript_partition_assignments)) {
	   box <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "transcript")
	   box_longer <- pivot_longer(box, cols = !transcript, names_to = "sample")
	  box_longer %>%filter(sample == "X131270_S158" | sample == "X131272_S159" | sample == "X131274_S160"| sample == "X131276_S161" | sample == "X131278_S162") %>%
	  mutate(cat = "asex") -> box_longer1
	  box_longer %>% filter(sample == "X131280_S163" | sample == "X131282_S164" | sample == "X131284_S165"| sample == "X131286_S166" | sample == "X131288_S167") %>%
	  mutate(cat = "sex") -> box_longer2
	  comb <- rbind(box_longer1, box_longer2)
	  comb$cat <-factor(comb$cat, levels = c("asex", "sex"))
	 g <- ggplot(comb, aes(x = cat, y = value)) +geom_boxplot(outlier.size = 0,outlier.shape = NA,alpha = 0.5,color = "#252A52",fill = "#252A52") +stat_smooth(aes(x = cat, y = value, group = "#252A52"), color = "#252A52",se = TRUE,method = "lm", formula = y~poly(x, splan)) +ggtitle(paste("cluster", i, "|", "contigs:", nrow(clusterlist[[i]]))) +scale_y_continuous(limits = c(-2, 2)) +
	      ylab("z-score") +
	      xlab(NULL) +
	      theme(legend.position = "none", panel.background = element_rect(colour = "darkgrey", size=1))
	    p[[i]] <- ggplotGrob(g)
	  }
	pdf(file = paste0("./cluster_boxplots_log2FC1.pdf"), width = 15, height = 3)
	grid.arrange(grobs = p, ncol = 2)
	dev.off()
	
![image](https://user-images.githubusercontent.com/34407101/216330407-570da456-1aa5-4f48-889f-6ef852eb2e9d.png)

	# GO import and heatmap
	library(clusterProfiler)
	library(org.Hs.eg.db)
	d1 <- read.table('cluster1.out')
	d2 <- read.table('cluster2.out')
	geneList1 <- d1[,2]
	geneList2 <- d2[,2]
	names(geneList1) <- as.character(d1[,1])
	names(geneList2) <- as.character(d2[,1])
	geneList1 <- sort(geneList1, decreasing = TRUE)
	geneList2 <- sort(geneList2, decreasing = TRUE)
	gene1 <- names(geneList1)[abs(geneList1) > 0]
	gene2 <- names(geneList2)[abs(geneList2) > 0]
	ego1 <- enrichGO(gene= gene1,universe= total_gene,OrgDb= org.Hs.eg.db,ont= "all",pAdjustMethod = "BH",pvalueCutoff=0.05,readable= TRUE)
	ego2 <- enrichGO(gene= gene2,universe= total_gene,OrgDb= org.Hs.eg.db,ont= "all",pAdjustMethod = "BH",pvalueCutoff=0.05,readable= TRUE)
	kk1 <- enrichKEGG(gene=gene1,organism= 'hsa',pvalueCutoff = 0.05)
	kk2 <- enrichKEGG(gene=gene2,organism= 'hsa',pvalueCutoff = 0.05)
	
	library(enrichplot)
	pdf("KEGG_enrichment_cluster1.pdf",width=10,height=22)
	dotplot(kk1, showCategory=500) + ggtitle("KEGG pathway analysis of Cluster1")
	dev.off()
	pdf("KEGG_enrichment_cluster2.pdf",width=10,height=22)
	dotplot(kk2, showCategory=500) + ggtitle("KEGG pathway analysis of Cluster2")
	dev.off() 
	
	pdf("GO_enrichment_cluster2.pdf",height=25)
	dotplot(ego2, split="ONTOLOGY",showCategory=300,  font.size = 7) + facet_grid(ONTOLOGY~., scale="free", space = "free")+ ggtitle("GO term of Cluster2")
	dev.off()                                  
	pdf("GO_enrichment_cluster1.pdf",height=25)
	dotplot(ego1, split="ONTOLOGY",showCategory=300,  font.size = 7) + facet_grid(ONTOLOGY~., scale="free", space = "free")+ ggtitle("GO term of Cluster1")
	dev.off()
	
	write.table(ego1, file = "GO_enrichment_cluster1.txt",sep = "\t", row.names = F,col.names = T)
	write.table(ego2, file = "GO_enrichment_cluster2.txt",sep = "\t", row.names = F,col.names = T)
	write.table(kk2, file = "KEGG_enrichment_cluster2.txt",sep = "\t", row.names = F,col.names = T)
	write.table(kk1, file = "KEGG_enrichment_cluster1.txt",sep = "\t", row.names = F,col.names = T)

	ggo1 <- groupGO(gene = gene1,OrgDb= org.Hs.eg.db,ont= "bp",level= 3)
	ggo2 <- groupGO(gene = gene2,OrgDb= org.Hs.eg.db,ont= "bp",level= 3)
	write.table(ggo1, file = "ggo_cluster1_level3.txt",sep = "\t", row.names = F,col.names = T)
	write.table(ggo2, file = "ggo_cluster2_level3.txt",sep = "\t", row.names = F,col.names = T)
	
	heat1 <- read.table('heatmap/ggo_cluster1_exp',col.name=c("genename","X131270_S158","X131272_S159","X131274_S160","X131276_S161","X131278_S162","X131280_S163","X131282_S164","X131284_S165","X131286_S166","X131288_S167"))
	heat1_exp<-heat1[,2:11]
	row.names(heat1_exp) <- AnnotationDbi::select(org.Hs.eg.db, keys=as.character(heat1$genename),columns=c("ENTREZID","GENENAME"), keytype="ENTREZID")$GENENAME
	pdf("GO_heatmap_cluster1.pdf")
	pheatmap(heat1_exp, scale="row", cluster_cols = F,  cluster_rows = F,  show_colnames = F, gaps_col = c(5),cellwidth = 10, cellheight = 10,main="sexual reproduction")
	dev.off()

	heat2 <- read.table('heatmap/ggo_cluster2_exp',col.name=c("genename","X131270_S158","X131272_S159","X131274_S160","X131276_S161","X131278_S162","X131280_S163","X131282_S164","X131284_S165","X131286_S166","X131288_S167"))
	heat2_exp<-heat2[,2:11]
	row.names(heat2_exp) <- AnnotationDbi::select(org.Hs.eg.db, keys=as.character(heat2$genename),columns=c("ENTREZID","GENENAME"), keytype="ENTREZID")$GENENAME
	pdf("GO_heatmap_cluster2.pdf")
	pheatmap(heat2_exp, scale="row", cluster_cols = F,  cluster_rows = F,  show_colnames = F, gaps_col = c(5),cellwidth = 10, cellheight = 10,main="sexual reproduction",fontsize_row=5)
	dev.off()
cluster1

![image](https://github.com/LittleFrogHill/Protist/assets/34407101/fefdd789-9689-42aa-b4b3-a4a61afd3cbf)

cluster2

![image](https://github.com/LittleFrogHill/Protist/assets/34407101/f73e1626-5e60-42c8-af1e-1085baeb4ee6)

supplementary add figures in pdf files
	GO_enrichment_cluster1.pdf
	GO_enrichment_cluster2.pdf
	KEGG_enrichment_cluster1.pdf
	KEGG_enrichment_cluster2.pdf

