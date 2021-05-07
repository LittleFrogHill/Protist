# Protist
Kenny Data

## 1.Trinity assembly
  LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH
  Trinity --seqType fq --max_memory 60G --CPU 20  \
  	--left ./A006200108_131270_S158_L003_R1_001_val_1.fq.gz,./A006200108_131272_S159_L003_R1_001_val_1.fq.gz,./A006200108_131274_S160_L003_R1_001_val_1.fq.gz,./A006200108_131276_S161_L003_R1_001_val_1.fq.gz,./A006200108_131278_S162_L003_R1_001_val_1.fq.gz,./A006200108_131280_S163_L003_R1_001_val_1.fq.gz,./A006200108_131282_S164_L003_R1_001_val_1.fq.gz,./A006200108_131284_S165_L003_R1_001_val_1.fq.gz,./A006200108_131286_S166_L003_R1_001_val_1.fq.gz,./A006200108_131288_S167_L003_R1_001_val_1.fq.gz \
	--right ./A006200108_131270_S158_L003_R2_001_val_2.fq.gz,./A006200108_131272_S159_L003_R2_001_val_2.fq.gz,./A006200108_131274_S160_L003_R2_001_val_2.fq.gz,./A006200108_131276_S161_L003_R2_001_val_2.fq.gz,./A006200108_131278_S162_L003_R2_001_val_2.fq.gz,./A006200108_131280_S163_L003_R2_001_val_2.fq.gz,./A006200108_131282_S164_L003_R2_001_val_2.fq.gz,./A006200108_131284_S165_L003_R2_001_val_2.fq.gz,./A006200108_131286_S166_L003_R2_001_val_2.fq.gz,./A006200108_131288_S167_L003_R2_001_val_2.fq.gz \
  --output ./trinity.out






## 1.Annotation with docker KOBAS

  docker run -it -v /home/shangao:/gpfs -v /home/shangao/Data/kobas_DB/seq_pep:/opt/kobas-3.0/seq_pep -v /home/shangao/Data/kobas_DB/sqlite3:/opt/kobas-3.0/sqlite3 -v /NVME2/Scratch/gaoshan/:/NVME2/Scratch/gaoshan/ kobas


  annotate.py -i ./DEGs.adj.gene.filter500.fa.transdecoder.pep -s hsa -t fasta:pro -o anno_seq.tsv -e 1e-5 -r 1 -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep

  identify.py -f anno_seq.tsv -o identify.out -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep
