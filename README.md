# Protist
Kenny Data

## 1.Annotation with docker KOBAS

  docker run -it -v /home/shangao:/gpfs -v /home/shangao/Data/kobas_DB/seq_pep:/opt/kobas-3.0/seq_pep -v /home/shangao/Data/kobas_DB/sqlite3:/opt/kobas-3.0/sqlite3 -v /NVME2/Scratch/gaoshan/:/NVME2/Scratch/gaoshan/ kobas


  annotate.py -i ./DEGs.adj.gene.filter500.fa.transdecoder.pep -s hsa -t fasta:pro -o anno_seq.tsv -e 1e-5 -r 1 -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep

  identify.py -f anno_seq.tsv -o identify.out -q /opt/kobas-3.0/sqlite3 -y /opt/kobas-3.0/seq_pep
