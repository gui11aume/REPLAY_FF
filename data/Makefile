all:
	wget ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.chr.gtf.gz
	zcat Mus_musculus.GRCm38.87.chr.gtf.gz | awk '{if ($$3 == "gene") print}' > GRCm38_genes.gtf
	wget http://www.replicationdomain.com/download_data.php?DataSetID=Int93235019&Build=mm10 -O RT_D3_EBM6_All.txt
	wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE93238&format=file&file=GSE93238%5Fgene%2Efpkm%2Etxt%2Egz" -O GSE93238_gene.fpkm.txt.gz

