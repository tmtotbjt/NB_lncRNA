## radimas genų, kurių promotorius/reguliacinis regionas
#turi CCTGCG motyvą (prie jo jungiasi CD44-ICD)
#ir sumapinimas prie Symbol ir ensID

# .fai pasidarymas
samtools faidx GRCh38.primary_assembly.genome.fa

# GTF → genePred → BED
gtfToGenePred gencode.v49.primary_assembly.annotation.gtf genes.genePred
genePredToBed genes.genePred genes.bed

# promotoriniu regionu pasidarymas
bedtools flank -i genes.bed -g GRCh38.primary_assembly.genome.fa.fai -l 1000 -r 100 -s > promoters.bed

# istraukimas promotoriniu seku
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name > promoters.fa

# kur randamas CCTGCG
grep -B 1 -E "CCTGCG" promoters.fa > motif_hits.fa

# transkriptu ID ištraukimas iš FASTA header
grep ">" motif_hits.fa | cut -d':' -f1 | sed 's/>//' | sort | uniq > candidate_transcripts.txt

# randamas genas transkriptui
awk '$3=="transcript" { match($0, /gene_id "([^"]+)";.*transcript_id "([^"]+)"/, arr); print arr[2] "\t" arr[1] }' \
    gencode.v49.primary_assembly.annotation.gtf | sort | uniq > transcript_to_gene.txt

# genai su motyvu promotoriuje (SYMBOL)
grep "gene_id" gencode.v49.primary_assembly.annotation.gtf | \
grep -F -f candidate_transcripts.txt | \
sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort | uniq > candidate_genes.txt

# genai su motybu promotoriuje (ENS ID)
grep -F -f candidate_transcripts.txt transcript_to_gene.txt | cut -f2 | sort | uniq > candidate_genes_ensembl.txt
