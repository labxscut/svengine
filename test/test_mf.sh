#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge -m example.meta example.fasta example.bed example.par example.fasta
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m deletion.meta example.fna example.bed example.par example.fna
samtools index deletion.lib1.bam
samtools index deletion.lib2.bam
plot_vcf.R deletion.bed -o deletion.lib1 -f bed deletion.lib1.bam
plot_vcf.R deletion.bed -o deletion.lib2 -f bed deletion.lib2.bam
#above gonna generate deletion.var, deletion.lib1.bam, deleltion.lib2.bam and everything prefixed with deletion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m inversion.meta example.fna example.bed example.par example.fna
samtools index inversion.lib1.bam
samtools index inversion.lib2.bam
plot_vcf.R inversion.bed -o inversion.lib1 -f bed inversion.lib1.bam
plot_vcf.R inversion.bed -o inversion.lib2 -f bed inversion.lib2.bam
#above gonna generate deletion.var, deletion.lib1.bam, deleltion.lib2.bam and everything prefixed with inversion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m insertion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m duplication.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m transposition.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m finsertion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 16  -m iduplication.meta example.fna example.bed example.par example.fna
swan_stat.R iduplication.lib1.bam
swan_stat.R iduplication.lib2.bam
swan_stat.R duplication.lib1.bam
swan_stat.R duplication.lib2.bam
swan_stat.R inversion.lib1.bam
swan_stat.R inversion.lib2.bam
swan_stat.R insertion.lib1.bam
swan_stat.R insertion.lib2.bam
swan_stat.R finsertion.lib1.bam
swan_stat.R finsertion.lib2.bam
swan_stat.R transposition.lib1.bam
swan_stat.R transposition.lib2.bam
swan_stat.R deletion.lib1.bam
swan_stat.R deletion.lib2.bam
samtools index insertion.lib1.bam
samtools index insertion.lib2.bam
plot_vcf.R insertion.bed -o insertion.lib1 -f bed insertion.lib1.bam
plot_vcf.R insertion.bed -o insertion.lib2 -f bed insertion.lib2.bam
samtools index duplication.lib1.bam
samtools index duplication.lib2.bam
plot_vcf.R duplication.bed -o duplication.lib1 -f bed duplication.lib1.bam
plot_vcf.R duplication.bed -o duplication.lib2 -f bed duplication.lib2.bam
samtools index transposition.lib1.bam
samtools index transposition.lib2.bam
plot_vcf.R transposition.bed -o transposition.lib1 -f bed transposition.lib1.bam
plot_vcf.R transposition.bed -o transposition.lib2 -f bed transposition.lib2.bam
samtools index iduplication.lib1.bam
samtools index iduplication.lib2.bam
plot_vcf.R iduplication.bed -o iduplication.lib1 -f bed iduplication.lib1.bam
plot_vcf.R iduplication.bed -o iduplication.lib2 -f bed iduplication.lib2.bam
samtools index finsertion.lib1.bam
samtools index finsertion.lib2.bam
plot_vcf.R finsertion.bed -o finsertion.lib1 -f bed finsertion.lib1.bam
plot_vcf.R finsertion.bed -o finsertion.lib2 -f bed finsertion.lib2.bam
