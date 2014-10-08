#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge -m example.meta example.fasta example.bed example.par example.fasta
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m deletion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m inversion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m insertion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m duplication.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m transposition.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m finsertion.meta example.fna example.bed example.par example.fna
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 64  -m iduplication.meta example.fna example.bed example.par example.fna
#swan_stat.R iduplication.lib1.bam
#swan_stat.R iduplication.lib2.bam
#swan_stat.R duplication.lib1.bam
#swan_stat.R duplication.lib2.bam
#swan_stat.R inversion.lib1.bam
#swan_stat.R inversion.lib2.bam
#swan_stat.R insertion.lib1.bam
#swan_stat.R insertion.lib2.bam
#swan_stat.R finsertion.lib1.bam
#swan_stat.R finsertion.lib2.bam
#swan_stat.R transposition.lib1.bam
#swan_stat.R transposition.lib2.bam
#swan_stat.R deletion.lib1.bam
#swan_stat.R deletion.lib2.bam
plot_vcf.R deletion.bed -o deletion.lib1 -f bed deletion.lib1.bam
plot_vcf.R deletion.bed -o deletion.lib2 -f bed deletion.lib2.bam
