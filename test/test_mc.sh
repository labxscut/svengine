#testing components
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DEL.meta -o DEL.mut example.fna example.gap.bed simple.par example.fna
samtools sort DEL.mut.lib1.bam DEL.mut
samtools index DEL.mut.bam
swan_plot.R -f bed -o DEL.mut DEL.mut.out.bed DEL.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DUP.meta -o DUP.mut example.fna example.gap.bed simple.par example.fna
samtools sort DUP.mut.lib1.bam DUP.mut
samtools index DUP.mut.bam
swan_plot.R -f bed -o DUP.mut DUP.mut.out.bed DUP.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DINS.meta -o DINS.mut example.fna example.gap.bed simple.par example.fna
samtools sort DINS.mut.lib1.bam DINS.mut
samtools index DINS.mut.bam
swan_plot.R -f bed -o DINS.mut DINS.mut.out.bed DINS.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m FINS.meta -o FINS.mut example.fna example.gap.bed simple.par example.fna
samtools sort FINS.mut.lib1.bam FINS.mut
samtools index FINS.mut.bam
swan_plot.R -f bed -o FINS.mut FINS.mut.out.bed FINS.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m INV.meta -o INV.mut example.fna example.gap.bed simple.par example.fna
samtools sort INV.mut.lib1.bam INV.mut
samtools index INV.mut.bam
swan_plot.R -f bed -o INV.mut INV.mut.out.bed INV.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m IDINS.meta -o IDINS.mut example.fna example.gap.bed simple.par example.fna
samtools sort IDINS.mut.lib1.bam IDINS.mut
samtools index IDINS.mut.bam
swan_plot.R -f bed -o IDINS.mut IDINS.mut.out.bed IDINS.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m TRA.meta -o TRA.mut example.fna example.gap.bed simple.par example.fna
samtools sort TRA.mut.lib1.bam TRA.mut
samtools index TRA.mut.bam
swan_plot.R -f bed -o TRA.mut TRA.mut.out.bed TRA.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m ITRA.meta -o ITRA.mut example.fna example.gap.bed simple.par example.fna
samtools sort ITRA.mut.lib1.bam ITRA.mut
samtools index ITRA.mut.bam
swan_plot.R -f bed -o ITRA.mut ITRA.mut.out.bed ITRA.mut.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m IDUP.meta -o IDUP.mut example.fna example.gap.bed simple.par example.fna
samtools sort IDUP.mut.lib1.bam IDUP.mut
samtools index IDUP.mut.bam
swan_plot.R -f bed -o IDUP.mut IDUP.mut.out.bed IDUP.mut.bam
