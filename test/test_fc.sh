#testing components
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DEL.meta -o DEL.fas example.gap.bed simple.par example.fna
samtools sort DEL.fas.lib1.bam DEL.fas
samtools index DEL.fas.bam
swan_plot.R -f bed -o DEL.fas DEL.fas.out.bed DEL.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DUP.meta -o DUP.fas example.gap.bed simple.par example.fna
samtools sort DUP.fas.lib1.bam DUP.fas
samtools index DUP.fas.bam
swan_plot.R -f bed -o DUP.fas DUP.fas.out.bed DUP.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DINS.meta -o DINS.fas example.gap.bed simple.par example.fna
samtools sort DINS.fas.lib1.bam DINS.fas
samtools index DINS.fas.bam
swan_plot.R -f bed -o DINS.fas DINS.fas.out.bed DINS.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m FINS.meta -o FINS.fas example.gap.bed simple.par example.fna
samtools sort FINS.fas.lib1.bam FINS.fas
samtools index FINS.fas.bam
swan_plot.R -f bed -o FINS.fas FINS.fas.out.bed FINS.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m INV.meta -o INV.fas example.gap.bed simple.par example.fna
samtools sort INV.fas.lib1.bam INV.fas
samtools index INV.fas.bam
swan_plot.R -f bed -o INV.fas INV.fas.out.bed INV.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m IDINS.meta -o IDINS.fas example.gap.bed simple.par example.fna
samtools sort IDINS.fas.lib1.bam IDINS.fas
samtools index IDINS.fas.bam
swan_plot.R -f bed -o IDINS.fas IDINS.fas.out.bed IDINS.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m TRA.meta -o TRA.fas example.gap.bed simple.par example.fna
samtools sort TRA.fas.lib1.bam TRA.fas
samtools index TRA.fas.bam
swan_plot.R -f bed -o TRA.fas TRA.fas.out.bed TRA.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m ITRA.meta -o ITRA.fas example.gap.bed simple.par example.fna
samtools sort ITRA.fas.lib1.bam ITRA.fas
samtools index ITRA.fas.bam
swan_plot.R -f bed -o ITRA.fas ITRA.fas.out.bed ITRA.fas.bam
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m IDUP.meta -o IDUP.fas example.gap.bed simple.par example.fna
samtools sort IDUP.fas.lib1.bam IDUP.fas
samtools index IDUP.fas.bam
swan_plot.R -f bed -o IDUP.fas IDUP.fas.out.bed IDUP.fas.bam
