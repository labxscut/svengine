.. contents:: Navigator

Manual
========

0. $SVENGINEPKG environmental variable
------------------------------------------------------------------------

$SVENGINEPKG should be set to the path where you unzip or cloned the svengine package.
Below is where you can locate the original script or example file within $SVENGINEPKG.

1. Short manuals for Python and R scripts are also available with '-h' option
-------------------------------------------------------------------------

$SVENGINEPKG/svengine/mf/mutforge.py. 

The command implements a parallelized algorithm that divides template genome into trunks of contigs, spikes structural variants into contigs, samples short reads from the altered contigs, and finally merges back the short-read sets into one and performs the alignment. 

.. code:: bash

  REPLACE_WITH_SVENGINE_MUTFORGE_HELP

$SVENGINEPKG/svengine/mf/tree2var.py. 

The tree2var command implements a procedure that determines variant fractions from an input phylogeny tree based on Equations (1) and (2) and a depth first search graph algorithm, and then substitutes these allele fractions in an input VAR file. (see SVEngine publication).

.. code:: bash

  REPLACE_WITH_SVENGINE_TREE2VAR_HELP

2. Short manuals for text based input are also found within the commented lines of the file
-------------------------------------------------------------------------------------------

$SVENGINEPKG/svengine/test/example.var.
The VAR format is intended for specifying exact information for individual variants, which includes variant id, parent id (if it is part of a complex event, e.g., the deletion part within a translocation), fraction, ploidy, chromosome, starting position, and the sequence length to be deleted and/or the sequence content to be inserted.

.. code:: bash

  REPLACE_WITH_SVENGINE_EXAMPLE_VARFILE

$SVENGINEPKG/svengine/test/example.meta.
The META format is intended for an easy higher-level control, allowing one to specify a desired meta distribution of variants, including variant type, and total number of events, size, allele fraction, and ploidy distributions per type. It is also possible to specify where and how to sample insertion sequence in the case of foreign insertion. For example, a user can easily request 100 deletions of size ranging from 100bp to 10kbp, of uniform allelic fraction distribution and fair coin Bernoulli distribution of homo- and heterozygosity by one line of test in the META file.  These parameters can be stipulated in a single line.

.. code:: bash

  REPLACE_WITH_SVENGINE_EXAMPLE_METAFILE

$SVENGINEPKG/svengine/test/example.par.

.. code:: bash

  REPLACE_WITH_SVENGINE_EXAMPLE_PARFILE

$SVENGINEPKG/svengine/test/egTree.newick.

The tree input to SVEngine is in standard NEWICK format  – a widely accepted format using parenthesis to encode nested tree structures. Each internal node is labelled by the population splitting variant and weighted by the conditional splitting fraction. Each leaf node is labelled by associated terminal genotype and weighted by the subpopulation fraction (optional). See http://evolution.genetics.washington.edu/phylip/newicktree.html. For instance, the NEWICK string for the example binary tree in SVEngine's publication is: 

.. code:: bash

  REPLACE_WITH_SVENGINE_EXAMPLE_TREENEWICK

$SVENGINEPKG/svengine/test/egTree.newick.

.. code:: bash

  REPLACE_WITH_SVENGINE_EXAMPLE_TREEVAR
