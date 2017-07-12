#!/bin/bash
export SVENGINEPKG=$HOME/setup/svengine                #where the svengine package is located
export SVENGINEPATH=$SVENGINEPKG/svengine           #where the python libs of svengine are located
export SCRIPTPATH=$SVENGINEPKG/scripts              #where the uitlity scripts svengine are located
export PYTHONPATH=$SVENGINEPATH
export PATH=$SCRIPTPATH:$PATH
cat Manual.master.rst | awk '/REPLACE_WITH_SVENGINE_EXAMPLE_TREEVAR/{system("cat ../test/egTree.var | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_EXAMPLE_TREENEWICK/{system("cat ../test/egTree.newick | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_EXAMPLE_PARFILE/{system("cat ../test/simple.par | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_EXAMPLE_METAFILE/{system("cat ../test/example.meta | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_EXAMPLE_VARFILE/{system("cat ../test/example.var | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_MUTFORGE_HELP/{system("python -m mf.mutforge -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_TREE2VAR_HELP/{system("python -m mf.tree2var -h | sed -e \"s/^/    /g\"");next}1' >Manual.rst
