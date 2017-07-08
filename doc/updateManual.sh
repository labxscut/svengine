#!/bin/bash
export SVENGINEPKG=$HOME/setup/svengine                #where the svengine package is located
export SVENGINEPATH=$SVENGINEPKG/svengine           #where the python libs of svengine are located
export SCRIPTPATH=$SVENGINEPKG/scripts              #where the uitlity scripts svengine are located
export PYTHONPATH=$SVENGINEPATH
export PATH=$SCRIPTPATH:$PATH
cat Manual.master.rst | awk '/REPLACE_WITH_SVENGINE_FASFORGE_HELP/{system("python -m mf.fasforge -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_MUTFORGE_HELP/{system("python -m mf.mutforge -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_SVENGINE_TREE2VAR_HELP/{system("python -m mf.tree2var -h | sed -e \"s/^/    /g\"");next}1' >Manual.rst
