# DONE:
#    this should allow user to specify a tree strucutre
# TODO:
#    generate a tree structure based on following parameters:
#    number of variants, cell populations will be that plus 1
#    once generated use topological sort to obtain an ordered cell groups, normal always first
#    conditional fractions of gaining each variants, a sequence of uniform
#    https://pythonhosted.org/DendroPy/primer/treesims.html

import pandas
import numpy
import dendropy
import argparse
import time
import os
import json

def main(args):
  varfile = args.varfile
  nwkfile = args.nwkfile
  plottree = args.plottree
  oname = varfile.name.rsplit(".")[0]
  oprefix = oname if not args.oprefix else args.oprefix

  #os.chdir("/home/lixia/setup/svengine/test")
  #egTree=dendropy.Tree.get(path="egTree.newick",schema="newick",suppress_internal_node_taxa=False)
  egTree=dendropy.Tree.get(path=nwkfile.name, schema="newick", suppress_internal_node_taxa=False)

  #for node in egTree.preorder_node_iter(): print node
  #for node in egTree.preorder_node_iter(): print node.edge.length

  if plottree: #plot tree in both stdout and a plot file 
    egTree.print_plot(show_internal_node_labels=True, plot_metric="depth")
    treePlotFileName = oprefix+".tree.plot"
    if os.path.isfile(treePlotFileName):
      quit("output treePlotFileName:" + treePlotFileName + "already exists, remove first!\n")
      print >>open(treePlotFileName, "w"), 
              egTree.as_ascii_plot(show_internal_node_labels=True, plot_metric="depth")

  #for node in egTree.preorder_node_iter():   #this gives topologically ordered 
    #if node.edge.length == None: print node;

  leafNodes = [ node for node in egTree.preorder_node_iter() if node.child_nodes() == [] ] #find leafs
  internalNodes = [ node for node in egTree.preorder_node_iter() if node.child_nodes() != [] ] #find internals
  internalLabels = [ node.taxon.label for node in egTree.preorder_node_iter() if node.child_nodes() != [] ] #internal labels
  internalFracs = [ node.edge.length for node in egTree.preorder_node_iter() if node.child_nodes() != [] ]  #conditional fractions
  #rootNode = egTree.find_node_with_taxon_label('V1')

  leafPath = dict()
  for leaf in leafNodes:
    #print leaf.taxon.label,
    pNode = leaf.parent_node
    cNode = leaf
    leafPath[ leaf.taxon.label ] = dict()
    while pNode != None:
      #print pNode.taxon.label,
      leafPath[ leaf.taxon.label ][ pNode.taxon.label ] = (pNode.child_nodes()[0] != cNode) + 1  
      #if first child, i..e non-carrier return 1, if second child, i.e. carrier,  return 2
      cNode = pNode 
      pNode = cNode.parent_node

  fV = pandas.DataFrame.from_dict(dict(zip(internalLabels, internalFracs)), orient='index')
  CVdf = pandas.DataFrame.from_dict(leafPath)
  SVdf = (((CVdf==2)*(numpy.multiply(CVdf-1,fV))+(CVdf==1)*numpy.multiply(2-CVdf,1-fV)))
  SVdf1 = SVdf.fillna(1)
  FC = SVdf1.product(axis=0)
  FV = numpy.sum((CVdf==2)*FC, axis=1)

  #so far fractions of variant have been computed
  #now to replace each line of var file by correct fractions

  varColNames = ["VID","MID","VARFREQ","VARHAP","VARCHR","VARPOS","VARDEL","VARDELSIZE","VARINS","VARINSSEQ"]
  varTable = pandas.read_table(varfile, comment='#', names=varColNames, sep='\t')
  varTable['VARFREQ'] = [ FV[v] for v in varTable['MID'] ]

  #output varTable
  treeVarFileName = oprefix+".tree.var"
  if os.path.isfile(treeVarFileName):
    quit("output treeVarFile:" + treeVarFileName + "already exists, remove first!\n")
  treeVarFile = open(treeVarFileName, "a")
  treeVarFile.write("#"+"\t".join(varColNames)+"\n")
  varTable.to_csv(treeVarFile, sep="\t", header=False, index=False, float_format="%.3f")
  treeVarFile.close()

  treeJsonFileName = oprefix+".tree.json"
  with open(treeJsonFileName, 'w') as outfile:
    json.dump({ 'FC':FC.to_dict(), 'FV':FV.to_dict(), 'fV':fV.to_dict() }, outfile)
    outfile.close()
  
  return #end of run

#use cases:
#  1. user provide a VAR file and a NEWICK file where MID and NEWICK file's internal ID matches

def run():
  parser = argparse.ArgumentParser(description="assign tree structures to VAR files")
  parser.add_argument('varfile', metavar='varfile', type=argparse.FileType('rU'), #example.bed
                        help="inputfile of SVEngine's VAR format, required")
  parser.add_argument('nwkfile', metavar='nwkfile', type=argparse.FileType('rU'), #example.par
                        help="inputfile of Newick format where the internal nodes are named and matches to MID used in varfile, required")
  parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #
                        help="prefix for output files, default is the same as varfile's prefix yet change suffix to .vtr")
  parser.add_argument('-p', '--plottree', dest='plottree', action='store_true', default=True,
                        help="generate a text plot of a evolution phylogenetic tree")

  args = parser.parse_args()
  main(args)

if __name__ == '__main__':
  starttime = time.time()
  run()
  print "total runtime: ", time.time()-starttime, " seconds"
