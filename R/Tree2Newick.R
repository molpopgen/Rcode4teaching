
parens=function(noden,tree,left,right,newick)
  {

    if( left[noden] == -1 )
      {
        newick[length(newick)+1] = noden
        newick[length(newick)+1] = ":"
        newick[length(newick)+1] = tree$times[tree$nodes[noden]]
      }
    else
      {
        newick[length(newick)+1] = "("
        newick = parens( left[noden], tree,left,right,newick )
        newick[length(newick)+1] = ","
        newick = parens( right[noden], tree,left,right,newick )
        if( tree$nodes[noden] == -1 )
          {
            newick[length(newick)+1] = ");"
          }
        else
          {
            time = tree$times[tree$nodes[noden]] - tree$times[noden]
            newick[length(newick)+1] = "):"
            newick[length(newick)+1] = time
          }
      }
    return(newick)
  }

tree2newick=function(tree)
  #returns coalescent tree in a newick format
  #to print this tree,use the ape package:
  #x=simWFtree(10)
  #x.tree=tree2newick(x)
  #cat(x.tree,file="temp")
  #library(ape)
  #x.tree=read.tree("temp")
  #plot.phylo(x.tree,direction="downwards",y.lim=c(0,x$tmrca),type="cladogram",show.tip.label=FALSE)
  {
    nsam =(length(tree$nodes)+1)/2
    left=array(dim=2*nsam-1,data=-1)
    right=array(dim=2*nsam-1,data=-1)
    print(length(left))
    for(i in 1:(length(left)-1))
      {
        if( left[ tree$nodes[i] ] == -1 )
          {
            left[ tree$nodes[i] ] = i
          }
        else
          {
            right[ tree$nodes[i] ] = i
          }
      }
    noden = 2*nsam - 1
    newick=array()
    newick=parens(noden,tree,left,right,newick)
    return(paste(newick[2:length(newick)],collapse=""))
  }
