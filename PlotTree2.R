
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
    #print(left)
    #print(right)
    #print(length(left))
    noden = 2*nsam - 1
    newick=array()
    newick=parens(noden,tree,left,right,newick)
    return(paste(newick[2:length(newick)],collapse=""))
  }
plot.coalTree = function(tree)
{
    nsam = (length(tree$nodes)+1)/2
  plot(c(0.5,nsam),c(0,2*nsam-1),ylim=c(0,tree$tmrca),pch="",xaxt="n",yaxt="n",xlab="",ylab="",main="An ugly tree!")
  xlabels=array(data=0,dim=nsam)
  ranks=rank(tree$nodes[1:nsam])

  k=1
  for(i in sort(unique(ranks)))
    {
      j = which(ranks == i)
      for(l in 1:length(j))
        {
          xlabels[k]=j[l]
          k=k+1
        }
    }
  for(i in 1:length(xlabels))
    {
      print(paste(i,xlabels[i],
                  tree$nodes[xlabels[i]],
                  tree$times[tree$nodes[xlabels[i]]]))
      if( i %% 2 != 0)
        {
          segments(i,0,i+0.5,tree$times[tree$nodes[xlabels[i]]])
        }
      else
        {
          segments(i,0,i-0.5,tree$times[tree$nodes[xlabels[i]]])
        }
      
    }

}
