#creates a plot good for figures in talks, etc.
prettyTree<-function(ntax=10,edge.width=10)
{
    plot(ladderize(chronogram(rtree(ntax)),right=FALSE),show.tip.label=FALSE,direction="upwards",edge.width=edge.width)
}
