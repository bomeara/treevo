## THIS IS A MINOR MODIFICATION OF GEIGER'S tworateTree FUNCTION -- CITE GEIGER IF YOU USE THIS
###TWORATE#####

tworateTreeModified<-function(phy, breakPoint, endRate, beginRate=1) 
{
	times<-branching.times(phy)     
	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=times[which(names(times)==phy$edge[i,1])] #gets tip to node length
		if((age-bl)<breakPoint) #identifies branches that are on the tip side of the break-point (i.e., young)
			phy$edge.length[i]<-(age-min(age, breakPoint))*beginRate+(min(age, breakPoint)-(age-bl))*endRate #If an edge is entirely to the right of the breakpoint it is simply multipled by f.  However, if the edge extends across the breakpoint, we want to leave the part to the left of the BP unchanged (this is the first part of this equation, which multiplies this part by beginRate).  We then want to multiply the part of the edge that is to the right of the BP by f (this is the second part of this line).
				else {
					phy$edge.length[i]<=bl*beginRate
				}
	}
		phy
}

