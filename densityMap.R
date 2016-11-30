# function plots posterior density of mapped states from stochastic mapping
# written by Liam J. Revell 2012, 2013, 2014, 2015
# minor patch on lines 154 and 177 where I added "subtitle=''" to add.color.bar so that the length of the scale is not printed. I found this confusing for the audience and not necessary for my purposes. [Balthasar Bickel, May 6, 2016]
# additional patch by adding a show.tip.label [Balthasar Bickel, Nov 6, 2016]

densityMap<-function(trees,res=100,fsize=NULL,ftype=NULL,lwd=3,check=FALSE,legend=NULL,
	outline=FALSE,type="phylogram",direction="rightwards",plot=TRUE,...){
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(states)) states<-list(...)$states
	else states<-NULL
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(length(lwd)==1) lwd<-rep(lwd,2)
	else if(length(lwd)>2) lwd<-lwd[1:2]
	tol<-1e-10
	if(!inherits(trees,"multiPhylo")&&inherits(trees,"phylo")) stop("trees not \"multiPhylo\" object; just use plotSimmap.")
	if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"multiPhylo\".")
	h<-sapply(unclass(trees),function(x) max(nodeHeights(x)))
	steps<-0:res/res*max(h)
	trees<-rescaleSimmap(trees,totalDepth=max(h))
	if(check){
		X<-matrix(FALSE,length(trees),length(trees))
		for(i in 1:length(trees)) X[i,]<-sapply(trees,all.equal.phylo,current=trees[[i]])
		if(!all(X)) stop("some of the trees don't match in topology or relative branch lengths")
	}
	tree<-trees[[1]]
	trees<-unclass(trees)
	if(is.null(states)) ss<-sort(unique(c(getStates(tree,"nodes"),getStates(tree,"tips"))))
	else ss<-states
	if(!all(ss==c("0","1"))){
		c1<-paste(sample(c(letters,LETTERS),6),collapse="")
		c2<-paste(sample(c(letters,LETTERS),6),collapse="")
		trees<-lapply(trees,mergeMappedStates,ss[1],c1)
		trees<-lapply(trees,mergeMappedStates,ss[2],c2)
		trees<-lapply(trees,mergeMappedStates,c1,"0")
		trees<-lapply(trees,mergeMappedStates,c2,"1")
	}	
	H<-nodeHeights(tree)
	message("sorry - this might take a while; please be patient")
	for(i in 1:nrow(tree$edge)){
		YY<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
			c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
		ZZ<-rep(0,nrow(YY))
		for(j in 1:length(trees)){
			XX<-matrix(0,length(trees[[j]]$maps[[i]]),2,dimnames=list(names(trees[[j]]$maps[[i]]),
				c("start","end")))
			XX[1,2]<-trees[[j]]$maps[[i]][1]
			if(length(trees[[j]]$maps[[i]])>1){
				for(k in 2:length(trees[[j]]$maps[[i]])){
					XX[k,1]<-XX[k-1,2]
					XX[k,2]<-XX[k,1]+trees[[j]]$maps[[i]][k]
				}
			}
			for(k in 1:nrow(YY)){
				lower<-which(XX[,1]<=YY[k,1]); lower<-lower[length(lower)]
				upper<-which(XX[,2]>=(YY[k,2]-tol))[1]; AA<-0
				names(lower)<-names(upper)<-NULL
				if(!all(XX==0)){
					for(l in lower:upper) 
						AA<-AA+(min(XX[l,2],YY[k,2])-max(XX[l,1],YY[k,1]))/(YY[k,2]-
							YY[k,1])*as.numeric(rownames(XX)[l])
				} else AA<-as.numeric(rownames(XX)[1])
				ZZ[k]<-ZZ[k]+AA/length(trees)
			}
		}
		tree$maps[[i]]<-YY[,2]-YY[,1]
		names(tree$maps[[i]])<-round(ZZ*1000)
	}
	cols<-rainbow(1001,start=0.7,end=0); names(cols)<-0:1000
	tree$mapped.edge<-makeMappedEdge(tree$edge,tree$maps)
	tree$mapped.edge<-tree$mapped.edge[,order(as.numeric(colnames(tree$mapped.edge)))]
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	x<-list(tree=tree,cols=cols,states=ss); class(x)<-"densityMap"
	if(plot) plot.densityMap(x,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,
		type=type,mar=mar,direction=direction,offset=offset,hold=hold)
	invisible(x)
}

## S3 plot method for objects of class "densityMap"
## also used internally by plot.contMap
## written by Liam J. Revell 2012, 2013, 2014, 2015, 2016

plot.densityMap<-function(x,show.tip.label=T, lend=0, ...){
	if(class(x)=="densityMap"){
		tree<-x$tree
		cols<-x$cols
	} else stop("x should be an object of class \"densityMap\"")
	H<-nodeHeights(tree)
	# get & set optional arguments
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-FALSE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-3
	if(length(lwd)==1) lwd<-rep(lwd,2)
	else if(length(lwd)>2) lwd<-lwd[1:2]
	if(hasArg(leg.txt)) leg.txt<-list(...)$leg.txt
	else leg.txt<-c("0",paste("PP(state=",x$states[2],")",sep=""),"1")
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(is.null(legend)) legend<-0.5*max(H)
	if(is.null(fsize)) fsize<-c(1,1)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(is.null(ftype)) ftype<-c("i","reg")
	if(length(ftype)==1) ftype<-c(ftype,"reg")
	if(!show.tip.label) ftype[1] <- 'off'
	# done optional arguments
	if(legend){
		if(legend>max(H)){ 
			message("legend scale cannot be longer than total tree length; resetting")
			legend<-0.5*max(H)
		}
	}
	if(hold) null<-dev.hold()
	if(type=="phylogram"){
		N<-length(tree$tip.label)
		if(legend&&is.null(ylim)) ylim<-c(1-0.12*(N-1),N)
		else if(is.null(ylim)) ylim<-NULL
		if(outline){
			par(col="transparent")
			plotTree(tree,fsize=fsize[1],lwd=lwd[1]+2,
				offset=offset+0.2*lwd[1]/3+0.2/3,ftype=ftype[1],xlim=xlim,
				ylim=ylim,mar=mar,direction=direction,hold=FALSE)
			par(col="black")
		}
		plotSimmap(tree,cols,pts=FALSE,lwd=lwd[1],fsize=fsize[1],mar=mar,ftype=ftype[1],add=outline,
			xlim=xlim,ylim=ylim,direction=direction,offset=offset,hold=FALSE, lend=lend)
		if(legend){
			ff<-function(dd){
				if(!("."%in%dd)) dig<-0
				else dig<-length(dd)-which(dd==".")
				dig
			}
			dig<-max(sapply(strsplit(leg.txt[c(1,3)],split=""),ff))
			add.color.bar(legend,cols,title=leg.txt[2],lims<-as.numeric(leg.txt[c(1,3)]),
				digits=dig,prompt=FALSE,x=if(direction=="leftwards") max(H)-legend else 0,
				y=1-0.08*(N-1),lwd=lwd[2], subtitle="",
				fsize=fsize[2],
				direction=if(!is.null(xlim)) if(xlim[2]<xlim[1]) "leftwards" else "rightwards" else "rightwards")
		}
	} else if(type=="fan"){
		if(outline){
			par(col="white")
			invisible(capture.output(plotTree(tree,type="fan",lwd=lwd[1]+2,
				mar=mar,fsize=fsize[1],
				ftype=ftype[1],xlim=xlim,ylim=ylim,hold=FALSE)))
			par(col="black")
		}
		invisible(capture.output(plotSimmap(tree,cols,lwd=lwd[1],
			mar=mar,fsize=fsize[1],add=outline,ftype=ftype[1],
			type="fan",xlim=xlim,ylim=ylim,hold=FALSE)))
		if(legend){
			ff<-function(dd){
				if(!("."%in%dd)) dig<-0
				else dig<-length(dd)-which(dd==".")
				dig
			}
			dig<-max(sapply(strsplit(leg.txt[c(1,3)],split=""),ff))
			add.color.bar(legend,cols,title=leg.txt[2],lims<-as.numeric(leg.txt[c(1,3)]),digits=dig,
				prompt=FALSE,x=0.9*par()$usr[1],y=0.9*par()$usr[3],lwd=lwd[2], subtitle="",
				fsize=fsize[2])
		}
	}
	if(hold) null<-dev.flush()
}

## S3 print method for object of class "densityMap"
## written by Liam J. Revell 2013

print.densityMap<-function(x,...){
	cat("Object of class \"densityMap\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with ",length(x$tree$tip.label)," tips and ",x$tree$Nnode," internal nodes.\n\n",sep=""))
	cat(paste("(2) The mapped posterior density of a discrete binary character with states (",x$states[1],", ",x$states[2],").\n\n",sep="")) 
}

## set new color map for object of class 'densityMap' or 'contMap'
## written by Liam J. Revell 2014

setMap<-function(x,...){
	if(hasArg(invert)) invert<-list(...)$invert
	else invert<-FALSE
	n<-length(x$cols)
	if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
	else x$cols[1:n]<-colorRampPalette(...)(n)
	x
}

## drop tips from an object of class 'densityMap'
## written by Liam J. Revell 2014, 2015

drop.tip.densityMap<-function(x,tip){
	if(inherits(x,"densityMap")){ 
		class(x)<-"contMap"
		x<-drop.tip.contMap(x,tip)
		class(x)<-"densityMap"
		return(x)
	} else cat("x should be an object of class \"densityMap\"\n")
}


makeMappedEdge<-function(edge,maps){
	st<-sort(unique(unlist(sapply(maps,function(x) names(x)))))
	mapped.edge<-matrix(0,nrow(edge),length(st))
	rownames(mapped.edge)<-apply(edge,1,function(x) paste(x,collapse=","))
	colnames(mapped.edge)<-st
	for(i in 1:length(maps)) 
		for(j in 1:length(maps[[i]])) 
			mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
	return(mapped.edge)
}


matchNodes<-function(tr1,tr2,method=c("descendants","distances"),...){
	require(phangorn)
	if(!inherits(tr1,"phylo")||!inherits(tr1,"phylo")) stop("tr1 & tr2 should both be objects of class \"phylo\".")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	method<-method[1]
	method<-matchType(method,c("descendants","distances"))
	if(method=="descendants"){
		# desc.tr1<-lapply(1:tr1$Nnode+length(tr1$tip),function(x) extract.clade(tr1,x)$tip.label)
		desc.tr1<-lapply(1:tr1$Nnode+length(tr1$tip),function(x) tr1$tip.label[unlist(Descendants(tr1, x, type="tips"))])
		names(desc.tr1)<-1:tr1$Nnode+length(tr1$tip)
		# desc.tr2<-lapply(1:tr2$Nnode+length(tr2$tip),function(x) extract.clade(tr2,x)$tip.label)
		desc.tr2<-lapply(1:tr2$Nnode+length(tr2$tip),function(x) tr2$tip.label[unlist(Descendants(tr2, x, type="tips"))])
		names(desc.tr2)<-1:tr2$Nnode+length(tr2$tip)
		Nodes<-matrix(NA,length(desc.tr1),2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:length(desc.tr1)){
			Nodes[i,1]<-as.numeric(names(desc.tr1)[i])
			for(j in 1:length(desc.tr2))
				if(all(desc.tr1[[i]]%in%desc.tr2[[j]])&&all(desc.tr2[[j]]%in%desc.tr1[[i]]))
					Nodes[i,2]<-as.numeric(names(desc.tr2)[j])
		}
	} else if(method=="distances"){
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-1e-6
		if(hasArg(corr)) corr<-list(...)$corr
		else corr<-FALSE
		if(corr) tr1$edge.length<-tr1$edge.length/max(nodeHeights(tr1))
		if(corr) tr2$edge.length<-tr2$edge.length/max(nodeHeights(tr2))
		D1<-dist.nodes(tr1)[1:length(tr1$tip),1:tr1$Nnode+length(tr1$tip)]
		D2<-dist.nodes(tr2)[1:length(tr2$tip),1:tr2$Nnode+length(tr2$tip)]
		rownames(D1)<-tr1$tip.label
		rownames(D2)<-tr2$tip.label
		common.tips<-intersect(tr1$tip.label,tr2$tip.label)
		D1<-D1[common.tips,]
		D2<-D2[common.tips,]
		Nodes<-matrix(NA,tr1$Nnode,2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:tr1$Nnode){
			if(corr) z<-apply(D2,2,function(X,y) cor(X,y),y=D1[,i])
			else z<-apply(D2,2,function(X,y) 1-sum(abs(X-y)),y=D1[,i])
			Nodes[i,1]<-as.numeric(colnames(D1)[i])
			if(any(z>=(1-tol))){
				a<-as.numeric(names(which(z>=(1-tol))))
				if(length(a)==1) Nodes[i,2]<-a
				else {
					Nodes[i,2]<-a[1]
					if(!quiet) warning("polytomy detected; some node matches may be arbitrary")
				}
			}
		}
	}
	return(Nodes)
}

matchType<-function(type,types){
	for(i in 1:length(types))
		if(all(strsplit(type,split="")[[1]]==strsplit(types[i],split="")[[1]][1:length(strsplit(type,split="")[[1]])]))
			type=types[i]
	return(type)
}
