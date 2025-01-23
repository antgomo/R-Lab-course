setwd("/media/IMID/Projects/FIS_Meth_analysis/Methylation_Paper/WB/Results/Quant")
library(igraph)
#Where spinmods is a list, with each element of the list a vector of Entrez IDs corresponding to the nodes
#of a community, pin is an iGraph graph object with nodes corrresponding to genes labelled by Entrez IDs,
#stat is a vector of node weights to be randomly permuted over the network and nrand is the number of
#random permutations to compute. In order to validate the T1D communities in the SKIN data set, it can be
#called as follows:




#However, the validation of a FEM module that you infer in some other dataset should not be too difficult.

#The way we had implemented this, if I recall correctly, was by comparing the weighted average statistic of DM for the module to the null 
#obtained by selecting statistics at random from the study and allocating them randomly to the nodes of the module. 
#The weights in the average do not change and would be +1 or -1 depending on the directionality of the observed statistic. 
#You would run say 1000 Monte-Carlo randomizations to get the null, from which you can then assess significance. Hope this helps. 

##############################3
##FUNCTIONS
#grap weighted according its methylation stats

statLabel = function(g, stat)
  # Label graph g with statistics stat (averaged on edges)
{
  A = get.adjacency(g); 
  A = A[names(stat),names(stat)]; # just to make sure
  astat = abs(stat);
  TMAX = max(astat);
  temp1 = apply(A, 1, function(v) return(v*astat))
  W = (temp1 + t(temp1))/(2*TMAX);
  print("Generating weighted graph");
  G = graph.adjacency(W, mode = "undirected", weighted=TRUE);
  V(G)$weight = astat;
  return(G);
}

Modularity = function(v, g)
  # Returns modularity score of subgraph induced by vertices v in iGraph g
{ 
  h = induced.subgraph(g, v);
  return(sum( E(h)$weight ))
}

nodeValidate = function(lv, g, tstats, nrandomizations = 1000)
  # Test elements of lv in g over node-randomizations of tstats 
{
  #A = get.adjacency(g);
  nrm.lv = list(); 
  nrp.v = vector();
  obs.v = vector();
  output = list();
  TMAX = max(tstats);
  
  # First restrict to g:
  mods = lapply(lv, function(v) {
    return( intersect(v,V(g)$name) )
  });
  print(sapply(mods, length));
  
  # Observed modularity values
  obs.v = as.vector(lapply(mods, function (j) { return(Modularity(j, g)) }));
  
  # Node-randomizations of modularity
  nmods = length(mods);
  nrm.lv = lapply(1:nmods, function (i) {
    j = mods[[i]];
    v = vector();
    h = induced.subgraph(g,j); # new iGraph uses induced.subgraph() not subgraph()
    B = get.adjacency(h, sparse=FALSE); # Do not want sparse object
    for (k in 1:nrandomizations) {
      print(paste("Testing significance: module",i,"of",nmods,"Randomization",k,"of",nrandomizations))
      atperm = sample(tstats, nrow(B) , replace=FALSE)
      temp1 = apply(B, 1, function(v) return(v*atperm))
      W = (temp1 + t(temp1))/(2*TMAX);
      v[k] = sum(W)/2; # W is weighted adj matrix (with diag=0) so every edge counted twice
      #v[k]<-sum(W)
    }
    return(v);
  }); names(nrm.lv) = names(mods);
  
  # Empirical p-values
  for (j in 1:length(mods)) {
    nrp.v[j] = length( which(nrm.lv[[j]] > obs.v[j]) )/nrandomizations;
  }; names(nrp.v) = names(mods);
  
  output[[1]] = nrp.v; output[[2]] = obs.v; output[[3]] = nrm.lv; 
  names(output) = c("fdr","Observed","Random");
  return(output);
}

# Code to assign colors to vectors
vect2color = function(v, palettev, breaks) {
  w = v; 
  for (i in 1:length(palettev)) { w[which( v >= breaks[i] & v < breaks[i+1] )] = palettev[i]; }
  return(w);
}

renderModule = function(eid, g, pval, stat) 
  # Renders module with ENTREZ IDs eid, on top of network g with pvalues pval and statistics stat:
{
  # Vertex palette:
  vertexPalette.v = maPalette(low="yellow", high= "blue", mid="white", k=49);
  vertexBreaks.v = seq(from = log10(0.05), to = -log10(0.05), 
                       by = -2*log10(0.05)/(-2+length(vertexPalette.v )))
  vertexBreaks.v = c(-(1+max(-log10(pval))), vertexBreaks.v, 1+max(-log10(pval)));
  length(vertexPalette.v); length(vertexBreaks.v);
  # Edge palette: grey to red via pink, red at top 10% of edges
  edgePalette.v = maPalette(low="lightgrey",high="red",mid="pink",k = 50);
  netPercentiles.v = seq(from=0,to=1,by=1/length(edgePalette.v))
  edgeBreaks.v = sapply(netPercentiles.v, function(x) quantile(E(g)$weight, x))
  edgeBreaks.v[1] = 0; edgeBreaks.v[length(edgeBreaks.v)] = 1;
  
  # Compute iGraph object
  h = induced.subgraph(g, eid);
  stat.v = stat[V(h)$name];
  pval.v = pval[V(h)$name];
  slpval.v = sign(stat.v)*-log10(pval.v); # signed, logged p-values
  
  par(mar=c(4,0,2,0))
  # Color edges between grey and red according to significance  
  E(h)$color = vect2color(E(h)$weight, edgePalette.v, edgeBreaks.v);
  # Color nodes blue to yellow according to hyper/hypo-methylation
  V(h)$color = vect2color(slpval.v, vertexPalette.v, vertexBreaks.v);
  
  if ( length(V(h)) <= 20 ) {  vl = unlist(entrez2symbol[V(h)$name])  } else { vl = "" }
  
  plot(h,layout = layout.fruchterman.reingold,  # same layout each time
       vertex.label = vl, 
       vertex.frame.color = "black", 
       vertex.label.dist = 1.1, 
       vertex.label.font = 3, vertex.label.color = "black", 
       #vertex.size = 15*13/length(V(h)),
       #edge.width = 160/length(V(h))
       vertex.size = if (length(V(h)) < 50) { 15 } else {15*13/length(V(h))},
       edge.width  = if (length(V(h)) < 50) {  6 } else { 160/length(V(h) )}
  );
}


## t test
myfunction <- function(x,y) {
  test = wilcox.test(x, y, paired=T)
  return(test$p.value)
}


####################################################################ANALYSIS
##load set to test P6
  
##original Data

##type of analysis

analysis<-"DAS28"


###

wk0.dir<-"wk0_FIS"

p6.dir<-"wk12_FIS"


          ##load wk0 EPimod
          
          load(paste0(wk0.dir,"/","EpiMod_wk0_DAS28.RData"))
          Spain.Epi<-EpiMod.o
          rm(EpiMod.o)
          gc()
          ##load P6 EPimod
          load(paste0(p6.dir,"/","Epimod_DAS28_FIS_wk12.RData"))
          
          ##now get module names to test
          
          
          lv<-lapply(Spain.Epi$topmod, function(x) {
            
            
            return(as.character(x[1]$EntrezID))
          } 
          
          )
          
          ####P6's graph to allocate randomly the stats
          
          g<-graph_from_adjacency_matrix(EpiMod.o$adj)
          
          ##methylation stats to allocate them randomly 
          load(paste0(p6.dir,"/","Stats_Epimod_DAS28_FIS_wk12.RData"))##antiTNF or QUANT depending analysis
          
          tstats<-intEpi.o$statM$t
          names(tstats)<-rownames(intEpi.o$statM)
          
          
          ##build graph
          
          g<-statLabel(g, tstats)
          
          # Node-randomizations of modularity
          
          Stats<-nodeValidate(lv,g, tstats,10000)##modules with vertex in EntrezID (FIS), graph to be test and meth stats (P6),number of randomizations
          
          stats.perm<-as.data.frame(Stats$fdr)
          ##now get correct direction
          
          ##for each component of the module, extract meth val and locate in its module
          
          stats.sp<-list()
          ##get the genes in P6 set
          stats.sp<-lapply(lv, function(x) {
            
            return(intEpi.o$statM[x,"t"])
            return(rownames(intEpi.o$statM[x]))} )
          
          #re-scale
          
          
          
          library(scales)
          stats.sp2<-lapply(stats.sp, function(x) {return(rescale(x, to = c(-2, 2)))} )
          
          #get original stats 
          
          
          ori.stats<-lapply(Spain.Epi$topmod, function(x)return(x$'stat(DNAm)'))
          
          ##re-scale
          ori.stats2<-lapply(ori.stats, function(x) return(rescale(x,to=c(-2,2))))
####dir to print results

	setwd("/media/IMID/Projects/FIS_Meth_analysis/Methylation_Paper/WB/Results/Quant/Replications_wk12")          
          ##Toni's option
          
          ##Sobre replicació moduls: de cada modul fer un vector dels estadístics dels gens (crec que era un tvalue: te significació i sentit)
          ##i comparar amb un paired t-test el vector del dataset 1 vs el 2. 
          #No hauria de donar significatiu per confirmar que estem veient el mateix tipus de mòdul.
          
          results.wilcox<-mapply(myfunction, x = ori.stats2, y = stats.sp2, SIMPLIFY = FALSE)
          
          results.wilcox<-as.data.frame(t(as.data.frame(results.wilcox)))
          
          
          Stats.r<-merge(stats.perm, results.wilcox, by="row.names")
          colnames(Stats.r)<-c("Module","pval_Permut","pval_Direction")
          
          write.csv(Stats.r,paste0("Validation_Modules_",analysis,"_",wk0.dir,"_",p6.dir,".csv"),row.names = F)
          
          #######################Perm graph
          
          pdf(paste0("Perm_Plot_validation_",analysis,"_",wk0.dir,"_",p6.dir,".pdf"))
          
          for( i in 1:length(Stats$Random)){
            
            xlims.0<-min(Stats$Random[[i]])
            xlims.1<-Stats$Observed[[i]]+5
            
            h<-hist(Stats$Random[[i]], col="lightblue",
                    xlab = "FDR estimate",
                    xlim = c(xlims.0,xlims.1),
                    main = paste0("Validation_ ",analysis,"_",wk0.dir,"_",p6.dir,"_",names(Stats$fdr)[i], sep="")
            )
            
            abline(v = Stats$Observed[[i]], col = "red", lwd = 4)
            abline(v = mean(Stats$Random[[i]]), col = "black", lwd = 4)
            
            
            xfit<-seq(min(Stats$Random[[i]]),max(Stats$Random[[i]]),length=40)
            yfit<-dnorm(xfit,mean=mean(Stats$Random[[i]]),sd=sd(Stats$Random[[i]]))
            yfit <- yfit*diff(h$mids[1:2])*length(Stats$Random[[i]])
            
            lines(xfit, yfit, col="black", type="l",lwd=2,lty=2) 
            arrows(mean(Stats$Random[[i]]),max(h$counts)-mean(h$counts),Stats$Observed[[i]],max(h$counts)-mean(h$counts),code=1)
            arrows(mean(Stats$Random[[i]]),max(h$counts)-mean(h$counts),Stats$Observed[[i]],max(h$counts)-mean(h$counts),code=2)
            
            
            quantiles <- quantile(Stats$Random[[i]], prob=0.95)
            abline(v = as.numeric(quantiles), col = "green", lwd = 4)
            
            
          }
          dev.off()
          
          
