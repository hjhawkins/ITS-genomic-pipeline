
##########################################################################
######################## BETA DIVERSITY (between-sample)##################
##########################################################################

#### 4 steps: ;  

# (1) transform; for beta, pre-processing is necessary

# (2) calculate distance: many options of distance measurement = ##bray## (= most common), Unifrac, wunifrac, dpcoa, jsd, manhattan, euclidean, Canberra, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, maximum, binary, minkowski, etc 

# (3) ordinate (many options); ordination methods in plyloseq are:for unconstrained -  "CCA", "DPCoA", "NMDS", "MDS", "PCoA", "DCA", "RDA","CAP", the most commom = ##PCoA and NMDS##.for constrained : "CCA", and specify the model (phyloseq~treatment, "CCA")

# (4) plot (per sample, per OTU, biplot, etc): types of plot : type = c("samples", "sites", "species", "taxa", "biplot", "split", "scree")


# Callahan

log_phyf2 = transform_sample_counts(phyf2, function(x) log(1 + x) )

#katie

total = median(sample_sums(phy))#find median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
M.std = transform_sample_counts(phy, standf)#apply to phyloseq object
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)

#with tree
M.stdtree = transform_sample_counts(phytree, standf)#apply to phyloseq object
M.ftree = filter_taxa(M.std,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)


# two different and widely-used ecological distances: Bray-Curtis and weighted UniFrac, are sensitive to differences in total counts and a deluge of rare taxa. Let's transform to relative abundance, and also apply those filtering criteria that we explored earlier. http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html

mdt = fast_melt(phy)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Define taxa to keep.
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
# Define new object with relative abundance
mpra = transform_sample_counts(mp, function(x) x / sum(x))
# Filter this new object
mpraf = prune_taxa(keepTaxa, mpra)



#########################################################
################## UNCONSTRAINED ORDINATION #############



#from Dan Knights video (vegan) (look at the code...)

d.euc <- dist(OTU) #euclidian
pc.euc <- cmdscale(d.euc, k=2) #Principal coordinates

CS <- cca(otu_table(phy)) # chi square. didn't work with cca(OTU)
PC <- CS$CA$u[,1:2]

BC <- vegdist(otu_table(phy)) #bray curtis = default for vegan
PCB <- cmdscale(BC, k=2)
plot(PCB[,1], PCB[,2], col = ..., cex = 3, pch = 16)


### PCoA and NMDS - BRAY-CURTIS distances ###

#D <- phyloseq::distance(logt, method = "bray") # distance () didn't work
#O <- ordinate(logt, "PCoA", distance = "D")
#p <- plot_ordination(logt, O, color="...", shape="...")


PCoA_bray_logphyf2 <- ordinate(log_phyf2, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
#evals <- PCoA_bray_logphyf2$values$Eigenvalues

p <- plot_ordination(log_phyf2, PCoA_bray_logphyf2, color="site", shape="treatm", title= "PCoA - Bray curtis distance - fungal populations") + coord_fixed(sqrt(evals[2]/evals[1])) + theme_bw()
p


p <- plot_ordination(log_phyf2, PCoA_bray_logphyf2, color="site", shape="FH", title= "PCoA - Bray curtis distance - fungal populations") + theme_bw()
p + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
p + stat_ellipse(type = "t")


set.seed(2) #For NMDS plots it's important to set a seed since the starting positions of samples in the alogrithm is random.select any start seed, in this case 2 (ensures reproducibility) 
NMDS_bray_logphyf2 <- ordinate(log_phyf2, "NMDS", "bray") #other arguments that might be added: k=2, trymax=100) # stress=0.06
evals <- NMDS_bray_logphyf2$values$Eigenvalues
p <- plot_ordination(log_phyf2, NMDS_bray_logphyf2, color="site", shape="FH", title= "NMDS - Bray curtis distance - fungal populations") + theme_bw()
p



GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100) # stress=0.06

color = c("FH")
shape = c("site")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")

MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"), 
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+ 
  theme_bw() + labs(color=color, shape=shape) + geom_point(size=5)
MDS.1
pdf(paste0(outDir,"/NMDS_tretment_Bray_Curtis.pdf"),8,5)
MDS.1
dev.off()


out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
out.pcoa.M.f <- ordinate(M.f, method = "PCoA", distance = "bray")

evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.M.f, type = "samples", color = c("site"), shape = c("FH"), title = "PCoA of bacterial Communities - BC distance"
) + labs(col = "Treatment (F and H)") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_bw()

# you can also add a facet_wrap(~Phylum, 3) 



### UNIFRAC distance (which takes OTU relatedness into account) ####

GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")
GP.ord.U

color = c("Treatm")
#shape = c("dog)
title=c("NMDS of 16S microbiome, Unifrac distance,k=2")
MDS = plot_ordination(M.f, GP.ord.U, color = color, title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
  theme_bw()+labs(color=color)+geom_point(size=5)
MDS.1

pdf(paste0(outDir,"/NMDS_treatment_Unifrac.pdf"),8,5)
MDS.1
dev.off()


### you can plot all the methods and compare them. http://joey711.github.io/phyloseq/plot_ordination-examples.html

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="site")
}, M.ftree, dist)

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=SampleType, shape=human, fill=SampleType))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
p

# If you want to replot a larger version of an individual plot, you can do by printing from the original plist from which pdataframe was made. Each element of plist is already a ggplot2 graphic. For example, printing the second element of the list.
plist[[2]] 



############ OTU plot and BIPLOT #############################


plot_ordination(M.f, out.pcoa.M.f, type = "taxa", color = c("Phylum"),shape = "site", title = "PCoA of bacterial Communities - BC distance") + labs(col = "Treatment (F and H)") +  theme_bw()

plot_ordination(M.f, out.pcoa.M.f, type = "split", color = c("Phylum"),shape = "site", title = "PCoA of bacterial Communities - BC distance") + labs(col = "Treatment (F and H)") +  theme_bw()

plot_ordination(M.f, out.pcoa.M.f, type = "split", color = c("site"), shape = c("FH"), title = "PCoA of bacterial Communities - BC distance"
) + labs(col = "Treatment (F and H)") + theme_bw()













