library(corrplot)

### For small p, show the matrix structure & sparsity
par(mfrow = c(2,3))
corrplot(Uns,
         is.corr = FALSE,
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         tl.srt=45, # text label rotation
         #addCoef.col="black", # colour of coefficients
         order="AOE", # ordering method
         #col=colorRampPalette(c("white","gray","black"))(10)
         mar=c(0,0,3,0),
         title = "Unstructured, random"
)

corrplot(AR1,
         is.corr = FALSE,
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         tl.srt=45, # text label rotation
         #addCoef.col="black", # colour of coefficients
         order="AOE", # ordering method
         mar=c(0,0,3,0),
         title = "Structured, AR(1)"
)

corrplot(Exp,
         is.corr = FALSE,
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         tl.srt=45, # text label rotation
         #addCoef.col="black", # colour of coefficients
         order="AOE", # ordering method
         mar=c(0,0,3,0),
         title = "Structured, Exponential"
)

corrplot(Blk,
         is.corr = FALSE,
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         tl.srt=45, # text label rotation
         #addCoef.col="black", # colour of coefficients
         order="AOE", # ordering method
         mar=c(0,0,3,0),
         title = "Structured, Block diagonal"
)

corrplot(Fct,
         is.corr = FALSE,
         method="shade", # visualisation method
         shade.col=NA, # colour of shade line
         tl.col="black", # colour of text label
         tl.srt=45, # text label rotation
         #addCoef.col="black", # colour of coefficients
         order="AOE", # ordering method
         mar=c(0,0,3,0),
         title = "Structured, Factor model"
)


### Uncertainty quantification
# 1. SBFIM
# 95\% posterior credible interval for selected entries of \Sigma
quantile(SBIFM[1,1,], c(0.025, 0.975))

# average estimated number of factors
quantile(numFacts)

