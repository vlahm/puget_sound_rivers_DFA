#run pca
env.pca<-prcomp(envdata, scale=TRUE) #uses correlation matrix (standardizes, so this
#is good for vars that are on different scales) - use princomp() if all vars are on
#same scale, as it works off the var-covar matrix
#chez it
summary(env.pca)
#eigenvalues (variance explained by each PC)
env.eig<-env.pca$sdev^2
pca.eigenval(env.pca) #another way to see all eigenvalues (from BIOSTATS.R)
#check significance of eigenvalues (of PCs)
screeplot(env.pca, bstick=TRUE) #lame way
#this says PC1 is nonsig, but ecologically, that much var explained is sig
ordi.monte(envdata,ord='pca',dim=5) #better way (limited to first 5 PCs for speed)
#check printout after the plots, p vals < 0.05 are sig
#check eigenvectors (variable loadings on each PC)
env.pca$rotation
#shows correlation of original variables with each PC - high abs vals show high contribution
#toward linear combination comprising each PC
pca.eigenvec(env.pca,dim=5,digits=3,cutoff=.3) #alternative way
#convert eigenvector coefficients to correlation coefficients
pca.structure(env.pca,envdata,dim=5,cutoff=.4)
#these are correlations between the original variables and the PC scores
#squaring these gives percentage of variance in each original var accounted
#for by each PC
#use these or the eigenvectors to generate an ecological interpretation of
#each PC
#get sample scores
env.scores<-env.pca$x
#standardized scores for each sample on each PC axis - they represent the position
#of each sample on each standardized PC axis
#visualize
biplot(env.pca) #easy way
ordiplot(env.pca, choices = c(1, 2), type="text", display='sites', xlab='PC 1 (27%)', ylab='PC 2 (18%)')
arrows(0,0,env.pca$rotation[,1]*5,env.pca$rotation[,2]*5,col='purple')
text(env.pca$rotation[,1]*5.2,env.pca$rotation[,2]*5.2,row.names(env.pca$rotation))
