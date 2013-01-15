# running fastpopsim and plotting results

# Set the work directory
setwd("~/")


# The objective here is to find how genetic variation changes over time 
# when the causal variants are epistatic and affect fitness directly.

# To this end, populations are simulated with genetic effects segregating,
# and we will assess:
# - How long the mutation remains polymorphic under selection
# - What the variance decomposition is (additive vs non-additive)
# - How different statistical tests perform in detecting them

# We can choose to vary a number of different conditions including:
# - Population size
# - Sample size for the statistical tests
# - Which GP maps to assess
# - How many QTLs in each population
# - What is the LD between observed SNPs and causal variants
# - What are the starting allele frequencies of the causal variants
# - How many generations to run each population
# - How many repeats of each population to run


#####################
# SETUP SIMULATIONS #
#####################


# These are the epistatic patterns that we wish to assess
# They are 3x3 genotype phenotype maps (i.e. for 2 biallelic SNPs)
enum <- rbind(
	c(0,0.25,0.5,0.25,0.5,0.75,0.5,0.75,1), # Additive
	c(1,1,1,1,1,1,1,1,0), # Canalisation 1
	c(1,1,1,1,0,0,1,0,0), # Canalisation 2
	c(0.0000, 0.0000, 0.5000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000), 
	c(0.0000, 1.0000, 0.0000, 1.0000, 0.2000, 0.2000, 0.0000, 0.8000, 0.2000),
	c(0.0000, 0.9987, 0.0018, 1.0000, 0.0398, 0.7887, 0.0020, 0.7913, 0.0908)
)

# To visualise the GP maps, try this code:
gen <- expand.grid(0:2,0:2)
patdat <- data.frame()
for(i in 1:nrow(enum))
{
	temp <- rep(i,9)
	patdat <- rbind(patdat,cbind(gen,temp,enum[i,]))
}
colnames(patdat) <- c("A","B","pattern","Fitness")
require(ggplot2)
ggplot(patdat,aes(x=factor(B),y=factor(A))) + geom_tile(aes(fill=Fitness)) + scale_fill_gradient(low="grey",high="black") + facet_wrap(~pattern) + scale_y_discrete("Locus A") + scale_x_discrete("Locus B")

# Set up the different conditions for the simulations.
# We can try different ranges of correlation between causal and observed SNPs (LD)
# We can try different epistatic genotype phenotype maps
# We can try different starting frequencies

# Which patterns to use (count from 0)
a <- 1:nrow(enum)-1
npat <- length(a)

# Which rsq values to use
rsq <- c(1,0.85,0.7)

# Choose the starting frequencies of the GP maps (e.g. 0.5 for locus A and locus B)
s2 <- as.matrix(expand.grid(0.5,0.5,a,rsq))

# Here we have defined 18 different conditions
# They will run 6 GP maps over 3 different LD levels.
conditions <- nrow(s2)


freq <- array(0,c(conditions,1,2))
freq[,1,1] <- s2[,1]
freq[,1,2] <- s2[,2]
rsq <- array(0,c(conditions,1,2))
rsq[,1,1] <- s2[,4]
rsq[,1,2] <- s2[,4]
patterns <- matrix(s2[,3],conditions,1)

# Choose the effect sizes for each condition
eff <- array(1,c(conditions,1))

# Choose the population size for each condition
N <- rep(1000,conditions)

# Choose the proportion of the population to be sampled for statistical testing
nprop <- rep(1,conditions)

# Choose the number of QTLs to simulate for each condition
nqtl <- rep(1,conditions)

# Choose the proportion of the variance to be genetic (starting value)
hsq <- rep(0.1,conditions)

# How many repeats for each condition
repeats <- rep(50,conditions)

# How many generations to run each simulation for
ngen <- rep(200,conditions)


increment <- rep(1,conditions)
seed <- 1:conditions
conditions <- s2[,3]+1


###################
# RUN SIMULATIONS #
###################


# The 'runsims' function will take all the above parameters and make a parameter file,
# and then run the analysis and store the results in a text file,
# then it will read the output file in for analysis

output_filename <- "test.txt"
parameter_filename <- "example_parameter.txt"

runsims <- function(enum,N,nprop,nqtl,hsq,rsq,pattern,eff,p,generations,repeats,filename,increment,seed,condition,config)
{
	if(file.exists(filename))
	{
		print("File exists")
		return(NULL);
	}
	command <- paste("./epiFit",config,sep=" ")
	
	for(i in 1:length(conditions))
	{
		print(i)
		# write enum
		con <- file(config,"wt")
		writeLines(as.character(nrow(enum)),con)
		for(j in 1:nrow(enum))
		{
			writeLines(paste(enum[j,],collapse=" "),con)
		}
		writeLines(paste("\n",N[i],nprop[i],nqtl[i],sep="\n"),con)
		for(j in 1:nqtl[i])
		{
			writeLines(paste(pattern[i,j],eff[i,j],p[i,j,1],p[i,j,2],rsq[i,j,1],rsq[i,j,2],sep=" "),con)
		}
		writeLines(paste(hsq[i],generations[i],repeats[i],filename,increment[i],seed[i],condition[i],sep="\n"),con)
		close(con)
		system(command)
	}
	print("Reading in data")
	dat <- read.table(filename,header=T)
	return(dat)
}

sdat <- runsims(enum,N,nprop,nqtl,hsq,rsq,patterns,eff,freq,ngen,repeats,output_filename,increment,seed,conditions,parameter_filename)
sdat$maf1 <- 1 - sdat$maf1
sdat$maf2 <- 1 - sdat$maf2

unlink(output_filename)
save(sdat,file="sims.RData")
#load("simspap.RData")



#######################
# ANALYSE SIMULATIONS #
#######################



# An example of visualising the output from the simulations
# How do allele frequencies change over time?
dat2 <- subset(dat,Statistic == "A1" | Statistic == "A2")
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("SNP 1","SNP 2")
d <- ggplot(dat2,aes(x=Generation, y=AF)) + geom_line(size=0.2,aes(group=factor(as.numeric(Repeat)*10000+as.numeric(Statistic)),colour=Statistic)) + facet_wrap(~Pattern)

# What is the change in the proportion of additive variance over time?
b <- dat$Var
totadd <- b[seq(1,length(b),8)]+b[seq(2,length(b),8)]
totgen <- b[seq(1,length(b),8)]+b[seq(2,length(b),8)]+b[seq(3,length(b),8)]+b[seq(4,length(b),8)]+b[seq(5,length(b),8)]+b[seq(6,length(b),8)]+b[seq(7,length(b),8)]+b[seq(8,length(b),8)]

dat2 <- subset(dat,Statistic == "A1")
dat2$Var <- totadd/totgen
names(dat2)[11] <- "PropAdditive"
d <- ggplot(dat2,aes(x=Generation, y=PropAdditive))
d + geom_line(size=0.2,aes(group=factor(as.numeric(Repeat)*10000))) + geom_smooth() + ylim(0,1) + facet_wrap(~Pattern)
d + geom_line(size=0.2,aes(group=factor(as.numeric(Repeat)*10000))) + ylim(0,1) + facet_wrap(~Pattern)

# How do variance components change?
dat2 <- subset(dat,Statistic == "A1" | Statistic == "A2" | Statistic == "D1")
totdom <- b[seq(3,length(b),8)]+b[seq(4,length(b),8)]
totepi <- b[seq(5,length(b),8)]+b[seq(6,length(b),8)]+b[seq(7,length(b),8)]+b[seq(8,length(b),8)]
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("Additive","Dominance","Epistatic")
dat2$Var[dat2$Statistic=="Additive"] <- totadd
dat2$Var[dat2$Statistic=="Dominance"] <- totdom
dat2$Var[dat2$Statistic=="Epistatic"] <- totepi
dat2$Var[dat2$Var==0] <- NA
d <- ggplot(dat2,aes(x=Generation, y=Var))
d + geom_line(aes(group=factor(as.numeric(Repeat)*10000+as.numeric(Statistic)),colour=Statistic)) + facet_wrap(Pattern~Statistic)


# Compare tests - 8df test, 4df tests, marginal tests
# 1D tests
dat2 <- subset(dat,Statistic == "A1" | Statistic == "A2" | Statistic == "D1" | Statistic == "D2")
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("A1", "A2", "AD1", "AD2")
dat2$Pval[dat2$Pval == 0] <- NA
d <- ggplot(dat2,aes(x=Generation, y=Pval))
d + geom_smooth(size=0.5,aes(colour = Statistic)) + scale_y_continuous(lim=c(0.0005,200),trans="log10") + facet_wrap(~Pattern) + scale_colour_brewer("Test")

# 2D tests
dat2 <- subset(dat,Statistic == "AD" | Statistic == "DA" | Statistic == "AA")
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("Marginal","Epistatic","Full")
dat2$Pval[dat2$Pval == 0] <- NA
d <- ggplot(dat2,aes(x=Generation, y=Pval))
d + geom_smooth(size=1,aes(colour = Statistic)) + scale_y_continuous(lim=c(0.0005,200),trans="log10") + facet_wrap(~Pattern) + scale_colour_brewer("Test")


	# // A1 - additive pval
	# // A2 - additive pval
	# // D1 - A+D pval
	# // D2 - A+D pval
	# // AA - marginal 2d pval
	# // AD - epistatic 2d pval
	# // DA - Full 2d pval


# Power
# 1D
dat2 <- subset(dat,Statistic == "A1" | Statistic == "A2" | Statistic == "D1" | Statistic == "D2")
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("A1", "A2", "AD1", "AD2")

a <- tapply(dat2$Sig,list(dat2$Pattern,dat2$Generation,dat2$Statistic),sum) / repeats[1] * 100
gennames <- as.numeric(dimnames(a)[[2]])
patnames <- dimnames(a)[[1]]
statnames <- dimnames(a)[[3]]
dat3 <- data.frame(expand.grid(patnames,gennames,statnames),c(a))
colnames(dat3) <- c("Pattern","Generation","Statistic","Power")
dat3$Power[dat3$Power == 0] <- 0.1
d <- ggplot(dat3,aes(x=Generation,y=Power))
d + geom_smooth(size=1,aes(colour=Statistic)) + facet_wrap(~Pattern) + scale_y_continuous(trans="log10") + scale_colour_brewer("Test")

# 2D
dat2 <- subset(dat,Statistic == "AA" | Statistic == "AD" | Statistic == "DA")
dat2$Statistic <- drop.levels(dat2$Statistic)
levels(dat2$Statistic) <- c("Marginal","Epistatic","Full")
a <- tapply(dat2$Sig,list(dat2$Pattern,dat2$Generation,dat2$Statistic),sum) / repeats[1] * 100
gennames <- as.numeric(dimnames(a)[[2]])
patnames <- dimnames(a)[[1]]
statnames <- dimnames(a)[[3]]
dat3 <- data.frame(expand.grid(patnames,gennames,statnames),c(a))
colnames(dat3) <- c("Pattern","Generation","Statistic","Power")
dat3$Power[dat3$Power == 0] <- 0.1
d <- ggplot(dat3,aes(x=Generation,y=Power))
d + geom_smooth(size=1,aes(colour=Statistic)) + facet_wrap(~Pattern) + scale_y_continuous(trans="log10") + scale_colour_brewer("Test")


# How much of the marginal additive variance is explained using the different methods?
b <- dat$Var
totadd <- b[seq(1,length(b),8)]+b[seq(2,length(b),8)]
tempdat <- subset(dat,Statistic=="A1" | Statistic=="A2")
temp <- tempdat$Sig * tempdat$Var
temp2 <- tempdat$Var
totadddetected <- temp[seq(1,length(temp),2)] + temp[seq(2,length(temp),2)]
tempdat <- subset(dat,Statistic=="D1" | Statistic=="D2")
temp <- tempdat$Sig * temp2
totadddomdetected <- temp[seq(1,length(temp),2)] + temp[seq(2,length(temp),2)]
totepidetected <- totadd * dat[dat$Statistic=="AD",]$Sig
totmargdetected <- totadd * dat[dat$Statistic=="AA",]$Sig
totfulldetected <- totadd * dat[dat$Statistic=="DA",]$Sig

dat3 <- subset(dat,Statistic == "A1" | Statistic == "A2" | Statistic == "D1" | Statistic == "D2" | Statistic == "AA" | Statistic == "AD")
dat3$Statistic <- drop.levels(dat3$Statistic)
levels(dat3$Statistic) <- c("Additive", "Marginal 1D", "Marginal 2D", "Epistatic", "Full","Total")
dat3[dat3$Statistic=="Additive",]$Var <- totadddetected
dat3[dat3$Statistic=="Marginal 1D",]$Var <- totadddomdetected
dat3[dat3$Statistic=="Full",]$Var <- totfulldetected
dat3[dat3$Statistic=="Marginal 2D",]$Var <- totmargdetected
dat3[dat3$Statistic=="Epistatic",]$Var <- totepidetected
dat3[dat3$Statistic=="Total",]$Var <- totadd


a <- tapply(dat3$Var,list(dat3$Pattern,dat3$Generation,dat3$Statistic),sum)

a[,,1] <- a[,,1] / a[,,6] * 100
a[,,2] <- a[,,2] / a[,,6] * 100
a[,,3] <- a[,,3] / a[,,6] * 100
a[,,4] <- a[,,4] / a[,,6] * 100
a[,,5] <- a[,,5] / a[,,6] * 100
a <- a[,,-6]

dat4 <- data.frame(expand.grid(dimnames(a)[[1]],as.numeric(dimnames(a)[[2]]),dimnames(a)[[3]]),c(a))
colnames(dat4) <- c("Pattern","Generation","Statistic","PercVar")
 
d <- ggplot(dat4,aes(x=Generation,y=PercVar))
d + geom_line(size=0.5,aes(colour=Statistic)) + facet_wrap(~Pattern) + scale_colour_brewer("Test") + scale_y_log10()

d <- ggplot(dat4,aes(x=Generation,y=PercVar))
d + geom_line(aes(colour=Statistic)) + facet_wrap(~Pattern) + scale_colour_brewer("Test")

b <- tapply(dat4$PercVar,list(dat4$Pattern,dat4$Statistic),function(x){sum(x,na.rm=T)})
dat5 <- data.frame(expand.grid(as.numeric(dimnames(b)[[1]]),dimnames(b)[[2]]),c(b))
names(dat5) <- c("Pattern","Test","VarExplained")
d <- ggplot(dat5,aes(Test,VarExplained))
d + geom_bar(aes(fill=Test)) + facet_wrap(~Pattern) + scale_colour_brewer("Test")



#########
# OTHER #
#########

# These are the patterns used in the supplementary materials

enum <- rbind(c(0,0,0,0,0,0,0,0,0), # NULL
c(0,0,0,0,0,0,0,0,1),
c(0,0,0,0,0,0,0,1,0),
c(0,0,0,0,0,0,0,1,1),
c(0,0,0,0,0,0,1,0,1),
c(0,0,0,0,0,0,1,1,1), # A+D
c(0,0,0,0,0,1,0,1,0),
c(0,0,0,0,0,1,0,1,1),
c(0,0,0,0,0,1,1,0,0),
c(0,0,0,0,0,1,1,0,1),
c(0,0,0,0,0,1,1,1,0),
c(0,0,0,0,0,1,1,1,1),
c(0,0,0,0,1,0,0,0,0),
c(0,0,0,0,1,0,0,0,1),
c(0,0,0,0,1,0,0,1,0),
c(0,0,0,0,1,0,0,1,1),
c(0,0,0,0,1,0,1,0,1),
c(0,0,0,0,1,0,1,1,1),
c(0,0,0,0,1,1,0,1,0),
c(0,0,0,0,1,1,0,1,1),
c(0,0,0,0,1,1,1,0,0),
c(0,0,0,0,1,1,1,0,1),
c(0,0,0,0,1,1,1,1,0),
c(0,0,0,1,0,1,0,0,0),
c(0,0,0,1,0,1,0,0,1),
c(0,0,0,1,0,1,0,1,0),
c(0,0,0,1,0,1,0,1,1),
c(0,0,0,1,0,1,1,0,1),
c(0,0,0,1,1,1,0,0,0), # D
c(0,0,0,1,1,1,0,0,1),
c(0,0,0,1,1,1,0,1,0),
c(0,0,0,1,1,1,0,1,1),
c(0,0,0,1,1,1,1,0,1),
c(0,0,1,0,0,0,1,0,0),
c(0,0,1,0,0,0,1,0,1),
c(0,0,1,0,0,0,1,1,0),
c(0,0,1,0,0,1,1,1,0),
c(0,0,1,0,1,0,1,0,0),
c(0,0,1,0,1,0,1,0,1),
c(0,0,1,0,1,0,1,1,0),
c(0,0,1,0,1,1,1,1,0),
c(0,0,1,1,0,0,0,0,1),
c(0,0,1,1,0,0,0,1,0),
c(0,0,1,1,0,0,0,1,1),
c(0,0,1,1,0,0,1,0,1),
c(0,0,1,1,0,1,0,1,0),
c(0,0,1,1,0,1,1,0,0),
c(0,0,1,1,1,0,0,0,1),
c(0,0,1,1,1,0,0,1,0),
c(0,1,0,1,0,1,0,1,0), # DxD
c(0,1,0,1,1,1,0,1,0),
c(1,0.5,0,0.5,0.5,0.5,0,0.5,1), # AxA
c(0,1,0,0.5,0.5,0.5,1,0,1), # AxD
c(1,1,1,0,0,0,1,1,1),
c(1,1,1,1,0,0,1,0,0),
c(0,0,0,0.5,0.5,0.5,1,1,1)) # A


enum <- abs(enum-1)

gen <- expand.grid(0:2,0:2)
patdat <- data.frame()
for(i in 1:nrow(enum))
{
	temp <- rep(i,9)-1
	patdat <- rbind(patdat,cbind(gen,temp,enum[i,]))
}
colnames(patdat) <- c("A","B","pattern","Fitness")

d <- ggplot(patdat,aes(x=factor(B),y=factor(A)))
d + geom_tile(aes(fill=Fitness)) + scale_fill_gradient(low="grey",high="black") + facet_wrap(~pattern)


