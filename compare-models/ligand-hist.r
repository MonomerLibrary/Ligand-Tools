
st=format(Sys.time(), '%d-%m-%Y')
pdf(paste("acedrg-ligands_",st, ".pdf", sep = ""))
x=rnorm(100,40,3)
y=rnorm(100,100,5)

par(mfrow=c(3,3), cex=0.5)

a = read.table("15-11-refmacs-deltas.table")

# hist(a$V6, breaks=30, col="lightblue", main="Overall Displacement Histogram",
#     xlab="Displacement (A)")

factors = factor(a$V3)
levels = levels(factors)
print(levels)

for (l in levels) {
    ligand_code = a$V2[a$V3==l][1]
    reso        = a$V4[a$V3==l][2]
    mt = l
    mt = paste(mt, ligand_code)
    mt = paste(mt, reso)
    mt = paste(mt,'\uc5')

    hist(a$V6[a$V3==l], col="lightblue", main=mt,
         xlim=c(0, 1.2), xlab="Displacement (A)")
    correl_table_fn <- paste("correlation-tables/results-", tolower(l), sep="")
    correl_table_fn <- paste(correl_table_fn, ".table", sep="")
    print(correl_table_fn)
    b <- read.table(correl_table_fn, header=TRUE)
    plot(b$rscc1, b$rscc2, xlim=c(0.3, 1.0), ylim=c(0.3, 1.0), xlab="correl", ylab="correl", main=mt)
    abline(a=0, b=1)
    plot(b$dist, b$b2, xlab="dist", ylab="b2", main=mt)
}

dev.off()
