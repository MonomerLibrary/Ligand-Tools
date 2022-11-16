
table_fn = "correlation-tables/results-8ers.table"
table_fn = "correlation-tables/results-8ewa.table"
table_fn = "correlation-tables/results-8dus.table"
table_fn = "correlation-tables/results-8du9.table"
table_fn = "correlation-tables/results-8bb5.table"
table_fn = "correlation-tables/results-8aqn.table"
table_fn = "correlation-tables/results-8duc.table"
table_fn = "correlation-tables/results-8a92.table"
table_fn = "correlation-tables/results-7fr6.table"

a = read.table(table_fn, header=TRUE)

plot(a$rscc1, a$rscc2, xlim=c(0.3, 1.0), ylim=c(0.3, 1.0), xlab="correl", ylab="correl")
abline(a=0, b=1)
