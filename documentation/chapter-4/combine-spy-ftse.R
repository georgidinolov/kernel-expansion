library("data.table")
ftse <- fread("FTSE.csv")
spy <- fread("SPY.csv")

ftse[, Date := as.Date(Date, "%Y-%m-%d")]
setkey(ftse, Date)

spy[, Date := as.Date(Date, "%Y-%m-%d")]
setkey(spy, Date)

all.data <- merge(x=spy, y=ftse, by="Date", all=FALSE, suffixes=c(".spy", ".ftse"))

plot(all.data[, Date], all.data[, log(Close.spy) - log(Open.spy)], type="l")
lines(all.data[, Date], all.data[, log(Close.ftse) - log(Open.ftse)], type="l", col="red")

for (i in seq(1,nrow(all.data))) {
    if (i==1) {
        write(x=paste0("rho=",0.0,";\n",
                       "sigma_x=",1.0,";\n",
                       "sigma_y=",1.0,";\n",
                       "x_0=",log(all.data[i,Open.spy]),";\n",
                       "y_0=",log(all.data[i,Open.ftse]),";\n",
                       "x_T=",log(all.data[i,Close.spy]),";\n",
                       "y_T=",log(all.data[i,Close.ftse]),";\n",
                       "t=",1,";\n",
                       "a=",log(all.data[i,Low.spy]),";\n",
                       "b=",log(all.data[i,High.spy]),";\n",
                       "c=",log(all.data[i,Low.ftse]),";\n",
                       "d=",log(all.data[i,High.ftse]),";"), file="spy-ftse.csv")
    } else {
        write(x=paste0("rho=",0.0,";\n",
                       "sigma_x=",1.0,";\n",
                       "sigma_y=",1.0,";\n",
                       "x_0=",log(all.data[i,Open.spy]),";\n",
                       "y_0=",log(all.data[i,Open.ftse]),";\n",
                       "x_T=",log(all.data[i,Close.spy]),";\n",
                       "y_T=",log(all.data[i,Close.ftse]),";\n",
                       "t=",1,";\n",
                       "a=",log(all.data[i,Low.spy]),";\n",
                       "b=",log(all.data[i,High.spy]),";\n",
                       "c=",log(all.data[i,Low.ftse]),";\n",
                       "d=",log(all.data[i,High.ftse]),";"),
              append=TRUE,
              file="spy-ftse.csv")
    }
}

