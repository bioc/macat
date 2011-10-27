library(macat)

demo.part(1)
loaddatapkg("stjudem")
data(stjude)

demo.part(2)
summary(stjude)

demo.part(3)              
stjude$geneName[1:10]
unique(stjude$labels)
table(stjude$labels)

demo.part(4)
sum(stjude$chromosome==1)

demo.part(5)
plotSliding(stjude, chrom=6, sample=3, kernel=rbf) 

demo.part(6)
evalknn6 <- evaluateParameters(stjude, class="T",
                               chromosome=6,kernel=kNN,
                               paramMultipliers=c(0.01,seq(0.25,2.0,0.25)))

demo.part(7)
plot(evalknn6)

demo.part(8)
evalres <- evalScoring(stjude, class="T", chromosome=6, nperms=1000,
                       kernel=kNN, kernelparams=evalknn6$best,
                       cross.validate = FALSE)

demo.part(9)
plot(evalres)

demo.part(10)
plot(evalres, output="html")

