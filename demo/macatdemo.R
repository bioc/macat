library(macat)

demo.part(1)
loaddatapkg("stjudem")
#data(stjude)

demo.part(2)
summary(stjude)

demo.part(3)              
stjude$geneName[1:10]
unique(stjude$labels)
table(stjude$labels)

demo.part(4)

length(stjude$geneName[stjude$chromosome==1])

demo.part(5)
plotSliding(stjude, chrom=6, sample=3, kernel=rbf) 



demo.part(6)
e1 = evalScoring(stjude, class="T", chromosome=6, nperms=1000)

demo.part(7)
e2 = evalScoring(stjude, class="T", chromosome=6, nperms=1000, kernel=kNN)


demo.part(8)
x11(width=12,height=12)
par(mfrow=c(2,1))
plot(e1, output="x11",new.device=FALSE)
mtext("rbf-kernel",3,line=0,font=2)
axis(1)
mtext("Coordinate",1,line=3,font=1)
plot(e2, output="x11",new.device=FALSE)
mtext("kNN-kernel",3,line=0,font=2)
axis(1)
mtext("Coordinate",1,line=3,font=1)

demo.part(9)
plot(e1, output="html")

