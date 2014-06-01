#------------------------------------------------------
#  plot of first four eigen vectors

if (Plot.png) png(file="eigenVectors.png")
par(mfrow=c(2,2))

siz <- 2*PC$rotation[,1]/max(abs(PC$rotation[,1]))
is.up <- siz > 0
siz <- abs(siz)
plot.phylo(rna.tree,
           cex=0.5,
           no.margin=T,
           show.tip.label=F,
           lab4ut="axial",
           type="unrooted",
           use.edge.length=F,
           edge.color=gray(0.65))
colorEdge(rna.tree,brk,200,sat=0.4)

nodelabels("",node.internal[is.up],
                cex=siz[is.up],
                pch=24,
                frame="none")

#  down pointing triangles for decreases
nodelabels("",node.internal[!is.up],
                cex=siz[!is.up],
                pch=25,
                frame="none")

siz <- 2*PC$rotation[,2]/max(abs(PC$rotation[,2]))
is.up <- siz > 0
siz <- abs(siz)
plot.phylo(rna.tree,
           cex=0.5,
           no.margin=T,
           show.tip.label=F,
           lab4ut="axial",
           type="unrooted",
           use.edge.length=F,
           edge.color=gray(0.65))
colorEdge(rna.tree,brk,200,sat=0.4)

nodelabels("",node.internal[is.up],
                cex=siz[is.up],
                pch=24,
                frame="none")

#  down pointing triangles for decreases
nodelabels("",node.internal[!is.up],
                cex=siz[!is.up],
                pch=25,
                frame="none")

siz <- 2*PC$rotation[,3]/max(abs(PC$rotation[,3]))
is.up <- siz > 0
siz <- abs(siz)
plot.phylo(rna.tree,
           cex=0.5,
           no.margin=T,
           show.tip.label=F,
           lab4ut="axial",
           type="unrooted",
           use.edge.length=F,
           edge.color=gray(0.65))
colorEdge(rna.tree,brk,200,sat=0.4)

nodelabels("",node.internal[is.up],
                cex=siz[is.up],
                pch=24,
                frame="none")

#  down pointing triangles for decreases
nodelabels("",node.internal[!is.up],
                cex=siz[!is.up],
                pch=25,
                frame="none")

siz <- 2*PC$rotation[,4]/max(abs(PC$rotation[,4]))
is.up <- siz > 0
siz <- abs(siz)
plot.phylo(rna.tree,
           cex=0.5,
           no.margin=T,
           show.tip.label=F,
           lab4ut="axial",
           type="unrooted",
           use.edge.length=F,
           edge.color=gray(0.65))
colorEdge(rna.tree,brk,200,sat=0.4)

nodelabels("",node.internal[is.up],
                cex=siz[is.up],
                pch=24,
                frame="none")

#  down pointing triangles for decreases
nodelabels("",node.internal[!is.up],
                cex=siz[!is.up],
                pch=25,
                frame="none")

if (Plot.png) graphics.off()



