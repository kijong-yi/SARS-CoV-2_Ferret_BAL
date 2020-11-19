# for figure 5A

```{run on bash after installing velocyto.py}
#!/bin/bash
for i in cellranger/old/*/outs; do
	id=$(basename $(dirname $i))
	crout=$(dirname $i)
	gtf="~/Projects/covid19/ferret/singlecell/ref/NC_045512.2_and_NC_004718.3_and_MusPutFur1.0/genes/genes2.gtf"

	velocyto run10x --samtools-threads 4 $crout $gtf &> $crout/velocyto.log

done
qfire Q=100 S=1
```


library(velocyto.R)
library(tidyverse)
library(Seurat)

gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

MP <- read_rds("data/SO.ferret.MP.rds")

ldat <- system(intern=T,"ls processed_data/velocity/*.loom") %>%
	lapply(read.loom.matrices)

emat.bak <- ldat %>% lapply("[[","spliced") %>% do.call(cbind,.)
nmat.bak <- ldat %>% lapply("[[","unspliced") %>% do.call(cbind,.)


par(mfrow=c(1,2))
hist(log10(rowSums(emat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
hist(log10(rowSums(nmat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')

emat.bak = emat.bak[,colnames(MP)]
nmat.bak = nmat.bak[,colnames(MP)]

cell.colors = gg_color_hue(length(unique(MP$Annotation)))[as.numeric(factor(MP$Annotation))]
names(cell.colors) = colnames(emat.bak)


emat_ <- filter.genes.by.cluster.expression(emat.bak, cell.colors, min.max.cluster.average = 0.08)
nmat_ <- filter.genes.by.cluster.expression(nmat.bak, cell.colors, min.max.cluster.average = 0.06)

# look at the resulting gene set
nrow(emat_)
nrow(nmat_)
length(intersect(rownames(emat_),rownames(nmat_)))
x <- intersect(rownames(emat_),rownames(nmat_))

emat.bak = emat.bak[x,]
nmat.bak = nmat.bak[x,]

rm(x,emat_,nmat_)

write_rds(emat.bak,"data/emat.bak.MP.Rds")
write_rds(nmat.bak,"data/nmat.bak.MP.Rds")


# manual knn pooling
library(velocyto.R)
library(tidyverse)
library(Seurat)
library(doMC)
registerDoMC(12)

pool <- function(emat, nmat, nn, k=20, n = 5000,thread=20){
	j = ncol(emat)
	i = nrow(emat)
	O.e = matrix(0,nrow=i,ncol=n)
	O.n = matrix(0,nrow=i,ncol=n)
	if(n==0){
		perm = 1:j
		n = j
	}else{
		perm = sort(sample(j, n, F))
	}
	D=-nn
	iter = split(1:n, sort(1:n%%thread)) # thread
	registerDoMC(thread)
	O.e = foreach(it = iter,.combine="cbind") %dopar% {
		O = matrix(0,nrow=i,ncol=length(it))
		for(m in 1:length(it)){
			r = rank(D[,perm[it[m]]],ties.method = "random") <= k
			O[,m] = Matrix::rowSums(emat[,r])
		}
		O
	}
	cat("calc emat done\n")
	O.n = foreach(it = iter,.combine="cbind") %dopar% {
		O = matrix(0,nrow=i,ncol=length(it))
		for(m in 1:length(it)){
			r = rank(D[,perm[it[m]]],ties.method = "random") <= k
			O[,m] = Matrix::rowSums(nmat[,r])
			
		}
		O
	}
	colnames(O.e) = colnames(emat)[perm]
	rownames(O.e) = rownames(emat)
	colnames(O.n) = colnames(nmat)[perm]
	rownames(O.n) = rownames(nmat)
	cat("calc nmat done\n")
	list(emat = O.e, nmat = O.n, perm=perm)
}


MP$ID_2 = str_replace(MP$ID_1, "-.*","")

MP$ncov_idx <- MP$ID_2 %in% c("Ctrl", "C2","C5")

MP.ncov <- MP[,MP$ncov_idx]
MP.ncov <- Seurat::FindNeighbors(MP.ncov,dims=1:11)

MP.ncov@misc$pooled <- pool(emat.bak[,colnames(MP)[MP$ncov_idx]],nmat.bak[,colnames(MP)[MP$ncov_idx]],nn=MP.ncov@graphs$SCT_snn,k=100,n=8000,thread=12)

MP.ncov@misc$rvel.qf  = gene.relative.velocity.estimates(MP.ncov@misc$pooled$emat, MP.ncov@misc$pooled$nmat, deltaT=1, kCells=1, fit.quantile=0.05, n.cores=10)

MP.ncov@misc$emb_sample <- Embeddings(MP.ncov, "umap")[MP.ncov@misc$pooled$perm,]

rownames(MP.ncov@misc$emb_sample) = colnames(MP.ncov@misc$pooled$emat)

write_rds(MP.ncov, "data/velocity.SO.MP.ncov.Rds")

show.velocity.on.embedding.cor2 = function (emb,
											vel, 
											n = 100,
											cell.colors = NULL, 
											corr.sigma = 0.05, 
											show.grid.flow = FALSE, 
											grid.n = 20, 
											grid.sd = NULL, 
											min.grid.cell.mass = 1, 
											min.arrow.size = NULL, 
											arrow.scale = 1,
											max.grid.arrow.length = NULL, 
											fixed.arrow.length = FALSE, 
											plot.grid.points = FALSE,
											scale = "log", 
											nPcs = NA,
											arrow.lwd = 1, xlab = "", ylab = "", n.cores = velocyto.R:::defaultNCores(), 
											do.par = T, show.cell = NULL, cell.border.alpha = 0.3, cc = NULL, 
											return.details = FALSE, expression.scaling = FALSE, ...) {
	randomize <- FALSE
	if (do.par) 
		par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2, 
																															0.65, 0), cex = 0.85)
	celcol <- "grey"
	if (is.null(show.cell)) {
		celcol <- cell.colors[rownames(emb)]
	}
	plot(emb, bg = celcol, pch = 21, col = velocyto.R::ac(1, alpha = cell.border.alpha), 
			 xlab = xlab, ylab = ylab)#, ...)
	em <- as.matrix(vel$current)
	ccells <- intersect(rownames(emb), colnames(em))
	em <- em[, ccells]
	emb <- emb[ccells, ]
	nd <- as.matrix(vel$deltaE[, ccells])
	cgenes <- intersect(rownames(em), rownames(nd))
	nd <- nd[cgenes, ]
	em <- em[cgenes, ]
	if (randomize) {
		nd <- t(apply(nd, 1, function(x) (rbinom(length(x), 1, 
																						 0.5) * 2 - 1) * abs(sample(x))))
	}
	if (is.null(cc)) {
		cat("delta projections ... ")
		if (scale == "log") {
			cat("log ")
			cc <- velocyto.R:::colDeltaCorLog10(em, (log10(abs(nd) + 1) * 
																							 	sign(nd)), nthreads = n.cores)
		}
		else if (scale == "sqrt") {
			cat("sqrt ")
			cc <- velocyto.R:::colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)), 
																				 nthreads = n.cores)
		}
		else if (scale == "rank") {
			cat("rank ")
			cc <- velocyto.R:::colDeltaCor((apply(em, 2, rank)), (apply(nd, 
																																	2, rank)), nthreads = n.cores)
		}
		else {
			cat("linear ")
			cc <- velocyto.R:::colDeltaCor(em, nd, nthreads = n.cores)
		}
		colnames(cc) <- rownames(cc) <- colnames(em)
		diag(cc) <- 0
	}
	cat("knn ... ")
	if (n > nrow(cc)) {
		n <- nrow(cc)
	}
	emb.knn <- velocyto.R:::balancedKNN(t(emb), k = n, maxl = nrow(emb), dist = "euclidean", 
																			n.threads = n.cores)
	diag(emb.knn) <- 1
	cat("transition probs ... ")
	tp <- exp(cc/corr.sigma) * emb.knn
	tp <- t(t(tp)/Matrix::colSums(tp))
	tp <- as(tp, "dgCMatrix")
	cat("done\n")
	if (!is.null(show.cell)) {
		i <- match(show.cell, rownames(emb))
		if (is.na(i)) 
			stop(paste("specified cell", i, "is not in the embedding"))
		points(emb, pch = 19, col = ac(val2col(tp[rownames(emb), 
																							show.cell], gradient.range.quantile = 1), alpha = 0.5))
		points(emb[show.cell, 1], emb[show.cell, 2], pch = 3, 
					 cex = 1, col = 1)
		di <- t(t(emb) - emb[i, ])
		di <- di/sqrt(Matrix::rowSums(di^2)) * arrow.scale
		di[i, ] <- 0
		dir <- Matrix::colSums(di * tp[, i])
		dic <- Matrix::colSums(di * (tp[, i] > 0)/sum(tp[, i] > 
																										0))
		dia <- dir - dic
		suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i], 
																												 2], emb[colnames(em)[i], 1] + dic[1], emb[colnames(em)[i], 
																												 																					2] + dic[2], length = 0.05, lwd = 1, col = "blue"))
		suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i], 
																												 2], emb[colnames(em)[i], 1] + dir[1], emb[colnames(em)[i], 
																												 																					2] + dir[2], length = 0.05, lwd = 1, col = "red"))
		suppressWarnings(arrows(emb[colnames(em)[i], 1] + dic[1], 
														emb[colnames(em)[i], 2] + dic[2], emb[colnames(em)[i], 
																																	1] + dir[1], emb[colnames(em)[i], 2] + dir[2], 
														length = 0.05, lwd = 1, lty = 1, col = "grey50"))
		suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i], 
																												 2], emb[colnames(em)[i], 1] + dia[1], emb[colnames(em)[i], 
																												 																					2] + dia[2], length = 0.05, lwd = 1, col = "black"))
	}else {
		cat("calculating arrows ... ")
		arsd <- data.frame(t(velocyto.R:::embArrows(emb, tp, arrow.scale, n.cores)))
		nan_idx = is.nan(arsd$X1) | is.nan(arsd$X2) | is.na(arsd$X1) | is.na(arsd$X2)
		arsd$X1[nan_idx] <- emb[nan_idx,1]
		arsd$X2[nan_idx] <- emb[nan_idx,2]
		rownames(arsd) <- rownames(emb)
		if (expression.scaling) {
			tpb <- tp > 0
			tpb <- t(t(tpb)/colSums(tpb))
			es <- as.matrix(em %*% tp) - as.matrix(em %*% as.matrix(tpb))
			pl <- pmin(1, pmax(0, apply(as.matrix(vel$deltaE[, 
																											 colnames(es)]) * es, 2, sum)/sqrt(colSums(es * 
																											 																						es))))
			arsd <- arsd * pl
		}
		ars <- data.frame(cbind(emb, emb + arsd))
		colnames(ars) <- c("x0", "y0", "x1", "y1")
		colnames(arsd) <- c("xd", "yd")
		rownames(ars) <- rownames(emb)
		cat("done\n")
		if (show.grid.flow) {
			cat("grid estimates ... ")
			rx <- range(c(range(ars$x0), range(ars$x1)))
			ry <- range(c(range(ars$y0), range(ars$y1)))
			gx <- seq(rx[1], rx[2], length.out = grid.n)
			gy <- seq(ry[1], ry[2], length.out = grid.n)
			if (is.null(grid.sd)) {
				grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2)/2
				cat("grid.sd=", grid.sd, " ")
			}
			if (is.null(min.arrow.size)) {
				min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] - 
																											gy[1])^2) * 0.01
				cat("min.arrow.size=", min.arrow.size, " ")
			}
			if (is.null(max.grid.arrow.length)) {
				max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx), 
																												length(gy)))^2)) * 0.25
				cat("max.grid.arrow.length=", max.grid.arrow.length, 
						" ")
			}
			garrows <- do.call(rbind, lapply(gx, function(x) {
				cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x - emb[, 1])^2)
				cw <- dnorm(cd, sd = grid.sd)
				gw <- Matrix::colSums(cw)
				cws <- pmax(1, Matrix::colSums(cw))
				gxd <- Matrix::colSums(cw * arsd$xd)/cws
				gyd <- Matrix::colSums(cw * arsd$yd)/cws
				al <- sqrt(gxd^2 + gyd^2)
				vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
				cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg], gy[vg] + 
								gyd[vg])
			}))
			colnames(garrows) <- c("x0", "y0", "x1", "y1")
			if (fixed.arrow.length) {
				suppressWarnings(arrows(garrows[, 1], garrows[, 2], garrows[, 3], garrows[, 4], length = 0.05, lwd = arrow.lwd))
			}else {
				alen <- pmin(max.grid.arrow.length, sqrt(((garrows[, 3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1, 2)]))^2 + ((garrows[, 4] - garrows[, 2]) * par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
				suppressWarnings(lapply(1:nrow(garrows), function(i) arrows(garrows[i, 1], garrows[i, 2], garrows[i, 3], garrows[i, 4], length = alen[i], lwd = arrow.lwd)))
			}
			if (plot.grid.points) {
				points(rep(gx, each = length(gy)), rep(gy, length(gx)), 
							 pch = ".", cex = 0.1, col = ac(1, alpha = 0.4))
			}
				
			cat("done\n")
			if (return.details) {
				cat("expression shifts .")
				scale.int <- switch(scale, log = 2, sqrt = 3, 
														1)
				if (!expression.scaling) {
					tpb <- tp > 0
					tpb <- t(t(tpb)/colSums(tpb))
					es <- as.matrix(em %*% tp) - as.matrix(em %*% 
																								 	as.matrix(tpb))
				}
				cat(".")
				gs <- do.call(cbind, parallel::mclapply(gx, function(x) {
					cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x - 
																											emb[, 1])^2)
					cw <- dnorm(cd, sd = grid.sd)
					gw <- Matrix::colSums(cw)
					cws <- pmax(1, Matrix::colSums(cw))
					cw <- t(t(cw)/cws)
					gxd <- Matrix::colSums(cw * arsd$xd)
					gyd <- Matrix::colSums(cw * arsd$yd)
					al <- sqrt(gxd^2 + gyd^2)
					vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
					if (any(vg)) {
						z <- es %*% cw[, vg]
					}else {
						NULL
					}
				}, mc.cores = n.cores, mc.preschedule = T))
				if (scale == "log") {
					nd <- (log10(abs(nd) + 1) * sign(nd))
				}else if (scale == "sqrt") {
					nd <- (sqrt(abs(nd)) * sign(nd))
				}
				cat(".")
				gv <- do.call(cbind, parallel::mclapply(gx, function(x) {
					cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x - 
																											emb[, 1])^2)
					cw <- dnorm(cd, sd = grid.sd)
					gw <- Matrix::colSums(cw)
					cws <- pmax(1, Matrix::colSums(cw))
					cw <- t(t(cw)/cws)
					gxd <- Matrix::colSums(cw * arsd$xd)
					gyd <- Matrix::colSums(cw * arsd$yd)
					al <- sqrt(gxd^2 + gyd^2)
					vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
					if (any(vg)) {
						z <- nd %*% cw[, vg]
					}else {
						NULL
					}
				}, mc.cores = n.cores, mc.preschedule = T))
				cat(". done\n")
				return(invisible(list(tp = tp, cc = cc, garrows = garrows, 
															arrows = as.matrix(ars), vel = nd, eshifts = es, 
															gvel = gv, geshifts = gs, scale = scale)))
			}
		}else {
			apply(ars, 1, function(x) {
				if (fixed.arrow.length) {
					suppressWarnings(arrows(x[1], x[2], x[3], x[4], 
																	length = 0.05, lwd = arrow.lwd))
				}
				else {
					ali <- sqrt(((x[3] - x[1]) * par("pin")[1]/diff(par("usr")[c(1, 
																																			 2)]))^2 + ((x[4] - x[2]) * par("pin")[2]/diff(par("usr")[c(3, 
																																			 																													 4)]))^2)
					suppressWarnings(arrows(x[1], x[2], x[3], x[4], 
																	length = min(0.05, ali), lwd = arrow.lwd))
				}
			})
		}
	}
	return(invisible(list(tp = tp, cc = cc, ars=ars)))
}




draw_arrows <- function(ars, fixed.arrow.length=F, arrow.lwd=1,minlen=0.2,maxlen=4,...){
	idx_clean <-ars %>% with(((x1-x0)^2+(y1-y0)^2)^0.5) %>% {.<maxlen & .>minlen}
	ars = ars[idx_clean,]
	apply(ars, 1, function(x) {
		if (fixed.arrow.length) {
			suppressWarnings(arrows(x[1], x[2], x[3], x[4], 
															length = 0.05, lwd = arrow.lwd,...))
		}
		else {
			ali <- sqrt(((x[3] - x[1]) * par("pin")[1]/diff(par("usr")[c(1, 
																																	 2)]))^2 + ((x[4] - x[2]) * par("pin")[2]/diff(par("usr")[c(3, 
																																	 																													 4)]))^2)
			suppressWarnings(arrows(x[1], x[2], x[3], x[4], 
															length = min(0.05, ali), lwd = arrow.lwd,...))
		}
	})
	invisible()
}




draw_arrows.grid <- function(garrows, arrows,grid.n, fixed.arrow.length=F, min.arrow.size = NULL,
														 max.grid.arrow.length = NULL,
														 arrow.lwd=1,maxlen1 = 4, minlen1 = 0.1,...){
	ars = as.data.frame(arrows)
	rx <- range(c(range(ars$x0), range(ars$x1)))
	ry <- range(c(range(ars$y0), range(ars$y1)))
	gx <- seq(rx[1], rx[2], length.out = grid.n)
	gy <- seq(ry[1], ry[2], length.out = grid.n)
	grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2)/2
	cat("grid.sd=", grid.sd, " ")
	if (is.null(min.arrow.size)) {
		min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] - 
																									gy[1])^2) * 0.01
		cat("min.arrow.size=", min.arrow.size, " ")
	}
	if (is.null(max.grid.arrow.length)) {
		max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx), 
																										length(gy)))^2)) * 0.25
		cat("max.grid.arrow.length=", max.grid.arrow.length, 
				" ")
	}
	cat("\n")
	print(head(garrows))
	alen <- pmin(max.grid.arrow.length, sqrt(((garrows[, 3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1, 2)]))^2 + ((garrows[, 4] - garrows[, 2]) * par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
	print("1")
	glen = ((garrows[,"x1"]-garrows[,"x0"])^2+(garrows[,"y1"]-garrows[,"y0"])^2)^0.5
	print("1")
	idx_clean <- glen<maxlen1 & glen>minlen1
	print("2")
	garrows = garrows[idx_clean,]
	print("3")
	alen = alen[idx_clean]
	print("4")
	suppressWarnings(lapply(1:nrow(garrows), function(i) arrows(garrows[i, 1], garrows[i, 2], garrows[i, 3], garrows[i, 4], length = alen[i], lwd = arrow.lwd,...)))
}



MP.ncov@misc$velocity = show.velocity.on.embedding.cor2(emb=MP.ncov@misc$emb_sample, vel = MP.ncov@misc$rvel.qf, scale="sqrt", return.details=T)



MP.ncov@misc$velocity.grid = show.velocity.on.embedding.cor2(
	emb=MP.ncov@misc$emb_sample, vel = MP.ncov@misc$rvel.qf, scale="sqrt", return.details=T,
	cc = MP.ncov@misc$velocity$cc, show.grid.flow = T, grid.n = 50)

# MP.ncov %>% write_rds("data/velocity.SO.MP.ncov.re.Rds")
# MP.ncov =    read_rds("data/velocity.SO.MP.ncov.re.Rds")

MP <- read_rds("data/SO.ferret.MP.rds")


table(MPcov$annot)
cellcol = c("#F8766D","#DE8C00","#B79F00","#7CAE00",
						"#00BA38","#00C08B","#00BFC4","#00B4F0",
						"#619CFF","#C77CFF","#F564E3","#9E9E9E")[as.numeric(factor(MP$annot))]

DimPlot2 <- function(obj, group.by=NULL, reduction="umap",alpha=0.3,
										 legend.x = "topleft",legend.y=NULL,pch=16,legend.bty="n",
										 legend.ncol=3, screen = F,
										 ...){
	gg_color_hue <- function(n) {
		hues = seq(15, 375, length = n + 1)
		hcl(h = hues, l = 65, c = 100)[1:n]
	}
	if(is.null(group.by)){
		v = Idents(obj)
	} else{
		v = obj@meta.data[[group.by]]
	}
	colpal = gg_color_hue(length(unique(v)))
	col = adjustcolor(colpal[as.numeric(factor(v))],alpha.f = alpha)
	plot(Embeddings(obj, reduction),pch=pch,col=col,bty="n",cex=1,xaxt="n",yaxt="n",...)
	if(screen){
		rect(-20,-20,20,20,col="#FFFFFF40")
	}
	legend(legend.x,legend.y,pch=pch,col=colpal,legend = levels(factor(v)),
				 bty = legend.bty,ncol = legend.ncol)
}

MP$disease[is.na(MP$disease)] = "Ctrl"




# pdf("figures/fig3.velocyto.MP.pdf",7,7)
pdf("final/Fig4b.velocyto.MP.pdf",7,7)
par(mar=c(1,1,1,1))
DimPlot2(MP,
				 group.by="Annotation", 
				 screen = T,
				 legend.x = NA,
				 legend.ncol = 2)
draw_arrows.grid(garrows = MP.ncov@misc$velocity.grid2$garrows,
								 arrows = MP.ncov@misc$velocity.grid2$arrows,
								 grid.n = 50,
								 col="black")
dev.off()


MP.ncov %>% write_rds("data/velocity.SO.MP.ncov.re.Rds")

# arrow length comparison

ncov.velocity <- MP.ncov@misc$velocity
cellcol = gg_color_hue(length(unique(MP$finalannot_nonum)))[as.numeric(factor(MP$finalannot_nonum))]

alens  <- ((ncov.velocity$ars$x0 - ncov.velocity$ars$x1)^2 + (ncov.velocity$ars$y0 - ncov.velocity$ars$y1)^2)^0.5
names(alens) = rownames(ncov.velocity$ars)
alens.df <- data.frame(cellid = names(alens),
					 cellgr = MP.ncov$finalannot_nonum[names(alens)],
					 expgr = MP.ncov$expgr[names(alens)],
					 alen = alens)

library(ggpubr)
ggboxplot(data = alens.df,x = "cellgr", y = "alen", fill = "cellgr",
				 orientation = c("vertical", "horizontal", "reverse")[2])

