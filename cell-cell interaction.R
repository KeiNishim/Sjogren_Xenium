0) Identities and quick sanity plots
Idents(xenium.obj) <- "TLS"
xenium.obj <- subset(xenium.obj, subset = !is.na(TLS))

Examples (adjust features/group.by as needed)
DotPlot(xenium.obj, features = c("TGFB1","SMAD2","SMAD3","TBX21","STAT1","IRF1"), group.by = "TLS")
DotPlot(xenium.obj, features = c("SMAD3"), group.by = "clusters")

1) Subset for analysis (example: IgA-related)
s.sub <- subset(xenium.obj, idents = c("IgA_TLS", "IgA_PC"))
Idents(s.sub) <- "clusters"

Optional quick plots
DotPlot(s.sub, features = c("CCR7","CXCR5","CXCR4","CXCR3"), group.by = "clusters")
DotPlot(s.sub, features = c("CXCL12","CCL19","CCL21","CXCL9"), group.by = "clusters")
VlnPlot(s.sub, features = c("IFNG_exposure","STAT1","IFNGR1"), pt.size = 0)

2) Normalize Xenium assay and prepare matrices
DefaultAssay(s.sub) <- "Xenium"
s.sub <- NormalizeData(s.sub)
data.input <- GetAssayData(s.sub, slot = "data", assay = "Xenium") # genes x cells
cells_keep <- colnames(data.input)

3) Coordinates: robust extractor
get_centroid_coords <- function(obj, offset_start = 0, offset_step = 0) {
md <- obj@meta.data
rn <- rownames(md)

candidate column names for x,y
cand_x <- c("x","X","pos_x","pxl_col","col","X_local")
cand_y <- c("y","Y","pos_y","pxl_row","row","Y_local")
cx <- cand_x[cand_x %in% colnames(md)][1]
cy <- cand_y[cand_y %in% colnames(md)][1]
if (!is.null(cx) && !is.null(cy)) {
coords <- as.matrix(md[, c(cx, cy)])
colnames(coords) <- c("x","y")
rownames(coords) <- rn
return(coords)
} else {
stop("Cannot find x/y columns in meta.data. Please add numeric columns 'x' and 'y' to s.sub@meta.data.")
}
}

centroid_coords <- get_centroid_coords(s.sub, offset_start = 10000, offset_step = 10000)

Align coordinates/meta/data rows
stopifnot(identical(rownames(centroid_coords), rownames(s.sub@meta.data)))
meta <- s.sub@meta.data
meta <- meta[cells_keep, , drop = FALSE]
spatial.locs <- centroid_coords[cells_keep, , drop = FALSE]
stopifnot(identical(colnames(data.input), rownames(meta)),
identical(rownames(spatial.locs), rownames(meta)))

4) Colors by clusters (consistent mapping)
cluster_names <- sort(unique(s.sub$clusters))
n_clusters <- length(cluster_names)
colors <- colorRampPalette(brewer.pal(12, "Paired"))(max(n_clusters, 3))
cluster_colors <- setNames(colors[seq_len(n_clusters)], cluster_names)

Highlight selected clusters if desired
hl <- c("FRC-like", "Immunofibro")
cluster_colors2 <- cluster_colors
others <- setdiff(names(cluster_colors2), hl)
cluster_colors2[others] <- lighten(cluster_colors2[others], amount = 0.7)
cluster_colors2[hl] <- darken(cluster_colors2[hl], amount = 0.3)

5) Spatial factors estimation from nearest-neighbor distances
k <- 20
dnn <- get.knn(spatial.locs, k = k)
min_distance <- min(dnn$nn.dist[dnn$nn.dist > 0])
conversion.factor <- 0.21 # pixels to micrometers (adjust to your instrument)
spot.size <- min_distance * conversion.factor
spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size / 2)

6) Create CellChat object (spatial mode)
cellchat <- createCellChat(
object = data.input,
meta = meta,
group.by = "clusters",
datatype = "spatial",
coordinates = spatial.locs,
spatial.factors = spatial.factors
)
cellchat@DB <- CellChatDB.human

7) Preprocessing and feature selection
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, min.cells = 1, thresh.pc = 0.001, thresh.fc = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)

Helper to inspect var.features
getFeature <- function(cellchat, feats, clusters = NULL, ignore.case = FALSE,
min_logFC = NULL, max_padj = NULL) {
df <- cellchat@var.features$features.info
m <- if (ignore.case) {
Reduce(|, lapply(feats, function(f) grepl(paste0("^", f, "$"), df$features, ignore.case = TRUE)))
} else {
df$features %in% feats
}
if (!is.null(clusters)) m <- m & df$clusters %in% clusters
if (!is.null(min_logFC)) m <- m & df$logFC >= min_logFC
if (!is.null(max_padj) && "pvalues.adj" %in% colnames(df)) m <- m & df$pvalues.adj <= max_padj
df[m, , drop = FALSE]
}

Examples:
getFeature(cellchat, c("CCR7","CCL21","CCL19","CXCL13","CXCL12","CXCL9","CXCR3","CCL5","CXCR5","CXCR4","TGFB1","TGFB2","TGFBR1","TGFBR2","IFNG","IFNGR1","ACVR1"))

8) Compute communication probability (tune distance params)
Tips:
- interaction.range/contact.range are in micrometers if you set spatial.factors$ratio accordingly.
- scale.distance scales the decay of distance; increase to penalize long distances more strongly.
Aim: average CCI ~ 1–2 for aggregate tissues (adjust scale.distance 0.1–2.0 accordingly).
cellchat <- computeCommunProb(
object = cellchat,
type = "truncatedMean",
trim = 0.1,
k.min = 1,
nboot = 50,
population.size = TRUE,
distance.use = TRUE,
interaction.range = 50,
contact.range = 20,
scale.distance = 0.2 # try 0.1–0.5; increase if CCIs are too large
)

9) Pathway probability, aggregate, and diagonal zeroing
cellchat <- computeCommunProbPathway(cellchat, thresh = 1)
cellchat <- aggregateNet(cellchat, thresh = 1)
cellchat <- filterCommunication(cellchat, min.cells = 1)

Set diagonal to zero for all pathway slices
for (i in 1:dim(cellchat@net$prob)[3]) {
diag(cellchat@net$prob[, , i]) <- 0
}

available_pathways <- cellchat@netP$pathways

10) Probability-based filtering and recompute pathway aggregation
df.signal <- subsetCommunication(cellchat, signaling = available_pathways)
thresh_prob <- quantile(df.signal$prob, 0.5) # median threshold (adjust)
cellchat_filtered <- cellchat
cellchat_filtered@net$prob[cellchat_filtered@net$prob < thresh_prob] <- 0

Recompute after masking low-prob edges
cellchat_filtered <- computeCommunProbPathway(cellchat_filtered, thresh = 1)
cellchat_filtered <- aggregateNet(cellchat_filtered, thresh = 1)

11) Mask LR edges by var.features (allowed clusters per gene)
mask_net_by_var_features <- function(cellchat, recompute_pathway = TRUE, verbose = TRUE) {
stopifnot(!is.null(cellchat@net$prob), !is.null(cellchat@net$pval))
prob <- cellchat@net$prob
pval <- cellchat@net$pval

senders <- rownames(prob)
receivers <- colnames(prob)
lrs <- dimnames(prob)[[3]]

vf <- cellchat@var.features
genes_keep <- toupper(vf$features)
fi <- vf$features.info
fi$features <- toupper(fi$features)
allowed_map <- split(as.character(fi$clusters), fi$features)

dfLR <- cellchat@LR$LRsig
need_cols <- c("ligand","receptor")
if (!all(need_cols %in% colnames(dfLR))) stop("LRsig must contain 'ligand' and 'receptor'.")
dfLR$ligand <- toupper(dfLR$ligand)
dfLR$receptor <- toupper(dfLR$receptor)
dfLR$lrkey <- paste(dfLR$ligand, gsub("\+", "", dfLR$receptor), sep = "")

m <- match(lrs, dfLR$lrkey)
if (verbose && any(is.na(m))) {
message("Warning: ", sum(is.na(m)), " LR names not matched in LRsig. Examples: ",
paste(head(lrs[is.na(m)]), collapse = ", "))
}
null2chr0 <- function(x) if (is.null(x)) character(0) else x

for (k in seq_along(lrs)) {
idx <- m[k]
if (is.na(idx)) next
lig_genes <- toupper(trimws(unlist(strsplit(dfLR$ligand[idx], "\+"))))
rec_genes <- toupper(trimws(unlist(strsplit(dfLR$receptor[idx], "\+"))))
all_genes <- unique(c(lig_genes, rec_genes))

1) LR containing genes not in var.features -> zero slice
if (!all(all_genes %in% genes_keep)) {
prob[, , k] <- 0; pval[, , k] <- 1; next
}

2) Allowed clusters intersection per side
allowed_send_sets <- lapply(lig_genes, function(g) null2chr0(allowed_map[[g]]))
allowed_recv_sets <- lapply(rec_genes, function(g) null2chr0(allowed_map[[g]]))
if (any(lengths(allowed_send_sets) == 0) || any(lengths(allowed_recv_sets) == 0)) {
prob[, , k] <- 0; pval[, , k] <- 1; next
}
allowed_send <- Reduce(intersect, allowed_send_sets)
allowed_recv <- Reduce(intersect, allowed_recv_sets)
if (length(allowed_send) == 0 || length(allowed_recv) == 0) {
prob[, , k] <- 0; pval[, , k] <- 1; next
}
not_send <- setdiff(senders, allowed_send)
not_recv <- setdiff(receivers, allowed_recv)
if (length(not_send) > 0) {
ii <- match(not_send, senders); prob[ii, , k] <- 0; pval[ii, , k] <- 1
}
if (length(not_recv) > 0) {
jj <- match(not_recv, receivers); prob[, jj, k] <- 0; pval[, jj, k] <- 1
}
}
cellchat@net$prob <- prob
cellchat@net$pval <- pval
cellchat@net$count <- apply(cellchat@net$pval <= 0.05, c(1, 2), sum)
cellchat@net$weight <- apply(cellchat@net$prob, c(1, 2), sum)

if (recompute_pathway) {
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
}
cellchat
}

cellchat_filtered2 <- mask_net_by_var_features(cellchat_filtered, recompute_pathway = TRUE, verbose = TRUE)

12) Visualizations (chemokine and TGFb)
available_pathways
netVisual_chord_gene(
cellchat_filtered,
targets.use = "Bcell",
scale = FALSE, lab.cex = 1,
legend.pos.x = 0, legend.pos.y = 50,
link.border = TRUE
)

netVisual_chord_gene(
cellchat_filtered2,
signaling = c("CXCL", "CCL", "TGFb"),

sources.use = c("Immunofibro", "FRC-like", "Fibroblast_CXCL12", "EC_AQP1"), # optional
targets.use = c("IgA_PC", "IgG_PC", "Bcell"),
scale = FALSE, lab.cex = 1,
legend.pos.x = 0, legend.pos.y = 50,
link.border = TRUE,
color.use = cluster_colors2
)

netVisual_chord_cell(
cellchat_filtered,
signaling = c("CXCL","CCL"),
color.use = cluster_colors,
link.border = TRUE
)
