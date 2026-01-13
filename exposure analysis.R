0) Pre-requisites: make sure Xenium assay is active and normalized
DefaultAssay(xenium.obj) <- "Xenium"
xenium.obj <- NormalizeData(xenium.obj)

1) Safe coordinate extractor (auto-detect x/y columns in meta.data)
get_coords <- function(obj, cells = NULL) {
md <- obj@meta.data
if (!is.null(cells)) {
cells <- intersect(cells, rownames(md))
md <- md[cells, , drop = FALSE]
}
cn <- colnames(md)
cand <- list(
c("x","y"), c("X","Y"),
c("Xenium_x","Xenium_y"),
c("spatial_x","spatial_y"),
c("pxl_col","pxl_row")
)
vars <- NULL
for (v in cand) if (all(v %in% cn)) { vars <- v; break }
if (is.null(vars)) stop("No coordinate columns found in meta.data. Please add numeric 'x' and 'y'.")
dat <- md[, vars, drop = FALSE]
coords <- as.matrix(dat); storage.mode(coords) <- "double"
rownames(coords) <- rownames(md); colnames(coords) <- c("x","y")
coords
}

2) Feature getter (returns numeric vector for a gene; NA if absent)
get_feature_vec <- function(obj, gene, cells = NULL) {
a <- DefaultAssay(obj)
m <- GetAssayData(obj, assay = a, slot = "data")
if (!gene %in% rownames(m)) return(setNames(rep(NA_real_, length(cells %||% colnames(m))), cells %||% colnames(m)))
v <- as.numeric(m[gene, cells %||% colnames(m)])
names(v) <- cells %||% colnames(m)
v
}

3) Gaussian-kernel exposure (receiver-to-sender weighted sum within radius R)
compute_exposure_gauss_knn <- function(obj, feature, recv_cells, send_cells,
R = 100, sigma = 20, k = 50, filter_zero = TRUE) {
coords <- get_coords(obj)
recv <- intersect(recv_cells, rownames(coords))
send <- intersect(send_cells, rownames(coords))
out <- setNames(rep(NA_real_, length(recv)), recv)
if (length(recv) == 0 || length(send) == 0) return(out)
expr_send <- get_feature_vec(obj, feature, cells = send)
expr_send[!is.finite(expr_send)] <- 0
if (filter_zero) {
keep <- which(expr_send > 0)
send <- send[keep]
expr_send <- expr_send[keep]
}
if (length(send) == 0) { out[] <- 0; return(out) }
Xs <- coords[send, , drop = FALSE]
Xr <- coords[recv, , drop = FALSE]
k_use <- max(1, min(k, nrow(Xs)))
nn <- FNN::get.knnx(data = Xs, query = Xr, k = k_use)
dist_mat <- nn$nn.dist
idx_mat <- nn$nn.index
W <- exp(-(dist_mat^2) / (2 * sigma^2)) # Gaussian weights
M <- (dist_mat <= R) # within radius mask
expr_mat <- expr_send[idx_mat]; dim(expr_mat) <- dim(dist_mat)
expo <- rowSums(W * M * expr_mat, na.rm = TRUE)
names(expo) <- recv
out[names(expo)] <- expo
out
}

4) Define sender/receiver sets
cells_all <- colnames(xenium.obj)

Receiver: IFNGR1+ if available, else all cells
ifngr1_vec <- get_feature_vec(xenium.obj, "IFNGR1", cells_all)
cells_receiver <- names(ifngr1_vec)[which(ifngr1_vec > 0)]
if (length(cells_receiver) == 0) cells_receiver <- cells_all

Sender: IFNG+ cells (generic; you can restrict to T/NK clusters if desired)
ifng_vec <- get_feature_vec(xenium.obj, "IFNG", cells_all)
cells_sender <- names(ifng_vec)[which(ifng_vec > 0)]
if (length(cells_sender) == 0) {
warning("No IFNG+ sender cells found; exposure will be zeros.")
cells_sender <- character(0)
}

5) Parameters (units must match your coordinates, e.g., pixels or micrometers)
R_expo <- 100 # interaction radius
sigma_exp <- 20 # Gaussian bandwidth
k_nn <- 50 # max neighbors per receiver

6) Compute per-sample to avoid cross-image distances (fallback: all together)
samples <- if ("sample" %in% colnames(xenium.obj@meta.data)) unique(xenium.obj$sample) else NA

Initialize output columns
xenium.obj$IFNG_exposure <- NA_real_
xenium.obj$IFNG_exposure_SCALED <- NA_real_

scale_minmax <- function(x) {
y <- rep(NA_real_, length(x))
ix <- is.finite(x)
if (sum(ix) >= 1) {
xmin <- min(x[ix]); xmax <- max(x[ix])
y[ix] <- if (xmax > xmin) (x[ix] - xmin) / (xmax - xmin) else 0
}
y
}

for (s in samples) {
recv_s <- cells_receiver
send_s <- cells_sender
if (!is.na(s)) {
recv_s <- intersect(recv_s, colnames(xenium.obj)[xenium.obj$sample == s])
send_s <- intersect(send_s, colnames(xenium.obj)[xenium.obj$sample == s])
}
if (length(recv_s) == 0) next
e <- compute_exposure_gauss_knn(
obj = xenium.obj, feature = "IFNG",
recv_cells = recv_s, send_cells = send_s,
R = R_expo, sigma = sigma_exp, k = k_nn, filter_zero = TRUE
)
xenium.obj@meta.data[names(e), "IFNG_exposure"] <- as.numeric(e)

Min–Max scaling within receivers of this sample
xv <- xenium.obj@meta.data[recv_s, "IFNG_exposure"]
xenium.obj@meta.data[recv_s, "IFNG_exposure_SCALED"] <- scale_minmax(as.numeric(xv))
}

7) Quick summary
cat(sprintf("Computed IFNG_exposure for %d receiver cells (samples: %s)\n",
sum(is.finite(xenium.obj$IFNG_exposure)),
paste(na.omit(samples), collapse = ", ")))
