library(FNN)
library(ggplot2)
library(dplyr)

=== 1. Prepare coordinate data (all FOVs, consolidated) ===
Collect centroid coordinates for all FOVs and apply offsets to avoid overlap
fov_names <- c("fov", paste0("fov.", 2:11))
offsets <- seq(10000, 110000, by = 10000)

centroid_list <- list()

for (i in seq_along(fov_names)) {
fov_name <- fov_names[i]
offset <- offsets[i]
fov_centroids <- s.sub[[fov_name]]$centroid
coords <- as.data.frame(fov_centroids@coords) + offset
rownames(coords) <- fov_centroids@cells
centroid_list[[fov_name]] <- coords
}

Combine all FOV coordinates
centroid_coords <- do.call(rbind, centroid_list)

Build cell-level metadata
cell_metadata <- data.frame(
cell_id = rownames(centroid_coords),
cluster = s.sub$clusters[rownames(centroid_coords)],
cluster_name = s.sub$cluster_name[rownames(centroid_coords)],
x = centroid_coords$x,
y = centroid_coords$y,
stringsAsFactors = FALSE
)

=== 2. Density scoring functions (realistic scale) ===
calculate_realistic_density_score <- function(cell_coords, radius = 50) {

Compute local neighbor density within a radius; normalize to a 0–5 scale
n_cells <- nrow(cell_coords)
if (n_cells < 3) {
return(list(density_score = 0, mean_neighbors = 0, theoretical_max = NA))
}
neighbor_counts <- numeric(n_cells)
for (i in 1:n_cells) {
distances <- sqrt((cell_coords$x - cell_coords$x[i])^2 +
(cell_coords$y - cell_coords$y[i])^2)
neighbor_counts[i] <- sum(distances <= radius & distances > 0)
}
mean_neighbors <- mean(neighbor_counts)

Normalize based on empirical expectation (tuned so 10–15 neighbors ~ score ~ 3)
theoretical_max <- 5
density_score <- mean_neighbors / theoretical_max
density_score <- min(density_score, 5)

return(list(
density_score = density_score,
mean_neighbors = mean_neighbors,
theoretical_max = theoretical_max
))
}

calculate_realistic_tls_score <- function(cell_data, tls_name, radius = 50) {

Score TLS based on B- and T-cell densities (realistic scale)
tls_cells <- cell_data[cell_data$cluster_name == tls_name, ]

total_cells <- nrow(tls_cells)
b_cells <- tls_cells[tls_cells$cluster == "Bcell", ]
t_cells <- tls_cells[tls_cells$cluster %in% c("CD4_Th", "Tfh", "Treg"), ]

b_count <- nrow(b_cells)
t_count <- nrow(t_cells)

if (total_cells < 10) {
return(list(
density_score = 0,
maturity_level = "Non_functional",
b_density_score = 0,
t_density_score = 0,
b_count = b_count,
t_count = t_count,
total_count = total_cells,
error = "insufficient_cells"
))
}

A. B-cell density
if (b_count >= 3) {
b_coords <- data.frame(x = b_cells$x, y = b_cells$y)
b_result <- calculate_realistic_density_score(b_coords, radius)
b_density_score <- b_result$density_score
b_mean_neighbors <- b_result$mean_neighbors
} else {
b_density_score <- 0
b_mean_neighbors <- 0
}

B. T-cell density
if (t_count >= 3) {
t_coords <- data.frame(x = t_cells$x, y = t_cells$y)
t_result <- calculate_realistic_density_score(t_coords, radius)
t_density_score <- t_result$density_score
t_mean_neighbors <- t_result$mean_neighbors
} else {
t_density_score <- 0
t_mean_neighbors <- 0
}

C. Overall density-based structure score (0–5 scale)
density_score <- (b_density_score + t_density_score) / 2

D. Maturity classification (realistic thresholds)
min_density <- min(b_density_score, t_density_score)
avg_density <- density_score

if (avg_density >= 1.5) {
maturity_level <- "Highly_Mature"
} else if (avg_density >= 1.0 && min_density >= 0.5) {
maturity_level <- "Mature_TLS"
} else if (avg_density >= 0.6 && min_density >= 0.3) {
maturity_level <- "Intermediate_TLS"
} else if (avg_density >= 0.3) {
maturity_level <- "Early_Stage"
} else {
maturity_level <- "Immature_TLS"
}

return(list(
density_score = density_score,
maturity_level = maturity_level,
b_density_score = b_density_score,
t_density_score = t_density_score,
b_count = b_count,
t_count = t_count,
total_count = total_cells,
error = ""
))
}

=== 3. TLS maturity evaluation (realistic scale) ===
unique_tls <- unique(cell_metadata$cluster_name)
unique_tls <- unique_tls[!is.na(unique_tls) & unique_tls != ""]

cat("Running T/B density-based TLS scoring (realistic scale)...\n")

tls_realistic_results <- list()
for (i in seq_along(unique_tls)) {
tls_name <- unique_tls[i]
cat(sprintf("Processing: %s (%d/%d)\n", tls_name, i, length(unique_tls)))
result <- calculate_realistic_tls_score(cell_metadata, tls_name, radius = 50)
tls_realistic_results[[tls_name]] <- result
}

Collect results
tls_scores_realistic <- sapply(tls_realistic_results, function(x) x$density_score)
tls_maturity_realistic <- sapply(tls_realistic_results, function(x) x$maturity_level)
b_density_scores_realistic <- sapply(tls_realistic_results, function(x) x$b_density_score)
t_density_scores_realistic <- sapply(tls_realistic_results, function(x) x$t_density_score)

tls_results_realistic <- data.frame(
TLS_name = names(tls_scores_realistic),
Structure_Score = tls_scores_realistic,
Maturity_Status = tls_maturity_realistic,
B_Density_Score = b_density_scores_realistic,
T_Density_Score = t_density_scores_realistic,
stringsAsFactors = FALSE
)

cat("=== T/B density-based structure scores (realistic): Highly Mature ≥ 1.5 ===\n")
print(tls_results_realistic[order(-tls_results_realistic$Structure_Score), ])

=== 4. Score distribution summary ===
cat("\n=== Score distribution (realistic) ===\n")
cat(sprintf("Min score: %.3f\n", min(tls_results_realistic$Structure_Score)))
cat(sprintf("Max score: %.3f\n", max(tls_results_realistic$Structure_Score)))
cat(sprintf("Mean score: %.3f\n", mean(tls_results_realistic$Structure_Score)))
cat(sprintf("Median: %.3f\n", median(tls_results_realistic$Structure_Score)))
cat(sprintf("SD: %.3f\n", sd(tls_results_realistic$Structure_Score)))

cat("\n=== Maturity distribution ===\n")
print(table(tls_results_realistic$Maturity_Status))

=== 5. Visualization (realistic scale) ===
p_realistic_hist <- ggplot(tls_results_realistic, aes(x = Structure_Score, fill = Maturity_Status)) +
geom_histogram(bins = 20, alpha = 0.7) +
theme_minimal() +
labs(title = "TLS Realistic Density Score Distribution",
subtitle = "Highly Mature: ≥1.5; Mature: ≥1.0; Intermediate: ≥0.6",
x = "Realistic Density Score", y = "Count") +
geom_vline(xintercept = c(0.6, 1.0, 1.5), linetype = "dashed", alpha = 0.7, color = "red") +
annotate("text", x = 1.5, y = Inf, label = "Highly Mature", vjust = 2, color = "red")

p_realistic_bar <- ggplot(tls_results_realistic, aes(x = reorder(TLS_name, Structure_Score),
y = Structure_Score, fill = Maturity_Status)) +
geom_col() +
coord_flip() +
theme_minimal() +
labs(title = "TLS Maturity by Realistic Density Scoring",
x = "TLS Name", y = "Realistic Density Score") +
theme(axis.text.y = element_text(size = 8)) +
geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", alpha = 0.7)

p_realistic_scatter <- ggplot(tls_results_realistic, aes(x = B_Density_Score, y = T_Density_Score,
color = Maturity_Status, size = Structure_Score)) +
geom_point(alpha = 0.8) +
geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", alpha = 0.5) +
geom_vline(xintercept = 1.5, linetype = "dashed", color = "red", alpha = 0.5) +
labs(title = "B-cell vs T-cell Density (Realistic Scale)",
x = "B-cell Density", y = "T-cell Density") +
theme_minimal()

print(p_realistic_hist)
print(p_realistic_bar)
print(p_realistic_scatter)

=== 6. FDC/Tfh scoring ===
calculate_fdc_tfh_score <- function(cell_data, tls_name) {
tls_cells <- cell_data[cell_data$cluster_name == tls_name, ]
total_cells <- nrow(tls_cells)

if (total_cells < 10) {
return(list(
fdc_tfh_score = 0,
fdc_count = 0,
tfh_count = 0,
fdc_ratio = 0,
tfh_ratio = 0,
fdc_score = 0,
tfh_score = 0,
error = "insufficient_cells"
))
}

fdc_cells <- tls_cells[tls_cells$cluster == "FDC", ]
tfh_cells <- tls_cells[tls_cells$cluster == "Tfh", ]

fdc_count <- nrow(fdc_cells)
tfh_count <- nrow(tfh_cells)

fdc_ratio <- fdc_count / total_cells
tfh_ratio <- tfh_count / total_cells

0–3 scale to align with density scores
fdc_score <- min(1, fdc_ratio * 20) # Full at ~5%
tfh_score <- min(1, tfh_ratio * 10) # Full at ~10%

fdc_tfh_score <- (fdc_score + tfh_score) / 2 * 3

return(list(
fdc_tfh_score = fdc_tfh_score,
fdc_count = fdc_count,
tfh_count = tfh_count,
fdc_ratio = fdc_ratio,
tfh_ratio = tfh_ratio,
fdc_score = fdc_score,
tfh_score = tfh_score,
error = ""
))
}

=== 7. Run FDC/Tfh scoring for all TLS ===
cat("Running FDC/Tfh scoring...\n")

fdc_tfh_realistic_results <- list()
for (i in seq_along(unique_tls)) {
tls_name <- unique_tls[i]
cat(sprintf("FDC/Tfh processing: %s (%d/%d)\n", tls_name, i, length(unique_tls)))
result <- calculate_fdc_tfh_score(cell_metadata, tls_name)
fdc_tfh_realistic_results[[tls_name]] <- result
}

fdc_tfh_scores_realistic <- sapply(fdc_tfh_realistic_results, function(x) x$fdc_tfh_score)
fdc_counts_realistic <- sapply(fdc_tfh_realistic_results, function(x) x$fdc_count)
tfh_counts_realistic <- sapply(fdc_tfh_realistic_results, function(x) x$tfh_count)
fdc_ratios_realistic <- sapply(fdc_tfh_realistic_results, function(x) x$fdc_ratio)
tfh_ratios_realistic <- sapply(fdc_tfh_realistic_results, function(x) x$tfh_ratio)

tls_results_realistic$FDC_Tfh_Score <- fdc_tfh_scores_realistic[tls_results_realistic$TLS_name]
tls_results_realistic$FDC_Count <- fdc_counts_realistic[tls_results_realistic$TLS_name]
tls_results_realistic$Tfh_Count <- tfh_counts_realistic[tls_results_realistic$TLS_name]
tls_results_realistic$FDC_Ratio <- fdc_ratios_realistic[tls_results_realistic$TLS_name]
tls_results_realistic$Tfh_Ratio <- tfh_ratios_realistic[tls_results_realistic$TLS_name]

cat("=== FDC/Tfh score results (realistic) ===\n")
print(tls_results_realistic[, c("TLS_name", "Structure_Score", "FDC_Tfh_Score",
"FDC_Count", "Tfh_Count", "FDC_Ratio", "Tfh_Ratio")])

=== 8. Combined score (realistic) ===
tls_results_realistic$Combined_Score <- (tls_results_realistic$Structure_Score +
tls_results_realistic$FDC_Tfh_Score) / 2

Combined functional classification (consistent labels)
tls_results_realistic$Combined_Maturity <- ifelse(
tls_results_realistic$Structure_Score >= 3 & tls_results_realistic$FDC_Tfh_Score >= 2, "Highly_Functional",
ifelse(tls_results_realistic$Structure_Score >= 2 & tls_results_realistic$FDC_Tfh_Score >= 1.2, "Functional",
ifelse(tls_results_realistic$Structure_Score >= 1 | tls_results_realistic$FDC_Tfh_Score >= 0.8, "Developing", "Not_Functional"))
)

cat("=== Combined score results (realistic) ===\n")
print(tls_results_realistic[order(-tls_results_realistic$Combined_Score),
c("TLS_name", "Structure_Score", "FDC_Tfh_Score",
"Combined_Score", "Combined_Maturity")])

=== 9. Map cluster_name to TLS ID (1–6) ===
tls_mapping <- data.frame(
cluster_name = character(),
tls_id = integer(),
stringsAsFactors = FALSE
)

correspondence_table <- table(s.sub$TLS, s.sub$cluster_name)
for (tls_id in 1:6) {
for (cluster_name in colnames(correspondence_table)) {
if (correspondence_table[as.character(tls_id), cluster_name] > 0) {
tls_mapping <- rbind(tls_mapping, data.frame(
cluster_name = cluster_name,
tls_id = tls_id,
stringsAsFactors = FALSE
))
}
}
}

cat("=== Mapping: cluster_name to TLS ID ===\n")
print(tls_mapping)

=== 10. Filter to TLS 1–6 for realistic analysis ===
tls_results_with_id_realistic <- merge(tls_results_realistic, tls_mapping,
by.x = "TLS_name", by.y = "cluster_name",
all.x = TRUE)

tls_results_filtered_realistic <- tls_results_with_id_realistic[
!is.na(tls_results_with_id_realistic$tls_id) &
tls_results_with_id_realistic$tls_id %in% 1:6, ]

cat("=== TLS 1–6 results (realistic) ===\n")
print(tls_results_filtered_realistic[order(tls_results_filtered_realistic$tls_id),
c("tls_id", "TLS_name", "Structure_Score",
"FDC_Tfh_Score", "Combined_Score", "Combined_Maturity")])

=== 11. Group statistics by TLS ID (realistic) ===
tls_group_stats_realistic <- tls_results_filtered_realistic %>%
group_by(tls_id) %>%
summarise(
n_clusters = n(),
mean_density_score = mean(Structure_Score),
mean_fdc_tfh_score = mean(FDC_Tfh_Score),
mean_combined_score = mean(Combined_Score),
sd_density_score = sd(Structure_Score),
min_density_score = min(Structure_Score),
max_density_score = max(Structure_Score),
highly_mature_count = sum(Maturity_Status == "Highly_Mature"),
mature_count = sum(Maturity_Status == "Mature_TLS"),
intermediate_count = sum(Maturity_Status == "Intermediate_TLS"),
highly_functional_count = sum(Combined_Maturity == "Highly_Functional"),
functional_count = sum(Combined_Maturity == "Functional"),
.groups = 'drop'
)

cat("=== TLS 1–6 group summary (realistic) ===\n")
print(tls_group_stats_realistic)

=== 12. Statistical tests (realistic) ===
if (length(unique(tls_results_filtered_realistic$tls_id)) > 2) {
anova_realistic <- aov(Structure_Score ~ factor(tls_id), data = tls_results_filtered_realistic)
cat("\n=== ANOVA: Structure score across TLS IDs (realistic) ===\n")
print(summary(anova_realistic))

if (summary(anova_realistic)[[1]][["Pr(>F)"]][1] < 0.1) {
tukey_realistic <- TukeyHSD(anova_realistic)
cat("\n=== Tukey HSD (realistic) ===\n")
print(tukey_realistic)
}

anova_fdc_tfh_realistic <- aov(FDC_Tfh_Score ~ factor(tls_id), data = tls_results_filtered_realistic)
cat("\n=== ANOVA: FDC/Tfh score across TLS IDs (realistic) ===\n")
print(summary(anova_fdc_tfh_realistic))

anova_combined_realistic <- aov(Combined_Score ~ factor(tls_id), data = tls_results_filtered_realistic)
cat("\n=== ANOVA: Combined score across TLS IDs (realistic) ===\n")
print(summary(anova_combined_realistic))

cor_realistic <- cor.test(tls_results_filtered_realistic$Structure_Score,
tls_results_filtered_realistic$FDC_Tfh_Score)
cat("\n=== Correlation: Structure vs FDC/Tfh (realistic) ===\n")
print(cor_realistic)
}

=== 13. Visualization: TLS 1–6 comparison (realistic) ===
if (require(gridExtra, quietly = TRUE)) {
use_gridarrange <- TRUE
} else {
use_gridarrange <- FALSE
cat("Note: gridExtra is not available. Plots will be printed individually.\n")
}

p_realistic_boxplot <- ggplot(tls_results_filtered_realistic,
aes(x = factor(tls_id), y = Structure_Score, fill = factor(tls_id))) +
geom_boxplot() +
geom_jitter(width = 0.2, alpha = 0.7) +
geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", alpha = 0.7) +
theme_minimal() +
labs(title = "Realistic Density Scores by TLS ID",
subtitle = "Red line: Highly Mature threshold (1.5)",
x = "TLS ID", y = "Realistic Density Score") +
theme(legend.position = "none")

p_realistic_boxplot2 <- ggplot(tls_results_filtered_realistic,
aes(x = factor(tls_id), y = FDC_Tfh_Score, fill = factor(tls_id))) +
geom_boxplot() +
geom_jitter(width = 0.2, alpha = 0.7) +
geom_hline(yintercept = 2, linetype = "dashed", color = "red", alpha = 0.7) +
theme_minimal() +
labs(title = "FDC/Tfh Scores by TLS ID",
subtitle = "Red line: Highly Functional threshold (2.0)",
x = "TLS ID", y = "FDC/Tfh Score") +
theme(legend.position = "none")

p_realistic_comparison <- ggplot(tls_results_filtered_realistic,
aes(x = Structure_Score, y = FDC_Tfh_Score,
color = factor(tls_id), size = Combined_Score)) +
geom_point(alpha = 0.8) +
geom_hline(yintercept = 2, linetype = "dashed", color = "red", alpha = 0.7) +
geom_vline(xintercept = 3, linetype = "dashed", color = "red", alpha = 0.7) +
labs(title = "Realistic Density vs FDC/Tfh Scores",
subtitle = "Red lines: Highly Functional thresholds",
x = "Realistic Density Score", y = "FDC/Tfh Score",
color = "TLS ID", size = "Combined Score") +
theme_minimal()

print(p_realistic_boxplot)
print(p_realistic_boxplot2)
print(p_realistic_comparison)
