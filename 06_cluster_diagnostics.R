# Extra clustering checks described in the response to reviewer 1

patients <- read.csv("final_analytic_dataset.csv", check.names = FALSE)
proteins <- readLines("mito_full_panel.txt")
stopifnot(all(proteins %in% names(patients)))

z <- scale(patients[, proteins])
z[is.na(z)] <- 0
pcs <- prcomp(z, center = FALSE, scale. = FALSE)$x[, 1:3]

set.seed(42)
fits <- vector("list", 7)
for (k in 2:8) {
  fits[[k - 1]] <- kmeans(pcs, centers = k, nstart = 100)
}
fit1 <- kmeans(pcs, centers = 1, nstart = 100)
results <- data.frame(
  k = 1:8,
  wcss = NA_real_,
  calinski_harabasz = NA_real_,
  davies_bouldin = NA_real_,
  mean_jaccard = NA_real_,
  jaccard_q025 = NA_real_,
  jaccard_q975 = NA_real_)
results$wcss[1] <- fit1$tot.withinss

# Elbow, Calinski-Harabasz and Davies-Bouldin
for (k in 2:8) {
  fit <- fits[[k - 1]]
  results$wcss[k] <- fit$tot.withinss
  results$calinski_harabasz[k] <-
    (fit$betweenss / (k - 1)) / (fit$tot.withinss / (nrow(pcs) - k))

  # Average distance from each patient to the center of their cluster
  radius <- numeric(k)
  for (j in 1:k) {
    cluster_points <- pcs[fit$cluster == j, , drop = FALSE]
    distances <- numeric(nrow(cluster_points))
    for (n in 1:nrow(cluster_points)) {
      distances[n] <- sqrt(sum((cluster_points[n, ] - fit$centers[j, ])^2))
    }
    radius[j] <- mean(distances)
  }

  largest_ratio <- numeric(k)
  for (j in 1:k) {
    ratios <- numeric()
    for (h in 1:k) {
      if (h != j) {
        center_distance <- sqrt(sum((fit$centers[j, ] - fit$centers[h, ])^2))
        ratios <- c(ratios, (radius[j] + radius[h]) / center_distance)
      }
    }
    largest_ratio[j] <- max(ratios)
  }
  results$davies_bouldin[k] <- mean(largest_ratio)
}

# Jaccard: repeat k-means on 80% of the patients, 200 times.
set.seed(42)
cluster_results <- data.frame()
for (k in 2:8) {
  original <- fits[[k - 1]]$cluster
  cluster_jaccard <- matrix(NA_real_, nrow = 200, ncol = k)
  stability <- numeric(200)

  for (b in 1:200) {
    selected <- sample.int(nrow(pcs), floor(0.8 * nrow(pcs)))
    refit <- kmeans(pcs[selected, , drop = FALSE], centers = k, nstart = 100)
    for (j in 1:k) {
      old_group <- original[selected] == j
      best <- 0
      for (h in 1:k) {
        new_group <- refit$cluster == h
        union <- sum(old_group | new_group)
        if (union > 0) {
          overlap <- sum(old_group & new_group) / union
          best <- max(best, overlap)
        }
      }
      cluster_jaccard[b, j] <- best
    }
    stability[b] <- mean(cluster_jaccard[b, ])
  }

  results$mean_jaccard[k] <- mean(stability)
  results$jaccard_q025[k] <- quantile(stability, 0.025)
  results$jaccard_q975[k] <- quantile(stability, 0.975)
  for (j in 1:k) {
    cluster_results <- rbind(cluster_results, data.frame(
      k = k, cluster = j, size = sum(original == j),
      mean_jaccard = mean(cluster_jaccard[, j])))
  }
}

write.csv(results, "cluster_diagnostics.csv", row.names = FALSE, na = "")
write.csv(cluster_results, "cluster_jaccard_by_cluster.csv", row.names = FALSE)
print(results, row.names = FALSE, digits = 4)

png("Figure_cluster_diagnostics.png", width = 1600, height = 1200, res = 180)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(results$k, results$wcss, type = "b", xlab = "Number of clusters",
     ylab = "Within-cluster sum of squares", main = "Elbow")
plot(2:8, results$calinski_harabasz[2:8], type = "b",
     xlab = "Number of clusters", ylab = "Calinski-Harabasz index")
plot(2:8, results$davies_bouldin[2:8], type = "b",
     xlab = "Number of clusters", ylab = "Davies-Bouldin index")
plot(2:8, results$mean_jaccard[2:8], type = "b", ylim = c(0, 1),
     xlab = "Number of clusters", ylab = "Mean Jaccard index")
dev.off()
