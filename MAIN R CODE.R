# Load necessary libraries
library(markovchain)
library(ggplot2)
library(gridExtra)

# Function to create random walk Markov chain
create_random_walk_chain <- function(n) {
  states <- as.character(1:n)
  P <- matrix(0, nrow = n, ncol = n)
  diag(P) <- 0.5
  for (i in 1:(n-1)) {
    P[i, i + 1] <- 0.5
    P[i + 1, i] <- 0.5
  }
  P[1, 2] <- 1
  P[n, n - 1] <- 1
  P <- P / rowSums(P) # Normalize rows to sum to 1
  return(new("markovchain", states = states, transitionMatrix = P))
}

# Function to reduce state space using simple clustering
reduce_state_space <- function(mc, k) {
  clusters <- cutree(hclust(dist(mc@transitionMatrix)), k = k)
  states <- unique(clusters)
  P_red <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      P_red[i, j] <- mean(mc@transitionMatrix[clusters == i, clusters == j])
    }
  }
  P_red <- P_red / rowSums(P_red) # Normalize rows to sum to 1
  return(list(mc = new("markovchain", states = as.character(states), transitionMatrix = P_red), clusters = clusters))
}

# Function to simulate a Markov chain and compute steady-state distribution
compute_steady_state <- function(mc) {
  steady <- steadyStates(mc)
  if (is.null(steady) || nrow(steady) == 0) {
    stop("Failed to compute steady-state distribution.")
  }
  return(as.numeric(steady))
}

# Plot transition matrices
plot_transition_matrix <- function(mc, title) {
  df <- as.data.frame(as.table(mc@transitionMatrix))
  colnames(df) <- c("From", "To", "Probability")
  p <- ggplot(df, aes(x = From, y = To, fill = Probability)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal(base_size = 10) +
    labs(title = title, x = "From State", y = "To State")
  return(p)
}

# Plot clusters overlaying the transition matrix
plot_clusters_on_transition_matrix <- function(mc, clusters, title) {
  df <- as.data.frame(as.table(mc@transitionMatrix))
  colnames(df) <- c("From", "To", "Probability")
  df$ClusterFrom <- as.factor(clusters[as.numeric(df$From)])
  df$ClusterTo <- as.factor(clusters[as.numeric(df$To)])
  p <- ggplot(df, aes(x = From, y = To, fill = Probability, color = ClusterFrom)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    scale_color_manual(values = rainbow(length(unique(clusters))), guide = "none") +
    theme_minimal(base_size = 10) +
    labs(title = title, x = "From State", y = "To State") +
    theme(legend.position = "bottom", legend.title = element_blank())
  return(p)
}

# Plot clusters alone
plot_clusters <- function(clusters, title) {
  df <- data.frame(State = 1:length(clusters), Cluster = as.factor(clusters))
  p <- ggplot(df, aes(x = State, y = Cluster, color = Cluster)) +
    geom_point(size = 4) +
    scale_color_manual(values = rainbow(length(unique(clusters)))) +
    theme_minimal(base_size = 10) +
    labs(title = title, x = "State", y = "Cluster") +
    theme(legend.position = "bottom", legend.title = element_blank())
  return(p)
}

# Plot comparison of original vs reduced state probabilities
plot_state_probabilities_comparison <- function(original_probs, reduced_probs, clusters) {
  original_df <- data.frame(State = 1:length(original_probs), Probability = original_probs, Type = "Original", Cluster = NA)
  
  reduced_probs_expanded <- sapply(1:length(clusters), function(x) reduced_probs[clusters[x]])
  reduced_df <- data.frame(State = 1:length(clusters), Probability = reduced_probs_expanded, Type = "Reduced", Cluster = as.factor(clusters))
  
  combined_df <- rbind(original_df, reduced_df)
  
  p <- ggplot() +
    geom_bar(data = original_df, aes(x = State, y = Probability, fill = Type), stat = "identity", color = "black", fill = "grey", linetype = "dashed", width = 0.6) +
    geom_bar(data = reduced_df, aes(x = State, y = Probability, fill = Cluster), stat = "identity", alpha = 0.5, width = 0.6) +
    scale_fill_manual(values = c("grey", rainbow(length(unique(clusters)))), guide = "legend") +
    theme_minimal(base_size = 10) +
    labs(title = "Comparison of Original vs Reduced State Probabilities", x = "State", y = "Probability", fill = "Model") +
    theme(legend.position = "bottom", legend.title = element_blank())
  return(p)
}

# Main simulation function
run_simulation <- function() {
  set.seed(123)
  
  # Create original random walk Markov chain
  n <- 20
  original_mc <- create_random_walk_chain(n)
  original_steady_state <- compute_steady_state(original_mc)
  
  # Reduce state space
  k <- 6
  reduced_result <- reduce_state_space(original_mc, k)
  reduced_mc <- reduced_result$mc
  clusters <- reduced_result$clusters
  reduced_steady_state <- compute_steady_state(reduced_mc)
  
  # Plot transition matrices
  p1 <- plot_transition_matrix(original_mc, "Original Transition Matrix")
  p2 <- plot_transition_matrix(reduced_mc, "Reduced Transition Matrix")
  
  # Plot steady-state distributions
  df1 <- data.frame(State = 1:length(original_steady_state), Probability = original_steady_state)
  
  cluster_assignment <- rep(NA, length(original_steady_state))
  for (i in 1:k) {
    cluster_assignment[clusters == i] <- i
  }
  df2 <- data.frame(State = 1:length(original_steady_state), Probability = reduced_steady_state[cluster_assignment], Cluster = as.factor(cluster_assignment))
  
  p3 <- ggplot(df1, aes(x = State, y = Probability)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 10) +
    labs(title = "Original Steady-State Distribution", x = "State", y = "Probability")
  
  p4 <- ggplot(df2, aes(x = State, y = Probability, fill = Cluster)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = rainbow(length(unique(clusters)))) +
    theme_minimal(base_size = 10) +
    labs(title = "Reduced Steady-State Distribution", x = "State", y = "Probability", fill = "Cluster")
  
  # Overlay clusters on transition matrix
  p5 <- plot_clusters_on_transition_matrix(original_mc, clusters, "Cluster Assignment on Original Transition Matrix")
  
  # Plot clusters alone
  p6 <- plot_clusters(clusters, "Cluster Assignment of States")
  
  # Plot comparison of original vs reduced state probabilities
  p7 <- plot_state_probabilities_comparison(original_steady_state, reduced_steady_state, clusters)
  
  # Save plots
  ggsave("original_transition_matrix.png", plot = p1, width = 8, height = 6)
  ggsave("reduced_transition_matrix.png", plot = p2, width = 8, height = 6)
  ggsave("original_steady_state_distribution.png", plot = p3, width = 8, height = 6)
  ggsave("reduced_steady_state_distribution.png", plot = p4, width = 8, height = 6)
  ggsave("clusters_on_original_transition_matrix.png", plot = p5, width = 8, height = 6)
  ggsave("cluster_assignment.png", plot = p6, width = 8, height = 6)
  ggsave("state_probabilities_comparison.png", plot = p7, width = 8, height = 6)
  
  # Display plots two at a time with relevancy
  grid.arrange(p1, p2, ncol = 2)
  grid.arrange(p3, p4, ncol = 2)
  grid.arrange(p5, p6, ncol = 2) # Cluster assignment and original transition matrix
  grid.arrange(p7, ncol = 1) # Comparison of original vs reduced state probabilities
}

# Run the simulation
run_simulation()
###########################

# Load necessary libraries
library(MASS)
library(markovchain)
library(ggplot2)
library(gridExtra)

# Function to create a transition matrix from the dataset
create_transition_matrix <- function(data, states) {
  n <- length(states)
  transition_matrix <- matrix(0, n, n)
  rownames(transition_matrix) <- states
  colnames(transition_matrix) <- states
  
  for (i in 1:(nrow(data)-1)) {
    from_state <- data[i]
    to_state <- data[i+1]
    transition_matrix[from_state, to_state] <- transition_matrix[from_state, to_state] + 1
  }
  
  transition_matrix <- transition_matrix / rowSums(transition_matrix)
  return(transition_matrix)
}

# Function to create a Markov chain from the transition matrix
create_markov_chain <- function(transition_matrix) {
  states <- rownames(transition_matrix)
  return(new("markovchain", states = states, transitionMatrix = transition_matrix))
}

# Function to reduce state space using simple clustering
reduce_state_space <- function(mc, k) {
  clusters <- cutree(hclust(dist(mc@transitionMatrix)), k = k)
  states <- unique(clusters)
  P_red <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      P_red[i, j] <- mean(mc@transitionMatrix[clusters == i, clusters == j])
    }
  }
  P_red <- P_red / rowSums(P_red) # Normalize rows to sum to 1
  return(list(mc = new("markovchain", states = as.character(states), transitionMatrix = P_red), clusters = clusters))
}

# Function to simulate a Markov chain and compute steady-state distribution
compute_steady_state <- function(mc) {
  steady <- steadyStates(mc)
  if (is.null(steady) || nrow(steady) == 0) {
    stop("Failed to compute steady-state distribution.")
  }
  return(as.numeric(steady))
}

# Function to compute total variation distance between two distributions
total_variation_distance <- function(p, q) {
  return(sum(abs(p - q)) / 2)
}

# Main simulation function with performance metrics
run_simulation <- function() {
  set.seed(123)
  
  # Load the dataset
  data("birthwt")
  birthwt <- birthwt[order(birthwt$id), ] # Ensure the data is ordered by id if needed
  
  # Create a state column combining race and smoke
  birthwt$state <- with(birthwt, paste(race, smoke, sep = "_"))
  states <- unique(birthwt$state)
  
  # Create original transition matrix and Markov chain
  original_transition_matrix <- create_transition_matrix(birthwt$state, states)
  original_mc <- create_markov_chain(original_transition_matrix)
  
  # Compute steady-state distribution for original Markov chain
  start_time <- Sys.time()
  original_steady_state <- compute_steady_state(original_mc)
  end_time <- Sys.time()
  original_computation_time <- end_time - start_time
  
  # Reduce state space
  k <- 6
  reduced_result <- reduce_state_space(original_mc, k)
  reduced_mc <- reduced_result$mc
  clusters <- reduced_result$clusters
  
  # Compute steady-state distribution for reduced Markov chain
  start_time <- Sys.time()
  reduced_steady_state <- compute_steady_state(reduced_mc)
  end_time <- Sys.time()
  reduced_computation_time <- end_time - start_time
  
  # Compute total variation distance between original and reduced steady-state distributions
  cluster_assignment <- rep(NA, length(original_steady_state))
  for (i in 1:k) {
    cluster_assignment[clusters == i] <- i
  }
  reduced_steady_state_expanded <- reduced_steady_state[cluster_assignment]
  steady_state_distance <- total_variation_distance(original_steady_state, reduced_steady_state_expanded)
  
  # Plot transition matrices
  p1 <- plot_transition_matrix(original_mc, "Original Transition Matrix")
  p2 <- plot_transition_matrix(reduced_mc, "Reduced Transition Matrix")
  
  # Plot steady-state distributions
  df1 <- data.frame(State = 1:length(original_steady_state), Probability = original_steady_state)
  df2 <- data.frame(State = 1:length(original_steady_state), Probability = reduced_steady_state_expanded, Cluster = as.factor(cluster_assignment))
  
  p3 <- ggplot(df1, aes(x = State, y = Probability)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 10) +
    labs(title = "Original Steady-State Distribution", x = "State", y = "Probability")
  
  p4 <- ggplot(df2, aes(x = State, y = Probability, fill = Cluster)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = rainbow(length(unique(clusters)))) +
    theme_minimal(base_size = 10) +
    labs(title = "Reduced Steady-State Distribution", x = "State", y = "Probability", fill = "Cluster")
  
  # Overlay clusters on transition matrix
  p5 <- plot_clusters_on_transition_matrix(original_mc, clusters, "Cluster Assignment on Original Transition Matrix")
  
  # Plot clusters alone
  p6 <- plot_clusters(clusters, "Cluster Assignment of States")
  
  # Plot comparison of original vs reduced state probabilities
  p7 <- plot_state_probabilities_comparison(original_steady_state, reduced_steady_state, clusters)
  
  # Save plots
  ggsave("original_transition_matrix.png", plot = p1, width = 8, height = 6)
  ggsave("reduced_transition_matrix.png", plot = p2, width = 8, height = 6)
  ggsave("original_steady_state_distribution.png", plot = p3, width = 8, height = 6)
  ggsave("reduced_steady_state_distribution.png", plot = p4, width = 8, height = 6)
  ggsave("clusters_on_original_transition_matrix.png", plot = p5, width = 8, height = 6)
  ggsave("cluster_assignment.png", plot = p6, width = 8, height = 6)
  ggsave("state_probabilities_comparison.png", plot = p7, width = 8, height = 6)
  
  # Display plots two at a time with relevancy
  grid.arrange(p1, p2, ncol = 2)
  grid.arrange(p3, p4, ncol = 2)
  grid.arrange(p5, p6, ncol = 2) # Cluster assignment and original transition matrix
  grid.arrange(p7, ncol = 1) # Comparison of original vs reduced state probabilities
  
  # Return performance metrics
  return(list(
    original_computation_time = original_computation_time,
    reduced_computation_time = reduced_computation_time,
    steady_state_distance = steady_state_distance
  ))
}

# Run the simulation and get performance metrics
performance_metrics <- run_simulation()
print(performance_metrics)
