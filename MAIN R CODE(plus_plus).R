# Load necessary libraries
library(MASS)
library(markovchain)
library(ggplot2)
library(gridExtra)

# Function to create a transition matrix from the dataset
create_transition_matrix <- function(data) {
  states <- unique(data)
  n <- length(states)
  transition_matrix <- matrix(0, n, n)
  rownames(transition_matrix) <- states
  colnames(transition_matrix) <- states
  
  for (i in 1:(length(data)-1)) {
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
  
  # Create a state column combining race and smoke
  birthwt$state <- with(birthwt, paste(race, smoke, sep = "_"))
  states <- unique(birthwt$state)
  
  # Create original transition matrix and Markov chain
  transition_matrix <- create_transition_matrix(birthwt$state)
  original_mc <- create_markov_chain(transition_matrix)
  
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
  
  # Plot transition matrix
  p1 <- plot_transition_matrix(reduced_mc, "Reduced Transition Matrix")
  
  # Plot steady-state distribution
  df <- data.frame(State = 1:length(reduced_steady_state), Probability = reduced_steady_state, Cluster = as.factor(clusters))
  
  p2 <- ggplot(df, aes(x = State, y = Probability, fill = Cluster)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = rainbow(length(unique(clusters)))) +
    theme_minimal(base_size = 10) +
    labs(title = "Reduced Steady-State Distribution", x = "State", y = "Probability", fill = "Cluster")
  
  # Overlay clusters on transition matrix
  p3 <- plot_clusters_on_transition_matrix(original_mc, clusters, "Cluster Assignment on Transition Matrix")
  
  # Plot clusters alone
  p4 <- plot_clusters(clusters, "Cluster Assignment of States")
  
  # Save plots
  ggsave("reduced_transition_matrix.png", plot = p1, width = 8, height = 6)
  ggsave("reduced_steady_state_distribution.png", plot = p2, width = 8, height = 6)
  ggsave("clusters_on_transition_matrix.png", plot = p3, width = 8, height = 6)
  ggsave("cluster_assignment.png", plot = p4, width = 8, height = 6)
  
  # Display plots two at a time with relevancy
  grid.arrange(p1, p2, ncol = 2)
  grid.arrange(p3, p4, ncol = 2)
  
  # Return performance metrics
  return(list(
    reduced_computation_time = reduced_computation_time
  ))
}

# Run the simulation and get performance metrics
performance_metrics <- run_simulation()
print(performance_metrics)
