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
  return(new("markovchain", states = as.character(states), transitionMatrix = P_red))
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
    theme_minimal() +
    labs(title = title, x = "From State", y = "To State")
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
  reduced_mc <- reduce_state_space(original_mc, k)
  reduced_steady_state <- compute_steady_state(reduced_mc)
  
  # Plot transition matrices
  p1 <- plot_transition_matrix(original_mc, "Original Transition Matrix")
  p2 <- plot_transition_matrix(reduced_mc, "Reduced Transition Matrix")
  
  # Plot steady-state distributions
  df1 <- data.frame(State = 1:length(original_steady_state), Probability = original_steady_state)
  df2 <- data.frame(State = 1:length(reduced_steady_state), Probability = reduced_steady_state)
  
  p3 <- ggplot(df1, aes(x = State, y = Probability)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = "Original Steady-State Distribution", x = "State", y = "Probability")
  
  p4 <- ggplot(df2, aes(x = State, y = Probability)) +
    geom_bar(stat = "identity", fill = "salmon") +
    theme_minimal() +
    labs(title = "Reduced Steady-State Distribution", x = "State", y = "Probability")
  
  # Save plots
  ggsave("original_transition_matrix.png", plot = p1, width = 8, height = 6)
  ggsave("reduced_transition_matrix.png", plot = p2, width = 8, height = 6)
  ggsave("original_steady_state_distribution.png", plot = p3, width = 8, height = 6)
  ggsave("reduced_steady_state_distribution.png", plot = p4, width = 8, height = 6)
  
  # Display plots
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

# Run the simulation
run_simulation()
