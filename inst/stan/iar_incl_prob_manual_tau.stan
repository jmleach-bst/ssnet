functions {
  real real_bernoulli_lpdf(vector x, vector theta, int J) {
    vector[J] xs; // define each term of sum in log "pdf"
    for (j in 1:J) xs[j] =  x[j] * log(theta[j]) +  (1 - x[j]) * log(1 - theta[j]); // evaluate each term of sum in log "pdf"
    return sum(xs); // up to constant, log "pdf" is this
  }
  real icar_normal_lpdf(vector psi, int J, int[] node1, int[] node2, real tau) {
    return ((J - 1.0) / 2) * log(tau) - 0.5 * tau * dot_self(psi[node1] - psi[node2])
      + normal_lpdf(sum(psi) | 0, 0.001 * J); // equivalent to mean(phi) ~ normal(0,0.001)
 }
}
data {
  int<lower = 0> J; // Change 'N' to J, for the number of sptatial locations.
  int<lower = 0> J_edges; // edge set
  int<lower = 1, upper = J> node1[J_edges];  // node1[i] ~ node2[i]
  int<lower = 1, upper = J> node2[J_edges];  // and node1[i] < node2[i]
  vector<lower = 0, upper = 1>[J] p; // binary outcomes (or conditional inclusion probs)
  real<lower = 0> tau; // manually set precision
}
parameters {
  vector[J] psi; // logit of prior probability of inclusion
}
transformed parameters{
  vector[J] theta = inv_logit(psi);
}
model {
  p ~ real_bernoulli(theta, J); // for conditional probabilities of inclusion
  psi ~ icar_normal_lpdf(J, node1, node2, tau);
}
