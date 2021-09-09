functions {
  real icar_normal_lpdf(vector psi, int J, int[] node1, int[] node2) {
    return -0.5 * dot_self(psi[node1] - psi[node2])
      + normal_lpdf(sum(psi) | 0, 0.001 * J); // equivalent to mean(phi) ~ normal(0,0.001)
 }
}
data {
  int<lower = 0> J; // Change 'N' to J, for the number of sptatial locations.
  int<lower = 0> J_edges; // edge set
  int<lower = 1, upper = J> node1[J_edges];  // node1[i] ~ node2[i]
  int<lower = 1, upper = J> node2[J_edges];  // and node1[i] < node2[i]
  vector<lower = 0, upper = 1>[J] p; // binary outcomes (or conditional inclusion probs)
}
transformed data{
  vector[J] psi = logit(p);
}
model {
  psi ~ icar_normal_lpdf(J, node1, node2);
}

