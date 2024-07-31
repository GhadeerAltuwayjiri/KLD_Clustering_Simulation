import numpy as np
from sklearn.mixture import GaussianMixture
import warnings

# Define the KLD function
def KLD(mu_p, sigma_p, mu_q, sigma_q):
    d = len(mu_p)
    KL1 = np.log(np.linalg.det(sigma_q)) - np.log(np.linalg.det(sigma_p)) - d
    KL2 = np.dot(np.dot((mu_p - mu_q).T, np.linalg.inv(sigma_q)), (mu_p - mu_q))
    KL3 = np.trace(np.dot(np.linalg.inv(sigma_q), sigma_p))
    KL = 0.5 * (KL1 + KL2 + KL3)
    return KL

# Define the cov_sim variations
cov_sim_values = [0, 0.1, 0.2, 0.5, 0.7]

# Initialize a list to hold the simulated data for each cov_sim value
simulated_data_lists = {}

# Example placeholders for mu and sigma
mu = np.random.rand(6, 9)
sigma = np.array([np.eye(9) for _ in range(6)])
pi = np.random.dirichlet(np.ones(6), size=1)[0]
G = 6
n = 100

# Iterate over the cov_sim values
for cov_sim in cov_sim_values:
    simulated_sample_list = []
    flow_sim_list = []
    MM_list = []
    Zhat_list1 = []
    Ztoy_list1 = []
    pi1 = []
    
    num_simulations = 5
    
    for j in range(num_simulations):
        np.random.seed(123 + j)
        
        sd_channel = np.std(mu, axis=1)
        
        mu_samp_list = []
        for k in range(5):
            mu_samp_list.append(mu + np.random.normal(scale=cov_sim * sd_channel[:, None], size=mu.shape))
        
        Ztoy = np.random.choice(range(G), size=n, p=pi)
        ytoy = np.array([np.random.multivariate_normal(mu_samp_list[j][:, z], sigma[z]) for z in Ztoy])
        
        try:
            flow_sim = GaussianMixture(n_components=G).fit(ytoy)
            Zhat = flow_sim.predict(ytoy)
            simulated_sample_list.append(ytoy)
            flow_sim_list.append(flow_sim)
            pi1.append(flow_sim.weights_)
            MM_list.append(ytoy)
            Zhat_list1.append(Zhat)
            Ztoy_list1.append(Ztoy)
        except Exception as e:
            warnings.warn(f"Clustering failed for simulation {j}: {e}")
    
    simulated_data_lists[f"cov_sim_{cov_sim}"] = {
        "simulated_samples": simulated_sample_list,
        "flow_sims": flow_sim_list,
        "MM": MM_list,
        "Zhat": Zhat_list1,
        "pi": pi1,
        "Ztoy": Ztoy_list1
    }
