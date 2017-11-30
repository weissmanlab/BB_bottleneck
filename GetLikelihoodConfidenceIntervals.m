function [MLE, CI_low, CI_high] = GetLikelihoodConfidenceIntervals(Nb, Nb_loglikelihoods)

locsOK = intersect(find(~isnan(Nb_loglikelihoods)), find(~isinf(Nb_loglikelihoods)));
max_loc = find(Nb_loglikelihoods == max(Nb_loglikelihoods(locsOK)));
MLE = Nb(max_loc);

max_likelihood_value = Nb_loglikelihoods(max_loc);

% 95% CI using likelihood ratio test
CI_likelihood_value = max_likelihood_value - 1.92;  
CI_likelihood_Nbs = find(Nb_loglikelihoods >= CI_likelihood_value); 
CI_likelihood_Nbs = intersect(find(Nb_loglikelihoods >= CI_likelihood_value), find(~isinf(Nb_loglikelihoods))); 

CI_low = Nb(min(CI_likelihood_Nbs)) - 1;
CI_high = Nb(max(CI_likelihood_Nbs)) + 1;