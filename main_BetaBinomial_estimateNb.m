function void = main_BetaBinomial_estimateNb(void)

clear all; close all; clc; format long g;

sim_data_results.infile = 'NGS_simulated_dataset';

sim_data_results.Nb_min = 5;
sim_data_results.Nb_max = 200;
sim_data_results.Nb_vals = sim_data_results.Nb_min:sim_data_results.Nb_max;

load(sim_data_results.infile);
outfile = strcat(sim_data_results.infile, '_betabinomial');

sim_data_results.log_likelihood_matrix = NaN*zeros(sim_data.n_variants, length(sim_data_results.Nb_vals)); % keeps log-likelihood of each var site for each Nb value

for site = 1:sim_data.n_variants
    [site sim_data.n_variants]
    sim_data_results.log_likelihood_matrix(site,:) = GetBetaBinomialLogLikelihoodAtSite(sim_data.donor_freqs_observed(site), sim_data.recipient_var_reads_observed(site), sim_data.recipient_total_reads(site), sim_data.recipient_freqs_observed(site), sim_data_results.Nb_vals, sim_data.var_calling_threshold);
    if mod(site, 50) == 0
        save(outfile, 'sim_data', 'sim_data_results')
        display('saved')
    end
end
sim_data_results.overall_log_likelihood = sum(sim_data_results.log_likelihood_matrix, 1);
[sim_data_results.Nb_MLE, sim_data_results.CI_low_Nb, sim_data_results.CI_high_Nb] = GetLikelihoodConfidenceIntervals(sim_data_results.Nb_vals, sim_data_results.overall_log_likelihood);

save(outfile, 'sim_data', 'sim_data_results');

plot(sim_data_results.Nb_vals, sim_data_results.overall_log_likelihood, 'g'); hold on;
xlabel('N_b'); ylabel('log-likelihood'); hold on; 
ax_vals = axis;
plot([sim_data_results.Nb_MLE, sim_data_results.Nb_MLE], [ax_vals(3), ax_vals(4)], 'g');
plot([sim_data_results.CI_low_Nb sim_data_results.CI_low_Nb], [ax_vals(3), ax_vals(4)], 'g-.');
plot([sim_data_results.CI_high_Nb sim_data_results.CI_high_Nb], [ax_vals(3), ax_vals(4)], 'g-.');
plot([sim_data.Nb, sim_data.Nb], [ax_vals(3), ax_vals(4)], 'k--');

