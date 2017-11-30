function log_likelihood_sitevector = GetBetaBinomialLogLikelihoodAtSite(nu_donor, var_reads_recipient, tot_reads_recipient, nu_recipient, Nb_vals, var_calling_threshold)

if nu_recipient < var_calling_threshold
    var_is_absent_from_recipient = 1;
else
    var_is_absent_from_recipient = 0;
end

cntr = 1;
for Nb = Nb_vals
    log_likelihood_sitevector(1,cntr) = GetBetaBinomialLogLikelihoodAtSiteForNb(nu_donor, var_reads_recipient, tot_reads_recipient, var_is_absent_from_recipient, var_calling_threshold, Nb);
    cntr = cntr + 1;
end
