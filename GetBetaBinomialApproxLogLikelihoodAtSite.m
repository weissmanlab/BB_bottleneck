function log_likelihood_sitevector = GetBetaBinomialApproxLogLikelihoodAtSite(nu_donor, nu_recipient, Nb_vals, var_calling_threshold)

if nu_recipient < var_calling_threshold
    var_is_absent_from_recipient = 1;
else
    var_is_absent_from_recipient = 0;
end

cntr = 1;
for Nb = Nb_vals
    log_likelihood_sitevector(1,cntr) = GetBetaBinomialApproxLogLikelihoodAtSiteForNb(nu_donor, nu_recipient, var_is_absent_from_recipient, var_calling_threshold, Nb);
    cntr = cntr + 1;
end
