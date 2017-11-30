function log_likelihood_site_forNb = GetBetaBinomialLogLikelihoodAtSiteForNb(nu_donor, var_reads_recipient, tot_reads_recipient, var_is_absent_from_recipient, var_calling_threshold, Nb)

prob_of_n_minor_var = binopdf(0:Nb, Nb, nu_donor);

if var_is_absent_from_recipient
    % just for the small number of variants code!!
    for n_minor_var = 0:Nb
        if n_minor_var == 0
           val = 1;
        elseif n_minor_var == Nb
            val = 0;
        else
            val = betacdf(var_calling_threshold, n_minor_var, (Nb-n_minor_var)); % does not incorporate read sampling stochasticity
        end
        %val = GetBetaBinomialCdf(floor(var_calling_threshold*tot_reads_recipient), tot_reads_recipient, n_minor_var, (Nb-n_minor_var));
        likelihood_site(1,n_minor_var + 1) = val;
    end
else
    for n_minor_var = 0:Nb
        %val = bbinopdf(var_reads_recipient, tot_reads_recipient, n_minor_var, (Nb-n_minor_var)); 
        val = GetBetaBinomialPdf(var_reads_recipient, tot_reads_recipient, n_minor_var, (Nb-n_minor_var));
        %if (isnan(val) || isinf(val))
        %    val = GetBetaBinomialPdf(var_reads_recipient, tot_reads_recipient, n_minor_var, (Nb-n_minor_var));
        %end
        likelihood_site(1,n_minor_var + 1) = val;
    end
end

likelihood_site_forNb = sum(prob_of_n_minor_var.*likelihood_site);
log_likelihood_site_forNb = log(likelihood_site_forNb);
