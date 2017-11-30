function log_likelihood_site_forNb = GetMethod3LogLikelihoodAtSiteForNbApprox(nu_donor, nu_recipient, var_is_absent_from_recipient, var_calling_threshold, Nb)

prob_of_n_minor_var = binopdf(0:Nb, Nb, nu_donor);

if var_is_absent_from_recipient
    for n_minor_var = 0:Nb
        if n_minor_var == 0
            likelihood_site(1,n_minor_var + 1) = 1;
        elseif n_minor_var == Nb
            likelihood_site(1,n_minor_var + 1) = 0;
        else
            likelihood_site(1,n_minor_var + 1) = betacdf(var_calling_threshold, n_minor_var, (Nb-n_minor_var)); % does not incorporate read sampling stochasticity
        end
    end
else
    for n_minor_var = 0:Nb
        if n_minor_var == 0
            likelihood_site(1,n_minor_var + 1) = 0;
        elseif n_minor_var == Nb
            likelihood_site(1,n_minor_var + 1) = 0;
        else
            likelihood_site(1,n_minor_var + 1) = betapdf(nu_recipient, n_minor_var, (Nb-n_minor_var)); % does not incorporate read sampling stochasticity
        end
    end
end

likelihood_site_forNb = sum(prob_of_n_minor_var.*likelihood_site);
log_likelihood_site_forNb = log(likelihood_site_forNb);
