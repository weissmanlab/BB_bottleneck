#  Put the long lists of donor and recipient frequencies in seperate text files that we then read from.  Same goes for n_variants and var_calling_threshold.
# First we read in the relevant files and set the parameters of the code
n_variants_dummy <- read.table("n_variants.txt")
n_variants <-  n_variants_dummy[1, 1]
# We read in and then save the number of variants
var_calling_threshold_dummy  <- read.table("var_calling_threshold.txt")
var_calling_threshold  <-  var_calling_threshold_dummy[1, 1]
# We read in and then save the variant calling threshold
donor_freqs_observed <- read.table("donor_freqs.txt")
# We read in and save the list of donor frequencies
recipient_total_reads <- matrix( 0, n_variants, 1)
# We read in and save the list of recipient frequencies
recipient_freqs_observed <- read.table("recipient_freqs.txt")
Nb_min <- 1
# Minimum bottleneck size we consider. 
Nb_max <- 200
# Maximum bottleneck size we consider.
############################################################.  LOAD IN DATA AND DECLARE ARRAYS TO BE FILLED
num_NB_values <- Nb_max -Nb_min + 1
likelihood_matrix <- matrix( 0, n_variants, num_NB_values)
#This matrix stores the likelihood values for each variant frequency and bottleneck size.  
log_likelihood_matrix <- matrix( 0, n_variants, num_NB_values)
#This matrix stores the log likelihood values for each variant frequency and bottleneck size.  We can sum the log likelihoods of all the variant alleles to get the total log likelihood. 
log_likelihood_function <- matrix( 0, Nb_max )
# create array of likelihoods for every variant and every Nb value
###################################################################################### BELOW THIS LINE CREATES LIKELIHOOD AND LOG LIKELIHOOD 
for (i in 1:n_variants) {for (j in 1:num_NB_values) {
  Nb_val <- (j - 1 + Nb_min)
    nu_donor <- donor_freqs_observed[i, 1]
    nu_recipient <- recipient_freqs_observed[i, 1]
    if (recipient_freqs_observed[i, 1] >= var_calling_threshold)
   	    { # implement variant calling threshold
for (k in 0:Nb_val){  
      likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
	(dbeta(nu_recipient, k, (Nb_val - k))*dbinom(k, size=Nb_val, prob= nu_donor)) 
	      }
	log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])  
	 }
 if (recipient_freqs_observed[i, 1] < var_calling_threshold)
   	    { # implement variant calling threshold
   
   likelihood_matrix[i, j] = 0
   log_likelihood_matrix[i,j] = 0
   for (k in 0:Nb_val){   likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
	             (pbeta(var_calling_threshold, k, (Nb_val - k))*dbinom(k, size=Nb_val, prob= nu_donor)) 
                      } 
log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])
            }
# Now we sum over log likelihoods of the variants at different loci to get the total log likelihood for each value of Nb
log_likelihood_function[ Nb_val] <- log_likelihood_function[ Nb_val] + log_likelihood_matrix[i,j]
# Shifting entries of log_likelihood function to ensure plot begins at proper point on x axis
}}
############################################################################################### ABOVE THIS LINE CALCULATES LIKELIHOOD AND LOG LIKELIHOOD 
###############################################################################################  BELOW THIS LINE DETERMINES PEAK LOG LIKELIHOOD AND CONFIDENCE INTERVALS
for (h in 1:(Nb_min )){  
		if(h< Nb_min)
		{log_likelihood_function[h] = - 999999999}	      # kludge for ensuring that these values less than Nb_min don't interfere with our search for the max of log likelihood in the interval of Nb_min to Nb_max
	  }
max_log_likelihood = which(log_likelihood_function == max(log_likelihood_function))  ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
max_val =  max(log_likelihood_function)
CI_height = max_val - 1.92  # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test
CI_index_lower = Nb_min
CI_index_upper = max_log_likelihood

for (h in 1:Nb_min){  
		if(h< Nb_min)
		{log_likelihood_function[h] = NA}	  #  Removing parameter values less than Nb_min from plot
	      }
## above loop just enforces our minimum bottleneck cutoff
for (h in Nb_min:max_log_likelihood){  
		test1 = (log_likelihood_function[CI_index_lower] - CI_height) * (log_likelihood_function[CI_index_lower] - CI_height)
		test2 = (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
        if( test2 < test1){  CI_index_lower = h  }  			
}

if(  (log_likelihood_function[CI_index_lower] - CI_height) > 0  ){CI_index_lower = CI_index_lower - 1   }  
# above loops use likelihood ratio test to find lower confidence interval
for (h in max_log_likelihood:Nb_max)
{  
	    test1 = (log_likelihood_function[CI_index_upper] - CI_height) * (log_likelihood_function[CI_index_upper] - CI_height)
	    test2 = (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
        if( test2 < test1  ){CI_index_upper = h   }  
}
if(  (log_likelihood_function[CI_index_upper] - CI_height) > 0  ){CI_index_upper = CI_index_upper + 1   }  
		  	# above loops use likelihood ratio test to find upper confidence interval
		  	##########################
##############################################################################################  ABOVE THIS LINE DETERMINES PEAK LOG LIKELIHOOD AND CONFIDENCE INTERVALS
 # Npw we plot the result
plot(log_likelihood_function)
abline(v = max_log_likelihood, col="black" )  # Draws a verticle line at Nb value for which log likelihood is maximized
abline(v = CI_index_lower, col="green" ) # confidence intervals
abline(v = CI_index_upper, col="green" )
print("Bottleneck size")
print(max_log_likelihood)
print("confidence interval left bound")
print(CI_index_lower)
print("confidence interval right bound")
print(CI_index_upper)