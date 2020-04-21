#########################################################
#install.packages("argparse")
library(argparse)
library(mdatools)
parser <- ArgumentParser()

parser$add_argument("--file", type="character", default= "no_input_return_error",
    help="file containing variant frequencies")
parser$add_argument("--plot_bool", type="logical", default= FALSE,
    help="determines whether pdf plot approx_plot.pdf is produced or not")
parser$add_argument("--var_calling_threshold", type="double", default= 0.03,
    help="variant calling threshold")
parser$add_argument("--Nb_min", type="integer", default= 1,
    help="Minimum bottleneck value considered")
parser$add_argument("--Nb_max", type="integer", default= 200,
    help="Maximum bottleneck value considered")
parser$add_argument("--confidence_level", type="double", default= .95,
    help="Confidence level (determines bounds of confidence interval)")
args <- parser$parse_args()
if (args$file == "no_input_return_error" ) { stop("file with lists of donor and recipient frequencies is a required argument.", call.=FALSE)}
 
 
plot_bool  <- args$plot_bool
var_calling_threshold  <- args$var_calling_threshold
Nb_min <- args$Nb_min # Minimum bottleneck size we consider. 
Nb_max <-  args$Nb_max
confidence_level <- args$confidence_level
donor_and_recip_freqs_observed <- read.table(args$file) 
original_row_count <- nrow(donor_and_recip_freqs_observed)

donor_and_recip_freqs_observed <- subset(donor_and_recip_freqs_observed, donor_and_recip_freqs_observed[, 1] >= var_calling_threshold)
new_row_count <- nrow(donor_and_recip_freqs_observed)
if(new_row_count != original_row_count )
{print("WARNING:  Rows of the input file with donor frequency less than variant calling threshold have been removed during analysis. ")}
#donor_freqs_observed <- read.table(args[1])
donor_freqs_observed <- as.data.frame(donor_and_recip_freqs_observed[, 1]) # We read in and save the list of donor frequencies
recipient_freqs_observed <- as.data.frame(donor_and_recip_freqs_observed[, 2]) # We read in and save the list of recipient frequencies
n_variants <- nrow(donor_and_recip_freqs_observed)#nrow(donor_freqs_observed) # number of variants 

#######################################################
#######################################################
generate_log_likelihood_function_approx <- function(donor_freqs_observed, recipient_freqs_observed, Nb_min, Nb_max, var_calling_threshold, confidence_level)
{

    num_NB_values <- Nb_max -Nb_min + 1
likelihood_matrix <- matrix( 0, n_variants, num_NB_values)#This matrix stores the likelihood values for each variant frequency and bottleneck size.  
log_likelihood_matrix <- matrix( 0, n_variants, num_NB_values)#This matrix stores the log likelihood values for each variant frequency and bottleneck size.  We can sum the log likelihoods of all the variant alleles to get the total log likelihood. 
log_likelihood_function <- matrix( 0, Nb_max )


for (i in 1:n_variants) {for (j in 1:num_NB_values) {
  Nb_val <- (j - 1 + Nb_min)
    nu_donor <- donor_freqs_observed[i, 1]
    nu_recipient <- recipient_freqs_observed[i, 1]
    if (recipient_freqs_observed[i, 1] >= var_calling_threshold)
        { # implement variant calling threshold
for (k in 0:Nb_val){  
      likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
  (dbeta(nu_recipient, k, (Nb_val - k))*dbinom(k, size=Nb_val, prob= nu_donor))
        }
  log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])  
   }
 if (recipient_freqs_observed[i, 1] < var_calling_threshold)
        { # implement variant calling threshold
   
   likelihood_matrix[i, j] = 0
   log_likelihood_matrix[i,j] = 0
   for (k in 0:Nb_val){   likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
               (pbeta(var_calling_threshold, k, (Nb_val - k))*dbinom(k,size=Nb_val, prob= nu_donor))
                      } 
log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])
            }
# Now we sum over log likelihoods of the variants at different loci to get the total log likelihood for each value of Nb
log_likelihood_function[ Nb_val] <- log_likelihood_function[ Nb_val] + log_likelihood_matrix[i,j]
# Shifting entries of log_likelihood function to ensure plot begins at proper point on x axis
}}


return(log_likelihood_function)

}
#################################################################
#################################################################
restrict_log_likelihood <- function(log_likelihood_function, Nb_min, Nb_max) # restricts log likelihood to the interval of interst
{
for (h in 1:(Nb_min )){  
    if(h< Nb_min)
    {log_likelihood_function[h] = - 999999999}        # kludge for ensuring that these values less than Nb_min don't interfere with our search for the max of log likelihood in the interval of Nb_min to Nb_max
    }

return(log_likelihood_function)
#print(erfinv(percent_confidence_interval)*sqrt(2))
}


return_bottleneck_size <- function(log_likelihood_function, Nb_min, Nb_max)
{ max_log_likelihood = which(log_likelihood_function == max(log_likelihood_function))

    return(max_log_likelihood)
}


return_CI_lower <- function(log_likelihood_function, Nb_min, Nb_max) ## returns lower bound of confidence interval
  { max_log_likelihood = which(log_likelihood_function == max(log_likelihood_function))  ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
max_val =  max(log_likelihood_function)  ## This is the maximum value of the log likelihood function, found when the index is our bottleneck estimate
CI_height = max_val - erfinv(confidence_level)*sqrt(2)  # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test


    CI_index_lower = Nb_min
  CI_index_upper = max_log_likelihood
for (h in 1:Nb_min){  
    if(h< Nb_min)
    {log_likelihood_function[h] = NA}   #  Removing parameter values less than Nb_min from plot
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
{       test1 = (log_likelihood_function[CI_index_upper] - CI_height) * (log_likelihood_function[CI_index_upper] - CI_height)
      test2 = (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
        if( test2 < test1  ){CI_index_upper = h   }  
}
if(  (log_likelihood_function[CI_index_upper] - CI_height) > 0  ){CI_index_upper = CI_index_upper + 1   }  
    
  return(CI_index_lower)

  }


return_CI_upper <- function(log_likelihood_function,  Nb_min, Nb_max) ## returns upper bound of confidence interval
  { max_log_likelihood = which(log_likelihood_function == max(log_likelihood_function))  ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
   max_val =  max(log_likelihood_function)  ## This is the maximum value of the log likelihood function, found when the index is our bottleneck estimate
   CI_height = max_val - erfinv(confidence_level)*sqrt(2)  # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test

    CI_index_lower = Nb_min
  CI_index_upper = max_log_likelihood
for (h in 1:Nb_min){  
    if(h< Nb_min)
    {log_likelihood_function[h] = NA}   #  Removing parameter values less than Nb_min from plot
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
{       test1 = (log_likelihood_function[CI_index_upper] - CI_height) * (log_likelihood_function[CI_index_upper] - CI_height)
      test2 = (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
        if( test2 < test1  ){CI_index_upper = h   }  
}
if(  (log_likelihood_function[CI_index_upper] - CI_height) > 0  ){CI_index_upper = CI_index_upper + 1   }  
    
  return(CI_index_upper)

  }


log_likelihood_function <- generate_log_likelihood_function_approx(donor_freqs_observed, recipient_freqs_observed, Nb_min, Nb_max, var_calling_threshold, confidence_level)
log_likelihood_function <- restrict_log_likelihood(log_likelihood_function, Nb_min, Nb_max)
bottleneck_size <- return_bottleneck_size(log_likelihood_function,  Nb_min, Nb_max)
CI_index_lower <- return_CI_lower(log_likelihood_function,  Nb_min, Nb_max)
CI_index_upper <- return_CI_upper(log_likelihood_function,  Nb_min, Nb_max)
        
        ##########################
##############################################################################################  ABOVE THIS LINE DETERMINES PEAK LOG LIKELIHOOD AND CONFIDENCE INTERVALS
 # Npw we plot the result
if(plot_bool == TRUE)
{pdf(file="approx_plot.pdf")
plot(log_likelihood_function)
abline(v = bottleneck_size, col="black" )  # Draws a verticle line at Nb value for which log likelihood is maximized
abline(v = CI_index_lower, col="green" ) # confidence intervals
abline(v = CI_index_upper, col="green" )
dev.off()
}

print("Bottleneck size")
print(bottleneck_size)
print("confidence interval left bound")
print(CI_index_lower)
print("confidence interval right bound")
print(CI_index_upper)

