#  Put the long lists of donor and recipient frequencies in seperate text files that we then read from.  Total number of reads and number of variant reads for each loci in recipient are also located in seperate text files.  Same goes for n_variants and var_calling_threshold.

# This code requires use of the R package rmutil 




library(rmutil)

 n_variants_dummy <- read.table("n_variants.txt")
n_variants <-  n_variants_dummy[1, 1]



 var_calling_threshold_dummy  <- read.table("var_calling_threshold.txt")
 
 var_calling_threshold  <-  var_calling_threshold_dummy[1, 1] 
 




donor_freqs_observed <- read.table("donor_freqs.txt")



recipient_total_reads <- read.table("recipient_total_reads.txt")



recipient_var_reads_observed <- read.table("recipient_var_reads.txt")


recipient_freqs_observed <- read.table("recipient_freqs.txt")

Nb_min <- 5

Nb_max <- 200

likelihood_matrix <- matrix( 0, n_variants, (Nb_max -Nb_min + 1))

log_likelihood_matrix <- matrix( 0, n_variants, (Nb_max -Nb_min + 1))

log_likelihood_function <- matrix( 0, (Nb_max ))

# create array of likelihoods for every variant and every Nb value
for (i in 1:n_variants) {for (j in 1:(Nb_max -Nb_min + 1)) {
  Nb_val <- (j - 1 + Nb_min)
    
  
  
           	    	nu_donor <- donor_freqs_observed[i, 1]
    #nu_recipient <- recipient_freqs_observed[i, 1]

   	     	variant_reads <- recipient_var_reads_observed[i, 1]
   
             total_reads <- recipient_total_reads[i, 1]
               
               


    if (variant_reads >= var_calling_threshold*total_reads)
   	    { # implement variant calling threshold

for (k in 0:Nb_val){  
	alpha <- k
	 	
	beta <- (Nb_val - k)
	
	if (alpha == 0)
   	    { alpha <- 0.00001 }

if (beta == 0)
   	    { beta <- 0.00001 }

	
	m <- alpha/(alpha + beta)
	
	s <- (alpha + beta)
	
	
	
	
	
	likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
	(dbetabinom( variant_reads, total_reads, m, s, log = FALSE)*dbinom(k, size=Nb_val, prob= nu_donor)) 
	
	      }
	
	 log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])  
	 }

    if (variant_reads < var_calling_threshold*total_reads)
   	    { # implement variant calling threshold

likelihood_matrix[i, j] = 0
   log_likelihood_matrix[i,j] = 0
   for (k in 0:Nb_val){  
	
	alpha <- k
	 	
	beta <- (Nb_val - k)	
	if (alpha == 0)
   	    { alpha <- 0.00001 }

	if (beta == 0)
   	    { beta <- 0.00001 }

	
	m <- alpha/(alpha + beta)
	
	s <- (alpha + beta)

	#likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
	#(pbeta(var_calling_threshold, k, (Nb_val - k))*dbinom(k, size=Nb_val, prob= nu_donor)) 

	
	
	likelihood_matrix[i, j] <- likelihood_matrix[i, j] + 
	(pbetabinom( floor(var_calling_threshold*total_reads), total_reads, m, s)*dbinom(k, size=Nb_val, prob= nu_donor)) 
	
	}

log_likelihood_matrix[i,j] = log(likelihood_matrix[i, j])



  }

# Now we sum over log likelihoods of the variants at different loci to get the total log likelihood for each value of Nb

log_likelihood_function[(j + Nb_min) -1] <- log_likelihood_function[(j + Nb_min) -1] + log_likelihood_matrix[i,j]




}}

for (h in 1:(Nb_min -1)){  
		
		log_likelihood_function[h] = - 999999999	  # trick for ensuring that these values less than Nbmin don't interfere with our search for the max of log likelihood in the interval of Nb_min to Nb_max

	      }

dummy = which(log_likelihood_function == max(log_likelihood_function))  ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized

max_val =  max(log_likelihood_function)

dummy2 = max(log_likelihood_function) - 1.92  # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test


dummy3 = Nb_min

dummy4 = dummy

for (h in 1:(Nb_min -1)){  
		
		log_likelihood_function[h] = NA	  #  Removing parameter values less than Nb_min from plot
	      }

for (h in Nb_min:dummy){  
		
		test1 = (log_likelihood_function[dummy3] - max_val + 1.92) * (log_likelihood_function[dummy3] - max_val + 1.92)
		
		test2 = (log_likelihood_function[h] - max_val + 1.92) * (log_likelihood_function[h] - max_val + 1.92)

		if( test2 < test1){  dummy3 = h  }  			
		
	      }

if(  (log_likelihood_function[dummy3] - max_val + 1.92) > 0  ){dummy3 = dummy3 - 1   }  



for (h in dummy:Nb_max){  
		
		
		test1 = (log_likelihood_function[dummy4] - max_val + 1.92) * (log_likelihood_function[dummy4] - max_val + 1.92)
		
		test2 = (log_likelihood_function[h] - max_val + 1.92) * (log_likelihood_function[h] - max_val + 1.92)

		  
		  if( test2 < test1  ){dummy4 = h   }  
		  	
		  	
		  	}
	      
	      
	      if(  (log_likelihood_function[dummy4] - max_val + 1.92) > 0  ){dummy4 = dummy4 + 1   }  
	      
	      #  The above loops find the edges of the confidence interval using the likelihood ratio test

	      
plot(log_likelihood_function)
abline(v = dummy, col="black" )  # Draws a verticle line at Nb value for which log likelihood is maximized

abline(v = dummy3, col="green" ) 
                                     # confidence intervals

abline(v = dummy4, col="green" )

print("Bottleneck size")
print(dummy)

print("confidence interval left bound")
print(dummy3)

print("confidence interval right bound")
print(dummy4)






