load('NGS_simulated_dataset.mat')
 
fid = fopen('recipient_total_reads.txt','wt');
 fprintf(fid,'%f\n',sim_data.recipient_total_reads);  % The format string is applied to each element of a
 fclose(fid);
 
 fid = fopen('recipient_var_reads.txt','wt');
 fprintf(fid,'%f\n',sim_data.recipient_var_reads_observed);  % The format string is applied to each element of a
 fclose(fid);
 
 fid = fopen('donor_freqs.txt','wt');
  fprintf(fid,'%f\n',sim_data.donor_freqs_observed);  % The format string is applied to each element of a
  fclose(fid);
  
   fid = fopen('recipient_freqs.txt','wt');
  fprintf(fid,'%f\n',sim_data.recipient_freqs_observed);  % The format string is applied to each element of a
  fclose(fid);
  
   fid = fopen('n_variants.txt','wt');
  fprintf(fid,'%f\n',sim_data.n_variants);  % The format string is applied to each element of a
  fclose(fid);
  
  fid = fopen('var_calling_threshold.txt','wt');
  fprintf(fid,'%f\n',sim_data.var_calling_threshold);  % The format string is applied to each element of a
  fclose(fid);