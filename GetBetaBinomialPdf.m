function return_value = GetBetaBinomialPdf(kP, nP, alphaP, betaP)
%x. n. a. b
%error('in betabinomialpdf - still need to update!!')

% kvk added this on 8/11/2016: this is the case because gamma(alphaP) and
% gamma(betaP) are both in the demoninator, and gamma(0) = Inf.
if alphaP == 0
    return_value = 0; return;
end
if betaP == 0
    return_value = 0; return;
end
% end of add

vector1_num = GetGammaVector(nP+1);
vector2_num = GetGammaVector(kP+alphaP);
vector3_num = GetGammaVector(nP - kP + betaP);
vector4_num = GetGammaVector(alphaP+betaP);
vector_num = [vector1_num vector2_num vector3_num vector4_num];

vector1_denom = GetGammaVector(kP+1);
vector2_denom = GetGammaVector(nP - kP+ 1);
vector3_denom = GetGammaVector(nP + alphaP + betaP);
vector4_denom = GetGammaVector(alphaP);
vector5_denom = GetGammaVector(betaP);
vector_denom = [vector1_denom vector2_denom vector3_denom vector4_denom vector5_denom];

standard_length = max(length(vector_num), length(vector_denom));
final_vector_num = ones(1, standard_length); final_vector_num(1:length(vector_num)) = vector_num; final_vector_num = sort(final_vector_num);
final_vector_denom = ones(1, standard_length); final_vector_denom(1:length(vector_denom)) = vector_denom; final_vector_denom = sort(final_vector_denom);
return_value = prod(final_vector_num./final_vector_denom);
