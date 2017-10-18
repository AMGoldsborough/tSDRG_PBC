function ee = eentropy(S)
%ee = eentropy(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entanglement (Von Neumann) Entropy
% function to calculate the entanglement (Von Neumann) entropy.
% ee = - \sum_{a}^{r} s_{a}^{2} log_{2} (s_{a}^{2}) 
% where S = diag(s_{1}, ... ,s_{r})
%
% input: matrix of singular values S
% output: the entanglement entropy
% 
% Andrew Goldsborough - 12/11/2012
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create array of singular values
Sarray = diag(S);
Sarray(Sarray==0) = [];

%calculate ee
ee = - sum(Sarray.*Sarray.*log2(Sarray.*Sarray));
