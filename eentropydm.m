function ee = eentropydm(dm)
%ee = eentropydm(dm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entanglement (Von Neumann) Entropy calculated from the density matrix
% function to calculate the entanglement (Von Neumann) entropy
% ee = - tr(dm * log2(dm))
%
% input: density matrix dm
% output: the entanglement entropy
%
% Andrew Goldsborough - 19/02/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make sure that the dm's are hermitean
dm = 0.5*(dm + dm');

[dmv,dmnrg] = eig(dm,'vector');

%check for negative eigs
if min(dmnrg) <= 0
    %check if numerically small and fix (for stability)
    if abs(min(dmnrg)) < 1E-15
        dmnrg(dmnrg<=0) = 1e-16;
    else
        error('dm is not positive definite');
    end
end

ee = -trace(dm*dmv*((diag(log2(dmnrg))+diag(log2(dmnrg))')/2)*dmv');
