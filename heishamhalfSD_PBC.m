function [W,d] = heishamhalfSD_PBC(L,J,Jz,h)
%[W,d] = heishamhalfSD_PBC(L,J,Jz,h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spin 1/2 Heisenberg model hamiltonian with PBCs
% input: chain length L, array of interaction strengths J, relative z strength Jz, array of on site strengths h
% output: MPOs for the Heisenberg Hamiltonian W, spin dimension d
% 
% Andrew Goldsborough 10/11/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = cell(L,1);

%spin-1/2 Heisenberg 
for i=1:(L)
    W{i} = zeros(2,2,5,5);
    W{i}(:,:,1,1) = [1 0;0 1];
    W{i}(:,:,1,2) = [0 0.5*J(i);0 0];
    W{i}(:,:,1,3) = [0 0;0.5*J(i) 0];
    W{i}(:,:,1,4) = [0.5*J(i)*Jz 0;0 -0.5*J(i)*Jz];
    W{i}(:,:,1,5) = [0.5*h(i) 0;0 -0.5*h(i)];
    W{i}(:,:,2,5) = [0 0;1 0];
    W{i}(:,:,3,5) = [0 1;0 0];
    W{i}(:,:,4,5) = [0.5 0;0 -0.5];
    W{i}(:,:,5,5) = [1 0;0 1];
end

%spin dimension for spin-1/2 is 2 (up, down)
d = 2;
