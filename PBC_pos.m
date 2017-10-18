%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = PBC_pos(x,L)
%gives position with PBCs
x = mod(x-1,L)+1;
end