function B = tsplit(A,leg,legsizes)
%B = tsplit(A,leg,legsizes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tsplit - split fused indices
% input: Tensor A, leg to be split 'leg', size of the resulting legs after the split 'legsizes'
% output: Tensor B with desired leg split
% 
% Andrew Goldsborough - 21/09/2012
%
% Note: check that the final order is as desired 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check that there are no zero sizes
zerosize = find(legsizes==0);
sizezerosize = size(zerosize);

if sizezerosize(2) ~= 0
    error('tsplit:zeroleg', 'legsize cannot contain a zero');
end

%get size of tensor
sizeA = size(A);
prodsizes = prod(legsizes);

%check that the product of legsizes equals the size of the leg to be split
if prodsizes ~= sizeA(leg)
    error('tsplit:sizeleg', 'product of legsizes must match the size of the leg to be split');
end

%flip the legsize vector to be in the same convention as MATLAB so can
%convert back properly at the end
legsizes = fliplr(legsizes);
%get final leg sizes
sizeB = cat(2,sizeA(1:(leg-1)),legsizes,sizeA((leg+1):end));

%reshape to new form
B = zeros(sizeB);
B = reshape(A,sizeB);

%permute to my index convention (left hand index slower)
noidx = size(sizeB);
nolegsizes = size(legsizes);
Bperm = [1:noidx(2)];
Bperm(leg:(leg+nolegsizes(2)-1)) = fliplr(Bperm(leg:(leg+nolegsizes(2)-1)));
B = permute(B,Bperm);
