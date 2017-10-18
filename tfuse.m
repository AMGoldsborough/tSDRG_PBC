function B = tfuse(A,Aid)
%B = tfuse(A,Aid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tfuse - fuse legs of a tensor
% input: Tensor A with desired final indices Aid where positive labels will be joined.
% output: Tensor B which is A with appropriate legs fused
% 
% Andrew Goldsborough 05/09/2012
% based on tcon function
% 
% example of Aid: [4,-1,-3,-2,4] will output a tensor with legs in the
% order of the absolute value, with the two legs labeled '4' fused.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check that there are no zero indices
zeroidxA = find(Aid==0);
sizezeroidxA = size(zeroidxA);

if sizezeroidxA(2) ~= 0        
    error('tfuse:zeroidx', 'zero index detected');
end

%get size of arrays
sizeAid = size(Aid);
sizeA = size(A);

%check that the size of A matches the size of Aid
if (size(sizeA,2) ~= sizeAid(2))
    error('tfuse:sizeAid','number of entries in Aid must match dimension of A');
end

%check that the unfused indices are not repeated
if (size(find(Aid<0),2) ~= 0) && (isequal(sort(Aid(find(Aid<0))), unique(Aid(find(Aid<0)))) ~= 1)
    fprintf('unfused index numbers must be unique\n');
    error('tfuse:Aidxunique','repeated unfused index detected');
end

%sort in increasing order keeping the ordering index array
[Aid,Aidoix] = sort(Aid);

%permute indices st the indices to be fused are on the right
A = permute(A,Aidoix);

%get size of new A
sizeA = size(A);

%need the number of negative indices
negidxA = find(Aid<0);
numnegA = size(negidxA);

%permute fused indices to match my index convention
A = permute(A,[1:numnegA(2),sizeAid(2):-1:(numnegA(2)+1)]);

%need the size of the fused indices (u)
sizeu = prod(sizeA((numnegA(2)+1):end));

%if not fully flatened
if (size(find(Aid<0),2) ~= 0)
    
    %fuse indices
    sizeA2 = cat(2,sizeA(1:numnegA(2)),sizeu);
    A2 = reshape(A,sizeA2);
    
    %permute to desired form
    [~,Bidoix] = sort(abs(unique(Aid)));
    B = permute(A2,Bidoix);
else
    %fuse indices
    sizeA2 = cat(2,sizeA(1:numnegA(2)),sizeu);
    A2 = reshape(A,sizeA2,1);
    
    B = A2;
end
