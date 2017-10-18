function [ee,n_a,chi] = TTNeeSVD_PBC(L,w,blocks,Jorder,tL,tR)
% [ee,n_a,chi] = TTNeeSVD_PBC(L,w,blocks,Jorder,tL,tR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entanglement entropy by SVD of the wavefunction
% input: length (L), isometries (w), vector showing which sites are in which
% block (blocks), order of contractions (Jorder), tensors below each leg (tL &
% tR)
% output: entanglement entropy (ee), number of bonds connecting the two
% blocks (n_A), product of leg sizes connecting the two blocks (chi)
% 
% Andrew Goldsborough 02/02/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate equivalent n_a and chi as if using density matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ver=1:2
    
    %Use original blocks for ver 1, inverse for ver 2
    if ver == 1
        blocks1 = blocks;
    end
    
    if ver == 2
        blocks1 = -blocks + 2;
    end
    
    %create a vector showing a 0(all B),1(A and B),2(all A) for the legs
    blockL = zeros(L-1,1);
    blockR = zeros(L-1,1);
    
    %fill block L and R
    for i=1:L-1
        blockL(i) = blocks1(Jorder(i));
        blockR(i) = blocks1(PBC_pos(Jorder(i)+1,size(blocks1,1)));
        
        if blockL(i) == 2 && blockR(i) == 2
            blocks1(Jorder(i)) = 2;
        elseif blockL(i) ~= 0 || blockR(i) ~= 0
            blocks1(Jorder(i)) = 1;
        elseif blockL(i) == 0 && blockR(i) == 0
            blocks1(Jorder(i)) = 0;
        end
        
        %remove the contracted index
        blocks1(PBC_pos(Jorder(i)+1,size(blocks1,1))) = [];
    end
    
    %need to contract down, keep a track of what is connected below
    tbelow = [tL(L-1),tR(L-1)];
    bbelow = [blockL(L-1),blockR(L-1)];
    
    %keep track of the size of the indices (# of elements)
    numele = [size(w{L-1,1},2),size(w{L-1,1},3)];
    
    %need to keep a track of the number of 1s in bbelow
    bones = size(find(bbelow==1),2);
    
    %start with top tensor, trace over any zeros
    if bbelow(1) == 0 && bbelow(2) ~= 0
        
        %update the index arrays etc
        bbelow(1) = [];
        tbelow(1) = [];
        numele(1) = [];
        
    elseif bbelow(2) == 0 && bbelow(1) ~= 0
        
        %update the index arrays etc
        bbelow(2) = [];
        tbelow(2) = [];
        numele(2) = [];
        
    elseif bbelow(1) == 0 && bbelow(2) == 0
        
        %update the index arrays etc
        bbelow = [];
        tbelow = [];
        numele = [];
        
    elseif bbelow(1) ~= 0 && bbelow(2) ~= 0
        
    end
    
    i=1;
    
    while bones ~= 0
        
        while tbelow(i) ~=0 && bbelow(i) ==1
            
            %check for zeros
            if blockL(tbelow(i)) == 0
                %update the index arrays etc
                numele(i) = size(w{tbelow(i),1},3);
                bbelow = [bbelow(1:(i-1)),blockR(tbelow(i)),bbelow((i+1):end)];
                bones = bones - 1 + size(find(blockR(tbelow(i))==1),2);
                tbelow = [tbelow(1:(i-1)),tR(tbelow(i)),tbelow((i+1):end)];
                
            elseif blockR(tbelow(i)) == 0
                %update the index arrays etc
                numele(i) = size(w{tbelow(i),1},2);
                bbelow = [bbelow(1:(i-1)),blockL(tbelow(i)),bbelow((i+1):end)];
                bones = bones - 1 + size(find(blockL(tbelow(i))==1),2);
                tbelow = [tbelow(1:(i-1)),tL(tbelow(i)),tbelow((i+1):end)];
                
            else
                %update the index arrays etc
                numele = [numele(1:(i-1)),size(w{tbelow(i),1},2),size(w{tbelow(i),1},3),numele((i+1):end)];
                bbelow = [bbelow(1:(i-1)),blockL(tbelow(i)),blockR(tbelow(i)),bbelow((i+1):end)];
                bones = bones - 1 + size(find([blockL(tbelow(i)),blockR(tbelow(i))]==1),2);
                tbelow = [tbelow(1:(i-1)),tL(tbelow(i)),tR(tbelow(i)),tbelow((i+1):end)];
                
            end
            
        end
        i = i+1;
        
    end
    
    chi(ver) = prod(numele);
    n_a(ver) = size(numele,2);
end

if chi(1) > chi(2)
    n_a = n_a(2);
    chi = chi(2);
else
    n_a = n_a(1);
    chi = chi(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate ee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a vector showing a 0(all B),1(A and B),2(all A) for the legs
blockL = zeros(L-1,1);
blockR = zeros(L-1,1);

%fill block L and R working up the tree
for i=1:L-1
    blockL(i) = blocks(Jorder(i));
    blockR(i) = blocks(PBC_pos(Jorder(i)+1,size(blocks,1)));
    
    if blockL(i) == 2 && blockR(i) == 2
        blocks(Jorder(i)) = 2;
    elseif blockL(i) ~= 0 || blockR(i) ~= 0
        blocks(Jorder(i)) = 1;
    elseif blockL(i) == 0 && blockR(i) == 0
        blocks(Jorder(i)) = 0;
    end
    
    %remove the contracted index
    blocks(PBC_pos(Jorder(i)+1,size(blocks,1))) = [];
end

%need to contract down, keep a track of what is connected below
tbelow = [tL(L-1),tR(L-1)];
bbelow = [blockL(L-1),blockR(L-1)];

%need to keep a track of the number of 1s in bbelow
bones = size(find(bbelow==1),2);

%start with psi = top tensor
psi = squeeze(w{L-1,1});

i=1;

while bones ~= 0
    
    while tbelow(i) ~=0 && bbelow(i) ~= 0 && bbelow(i) ~=2
               
        %contract psi with the tensor below
        psi = ncon({psi,w{tbelow(i),1}},{-[1:(i-1),-1,i+2:(size(size(psi),2)+1)],-[-1,i,i+1]});
        
        %update the index arrays etc
        bbelow = [bbelow(1:(i-1)),blockL(tbelow(i)),blockR(tbelow(i)),bbelow((i+1):end)];
        bones = bones - 1 + size(find([blockL(tbelow(i)),blockR(tbelow(i))]==1),2);
        tbelow = [tbelow(1:(i-1)),tL(tbelow(i)),tR(tbelow(i)),tbelow((i+1):end)];
    end
    i = i+1;
    
end

%create indices for which parts are A(2) and which are B(0)
svd_v = bbelow;
svd_v(bbelow==0) = 1;
svd_v(bbelow==2) = -(2:(sum(bbelow==2)+1));
eeM = tfuse(tfuse(psi,svd_v),[-1,2*ones(1,sum(bbelow==2))]);
[~,S,~] = svd(eeM,'econ');
ee = eentropy(S);

end