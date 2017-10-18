function [ee,n_a,chi] = TTNeeDM_PBC(L,w,blocks,Jorder,tL,tR)
%[ee,n_a,chi] = TTNeeDM2_PBC(L,w,blocks,Jorder,tL,tR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entanglement entropy by reduced density matrix
%
% input: length (L), isometries (w), vector showing which sites are in which
% block (blocks), order of contractions (Jorder), tensors below each leg (tL &
% tR)
% output: entanglement entropy (ee), number of bonds connecting the two
% blocks (n_A), product of leg sizes connecting the two blocks (chi)
%
% Andrew Goldsborough - 02/02/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate n_a and chi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find which way is best (trace over A or B)
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
    
    %rename A <-> B
    blocks = -blocks + 2;
    
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

%fill block L and R
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

%start with top tensor, trace over any zeros
if bbelow(1) == 0 && bbelow(2) ~= 0
    rdm = ncon({w{L-1,1},conj(w{L-1,1})},{[1,2,-1],[1,2,-2]});
    
    %update the index arrays etc
    bbelow(1) = [];
    tbelow(1) = [];
    
elseif bbelow(2) == 0 && bbelow(1) ~= 0
    rdm = ncon({w{L-1,1},conj(w{L-1,1})},{[1,-1,2],[1,-2,2]});
    
    %update the index arrays etc
    bbelow(2) = [];
    tbelow(2) = [];

elseif bbelow(1) == 0 && bbelow(2) == 0
    rdm = ncon({w{L-1,1},conj(w{L-1,1})},{[1,2,3],[1,2,3]});
    
    %update the index arrays etc
    bbelow = [];
    tbelow = [];
    
elseif bbelow(1) ~= 0 && bbelow(2) ~= 0
    rdm = ncon({w{L-1,1},conj(w{L-1,1})},{[1,-1,-2],[1,-3,-4]});
end

complete = 1;

i=1;

while bones ~= 0
    
    while tbelow(i) ~=0 && bbelow(i) ==1
        
        %find the number of downward legs
        nolegsdown = size(size(rdm),2)/2;
        
        %contract rdm with the tensor below
        contract1 = -[1:(i-1),-1,i+2:(2*nolegsdown+1)];
        contract2 = -[-1,i,i+1];
        rdm = ncon({rdm,w{tbelow(i),1}},{contract1,contract2});
        
        %check for zeros
        if blockL(tbelow(i)) == 0
            contract1 = -[1:(i-1),-1,i:(nolegsdown+i-1),-2,(nolegsdown+1+i):(2*nolegsdown)];
            contract2 = -[-2,-1,(nolegsdown+i)];
            
            %perform contraction
            rdm = ncon({rdm,conj(w{tbelow(i),1})},{contract1,contract2});
            
            %update the index arrays etc
            bbelow = [bbelow(1:(i-1)),blockR(tbelow(i)),bbelow((i+1):end)];
            bones = bones - 1 + size(find(blockR(tbelow(i))==1),2);
            tbelow = [tbelow(1:(i-1)),tR(tbelow(i)),tbelow((i+1):end)];
        elseif blockR(tbelow(i)) == 0
            contract1 = -[1:i,-1,i+1:(nolegsdown+i-1),-2,(nolegsdown+1+i):(2*nolegsdown)];
            contract2 = -[-2,(nolegsdown+i),-1];
            
            %perform contraction
            rdm = ncon({rdm,conj(w{tbelow(i),1})},{contract1,contract2});
            
            %update the index arrays etc
            bbelow = [bbelow(1:(i-1)),blockL(tbelow(i)),bbelow((i+1):end)];
            bones = bones - 1 + size(find(blockL(tbelow(i))==1),2);
            tbelow = [tbelow(1:(i-1)),tL(tbelow(i)),tbelow((i+1):end)];
        else
            contract1 = -[1:(nolegsdown+i),-1,(nolegsdown+3+i):(2*nolegsdown+2)];
            contract2 = -[-1,(nolegsdown+1+i),(nolegsdown+2+i)];
            
            %perform contraction
            rdm = ncon({rdm,conj(w{tbelow(i),1})},{contract1,contract2});
            
            %update the index arrays etc
            bbelow = [bbelow(1:(i-1)),blockL(tbelow(i)),blockR(tbelow(i)),bbelow((i+1):end)];
            bones = bones - 1 + size(find([blockL(tbelow(i)),blockR(tbelow(i))]==1),2);
            tbelow = [tbelow(1:(i-1)),tL(tbelow(i)),tR(tbelow(i)),tbelow((i+1):end)];
        end
    end
    i = i+1;
end

if complete == 1
    
    %fuse to get a matrix
    rdm = tfuse(rdm,[ones(1,size(bbelow,2)),-(2:size(bbelow,2)+1)]);
    rdm = tfuse(rdm,[-1,2*ones(1,size(bbelow,2))]);
    
    ee = eentropydm(rdm);
end
