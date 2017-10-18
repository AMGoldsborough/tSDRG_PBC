function tSDRG_PBC(L,Jstr,Jdis,Jz,chi,Pdist,Jseed)
%tSDRG_PBC(L,Jstr,Jdis,Jz,chi,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tSDRG with Periodic Boundary Conditions
% based on 10.1103/PhysRevB.89.214203 by Goldsborough and Roemer
%
% Andrew Goldsborough - 10/11/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%inputs
% L = 30;         %chain length
% Jz = 1;         %anisotropy
% Jstr = 1;       %overall J strength
% Jdis = 2;       %disorder strength
% Jseed = 1;      %seed for rng, 0 => shuffle
% chi = 4;        %max chi

%coupling distribution
% Pdist = 4;
%0 => manual
%1 => 2 theta(K-1/2)
%2 => 1
%3 => uniform around Jstr normalised by Jstr
%4 => uniform around Jstr un-normalised Jdis
%5 => box distribution of Hikihara AF (10.1103/PhysRevB.60.12116)
%6 => Laflorencie's infinite disorder distribution (10.1103/PhysRevB.72.140408)

%when compiled the command line inputs are strings, convert to numbers
if ischar(L)==1
  L = str2double(L);
end
if ischar(Jstr)==1
  Jstr = str2double(Jstr);
end
if ischar(Jdis)==1
  Jdis = str2double(Jdis);
end
if ischar(Jz)==1
  Jz = str2double(Jz);
end
if ischar(chi)==1
  chi = str2double(chi);
end
if ischar(Pdist)==1
  Pdist = str2double(Pdist);
end
if ischar(Jseed)==1
  Jseed = str2double(Jseed);
end

%turn off warnings
warning('off','MATLAB:eigs:SigmaChangedToSA');
warning('off','ncon:suboptimalsequence');

%storage for tensors
w = cell(L-1,1);

%create array for the order
Jorder = zeros((L-1),1);

%create an array for the tensor on the left and right leg of each
tR = zeros((L-1),1);
tL = zeros((L-1),1);

%create an array to find tL and tR (curt = current tensor below)
curt = zeros(L,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate couplings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set seed, set by clock if zero
if Jseed > 0
    rng(Jseed);
else
    rng('shuffle');
end

%set the probability distribution
if Pdist==0
    %custom set
    J = [0.1,Jdis,0.1,1,Jdis,1,0.1,Jdis,0.1,0.01];
elseif Pdist==1
    %P(K) = 2 theta(K-1/2)
    J = zeros(1,L) + Jstr*random('unif',0.5,1,[1,L]);
elseif Pdist==2
    %P(K) = 1
    J = zeros(1,L) + Jstr*(rand(L,1));
elseif Pdist==3
    %uniform around Jstr normalised by Jstr
    J = zeros(1,L) + Jstr + Jstr*Jdis*(rand(1,L) - 0.5);
elseif Pdist==4
    %uniform around Jstr un-normalised Jdis
    J = zeros(1,L) + Jstr + Jdis*(rand(1,L) - 0.5);
elseif Pdist==5
    %box distribution of Hikihara AF
    J = zeros(1,L) + Jdis*random('unif',0,1,[1,L]);
elseif Pdist==6
    %Laflorencie's infinite disorder distribution
    J = rand(1,L).^Jdis;
end

rng('shuffle');
h = zeros(L,1);

%print J to file
fprintf('printing interaction strengths\n');
fname = strcat('./J/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_J.txt');
fidJ = fopen(fname, 'w');
for i = 1:L
    fprintf(fidJ,'%.15e\n',J(i));
end
fclose(fidJ);

%normalise J to 1 (for stability)
J_norm = max(J);
J = J./J_norm;

%import hamiltonian MPOs "W" and spin dimension d
[W,~] = heishamhalfSD_PBC(L,J,Jz,h);

%set J to be the gap
if Jz > 1
    %gap = 1
elseif Jz >= 0
    J = J.*((1+Jz)/2);
else
    error('only positive Jz currently implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(L-2)
    
    %find the maximum J
    [~,Jmaxpos] = max(abs(J));
    
    %store Jmaxpos
    Jorder(i) = Jmaxpos;
    
    %save tL and tR
    tL(i,1) = curt(Jmaxpos);
    tR(i,1) = curt(PBC_pos(Jmaxpos+1,L-(i-1)));
    
    %update curt
    curt(Jmaxpos) = i;
    curt(PBC_pos(Jmaxpos+1,L-(i-1))) = [];
    
    %block Jmaxpos and Jmaxpos+1
    block = ncon({W{Jmaxpos}(:,:,1,:),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]});
    
    %create a tensor for course graining ansatz
    leg2 = size(block,1);
    leg3 = size(block,3);
    
    %fuse to create matrix
    block = tfuse(block,[1,-2,1,-3]);
    block = tfuse(block,[-1,2,2]);
    
    %force hermiticity
    block = 0.5*(block+block');
    
    sizeblock = size(block,1);
    
    [v,nrg] = eig(block);
    
    if chi < (sizeblock)
        %take only full blocks
        nrgvec = sort(real(diag(nrg)),'ascend');
        diffnrg = diff(nrgvec);
        diffloc = find(diffnrg>=1E-12);
        ltm = find(diffloc<=chi);
        blockend = diffloc(ltm(end));
        nrg = nrg(1:blockend,1:blockend);
        v = v(:,1:blockend);
    end
    
    %split set of vectors into isometry
    w{i,1} = tsplit((v'),2,[leg2,leg3]);
    
    %find sizes of matrices
    sizenew = size(nrg,1);

    %contract isometries with full two site Hamiltonian
    %faster to do each element individually than explicitly contract
    blockfull = zeros(sizenew,sizenew,5,5);
    blockfull(:,:,1,1) = eye(sizenew);
    blockfull(:,:,1,2) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,1,2)))*v;
    blockfull(:,:,1,3) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,1,3)))*v;
    blockfull(:,:,1,4) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,1,4)))*v;
    blockfull(:,:,1,5) = nrg;
    blockfull(:,:,2,5) = v'*(kron(W{Jmaxpos}(:,:,2,5),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,5,5)))*v;
    blockfull(:,:,3,5) = v'*(kron(W{Jmaxpos}(:,:,3,5),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,5,5)))*v;
    blockfull(:,:,4,5) = v'*(kron(W{Jmaxpos}(:,:,4,5),W{PBC_pos(Jmaxpos+1,L-(i-1))}(:,:,5,5)))*v;
    blockfull(:,:,5,5) = eye(sizenew);
    
    %set this new block as ham for site Jmaxpos
    W{Jmaxpos} = blockfull;
    
    %remove combined site
    W(PBC_pos(Jmaxpos+1,L-(i-1))) = [];
    J(Jmaxpos) = [];

    %need to get the new gaps for each side of the block
    if i ~= L-2
        %LHS of block
        
        % H' = HB + HC + HB
        if Jmaxpos == L-(i-1)
            %PBC term
            Hc = ncon({W{PBC_pos(Jmaxpos-2,L-(i-1)-1)}(:,:,1,:),W{PBC_pos(Jmaxpos-1,L-(i-1)-1)}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]});
        else
            Hc = ncon({W{PBC_pos(Jmaxpos-1,L-(i-1)-1)}(:,:,1,:),W{PBC_pos(Jmaxpos,L-(i-1)-1)}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]});
        end
        
        %fuse to create matrix
        Hc = tfuse(Hc,[1,-2,1,-3]);
        Hc = tfuse(Hc,[-1,2,2]);
        
        %force hermiticity
        Hc = 0.5*(Hc+Hc');
        gap = eig(Hc);
        gap = sort(gap,'ascend');
        diffgap = diff(gap);
        
        %heighest gap in chi
        if size(gap,1) <= chi
            loc = find(abs(diffgap)>=1E-12);
        else
            loc = find(abs(diffgap(1:chi))>=1E-12);
        end
        
        if size(loc,1) == 0
            %no gap found, exit
            error('SDRGgaptestSU2:nogap','no gap detected');
        else
            gaploc = loc(end);
            
            J(PBC_pos(Jmaxpos-1,L-(i-1)-1)) = diffgap(gaploc);
        end
        
        %RHS of block
        
        % H' = HB + HC + HB
        if Jmaxpos == L-(i-1)
            %PBC term
            Hc = ncon({W{PBC_pos(Jmaxpos-1,L-(i-1)-1)}(:,:,1,:),W{PBC_pos(Jmaxpos,L-(i-1)-1)}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]});
        else
            Hc = ncon({W{PBC_pos(Jmaxpos,L-(i-1)-1)}(:,:,1,:),W{PBC_pos(Jmaxpos+1,L-(i-1)-1)}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]});
        end
        
        %fuse to create matrix
        Hc = tfuse(Hc,[1,-2,1,-3]);
        Hc = tfuse(Hc,[-1,2,2]);
        
        %force hermiticity
        Hc = 0.5*(Hc+Hc');
        gap = eig(Hc);
        gap = sort(gap,'ascend');
        diffgap = diff(gap);
        
        %heighest gap in chi
        if size(gap,1) <= chi
            loc = find(abs(diffgap)>=1E-12);
        else
            loc = find(abs(diffgap(1:chi))>=1E-12);
        end
        
        if size(loc,1) == 0
            %no gap found, exit
            error('SDRGgaptestSU2:nogap','no gap detected');
        else
            gaploc = loc(end);
            
            J(PBC_pos(Jmaxpos+1,L-(i-1)-1)) = diffgap(gaploc);
        end
    end
end

%final step builds the full system block
i = L-1;
Jmaxpos = 1;
Jorder(i) = Jmaxpos;

%save tL and tR
tL(i,1) = curt(Jmaxpos);
tR(i,1) = curt(Jmaxpos+1);

%block Jmaxpos and Jmaxpos+1
block = ncon({W{Jmaxpos}(:,:,1,:),W{Jmaxpos+1}(:,:,:,end)},{[-1,-2,-5,1],[-3,-4,1,-6]})...
    + permute(ncon({W{Jmaxpos+1}(:,:,1,2:end-1),W{Jmaxpos}(:,:,2:end-1,end)},{[-1,-2,-5,1],[-3,-4,1,-6]}),[3,4,1,2]);

leg2 = size(block,1);
leg3 = size(block,3);

block = tfuse(block,[1,-2,1,-3]);
block = tfuse(block,[-1,2,2]);

%force hermiticity
block = 0.5*(block+block');
[v,nrg] = eig(block);

%take ground state as the top tensor
v = v(:,1);

%split into top tensor
w{i,1} = tsplit(v',2,[leg2,leg3]);
fprintf('gs energy is %.15f\n',nrg(1,1)*J_norm);

%print energy to file
fprintf('printing energy\n');
fname = strcat('./energy/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_energy_tSDRG_PBC.txt');
fidenergy = fopen(fname, 'w');
fprintf(fidenergy,'%.15e\n',nrg(1,1)*J_norm);
fclose(fidenergy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sz correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open file to print to
fprintf('printing Sz correlation functions\n');
fname = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_tSDRG_PBC.txt');
fidcorr = fopen(fname, 'w');

%make operators (SzSz)
Oi = 0.5*[1 0;0 -1];
Oj = 0.5*[1 0;0 -1];

for Si = 1:(L-1)
    for Sj = (Si+1):L
        [corr,~] = TTNcorr_PBC(L,w,Jorder,Oi,Si,Oj,Sj);
        fprintf(fidcorr,'%d %d %.15e\n',Si,Sj,corr);
    end
end
fclose(fidcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SpSm correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open file to print to
fprintf('printing SpSm correlation functions\n');
fname = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_SpSmcorr_tSDRG_PBC.txt');
fidSpSmcorr = fopen(fname, 'w');

%make operators (SpSm)
Oi = [0 1;0 0];
Oj = [0 0;1 0];

for Si = 1:(L-1)
    for Sj = (Si+1):L
        [corr,~] = TTNcorr_PBC(L,w,Jorder,Oi,Si,Oj,Sj);
        fprintf(fidSpSmcorr,'%d %d %.15e\n',Si,Sj,corr);
    end
end
fclose(fidSpSmcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SmSp correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open file to print to
fprintf('printing SmSp correlation functions\n');
fname = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_SmSpcorr_tSDRG_PBC.txt');
fidSmSpcorr = fopen(fname, 'w');

%make operators (SmSp)
Oi = [0 0;1 0];
Oj = [0 1;0 0];

for Si = 1:(L-1)
    for Sj = (Si+1):L
        [corr,~] = TTNcorr_PBC(L,w,Jorder,Oi,Si,Oj,Sj);
        fprintf(fidSmSpcorr,'%d %d %.15e\n',Si,Sj,corr);
    end
end
fclose(fidSmSpcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%entanglement entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if L <= 30
    
    %print to file
    fprintf('printing ee\n');
    fname = strcat('./ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_ee_tSDRG_PBC.txt');
    fidee = fopen(fname, 'w');
    
    %loop over the bipartition positions (splitting chain into two)
    % AAAA/BB/AAAA
    %output format [spos ee n_a chi]
    %where spos = size of central block, ee = entanglement entropy, n_a =
    %number of legs connecting the two blocks, chi = product of the leg sizes

    for i=1:L-1
        for j=i+1:L
            %create indices for which parts are A(2) and which are B(0)
            blocks = [zeros(i,1);2*ones(j-i,1);zeros(L-j,1)];
            
            try
                %ee from SVD
                [ee,n_a,chi] = TTNeeSVD_PBC(L,w,blocks,Jorder,tL,tR);
                fprintf(fidee,'%d %d %.15e %d %d\n',i,j,ee,n_a,chi);
            catch
                try
                    %ee from density matrix
                    [ee,n_a,chi] = TTNeeDM_PBC(L,w,blocks,Jorder,tL,tR);
                    fprintf(fidee,'%d %d %.15e %d %d\n',i,j,ee,n_a,chi);
                catch 
                    fprintf('ee failed: %d %d\n',i,j);
                    fprintf(fidee,'%d %d %.15e %d %d\n',i,j,NaN,NaN,NaN);
                    continue
                end
            end
        end
    end
    fclose(fidee);
end

toc
end
