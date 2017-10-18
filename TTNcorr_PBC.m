function [corr,tnum] = TTNcorr_PBC(L,w,Jorder,op1,site1,op2,site2)
%[corr,tnum] = TTNcorr_PBC(L,w,Jorder,op1,site1,op2,site2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-Point Correlation function for PBC tree tensor network
% input: Length L, Tree tensor network w, order in which the network is
% built Jorder, the two operators, the sites under investigation site1 < site2
% output: correlator corr, number of tensors in the geodesic tnum
% 
% Andrew Goldsborough 08/12/2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check that site1 is left of site2
if (site1 >= site2)
    fprintf('site 1 has to be to the left of site 2\n');
    
    error('TTNSzcorr:sitepos', 'incompatible site definitions');
end

%create vector showing the active sites
vec = [zeros(site1-1,1);1;zeros(site2-site1-1,1);1;zeros(L-site2,1)];

%create the operators as cells with zeros for the empty sites
coro = cell(L,1);
coro{site1} = op1;
coro{site2} = op2;

%count how many tensors have been used
tleft = 0;
tright = 0;

%calculate correlator
for i=1:L-1
    
    pos = Jorder(i);
    Lcur = L-(i-1);
        
    %check to see whether tensor is used
    if vec(pos) == 0 && vec(PBC_pos(pos+1,Lcur)) == 0
        %not used, remove site
        coro(PBC_pos(pos+1,Lcur)) = [];
        vec(PBC_pos(pos+1,Lcur)) = [];
        
    elseif vec(PBC_pos(pos+1,Lcur)) == 0
        %left leg used
        coro{pos} = ncon({coro{pos},w{i},conj(w{i})},{[1,2],[-1,1,3],[-2,2,3]});
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tleft = tleft + 1;
        
        %remove site
        coro(PBC_pos(pos+1,Lcur)) = [];
        vec(PBC_pos(pos+1,Lcur)) = [];
        
    elseif vec(pos) == 0
        %right leg used
        coro{pos} = ncon({coro{PBC_pos(pos+1,Lcur)},w{i},conj(w{i})},{[1,2],[-1,3,1],[-2,3,2]});
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tright = tright + 1;
        
        %remove site
        coro(PBC_pos(pos+1,Lcur)) = [];
        vec(PBC_pos(pos+1,Lcur)) = [];
        
    else
        %both legs used
        coro{pos} = ncon({coro{pos},coro{PBC_pos(pos+1,Lcur)},w{i},conj(w{i})},{[1,3],[2,4],[-1,1,2],[-2,3,4]});
        
        %propagate vec
        vec(pos) = 1;
        
        %output the number of tensors before the path is connected
        tnum = tleft + tright + 1;
        
        %remove site
        coro(PBC_pos(pos+1,Lcur)) = [];
        vec(PBC_pos(pos+1,Lcur)) = [];
    end 
end

%in the end corr is the value of the 1 element cell
corr = coro{1};