function make_Spcorr_tSDRG(L,Jstr,Jdis,Jz,chi,Pdist,Jseed)
%make_Spcorr_tSDRG(L,Jstr,Jdis,Jz,chi,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_Spcorr_tSDRG
% function to make spin corr from components
% S.S = SzSz + 0.5*(SpSm + SmSp)
%
% Andrew Goldsborough - 08/12/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open files to read in data
fnameSz = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_tSDRG_PBC.txt');
Spcorr = importdata(fnameSz);

%check length
if size(Spcorr,1) == 0.5*L*(L-1)
    
    fnameSpSm = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_SpSmcorr_tSDRG_PBC.txt');
    corr = importdata(fnameSpSm);
    Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);
    
    fnameSmSp = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_SmSpcorr_tSDRG_PBC.txt');
    corr = importdata(fnameSmSp);
    Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);
    
    %open files to write to
    fnameSp = strcat('./Spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_Spcorr_tSDRG_PBC.txt');
    fprintf(strcat(fnameSp,'\n'));
    fidSpcorr = fopen(fnameSp, 'w');
    
    %print to file
    for i=1:size(Spcorr,1)
        fprintf(fidSpcorr,'%d %d %.15e\n',Spcorr(i,1),Spcorr(i,2),Spcorr(i,3));
    end
    
    %close file
    fclose(fidSpcorr);
    
else
    fprintf(fnameSz);
    fprintf(' : file not full\n');
end

