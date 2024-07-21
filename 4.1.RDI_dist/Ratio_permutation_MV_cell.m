Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

% parpool('local',8); % Change 4 to the number of workers you want
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RegionList = {'CA1','SUB'};
for reg=1:2
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];
% thisRegion0 = 'CA1';
% thisRegion = 'CA1';
% thisRegion2 = 'CA1_field';

% thisRegion0 = 'SUB';
% thisRegion = 'SUB';
% thisRegion2 = 'SUB_field';


% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
ReactTable_A = readtable([ROOT.Save '\ReactTable_' thisRegion0 '_' thisRegion '.xlsx']);
ReactTable_B = readtable([ROOT.Save '\ReactTable_' thisRegion0 '_' thisRegion2 '.xlsx']);


UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=2000;
%%
% UnitsTable_A =UnitsTable;
% sz = size(UnitsTable_A,2);
% sa=29; sb=31;
% col_name_A = string(UnitsTable_A.Properties.VariableNames(sa:sz));
% UnitsTable_B(1,sb:sz+2) = array2table(nan(1,sz-sa+1));
% UnitsTable_B = renamevars(UnitsTable_B, UnitsTable_B.Properties.VariableNames(sb:sz+2),col_name_A);
% for uid=1:size(UnitsTable_B,1)
%     tar = find(strncmp(UnitsTable_A.ID,UnitsTable_B.ID(uid),12));
%     if isempty(tar)
%         UnitsTable_B(uid,sb)={nan};
%         continue;
%     end
%     UnitsTable_B(uid,sb:sz+2) = UnitsTable_A(tar,sa:sz);
%     UnitsTable_B.NumField(uid) = str2double(UnitsTable_B.ID{uid}(end));
% end
% %%
% UnitsTable_B(isnan(UnitsTable_B.NumField),:)=[];
%    writetable(UnitsTable_B,[ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx'],'writemode','replacefile');
UnitsTable = UnitsTable_B;

%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for rid=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{rid};
    thisReact_A = ReactTable_A(strcmp(ReactTable_A.RippleID,RipID),:);
    thisReact_A = unique(thisReact_A(:,1:2));

    thisReact_B = ReactTable_B(strcmp(ReactTable_B.RippleID,RipID),:);
    thisReact_B = unique(thisReact_B(:,1:2));

    RipplesTable_p.nPCs(rid) = size(thisReact_A,1);
%     if RipplesTable_p.nPCs(rid)<3, continue; end
    RipplesTable_p.nFields(rid) = size(thisReact_B,1);

    thisUnits=table;
    for u=1:size(thisReact_B,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact_B.UnitID(u)),UnitsTable.ID)),:)];
    end


%  RipplesTable_p.nFields(rid) = size(thisUnits,1);
 
    if size(thisUnits,1)<3

        RipplesTable_p.pRDI_L(rid)=nan;
        RipplesTable_p.pRDI_R(rid)=nan;
        RipplesTable_p.pRDI_C(rid)=nan;

        RipplesTable_p.pRDI_L_UV(rid)=nan;
        RipplesTable_p.pRDI_R_UV(rid)=nan;
        RipplesTable_p.pRDI_C_UV(rid)=nan;
        continue;
    end



     thisPool_A = UnitsTable_A;
     thisPool_B = UnitsTable_B;
    


%     thisUnits_Sig_L = thisUnits(thisUnits.Selectivity_LScene>=1  & abs(thisUnits.RDI_LScene)>0.1,:);
%     thisUnits_Sig_R = thisUnits(thisUnits.Selectivity_RScene>=1 & abs(thisUnits.RDI_RScene)>0.1,:);
%     thisUnits_Sig_C = thisUnits(thisUnits.Selectivity_LR>=1  & abs(thisUnits.RDI_LR)>0.1,:);
% 
%     thisPool_A_Sig_L =  thisPool_A(thisPool_A.Selectivity_LScene>=1 ,:);
%     thisPool_A_Sig_R = thisPool_A(thisPool_A.Selectivity_RScene>=1 ,:);
%     thisPool_A_Sig_C =  thisPool_A(thisPool_A.Selectivity_LR>=1 ,:);


    thisUnits_Sig_L = thisUnits;
    thisUnits_Sig_R = thisUnits;
    thisUnits_Sig_C = thisUnits;

    thisPool_A_Sig_L =  thisPool_A;
    thisPool_A_Sig_R = thisPool_A;
    thisPool_A_Sig_C =  thisPool_A;




    Ratio_L.dist=table;
    Ratio_R.dist=table;
    Ratio_C.dist=table;

    Ratio_L.dist=table;
    Ratio_R.dist=table;
    Ratio_C.dist=table;


    temp1=nan; temp2=nan; temp3=nan; temp4=nan; temp5=nan; temp6=nan;
  

         tic
    try
        parfor i=1:randN


            [x1,x2] = rand_samp_ratio(thisPool_A_Sig_L,thisPool_B,thisUnits_Sig_L,'RDI_LScene');
            temp1(i) = x1; temp2(i) = x2;
            [x1,x2] = rand_samp_ratio(thisPool_A_Sig_R,thisPool_B,thisUnits_Sig_R,'RDI_RScene');
            temp3(i) = x1; temp4(i) = x2;
            [x1,x2] = rand_samp_ratio(thisPool_A_Sig_C,thisPool_B,thisUnits_Sig_C,'RDI_LR');
            temp5(i) = x1; temp6(i) = x2;
        end
    catch
        disp([RipID ' perm fail'])
    end
         toc



Ratio_L.dist.bias = temp1; Ratio_L.dist.posneg = temp2;
Ratio_R.dist.bias = temp3; Ratio_R.dist.posneg = temp4;
Ratio_C.dist.bias = temp5; Ratio_C.dist.posneg = temp6;
%%
[x1,x2] = samp_ratio(thisUnits_Sig_L,'RDI_LScene');
Ratio_L.act_bias = x1; Ratio_L.act_posneg = x2;
[x1,x2] = samp_ratio(thisUnits_Sig_R,'RDI_RScene');
Ratio_R.act_bias = x1; Ratio_R.act_posneg = x2;
[x1,x2] = samp_ratio(thisUnits_Sig_C,'RDI_LR');
Ratio_C.act_bias = x1; Ratio_C.act_posneg = x2;
%%

    Ratio_L.p_bias = min(sum(Ratio_L.act_bias>Ratio_L.dist.bias),sum(Ratio_L.act_bias<Ratio_L.dist.bias))/sum(~isnan(Ratio_L.dist.bias));
    Ratio_R.p_bias = min(sum(Ratio_R.act_bias>Ratio_R.dist.bias),sum(Ratio_R.act_bias<Ratio_R.dist.bias))/sum(~isnan(Ratio_R.dist.bias));
    Ratio_C.p_bias = min(sum(Ratio_C.act_bias>Ratio_C.dist.bias),sum(Ratio_C.act_bias<Ratio_C.dist.bias))/sum(~isnan(Ratio_C.dist.bias));
    RipplesTable_p.pBias_L_UV(rid) = Ratio_L.p_bias;
    RipplesTable_p.pBias_R_UV(rid) = Ratio_R.p_bias;
    RipplesTable_p.pBias_C_UV(rid) = Ratio_C.p_bias;
    RipplesTable_p.mBias_L_UV(rid) = Ratio_L.act_bias;
    RipplesTable_p.mBias_R_UV(rid) = Ratio_R.act_bias;
    RipplesTable_p.mBias_C_UV(rid) = Ratio_C.act_bias;


        Ratio_L.p_posneg = min(sum(Ratio_L.act_posneg>Ratio_L.dist.posneg),sum(Ratio_L.act_posneg<Ratio_L.dist.posneg))/sum(~isnan(Ratio_L.dist.posneg));
    Ratio_R.p_posneg = min(sum(Ratio_R.act_posneg>Ratio_R.dist.posneg),sum(Ratio_R.act_posneg<Ratio_R.dist.posneg))/sum(~isnan(Ratio_R.dist.posneg));
    Ratio_C.p_posneg = min(sum(Ratio_C.act_posneg>Ratio_C.dist.posneg),sum(Ratio_C.act_posneg<Ratio_C.dist.posneg))/sum(~isnan(Ratio_C.dist.posneg));
    RipplesTable_p.pPosneg_L_UV(rid) = Ratio_L.p_posneg;
    RipplesTable_p.pPosneg_R_UV(rid) = Ratio_R.p_posneg;
    RipplesTable_p.pPosneg_C_UV(rid) = Ratio_C.p_posneg;
        RipplesTable_p.mPosneg_L_UV(rid) = Ratio_L.act_posneg;
    RipplesTable_p.mPosneg_R_UV(rid) = Ratio_R.act_posneg;
    RipplesTable_p.mPosneg_C_UV(rid) = Ratio_C.act_posneg;

    if exist([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']])
    save([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']],'thisUnits','Ratio_L','Ratio_R','Ratio_C','-append')
    else
        save([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']],'thisUnits','Ratio_L','Ratio_R','Ratio_C')
    end
    disp([RipID ' is finished!'])




    %     if size(thisUnits_UV_SC,1)<5, RipplesTable_p.pRDI_L_UV_scuv(r)=nan; RipplesTable_p.pRDI_R_UV_scuv(r)=nan; RipplesTable_p.pRDI_C_UV_scuv(r)=nan; end

end
% RipplesTable_p = RipplesTable_p(RipplesTable_p.nFields>=5,:);
suff = '_forAnalysis';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion0 suff '.xlsx'],'writemode','replacefile')
delete(gcp('nocreate'));

RipplesTable_p.pRDI_M_UV = nanmin([RipplesTable_p.pRDI_L_UV,RipplesTable_p.pRDI_R_UV,RipplesTable_p.pRDI_C_UV],[],2);
sum(RipplesTable_p.pRDI_M_UV<0.05)
sum(RipplesTable_p.DecodingP_all<0.05)
end
%%


%%



function [x1,x2] = rand_samp(thisPool_A,thisPool_B,thisUnits,rdi)
samp = datasample(thisPool_A,size(thisUnits,1),1);
if ~isempty(samp)
    samp2=table;
    for j=1:size(thisPool_B,1)
        if max(strcmp(thisPool_B.ID{j}(1:12),samp.ID)) & abs(thisPool_B.(rdi)(j))>0.1
            samp2=[samp2;thisPool_B(j,:)];
        end
    end
    x1 = nanmean(samp2.(rdi)); x2 = nanmedian(samp2.(rdi));
else
    x1 = nan; x2=nan;
end
end

function [x1] = rand_samp_snr(thisPool_A,thisPool_B,thisUnits,rdi)
samp = datasample(thisPool_A,size(thisUnits,1),1);
if ~isempty(samp)
    samp2=table;
    for j=1:size(thisPool_B,1)
        if max(strcmp(thisPool_B.ID{j}(1:12),samp.ID)) & abs(thisPool_B.(rdi)(j))>0.1
            samp2=[samp2;thisPool_B(j,:)];
        end
    end
    x1 = cal_snr(samp2.(rdi));
else
    x1 = nan; 
end
end

function [x1,x2] = rand_samp_ratio(thisPool_A,thisPool_B,thisUnits,rdi)
samp = datasample(thisPool_A,size(thisUnits,1),1);
if ~isempty(samp)
    samp2=table;
    for j=1:size(thisPool_B,1)
        if max(strcmp(thisPool_B.ID{j}(1:12),samp.ID)) & abs(thisPool_B.(rdi)(j))>0.1
            samp2=[samp2;thisPool_B(j,:)];
        end
    end
    if ~isempty(samp2)
    pos = sum(samp2.(rdi)>=0.1); neg = sum(samp2.(rdi)<=-0.1);
    x1 = (pos-neg)/(pos+neg); x2 = (pos-neg);
    else
        x1 = nan; x2=nan;
    end
else
    x1 = nan; x2=nan;
end
end

function [x1,x2] = samp_ratio(samp2,rdi)

    if ~isempty(samp2)
    pos = sum(samp2.(rdi)>=0.1); neg = sum(samp2.(rdi)<=-0.1);
    x1 = (pos-neg)/(pos+neg); x2 = (pos-neg);
    else
        x1 = nan; x2=nan;
    end

end