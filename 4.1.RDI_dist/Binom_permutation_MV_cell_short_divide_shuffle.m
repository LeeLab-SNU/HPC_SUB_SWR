% function [RipplesTable_A, prop] = Binom_permutation_shuffle(ROOT,RipplesTable,ReactTable)

addpath('D:\HPC-SWR project\Analysis Program')
Initial_SWRFilter_common;
warning off
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


% thisRegion0 = 'CA1';
% thisRegion = 'CA1';
% thisRegion2 = 'CA1_field';
% 
thisRegion0 = 'SUB';
thisRegion = 'SUB';
thisRegion2 = 'SUB_field';


% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
ReactTableA = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactTableB = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);


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

UnitsTable = UnitsTable_B;
UnitsTableA = UnitsTable_A;

%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
ratio_shuffle=[];
n_rep=2000;
parfor nr=1:n_rep
    shuffle_L=[]; shuffle_R=[]; shuffle_C=[];
    for rid=1:size(RipplesTable_p,1)
        RipID = RipplesTable_p.ID{rid};
        thisReact_A = ReactTableA(strcmp(ReactTableA.RippleID,RipID),:);
        thisReact_A = unique(thisReact_A(:,1:2));

        thisReact_B = ReactTableB(strcmp(ReactTableB.RippleID,RipID),:);
        thisReact_B = unique(thisReact_B(:,1:2));

        for uid=1:size(thisReact_A,1)
            try
                idx=find(strcmp(thisReact_A.UnitID{uid},UnitsTable_A.ID));
                thisReact_A.sel_L(uid)=UnitsTable_A.Selectivity_LScene(idx);
                thisReact_A.sel_R(uid)=UnitsTable_A.Selectivity_RScene(idx);
                thisReact_A.sel_C(uid)=UnitsTable_A.Selectivity_LR(idx);
            catch
            end
        end

        for uid=1:size(thisReact_B,1)
            try
                idx=find(strcmp(thisReact_B.UnitID{uid},UnitsTable_B.ID));
                thisReact_B.sel_L(uid)=UnitsTable_B.RDI_LScene(idx);
                thisReact_B.sel_R(uid)=UnitsTable_B.RDI_RScene(idx);
                thisReact_B.sel_C(uid)=UnitsTable_B.RDI_LR(idx);
            catch
            end
        end


        try
            %     if exist([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']])

            %         load([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
            %     disp([RipID ' is finished!'])

            thisUnits=table;
            for u=1:size(thisReact_B,1)
                thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact_B.UnitID(u)),UnitsTable.ID)),:)];
            end

            thisUnitsA=table;
            for u=1:size(thisReact_A,1)
                thisUnitsA =[thisUnitsA;UnitsTableA(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UnitsTableA.ID)),:)];
            end

            %         RipplesTable_p.mRDI_L_UV(rid) = RDI_L.act_mean;
            %         RipplesTable_p.mRDI_R_UV(rid) = RDI_R.act_mean;
            %         RipplesTable_p.mRDI_C_UV(rid) = RDI_C.act_mean;
            %%

            [m,p1,p2,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisReact_A, UnitsTable_B,'RDI_LScene');
            p1=0;p2=0;
            shuffle_L(rid,1) = min([b1,b2]);

            [m,p1,p2,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisReact_A, UnitsTable_B,'RDI_RScene');
            p1=0;p2=0;
            shuffle_R(rid,1) = min([b1,b2]);

            [m,p1,p2,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisReact_A, UnitsTable_B,'RDI_LR');
            p1=0;p2=0;
            shuffle_C(rid,1) = min([b1,b2]);

        catch
            shuffle_L(rid,1) = nan;
            shuffle_R(rid,1) = nan;
            shuffle_C(rid,1) = nan;
        end

    end

    shuffle_M = min([shuffle_L,shuffle_R,shuffle_C],[],2);


    ratio_shuffle(nr,:) = [sum(shuffle_L<0.05) sum(shuffle_R<0.05) sum(shuffle_C<0.05) sum(shuffle_M<0.05)] ./size(shuffle_M,1);


end


ratio_shuffle;

save([ROOT.Save '\FieldRandomize_' thisRegion '.mat'],"ratio_shuffle")
figure
histogram(ratio_shuffle(:,4))
%%
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==20) =25;
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==-20) =-25;
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


%%
function [m,p1,p2,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisReact_A, UnitsTable_B,var)

thisUnits_p=table; thisUnits_n=table; thisUnits_r =table;
for u=1:size(thisReact_A,1)
    temp =UnitsTable_B(find(strncmp(thisReact_A.UnitID(u),UnitsTable_B.ID,12)),:);
    if ~isempty(temp)
        if ~isnan(max(temp.(var)))
            [~,m] = max(temp.(var));
            thisUnits_p =[thisUnits_p;temp(m,:)];

            [~,m] = min(temp.(var));
            thisUnits_n =[thisUnits_n;temp(m,:)];

            r = randi([1 size(temp,1)],1,1);
            thisUnits_r = [thisUnits_r;temp(r,:)];
        else
            thisUnits_r =[thisUnits_r;temp(1,:)];
            thisUnits_p =[thisUnits_p;temp(1,:)];
            thisUnits_n =[thisUnits_n;temp(1,:)];
        end
    end

end

m = median(UnitsTable_B.(var)(abs(UnitsTable_B.(var))>=0.1));
p1 = sum(UnitsTable_B.(var)>=0.1)/sum(abs(UnitsTable_B.(var))>=0.1); p2 = sum(UnitsTable_B.(var)<=-0.1)/sum(abs(UnitsTable_B.(var))>=0.1);

if sum(abs(thisUnits.(var))>=0.1)<1
    b1 = nan; b2=nan;
else
    b1=myBinomTest(sum(thisUnits_r.(var)>=0.1),sum(abs(thisUnits_r.(var))>=0.1),p1,'one');
    b2= myBinomTest(sum(thisUnits_r.(var)<=-0.1),sum(abs(thisUnits_r.(var))>=0.1),p2,'one');


end
end