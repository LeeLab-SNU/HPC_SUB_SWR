Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

% parpool('local',8); % Change 4 to the number of workers you want
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
%% 
thisRegion2 = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_field_RDIs_UV.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);

UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=5000;

RPR_table_L=nan; 
RPR_table_R=nan;
RPR_table_C=nan;
UnitsTable = UnitsTable_A;
%%
RipplesTable_p = RipplesTable;
for rid=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{rid};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable.ID)),:)];
    end


    if size(thisUnits,1)==0

        RipplesTable_p.pRDI_L(rid)=nan;
        RipplesTable_p.pRDI_R(rid)=nan;
        RipplesTable_p.pRDI_C(rid)=nan;

        RipplesTable_p.pRDI_L_scuv(rid)=nan;
        RipplesTable_p.pRDI_R_scuv(rid)=nan;
        RipplesTable_p.pRDI_C_scuv(rid)=nan;
        continue;
    end


    RipplesTable_p.nFields(rid) = size(thisUnits,1);
    if RipplesTable_p.nFields(rid)<5, continue; end




    thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(rid) & UnitsTable.session==RipplesTable_p.session(rid),:);
    thisUnits = thisPool;

    RPR_table_L(rid,1) = sum(thisUnits.Selectivity_LScene==0) / size(thisUnits,1);
    RPR_table_L(rid,2) = sum(thisUnits.Selectivity_LScene==0.5) / size(thisUnits,1);
    RPR_table_L(rid,3) = sum(thisUnits.Selectivity_LScene==1) / size(thisUnits,1);
    RPR_table_L(rid,4) = sum(thisUnits.Selectivity_LScene==2) / size(thisUnits,1);
    RPR_table_L(rid,5) = sum(thisUnits.Selectivity_LScene==3) / size(thisUnits,1);
    RPR_table_L(rid,6) = sum(thisUnits.Selectivity_LScene==4) / size(thisUnits,1);

        RPR_table_R(rid,1) = sum(thisUnits.Selectivity_RScene==0) / size(thisUnits,1);
    RPR_table_R(rid,2) = sum(thisUnits.Selectivity_RScene==0.5) / size(thisUnits,1);
    RPR_table_R(rid,3) = sum(thisUnits.Selectivity_RScene==1) / size(thisUnits,1);
    RPR_table_R(rid,4) = sum(thisUnits.Selectivity_RScene==2) / size(thisUnits,1);
    RPR_table_R(rid,5) = sum(thisUnits.Selectivity_RScene==3) / size(thisUnits,1);
    RPR_table_R(rid,6) = sum(thisUnits.Selectivity_RScene==4) / size(thisUnits,1);

        RPR_table_C(rid,1) = sum(thisUnits.Selectivity_LR==0) / size(thisUnits,1);
    RPR_table_C(rid,2) = sum(thisUnits.Selectivity_LR==0.5) / size(thisUnits,1);
    RPR_table_C(rid,3) = sum(thisUnits.Selectivity_LR==1) / size(thisUnits,1);
    RPR_table_C(rid,4) = sum(thisUnits.Selectivity_LR==2) / size(thisUnits,1);
    RPR_table_C(rid,5) = sum(thisUnits.Selectivity_LR==3) / size(thisUnits,1);
    RPR_table_C(rid,6) = sum(thisUnits.Selectivity_LR==4) / size(thisUnits,1);



end

%%
c = [mean(RPR_table_L(:,1)),mean(RPR_table_L(:,2)),mean(RPR_table_L(:,3)),...
    mean(RPR_table_L(:,4)),mean(RPR_table_L(:,5)),mean(RPR_table_L(:,6))];
c=c/sum(c);

figure
pie(c)