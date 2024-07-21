% Add_RDI_only

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R1'];
ROOT.Unit = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'SUB_field';

UnitsTable_all =table;

Cluster_List = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
exper = {'LSM'};
% exper2 = {'JS','LSM'};
% exper = {'SEB'};
TargRegion = 'SUB_field';
% gcp = parpool(4);

for cid=1:size(Cluster_List,1)
    id = Cluster_List.ID{cid};
    if ~ismember(Cluster_List.experimenter{cid},exper), continue; end

    thisSID = id(1:6);

    Recording_region_TT = Recording_region({thisSID},:);

    TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));

    if strcmp(thisRegion,'SUB_field') | strcmp(thisRegion,'CA1_field')
        Spike = LoadSpikeData_field(ROOT, thisSID, TargetTT,Params.cellfindn);
    else
        Spike = LoadSpikeData(ROOT, thisSID, TargetTT,Params.cellfindn);
    end
    load([ROOT.Behav '\' thisSID '.mat'])

    UnitsTable = MkUnitTable(ROOT,Behav,Spike,id,thisRegion,TargetTT,Params);
    CalRDI_forEachBin(id,ROOT,Behav,Spike)


end
%%
