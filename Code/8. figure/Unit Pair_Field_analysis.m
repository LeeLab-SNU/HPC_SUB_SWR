%% Unit Pair_Field_analysis
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4_SUB_refCA1'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R25_ca1'];
ROOT.Unit1 = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];


CList = [ [207 8 23]/255;[23 84 181]/255];


RegionList = {'SUB','CA1'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%%
UnitPairF.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);
UnitPairF.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);

thisRegion = 'SUB';
UP = UnitPairF.(thisRegion); UT = UnitsTable.(thisRegion); UTF = UnitsTable_field.(thisRegion); FR = FRMaps.(thisRegion); FRF = FRMaps.([thisRegion '_field']);

for pid = 1:size(UP,1)
    ID1 = UP.UID1{pid}; id1 = find(strcmp(ID1,UT.ID));
    ID2 = UP.UID2{pid}; id2 = find(strcmp(ID2,UT.ID));

    FID = find(strcmp(UTF.ID,ID1));
    Norm_FRf_1 = FRF(1,:,FID) ./ max(FRF(1,:,FID));

    FID = find(strcmp(UTF.ID,ID2));
    Norm_FRf_2 = FRF(1,:,FID) ./ max(FRF(1,:,FID));

    F = load([ROOT.Unit1 '\' ID1 '.mat']); F= F.RDIs_field;
    RDI_L1 = [interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
    RDI_R1 = [interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
    RDI_C1 = [interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];

    F = load([ROOT.Unit1 '\' ID2 '.mat']); F= F.RDIs_field;
    RDI_L2 = [interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
    RDI_R2 = [interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
    RDI_C2 = [interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];

    RDI_L1(isnan(RDI_L1))=0; RDI_L2(isnan(RDI_L2))=0;
    RDI_R1(isnan(RDI_R1))=0; RDI_R2(isnan(RDI_R2))=0;
    RDI_C1(isnan(RDI_C1))=0; RDI_C2(isnan(RDI_C2))=0;

    RDI_L1m = nanmean(RDI_L1,1); RDI_R1m = nanmean(RDI_R1,1);  RDI_C1m = nanmean(RDI_C1,1);
    RDI_L2m = nanmean(RDI_L2,1); RDI_R2m = nanmean(RDI_R2,1);  RDI_C2m = nanmean(RDI_C2,1);

    UnitPairF.(thisRegion).Sp(pid) = corr(Norm_FRf_1',Norm_FRf_2');
            UnitPairF.(thisRegion).Nsp_L(pid) = corr(RDI_L1m',RDI_L2m');
            UnitPairF.(thisRegion).Nsp_R(pid) = corr(RDI_R1m',RDI_R2m');
      UnitPairF.(thisRegion).Nsp_C(pid) = corr(RDI_C1m',RDI_C2m');

end

    UFca1 = UnitPairF.CA1;
    UFsub = UnitPairF.SUB;

UFsub.crp = UFsub.pq./UFsub.un; UFsub.crp(UFsub.q==0 | UFsub.q==0)=nan;
UFca1.crp = UFca1.pq./UFca1.un; UFca1.crp(UFca1.q==0 | UFca1.q==0)=nan;

figure;
scatter(UFsub. ,UFsub.crp )