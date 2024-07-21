Initial_SWRFilter_common;
warning off

ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R23_sub'];
ROOT.Unit1 = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RegionList = {'SUB','CA1'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB_field.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1_field.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%%
for Reg = 2:numel(RegionList)
    thisRegion = RegionList{Reg};
    % thisRegion = 'SUB';
    cl=Reg;

    UTf = UnitsTable.(thisRegion);
    UTf = UnitsTable_field.(thisRegion);
    RT = RipplesTable.(thisRegion);
    ReT = ReactTable.(thisRegion);
    FRs  = FRMaps.(thisRegion);
    FRf = FRMaps.([thisRegion '_field']);

    UnitPairF.(thisRegion)=table; u=1;
    for uid=1:size(UTf,1)
        p0 = sum(RT.rat==UTf.rat(uid) & RT.session==UTf.session(uid));
        thisUnits = UTf(UTf.rat==UTf.rat(uid) & UTf.session==UTf.session(uid),:);
        thisRips = unique(ReT(find(cellfun(Params.cellfind(UTf.ID{uid}),ReT.UnitID)),[1 2]));
        p1 = size(thisRips,1);

        for uid2=1:size(thisUnits,1)
            if strncmp(thisUnits.ID{uid2},UTf.ID{uid},14), continue; end
            
            if uid~=1
                if(find(strcmp(UTf.ID{uid},UnitPairF.(thisRegion).UID2) & strcmp(thisUnits.ID{uid2},UnitPairF.(thisRegion).UID1)))
                    continue;
                end
            end
            thisRips2=unique(ReT(find(cellfun(Params.cellfind(thisUnits.ID{uid2}),ReT.UnitID)),[1 2]));
            CoActRips = intersect(thisRips.RippleID,thisRips2.RippleID);

            q1 = size(thisRips2,1);
            pq = size(CoActRips,1);

            UnitPairF.(thisRegion).UID1{u}=UTf.ID{uid};
            UnitPairF.(thisRegion).UID2{u}=thisUnits.ID{uid2};
            UnitPairF.(thisRegion).SameUnit{u} = strncmp(thisUnits.ID{uid2},UTf.ID{uid},12);
            UnitPairF.(thisRegion).p(u)=p1;
            UnitPairF.(thisRegion).q(u)=q1;
            UnitPairF.(thisRegion).pq(u)=pq;
            UnitPairF.(thisRegion).un(u) = size(union(thisRips.RippleID,thisRips2.RippleID),1);
            UnitPairF.(thisRegion).p0(u)=p0;

      
ID1 = UnitPairF.(thisRegion).UID1{u};
ID2 = UnitPairF.(thisRegion).UID2{u};


            FID =  find(strcmp(UTf.ID,UTf.ID{uid}));
    Norm_FRf_1 = FRf(1,:,FID) ./ max(FRf(1,:,FID));

    FID = find(strcmp(UTf.ID,thisUnits.ID{uid2}));
    Norm_FRf_2 = FRf(1,:,FID) ./ max(FRf(1,:,FID));

    if ~(exist([ROOT.Unit1 '\' ID1 '.mat']) & exist([ROOT.Unit1 '\' ID2 '.mat']))
            UnitPairF.(thisRegion).Sp(u) = nan;
            UnitPairF.(thisRegion).Nsp_L(u) =  nan;
            UnitPairF.(thisRegion).Nsp_R(u) =  nan;
      UnitPairF.(thisRegion).Nsp_C(u) =  nan;
        continue;
    end
    F1 = load([ROOT.Unit1 '\' ID1 '.mat']);
    F2 = load([ROOT.Unit1 '\' ID2 '.mat']);
    if ~(isfield(F1,'thisFieldMap1D_trial') & isfield(F2,'thisFieldMap1D_trial'))

    UnitPairF.(thisRegion).Sp(u) = nan;
            UnitPairF.(thisRegion).Nsp_L(u) =  nan;
            UnitPairF.(thisRegion).Nsp_R(u) =  nan;
      UnitPairF.(thisRegion).Nsp_C(u) =  nan;
        continue;
    end

     FR= F1.RDIs_field; 
     if size(FR,1)<3, FR(3,:)=nan; end
    RDI_L1 = [interp1(1:size(FR,2),FR(1,:),linspace(1,size(FR,2),40))];
    RDI_R1 = [interp1(1:size(FR,2),FR(2,:),linspace(1,size(FR,2),40))];
    RDI_C1 = [interp1(1:size(FR,2),FR(3,:),linspace(1,size(FR,2),40))];

    FR= F2.RDIs_field;
    if size(FR,1)<3, FR(3,:)=nan; end
    RDI_L2 = [interp1(1:size(FR,2),FR(1,:),linspace(1,size(FR,2),40))];
    RDI_R2 = [interp1(1:size(FR,2),FR(2,:),linspace(1,size(FR,2),40))];
    RDI_C2 = [interp1(1:size(FR,2),FR(3,:),linspace(1,size(FR,2),40))];

    RDI_L1(isnan(RDI_L1))=0; RDI_L2(isnan(RDI_L2))=0;
    RDI_R1(isnan(RDI_R1))=0; RDI_R2(isnan(RDI_R2))=0;
    RDI_C1(isnan(RDI_C1))=0; RDI_C2(isnan(RDI_C2))=0;

    RDI_L1m = nanmean(RDI_L1,1); RDI_R1m = nanmean(RDI_R1,1);  RDI_C1m = nanmean(RDI_C1,1);
    RDI_L2m = nanmean(RDI_L2,1); RDI_R2m = nanmean(RDI_R2,1);  RDI_C2m = nanmean(RDI_C2,1);

 Norm_FRf_1(isnan(Norm_FRf_1))=0;
Norm_FRf_2(isnan(Norm_FRf_2))=0;

    UnitPairF.(thisRegion).Sp(u) = corr(Norm_FRf_1',Norm_FRf_2');
            UnitPairF.(thisRegion).Nsp_L(u) = corr(RDI_L1m',RDI_L2m');
            UnitPairF.(thisRegion).Nsp_R(u) = corr(RDI_R1m',RDI_R2m');
      UnitPairF.(thisRegion).Nsp_C(u) = corr(RDI_C1m',RDI_C2m');

            UnitPairF.(thisRegion).L1(u)=UTf.RDI_LScene(uid);
            UnitPairF.(thisRegion).R1(u)=UTf.RDI_RScene(uid);
            UnitPairF.(thisRegion).C1(u)=UTf.RDI_LR(uid);

            UnitPairF.(thisRegion).L2(u)=thisUnits.RDI_LScene(uid2);
            UnitPairF.(thisRegion).R2(u)=thisUnits.RDI_RScene(uid2);
            UnitPairF.(thisRegion).C2(u)=thisUnits.RDI_LR(uid2);

            %%
            UnitPairF.(thisRegion).Lt1(u) = UTf.Selectivity_LScene(uid) * max(UnitPairF.(thisRegion).L1(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).L1(u)))]);
            UnitPairF.(thisRegion).Rt1(u) = UTf.Selectivity_RScene(uid) * max(UnitPairF.(thisRegion).R1(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).R1(u)))]);
            UnitPairF.(thisRegion).Ct1(u) = UTf.Selectivity_LR(uid) * max(UnitPairF.(thisRegion).C1(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).C1(u)))]);

            UnitPairF.(thisRegion).Lt2(u) = thisUnits.Selectivity_LScene(uid2) * max(UnitPairF.(thisRegion).L2(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).L2(u)))]);
            UnitPairF.(thisRegion).Rt2(u) = thisUnits.Selectivity_RScene(uid2) * max(UnitPairF.(thisRegion).R2(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).R2(u)))]);
            UnitPairF.(thisRegion).Ct2(u) = thisUnits.Selectivity_LR(uid2) * max(UnitPairF.(thisRegion).C2(u)) /...
                max([10^(-300),abs(max(UnitPairF.(thisRegion).C2(u)))]);


            UnitPairF.(thisRegion).region(u) = cl;
            %              UnitPair.NumField1(u) = UT.NumField(uid);
            %              UnitPair.NumField2(u) = thisUnits.NumField(uid2);

            %         UnitPair.NumField1(u) = sum(UT.AvgFR(uid)==UT.AvgFR);
            %          UnitPair.NumField2(u) = sum(thisUnits.AvgFR(uid2)==thisUnits.AvgFR);

            u=u+1;
        end
    end
                UnitPairF.(thisRegion).crp =  UnitPairF.(thisRegion).pq./UnitPairF.(thisRegion).un; 
UnitPairF.(thisRegion).crp(UnitPairF.(thisRegion).q==0 | UnitPairF.(thisRegion).q==0)=nan;
end
% writetable(UnitPair,[ROOT.Save '\UnitPair_' thisRegion '.xlsx'],'writemode','replacefile');

%%

