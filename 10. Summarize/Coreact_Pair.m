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
RegionList = {'CA1','SUB'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);
UnitPair_field.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);
UnitPair_field.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%%
for Reg = 1:numel(RegionList)
    thisRegion = RegionList{Reg};
    % thisRegion = 'SUB';
    cl=Reg;

    UT = UnitsTable.(thisRegion);
    UTf = UnitsTable_field.(thisRegion);
    RT = RipplesTable.(thisRegion);
    ReT = ReactTable.(thisRegion);
    FR  = FRMaps.(thisRegion);
    FRf = FRMaps.([thisRegion '_field']);
    UPF = UnitPair_field.(thisRegion);

    UnitPair=table; u=1;
    for uid=1:size(UT,1)
        p0 = sum(RT.rat==UT.rat(uid) & RT.session==UT.session(uid));
        thisUnits = UT(UT.rat==UT.rat(uid) & UT.session==UT.session(uid),:);
        thisRips = unique(ReT(find(cellfun(Params.cellfind(UT.ID{uid}),ReT.UnitID)),[1 2]));
        p1 = size(thisRips,1);

        FR_u1 = FR(1,:,uid);
        FR_u1(isnan(FR_u1)) = 0;
        FR_u1 = FR_u1 ./ max(FR_u1);

        for uid2=1:size(thisUnits,1)
            if strncmp(thisUnits.ID{uid2},UT.ID{uid},12), continue; end

            if uid~=1
                if(find(strcmp(UT.ID{uid},UnitPair.UID2) & strcmp(thisUnits.ID{uid2},UnitPair.UID1)))
                    continue;
                end
            end
            thisRips2=unique(ReT(find(cellfun(Params.cellfind(thisUnits.ID{uid2}),ReT.UnitID)),[1 2]));
            CoActRips = intersect(thisRips.RippleID,thisRips2.RippleID);

            uid1_2 = find(strcmp(thisUnits.ID{uid2},UT.ID));

            FR_u2 = FR(1,:,uid1_2);
            FR_u2(isnan(FR_u2)) = 0;
            FR_u2 = FR_u2 ./ max(FR_u2);

            q1 = size(thisRips2,1);
            pq = size(CoActRips,1);

            UnitPair.UID1{u}=UT.ID{uid};
            UnitPair.UID2{u}=thisUnits.ID{uid2};
            UnitPair.p(u)=p1;
            UnitPair.q(u)=q1;
            UnitPair.pq(u)=pq;
            UnitPair.un(u) = size(union(thisRips.RippleID,thisRips2.RippleID),1);
            UnitPair.p0(u)=p0;
            if p1==0 || q1==0
                UnitPair.crp(u)= nan;
            else
                UnitPair.crp(u) = pq/UnitPair.un(u);
            end

            %%
            thisField_1 = UTf(find(strncmp(UTf.ID,UT.ID{uid},12)),:);
            thisField_2 = UTf(find(strncmp(UTf.ID,thisUnits.ID{uid2},12)),:);
            if isempty(thisField_1) | isempty(thisField_2)
                UnitPair.region(u)=1;
                continue;
            end

            Norm_FRf_1 = [];
            for f = 1:size(thisField_1,1)
                FID = find(strcmp(UTf.ID,thisField_1.ID{f}));
                Norm_FRf_1(f,:) = FRf(1,:,FID) ./ max(FRf(1,:,FID));
            end
            Norm_FRf_1 = nanmean(Norm_FRf_1,1); Norm_FRf_1(isnan(Norm_FRf_1))=0;

            Norm_FRf_2 = [];
            for f = 1:size(thisField_2,1)
                FID = find(strcmp(UTf.ID,thisField_2.ID{f}));
                Norm_FRf_2(f,:) = FRf(1,:,FID) ./ max(FRf(1,:,FID));
            end
            Norm_FRf_2 = nanmean(Norm_FRf_2,1); Norm_FRf_2(isnan(Norm_FRf_2))=0;




            RDI_L1 = []; RDI_R1 = []; RDI_C1 = [];
            for f = 1:size(thisField_1,1)
                F = load([ROOT.Unit1 '\' thisField_1.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
                if size(F,1)<3, F(3,:)=nan; end
                RDI_L1 = [RDI_L1;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
                RDI_R1 = [RDI_R1; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
                RDI_C1 = [RDI_C1; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
            end

            %
            %             [~,o] = sort(thisField_1.RDI_LScene); RDI_L1 = RDI_L1(o,:);
            %             [~,o] = sort(thisField_1.RDI_RScene); RDI_R1 = RDI_R1(o,:);
            %             [~,o] = sort(thisField_1.RDI_LR); RDI_C1 = RDI_C1(o,:);


            RDI_L2 = []; RDI_R2 = []; RDI_C2 = [];
            for f = 1:size(thisField_2,1)
                F = load([ROOT.Unit1 '\' thisField_2.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
                if size(F,1)<3, F(3,:)=nan; end
                RDI_L2 = [RDI_L2;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
                RDI_R2 = [RDI_R2; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
                RDI_C2 = [RDI_C2; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
            end

            %             [~,o] = sort(thisField_2.RDI_LScene); RDI_L2 = RDI_L2(o,:);
            %             [~,o] = sort(thisField_2.RDI_RScene); RDI_R2 = RDI_R2(o,:);
            %             [~,o] = sort(thisField_2.RDI_LR); RDI_C2 = RDI_C2(o,:);

            %             s1 = size(RDI_L1,1); s2 = size(RDI_L2,1);
            %             if s1<s2,RDI_L1(s1+1:s2,:)=nan; elseif s1>s2,RDI_L2(s2+1:s1,:)=nan; end
            RDI_L1(isnan(RDI_L1))=0; RDI_L2(isnan(RDI_L2))=0;
            %             RDI_L1(abs(thisField_1.RDI_LScene)<0.1,:)=[]; if(isempty(RDI_L1)), RDI_L1=nan(1,40); end
            %             RDI_L2(abs(thisField_2.RDI_LScene)<0.1,:)=[]; if(isempty(RDI_L2)), RDI_L2=nan(1,40); end

            %              s1 = size(RDI_R1,1); s2 = size(RDI_R2,1);
            %             if s1<s2,RDI_R1(s1+1:s2,:)=nan; elseif s1>s2,RDI_R2(s2+1:s1,:)=nan; end
            RDI_R1(isnan(RDI_R1))=0; RDI_R2(isnan(RDI_R2))=0;
            %             RDI_R1(abs(thisField_1.RDI_RScene)<0.1,:)=[]; if(isempty(RDI_R1)), RDI_R1=nan(1,40); end
            %             RDI_R2(abs(thisField_2.RDI_RScene)<0.1,:)=[]; if(isempty(RDI_R2)), RDI_R2=nan(1,40); end

            %              s1 = size(RDI_C1,1); s2 = size(RDI_C2,1);
            %             if s1<s2,RDI_C1(s1+1:s2,:)=nan; elseif s1>s2,RDI_C2(s2+1:s1,:)=nan; end
            RDI_C1(isnan(RDI_C1))=0; RDI_C2(isnan(RDI_C2))=0;
            %             RDI_C1(abs(thisField_1.RDI_LR)<0.1,:)=[]; if(isempty(RDI_C1)), RDI_C1=nan(1,40); end
            %             RDI_C2(abs(thisField_2.RDI_LR)<0.1,:)=[]; if(isempty(RDI_C2)), RDI_C2=nan(1,40); end

            RDI_L1m = nanmean(RDI_L1,1); RDI_R1m = nanmean(RDI_R1,1);  RDI_C1m = nanmean(RDI_C1,1);
            RDI_L2m = nanmean(RDI_L2,1); RDI_R2m = nanmean(RDI_R2,1);  RDI_C2m = nanmean(RDI_C2,1);

            %             UnitPair.Sp(u) = corr(FR_u1',FR_u2');
            %             UnitPair.Sp(u) = corr(Norm_FRf_1',Norm_FRf_2');
            %
            %             UnitPair. Nsp_L(u) = corr(RDI_L1m',RDI_L2m');
            %             UnitPair. Nsp_R(u) = corr(RDI_R1m',RDI_R2m');
            %             UnitPair. Nsp_C(u) = corr(RDI_C1m',RDI_C2m');

            L1={thisField_1.RDI_LScene};
            R1={thisField_1.RDI_RScene};
            C1={thisField_1.RDI_LR};

            L2={thisField_2.RDI_LScene};
            R2={thisField_2.RDI_RScene};
            C2={thisField_2.RDI_LR};

            %%

            thisPair = table; f=1;
            for f1 = 1: size(thisField_1)
                for f2 = 1: size(thisField_2)
                    fi = find((strcmp(UPF.UID1,thisField_1.ID{f1}) & strcmp(UPF.UID2,thisField_2.ID{f2})) |...
                        (strcmp(UPF.UID2,thisField_1.ID{f1}) & strcmp(UPF.UID1,thisField_2.ID{f2})));
                    if ~isempty(fi)
                        thisPair(f,:) = UPF(fi,:);
                        f=f+1;
                    end
                end
            end
            thisPair(thisPair.SameUnit,:)=[];

            UnitPair.SpCorr_Max(u) = nanmax(thisPair.Sp);
            UnitPair.SpCorr_min(u) = nanmin(thisPair.Sp);
            UnitPair.SpCorr_Avg(u) = nanmean(thisPair.Sp);

            UnitPair.NSpCorr_L_Max(u) = nanmax(thisPair.Nsp_L);
            UnitPair.NSpCorr_L_min(u) = nanmin(thisPair.Nsp_L);
            UnitPair.NSpCorr_L_Avg(u) = nanmean(thisPair.Nsp_L);

            UnitPair.NSpCorr_R_Max(u) = nanmax(thisPair.Nsp_R);
            UnitPair.NSpCorr_R_min(u) = nanmin(thisPair.Nsp_R);
            UnitPair.NSpCorr_R_Avg(u) = nanmean(thisPair.Nsp_R);

            UnitPair.NSpCorr_C_Max(u) = nanmax(thisPair.Nsp_C);
            UnitPair.NSpCorr_C_min(u) = nanmin(thisPair.Nsp_C);
            UnitPair.NSpCorr_C_Avg(u) = nanmean(thisPair.Nsp_C);

            UnitPair.numFields_1(u) = size(thisField_1,1);
            UnitPair.numFields_2(u) = size(thisField_2,1);


            UnitPair.NSpCorr_M_Max(u) = nanmax([UnitPair.NSpCorr_L_Max(u),UnitPair.NSpCorr_R_Max(u),UnitPair.NSpCorr_C_Max(u)]);
            %%
            UT.Selectivity_LScene(uid) = getSel(L1{1});
            UT.Selectivity_RScene(uid) = getSel(R1{1});
            UT.Selectivity_LR(uid) = getSel(C1{1});

            thisUnits.Selectivity_LScene(uid2) = getSel(L2{1});
            thisUnits.Selectivity_RScene(uid2) = getSel(R2{1});
            thisUnits.Selectivity_LR(uid2) = getSel(C2{1});
            %%
            UnitPair.Lt1(u) = UT.Selectivity_LScene(uid) * max(L1{1}) / max([10^(-300),abs(max(L1{1}))]);
            UnitPair.Rt1(u) = UT.Selectivity_RScene(uid) * max(R1{1}) / max([10^(-300),abs(max(R1{1}))]);
            UnitPair.Ct1(u) = UT.Selectivity_LR(uid) * max(C1{1}) / max([10^(-300),abs(max(C1{1}))]);

            UnitPair.Lt2(u) = thisUnits.Selectivity_LScene(uid2) * max(L2{1}) / max([10^(-300),abs(max(L2{1}))]);
            UnitPair.Rt2(u) = thisUnits.Selectivity_RScene(uid2) * max(R2{1}) / max([10^(-300),abs(max(R2{1}))]);
            UnitPair.Ct2(u) = thisUnits.Selectivity_LR(uid2) * max(C2{1}) / max([10^(-300),abs(max(C2{1}))]);

            UnitPair.Lm1(u)=mean(thisField_1.RDI_LScene(abs(thisField_1.RDI_LScene)>0));
            UnitPair.Rm1(u)=mean(thisField_1.RDI_RScene(abs(thisField_1.RDI_RScene)>0));
            UnitPair.Cm1(u)=mean(thisField_1.RDI_LR(abs(thisField_1.RDI_LR)>0));

            UnitPair.Lm2(u)=mean(thisField_2.RDI_LScene(abs(thisField_2.RDI_LScene)>0));
            UnitPair.Rm2(u)=mean(thisField_2.RDI_RScene(abs(thisField_2.RDI_RScene)>0));
            UnitPair.Cm2(u)=mean(thisField_2.RDI_LR(abs(thisField_2.RDI_LR)>0));


            UnitPair.region(u) = cl;
            %              UnitPair.NumField1(u) = UT.NumField(uid);
            %              UnitPair.NumField2(u) = thisUnits.NumField(uid2);

            %         UnitPair.NumField1(u) = sum(UT.AvgFR(uid)==UT.AvgFR);
            %          UnitPair.NumField2(u) = sum(thisUnits.AvgFR(uid2)==thisUnits.AvgFR);

            if u==1032
                1;
            end

            u=u+1;
        end
    end
    save([ROOT.Save '\UnitPair_' thisRegion '.mat'],'UnitPair');
    writetable(UnitPair,[ROOT.Save '\UnitPair_' thisRegion '.xlsx'],'writemode','replacefile');
end


%%
UnitPairT=struct;
UnitPairT.CA1 = load([ROOT.Save '\UnitPair_CA1.mat']); UnitPairT.CA1= UnitPairT.CA1.UnitPair;
UnitPairT.SUB = load([ROOT.Save '\UnitPair_SUB.mat']); UnitPairT.SUB= UnitPairT.SUB.UnitPair;
%
% writetable(UnitPairT.CA1,[ROOT.Save '\UnitPair_CA1.xlsx'],'writemode','replacefile')
% writetable(UnitPairT.SUB,[ROOT.Save '\UnitPair_SUB.xlsx'],'writemode','replacefile')
%%
UnitPair = UnitPairT.SUB;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
CR_SUB.SF = UnitPair.CoactP((abs(UnitPair.Lt1)==1 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==1 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==1 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_SUB.SMF = UnitPair.CoactP((abs(UnitPair.Lt1)==2 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==2 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==2 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_SUB.MF = UnitPair.CoactP((abs(UnitPair.Lt1)==3 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==3 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==3 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));

UnitPair = UnitPairT.CA1;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
CR_CA1.SF = UnitPair.CoactP((abs(UnitPair.Lt1)==1 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==1 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==1 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_CA1.SMF = UnitPair.CoactP((abs(UnitPair.Lt1)==2 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==2 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==2 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_CA1.MF = UnitPair.CoactP((abs(UnitPair.Lt1)==3 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==3 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==3 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));


data = [nanmean([CR_SUB.SF;CR_CA1.SF]) nanmean(CR_SUB.SF) nanmean(CR_CA1.SF);...
    nanmean([CR_SUB.SMF;CR_SUB.MF;CR_CA1.SMF;CR_CA1.MF]) nanmean([CR_SUB.SMF;CR_SUB.MF]) nanmean([CR_CA1.SMF;CR_CA1.MF])];

figure;
bar(data')
xticklabels({'All','SUB','CA1'})
legend({'SF-SF','MF-MF'},'location','best')
ylabel('Co-Reactivation Probability')
set(gca,'fontsize',12,'fontweight','b')
[h,p] = ttest2([CR_CA1.SF],[CR_CA1.SMF;CR_CA1.MF])
%%  SF-SF, MF-MF coreactivation prob.
figure;
UnitPair = UnitPairT.SUB;
% UnitPair(UnitPair.Lt1==1 | UnitPair.Lt1==0 | UnitPair.Lt2==1 | UnitPair.Lt2==0,:) = [];
UnitPair(UnitPair.Lt1>1 | UnitPair.Lt1==0.5 | UnitPair.Lt2>1 | UnitPair.Lt2==0.5,:) = [];
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
d_List = {'L','R','C'};

for d0 = 1:numel(d_List)
    d = d_List{d0};
    Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
    Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
    Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
    Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
    Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
        UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
    Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
    Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
    Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;
end

data_sub = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,4,2)
bar(data_sub')
ylabel('Co-Reactivation Probability'); ylim([0 0.18])
set(gca,'fontsize',12,'fontweight','b')
title('SUB')
xticklabels({'Left','Right','Choice'})


UnitPair = UnitPairT.CA1;
% UnitPair(UnitPair.Lt1==1 | UnitPair.Lt1==0 | UnitPair.Lt2==1 | UnitPair.Lt2==0,:) = [];
UnitPair(UnitPair.Lt1>1 | UnitPair.Lt1==0.5 | UnitPair.Lt2>1 | UnitPair.Lt2==0.5,:) = [];
UnitPair.CoactP = UnitPair.pq./UnitPair.un;

for d0 = 1:numel(d_List)
    d = d_List{d0};
    Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
    Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
    Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
    Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
    Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
        UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
    Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
    Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
    Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;
end

data_ca1 = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,2,2)
bar(data_ca1')
ylabel('Co-Reactivation Probability'); ylim([0 0.18])
set(gca,'fontsize',12,'fontweight','b')
title('CA1')
xticklabels({'Left','Right','Choice'})
legend({'Non-Non','Non-Remap','Remap-Remap (diff.)','Remap-Remap (same)'},'location','eastoutside','Orientation','vertical')
% data_ca1 = [nanmean([Coreacts.L1; Coreacts.R1;Coreacts.C1]);...
%     nanmean([Coreacts.L2; Coreacts.R2;Coreacts.C2]);...
%     nanmean([Coreacts.L3; Coreacts.R3;Coreacts.C3]);...
%     nanmean([Coreacts.L4; Coreacts.R4;Coreacts.C4])];

[h,p] = ttest2(Coreacts.C4,Coreacts.C3)

%%
data = [data_sub(:,1)' ; data_ca1(:,1)'];
figure;
bar(data)
ylabel('Co-Reactivation Probability');
set(gca,'fontsize',12,'fontweight','b')
title('SF')
xticklabels({'SUB','CA1'})
legend({'Non-Non','Non-Remap','Remap-Remap (diff.)','Remap-Remap (same)'},'location','eastoutside','Orientation','vertical')

%% CoReactivation scatterplot

d='R';
w = 0.5; n=2;
for r = 1:numel(RegionList)
    thisRegion = RegionList{r};
    figure; sgtitle(thisRegion)
    for i=1:n
        UP = UnitPairT.(thisRegion);

        UP.CoactP = UP.pq./UP.un;
        UP(UP.p==0,:)=[]; UP(UP.q==0,:)=[];
        % UP(UP.(['Nsp_' d])<0,:)=[];
        UP = UP((UP.Sp)>1-w*i & (UP.Sp)<=(1+w)-w*i,:);
        % x = UP.(['Nsp_' d]); s=UP.Sp+1; y=UP.CoactP;
        x = min([UP.Nsp_L,UP.Nsp_R,UP.Nsp_C],[],2); s=UP.Sp+1; y=UP.CoactP;
        %  x(y>0.25)=[]; y(y>0.25)=[];
        % id1 = ((UP.([d 't1'])>1|UP.([d 't1'])==0.5)& UP.([d 't1'])<4) & ((UP.([d 't2'])>1|UP.([d 't2'])==0.5)& UP.([d 't1'])<4);
        % id2 = (UP.([d 't1'])==1|UP.([d 't1'])==0) & (UP.([d 't2'])==1|UP.([d 't2'])==0);

        ax = subplot(1,n,i);
        scatter(x,y,40,CList(r,:),'filled')
        hold on
        mdl = fitlm(x,y);
        h=plot(mdl);
        delete(h(1))
        set(h(2),'LineWidth',3,'color','k')
        set(h(3),'LineWidth',2,'color','k')
        set(h(4),'LineWidth',2,'color','k')
        legend off

        hold on
        % y(isnan(x))=[]; x(isnan(x))=[];

        % scatter3(x(id2),s(id2),y(id2),20,'b')
        % xlabel(['Non-Spatial correlation (' d ')']);
        xlabel(['Non-Spatial correlation']);
        ylabel('Co-Reactivation Probability');
        % zlabel('Co-Reactivation Probability'); zlim([0 1])
        title([ num2str(1-i*w) ' < Spatial Corr. < ' num2str((1+w)-i*w)])
        set(gca,'fontsize',12,'fontweight','b')
        %  lsline(ax)

        x(isnan(y))=[]; y(isnan(y))=[];
        y(isnan(x))=[]; x(isnan(x))=[];


        hold on


        [p1,p2] = corr(x,y);
        text(.3,.6,['corr = ' jjnum2str(p1,3)],'color','k')
        text(.3,.57,['p = ' num2str(p2)],'color','k')
        xlim([-1 1]); ylim([0 .7])



    end
end

%%
UP = UnitPairT.SUB;

n=0.62876;
[val,idx]=min(abs(UP.(['Nsp_' d])-n));
minVal=UP.Nsp_L(idx)

UP.UID1{idx}
UP.UID2{idx}
%% Sp, NSp boxplot
figure
subplot(2,2,1)
x1 = UnitPair.CA1.Sp;
x2 = UnitPair.SUB.Sp;
x = [x1;x2];
g = [repmat({['CA1 (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['SUB (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('Spatial Correlation')
xlabel('')
ylim([-1 1])
[~,p] = ttest2(x1,x2)


subplot(2,2,1)
x1 = nanmax([UnitPair.CA1.Nsp_L,UnitPair.CA1.Nsp_R,UnitPair.CA1.Nsp_C],[],2);
x2 = nanmax([UnitPair.SUB.Nsp_L,UnitPair.SUB.Nsp_R,UnitPair.SUB.Nsp_C],[],2);
x = [x1;x2];
g = [repmat({['CA1 (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['SUB (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('Non-Spatial Correlation')
xlabel('')
ylim([-1 1])
[~,p] = ttest2(x1,x2)

subplot(2,2,2)
x1 = UnitPair.CA1.Nsp_L;
x2 = UnitPair.SUB.Nsp_L;
x = [x1;x2];
g = [repmat({['CA1 (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['SUB (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('Non-Spatial Correlation (L)')
xlabel('')
ylim([-1 1])
[~,p] = ttest2(x1,x2)

subplot(2,2,3)
x1 = UnitPair.CA1.Nsp_R;
x2 = UnitPair.SUB.Nsp_R;
x = [x1;x2];
g = [repmat({['CA1 (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['SUB (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('Non-Spatial Correlation (R)')
xlabel('')
ylim([-1 1])
[~,p] = ttest2(x1,x2)

subplot(2,2,4)
x1 = UnitPair.CA1.Nsp_C;
x2 = UnitPair.SUB.Nsp_C;
x = [x1;x2];
g = [repmat({['CA1 (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['SUB (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('Non-Spatial Correlation (C)')
xlabel('')
ylim([-1 1])
[~,p] = ttest2(x1,x2)

%%
UP = UnitPairT.SUB;
UP = UP((UP.Lt1>1 | UP.Lt1==0.5) & (UP.Lt2>1 | UP.Lt2==0.5),:);

% UP = UP((UP.Lt1==1 | UP.Lt1==0) & (UP.Lt2==1 | UP.Lt2==0),:);
UP.CoactP = UP.pq ./ UP.un;
x2 = 1-abs(UP.([d 'm1'])-UP.([d 'm2'])); x1=abs(UP.Sp); y=UP.CoactP;
X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X)    % Removes NaN data

scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):.01:max(x1);
x2fit = min(x2):.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Sp Sim.')
ylabel('NonSp Sim.')
zlabel('CoReact Prob')
view(50,10)
hold off

%%
UnitPairT.SUB.CoactP = UnitPairT.SUB.pq ./ UnitPairT.SUB.un;
UnitPairT.CA1.CoactP = UnitPairT.CA1.pq ./ UnitPairT.CA1.un;
figure
x1 = UnitsTable.CA1.FieldSize_TP(UnitsTable.CA1.NumField==1);
x2 = UnitsTable.CA1.FieldSize_TP(UnitsTable.CA1.NumField>1);
x3 = UnitPairT.SUB.CoactP(UnitsTable.SUB.NumField==1);
x4 = UnitPairT.SUB.CoactP(UnitsTable.SUB.NumField>1);
x = [x1;x2;x3;x4].*2;
g = [repmat({['CA1-SF (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1-MF (n=' num2str(size(x2,1)) ')']},size(x2,1),1);...
    repmat({['SUB-SF (n=' num2str(size(x3,1)) ')']},size(x3,1),1);repmat({['SUB-MF (n=' num2str(size(x4,1)) ')']},size(x4,1),1)];
boxplot(x,g)
%%
function s = getSel(rdis)
rdisig = rdis(abs(rdis)>0.1);
if length(rdis)>1
    if length(rdisig)==1
        s=2;
    elseif length(rdisig)>1
        if min(rdisig>0) | min(rdisig<0)
            s=3;
        else
            s=4;
        end
    else
        s=0.5;
    end

else
    if length(rdisig)>0
        s=1;
    else
        s=0;
    end

end
end


%%
function line_w_shade(x1,y)
a=0.05;
X = [ones(size(x1)) x1];
[b,bint,~,~,stats] = regress(y,X) ;
xval = min(x1)-a:0.01:max(x1)+a;
yhat = b(1)+b(2)*xval;
ylow = bint(1,1)+bint(2,1)*xval;
yupp = bint(1,2)+bint(2,2)*xval;


tbl = table(x1,y,'VariableNames',{'NonSpatialCorr','CoReact'});
mdl = fitlm(x1,y);
% xpt = find(yupp>=yhat,1,'first');
%
% temp = ylow(1:xpt-1);
% ylow = [yupp(1:xpt-1) ylow(xpt:end)];
% yupp = [temp yupp(xpt:end)];
% plot(x1,y,'ks', 'LineWidth', 3, 'MarkerSize', 2);
% p5.Color(4) = 0.5;
hold on;
p6=plot(xval,yhat,'r','linewidth',3);
p6.Color(4) = 0.5;

% scatter(lin(:,1),lin(:,2))
% hold off



hold on
patch([xval fliplr(xval)], [ylow fliplr(yupp)], 'r','EdgeColor','none','facealpha',0.2)


end
