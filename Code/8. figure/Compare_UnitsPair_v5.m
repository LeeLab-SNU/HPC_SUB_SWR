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
UnitPair.SUB = readtable([ROOT.Save '\UnitPair_SUB.xlsx']);
UnitPair_field.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);
UnitPair.CA1= readtable([ROOT.Save '\UnitPair_CA1.xlsx']);
UnitPair_field.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%%


rb = redblue(256);
thisRegion = 'CA1';

UP = UnitPair.(thisRegion); 
UPF = UnitPair_field.(thisRegion);
UT = UnitsTable.(thisRegion); UTF = UnitsTable_field.(thisRegion);
FR = FRMaps.(thisRegion); FRF = FRMaps.([thisRegion '_field']);

for pid = 1:size(UP,1)
    ID1 = UP.UID1{pid}; id1 = find(strcmp(ID1,UT.ID));
    ID2 = UP.UID2{pid}; id2 = find(strcmp(ID2,UT.ID));

    Clist_f = {'#FF6A68','#00B455','#00A2FF','#FA7344','#00B7A1'};
    Maps_f = FRMaps.SUB_field;


    thisField_1 = UTF(find(strncmp(UTF.ID,UP.UID1{pid},12)),:);
    thisField_2 = UTF(find(strncmp(UTF.ID,UP.UID2{pid},12)),:);

            Norm_FRf_1 = [];
            for f = 1:size(thisField_1,1)
                FID = find(strcmp(UTF.ID,thisField_1.ID{f}));
                Norm_FRf_1(f,:) = FRF(1,:,FID) ./ max(FRF(1,:,FID));
            end
%             Norm_FRf_1 = nanmean(Norm_FRf_1,1); Norm_FRf_1(isnan(Norm_FRf_1))=0;

            Norm_FRf_2 = [];
            for f = 1:size(thisField_2,1)
                FID = find(strcmp(UTF.ID,thisField_2.ID{f}));
                Norm_FRf_2(f,:) = FRF(1,:,FID) ./ max(FRF(1,:,FID));
            end
%             Norm_FRf_2 = nanmean(Norm_FRf_2,1); Norm_FRf_2(isnan(Norm_FRf_2))=0;

 RDI_L1 = []; RDI_R1 = []; RDI_C1 = [];
    for f = 1:size(thisField_1,1)
        F = load([ROOT.Unit1 '\' thisField_1.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
        if size(F,1)<3, F(3,:)=nan; end
        RDI_L1 = [RDI_L1;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
        RDI_R1 = [RDI_R1; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
        RDI_C1 = [RDI_C1; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
    end

    RDI_L2 = []; RDI_R2 = []; RDI_C2 = [];
    for f = 1:size(thisField_2,1)
        F = load([ROOT.Unit1 '\' thisField_2.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
        if size(F,1)<3, F(3,:)=nan; end
        RDI_L2 = [RDI_L2;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
        RDI_R2 = [RDI_R2; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
        RDI_C2 = [RDI_C2; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
    end

    RDI_L1(isnan(RDI_L1))=0; RDI_L2(isnan(RDI_L2))=0;
    %             RDI_L1(abs(thisField_1.RDI_LScene)<0.1,:)=[]; if(isempty(RDI_L1)), RDI_L1=nan(1,40); end
    %             RDI_L2(abs(thisField_2.RDI_LScene)<0.1,:)=[]; if(isempty(RDI_L2)), RDI_L2=nan(1,40); end

    RDI_R1(isnan(RDI_R1))=0; RDI_R2(isnan(RDI_R2))=0;
    %             RDI_R1(abs(thisField_1.RDI_RScene)<0.1,:)=[]; if(isempty(RDI_R1)), RDI_R1=nan(1,40); end
    %             RDI_R2(abs(thisField_2.RDI_RScene)<0.1,:)=[]; if(isempty(RDI_R2)), RDI_R2=nan(1,40); end

    RDI_C1(isnan(RDI_C1))=0; RDI_C2(isnan(RDI_C2))=0;
    %             RDI_C1(abs(thisField_1.RDI_LR)<0.1,:)=[]; if(isempty(RDI_C1)), RDI_C1=nan(1,40); end
    %             RDI_C2(abs(thisField_2.RDI_LR)<0.1,:)=[]; if(isempty(RDI_C2)), RDI_C2=nan(1,40); end


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
%     mdl = corr(RDI_C1m' ,RDI_C2m' )
    %%
    figure('position',[135,80,1041,885])


    subplot('position',[.1 .8 .3 .08])
    FL = FR(1,:,id1)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc((Norm_FRf_1)); axis off
    idf1 = find(strncmp(UT.ID{id1},UTF.ID,12));
    title(['Cell 1 (' UT.ID{id1} ')'],'FontSize',15)



    subplot('position',[.5 .8 .3 .08])
    FL = FR(1,:,id2)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc(Norm_FRf_2); axis off
    idf2 = find(strncmp(UT.ID{id2},UTF.ID,12));
        title(['Cell 2 (' UT.ID{id2} ')'],'FontSize',15)


    colormap jet


    subplot('position',[.85 .85 .1 .05])
    text(0,0,sprintf(['Sp.Corr \n Max ' jjnum2str((UP.SpCorr_Max(pid)),3)... 
        '\n min ' jjnum2str((UP.SpCorr_min(pid)),3) '\n Avg ' jjnum2str((UP.SpCorr_Avg(pid)),3)]),'fontsize',15,'FontWeight','b')
    axis off
    %%
    subplot('position',[.1 .45 .25 .25])
    hold on
scatter(thisPair.Sp,thisPair.Nsp_L,40,'k','filled')
line([-1 1], [-1 1])
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial corr'); ylabel('In-field scene(L) selectivity corr')
    set(gca,'FontSize',12,'FontWeight','b')

    subplot('position',[.5 .45 .25 .25])
    hold on
scatter(thisPair.Sp,thisPair.Nsp_R,40,'k','filled')
line([-1 1], [-1 1])
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial corr'); ylabel('In-field scene(R) selectivity corr')
    set(gca,'FontSize',12,'FontWeight','b')

  subplot('position',[.1 .1 .25 .25])
    hold on
scatter(thisPair.Sp,thisPair.Nsp_C,40,'k','filled')
line([-1 1], [-1 1])
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial corr'); ylabel('In-field choice selectivity corr')
    set(gca,'FontSize',12,'FontWeight','b')

    subplot('position',[.5 .1 .25 .25])
    hold on
    thisPair.Nsp_M = max([thisPair.Nsp_L thisPair.Nsp_R thisPair.Nsp_C],[],2);
    [~,id1] = max(thisPair.Sp);   [~,id2] = max(thisPair.Nsp_M); 
scatter(thisPair.Sp,thisPair.Nsp_M,40,'k','filled')
scatter(thisPair.Sp(id1),thisPair.Nsp_M(id1),40,'r','filled')
scatter(thisPair.Sp(id2),thisPair.Nsp_M(id2),40,'r','filled')
line([-1 1], [-1 1])
line([-1 thisPair.Sp(id1)], [-1 thisPair.Nsp_M(id1)],'color','r','linestyle',':')
line([-1 thisPair.Sp(id2)], [-1 thisPair.Nsp_M(id2)],'color','r','linestyle',':')
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial corr'); ylabel('In-field max. selectivity corr')
    set(gca,'FontSize',12,'FontWeight','b')
    %% spatial-nonspatial heterogeneity
subplot('position',[.76 .15 .05 .05])
v0 = [2 2]; v1=[thisPair.Sp(id1) thisPair.Nsp_M(id1)]+1; v2=[thisPair.Sp(id2) thisPair.Nsp_M(id2)]+1;

% 벡터의 크기 계산
norm_v0 = norm(v0);
norm_v1 = norm(v1);
norm_v2 = norm(v2);

% 코사인 각도 계산
r1 = acos(dot(v0, v1) / (norm_v0 * norm_v1));
d1 = rad2deg(r1);

r2 = acos(dot(v2, v0) / (norm_v0 * norm_v2));
d2 = rad2deg(r2);

r3 = acos(dot(v2, v1) / (norm_v1 * norm_v2));
d3 = rad2deg(r3);

het = sin(r1-pi/4) * sin(r2-pi/4);

text(0,2,['θ(Sp) = ' jjnum2str(d1,2) '°'])
text(0,1,['θ(Nsp) = ' jjnum2str(d2,2) '°'])
text(0,0,sprintf(['spatial-nonspatial heterogeneity strength \n = ' jjnum2str(het,3)]))
text(0,-1,sprintf(['θ(Sp-Nsp) ' jjnum2str(d3,2)]))

axis off

UnitPair.(thisRegion).deg_sp(pid) = d1;
UnitPair.(thisRegion).deg_sp(pid) = d2;
UnitPair.(thisRegion).deg_sp_nsp(pid) = d3;
UnitPair.(thisRegion).hetero(pid) = het;
    %%
    ROOT.Fig = [ROOT.Mother '\Processed Data\units_mat\ProfilingSheet\U11_ca1'];

    if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end
    saveas(gca,[ROOT.Fig '\'  UP.UID1{pid} '_' UP.UID2{pid} '.svg'])
    saveas(gca,[ROOT.Fig '\'  UP.UID1{pid} '_' UP.UID2{pid} '.png'])
    close all
end
%%
function settick(Dv,mx)
xlabel('position')
xticks([0 Dv 48])
xticklabels({'StBox','Dv','Fd'})
xlim([0 48])
ylabel('firing rate (Hz)')
ylim([0 mx*1.1])
line([Dv Dv],[0 mx*1.1],'linestyle','--','color','k')
end