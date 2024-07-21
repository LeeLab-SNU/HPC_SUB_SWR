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

UnitPair.(thisRegion).SpCorr_Max(pid) = nanmax(thisPair.Sp);
UnitPair.(thisRegion).SpCorr_min(pid) = nanmin(thisPair.Sp);
UnitPair.(thisRegion).SpCorr_Avg(pid) = nanmean(thisPair.Sp);

UnitPair.(thisRegion).NSpCorr_L_Max(pid) = nanmax(thisPair.Nsp_L);
UnitPair.(thisRegion).NSpCorr_L_min(pid) = nanmin(thisPair.Nsp_L);
UnitPair.(thisRegion).NSpCorr_L_Avg(pid) = nanmean(thisPair.Nsp_L);

UnitPair.(thisRegion).NSpCorr_R_Max(pid) = nanmax(thisPair.Nsp_R);
UnitPair.(thisRegion).NSpCorr_R_min(pid) = nanmin(thisPair.Nsp_R);
UnitPair.(thisRegion).NSpCorr_R_Avg(pid) = nanmean(thisPair.Nsp_R);

UnitPair.(thisRegion).NSpCorr_C_Max(pid) = nanmax(thisPair.Nsp_C);
UnitPair.(thisRegion).NSpCorr_C_min(pid) = nanmin(thisPair.Nsp_C);
UnitPair.(thisRegion).NSpCorr_C_Avg(pid) = nanmean(thisPair.Nsp_C);
end
writetable(UnitPair.(thisRegion),[ROOT.Save '\UnitPair_' (thisRegion) '.xlsx'],'writemode','replacefile')
%%



pid = 12;
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



%     mdl = corr(RDI_C1m' ,RDI_C2m' )
    %%
    figure('position',[135,80,1041,885])


    subplot('position',[.1 .65 .3 .08])
    FL = FR(1,:,id1)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc((Norm_FRf_1)); axis off
    idf1 = find(strncmp(UT.ID{id1},UTF.ID,12));




    subplot('position',[.5 .65 .3 .08])
    FL = FR(1,:,id2)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc(Norm_FRf_2); axis off
    
    idf2 = find(strncmp(UT.ID{id2},UTF.ID,12));


    colormap jet


    subplot('position',[.85 .7 .1 .05])
    text(0,0,sprintf(['Sp.Corr \n Max ' jjnum2str((UP.SpCorr_Max(pid)),3)... 
        '\n min ' jjnum2str((UP.SpCorr_min(pid)),3) '\n Avg ' jjnum2str((UP.SpCorr_Avg(pid)),3)]),'fontsize',15,'FontWeight','b')
    axis off
    %%
    subplot('position',[.1 .8 .3 .14])
    hold on
    mx=max(FR(1,:,[id1]),[],'all');
    plot(smooth(FR(1,:,id1)')','linewidth',4,'color','k')
    for f = 1:size(idf1,1)
        plot(smooth(Maps_f(1,:,idf1(f))),'linewidth',4,'color',hex2rgb(Clist_f{f}))
        xlim([0 30]);
        ylim([0 mx*1.1])

        settick(31,mx*1.1)
        xlabel('')
        title(['Cell 1 (' UT.ID{id1} ')'],'FontSize',15)
    end
    % title('Ratemap for Left scene trials')
    set(gca,'FontSize',12,'FontWeight','b')

    subplot('position',[.5 .8 .3 .14])
    hold on
    mx=max(FR(1,:,[id2]),[],'all');
    thisField_2 = UTF(find(strncmp(UTF.ID,UT.ID{id2},12)),:);
    plot(smooth(FR(1,:,id2)')','linewidth',4,'color','k')
    for f = 1:size(idf2,1)
        plot(smooth(Maps_f(1,:,idf2(f))),'linewidth',4,'color',hex2rgb(Clist_f{f}))
        %     p = plot(Maps_f(3,:,idf2(f)),'linewidth',4,'color',hex2rgb(Clist_f{f})); p.Color(4)=0.2;
        %     plot(Maps_f(3,:,idf2(f)),'linewidth',4,'color',hex2rgb(Clist_f{f}),'linestyle',':')
        %     if abs(UTF.RDI_LScene(idf2(f)))>0.1, c=hex2rgb(Clist_f{f}); else, c='k'; end
        xlim([0 30]);
        ylim([0 mx*1.1])

        settick(31,mx*1.1)
        xlabel('')
        title(['Cell 2 (' UT.ID{id2} ')'],'FontSize',15)
    end
    % title('Ratemap for Left scene trials')
    set(gca,'FontSize',12,'FontWeight','b')
    %%

    ax2 = subplot('position',[.1 .45 .3 .12]);
    imagesc((RDI_L1)); axis off
    colormap(ax2,rb); caxis([-1 1])
    title('In-field RDI mean (Left scene)','fontsize',12)

    ax3 = subplot('position',[.5 .45 .3 .12]);
    imagesc((RDI_L2)); axis off
    colormap(ax3,rb); caxis([-1 1])
    title('In-field RDI mean (Left scene)','fontsize',12)


    ax4 = subplot('position',[.1 .25 .3 .12]);
    imagesc((RDI_R1)); axis off
    colormap(ax4,rb); caxis([-1 1])
    title('In-field RDI mean (Right scene)','fontsize',12)

    ax5 = subplot('position',[.5 .25 .3 .12]);
    imagesc(RDI_R2); axis off
    colormap(ax5,rb); caxis([-1 1])
    title('In-field RDI mean (Right scene)','fontsize',12)


    ax6 = subplot('position',[.1 .05 .3 .12]);
    imagesc(RDI_C1); axis off
    colormap(ax6,rb); caxis([-1 1])
    title('In-field RDI mean (Choice)','fontsize',12)

    ax7 = subplot('position',[.5 .05 .3 .12]);
    imagesc(RDI_C2); axis off
    colormap(ax7,rb); caxis([-1 1])
    title('In-field RDI mean (Choice)','fontsize',12)


    subplot('position',[.85 .5 .1 .05])
    text(0,0,sprintf(['NSp.Corr (L) \n Max ' jjnum2str((UP.NSpCorr_L_Max(pid)),3)...
        '\n min ' jjnum2str((UP.NSpCorr_L_min(pid)),3) '\n Avg ' jjnum2str((UP.NSpCorr_L_Avg(pid)),3)]),'fontsize',15,'FontWeight','b')
    axis off

       subplot('position',[.85 .3 .1 .05])
    text(0,0,sprintf(['NSp.Corr (R) \n Max ' jjnum2str((UP.NSpCorr_R_Max(pid)),3)...
        '\n min ' jjnum2str((UP.NSpCorr_R_min(pid)),3) '\n Avg ' jjnum2str((UP.NSpCorr_R_Avg(pid)),3)]),'fontsize',15,'FontWeight','b')
    axis off

       subplot('position',[.85 .1 .1 .05])
    text(0,0,sprintf(['NSp.Corr (C) \n Max ' jjnum2str((UP.NSpCorr_C_Max(pid)),3)...
        '\n min ' jjnum2str((UP.NSpCorr_C_min(pid)),3) '\n Avg ' jjnum2str((UP.NSpCorr_C_Avg(pid)),3)]),'fontsize',15,'FontWeight','b')
    axis off



    %%

    subplot('Position',[.85 .8 .14 .15])

    bar([1 2 3 4],[UP.p(pid) UP.q(pid) UP.un(pid) UP.pq(pid)]./ UP.p0(pid))
    xticklabels({'Cell 1', 'Cell 2', '1 or 2', '1 and 2'})
    ylabel('Ripple participation rate')
    title(['Co-React. Prob. = ' jjnum2str(UP.pq(pid)/UP.un(pid),3)])
%     set(gca,'fontsize',12,'FontWeight','b')

    ylim([0 .8])

    
    %%
    ROOT.Fig = [ROOT.Mother '\Processed Data\units_mat\ProfilingSheet\U9_ca1'];

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