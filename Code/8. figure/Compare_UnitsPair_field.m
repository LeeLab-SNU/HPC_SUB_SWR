UnitPair=UnitPair_field;

thisRegion = 'SUB';
pid = 12;

rb = redblue(256);

UP = UnitPair.(thisRegion); UT = UnitsTable.(thisRegion); UTF = UnitsTable_field.(thisRegion);
FR = FRMaps.(thisRegion); FRF = FRMaps.([thisRegion '_field']);
for pid = 1:size(UP,1)
    ID1 = UP.UID1{pid}; id1 = find(strcmp(ID1,UTF.ID));
    ID2 = UP.UID2{pid}; id2 = find(strcmp(ID2,UTF.ID));

    Clist_f = {'#FF6A68','#00B455','#00A2FF','#FA7344','#00B7A1'};
    Maps_f = FRMaps.SUB_field;


    thisField_1 = UTF(find(strncmp(UTF.ID,UP.UID1{pid},12)),:);
    thisField_2 = UTF(find(strncmp(UTF.ID,UP.UID2{pid},12)),:);
    %%
            Norm_FRf_1 = [];

                Norm_FRf_1 = FRF(1,:,id1) ./ max(FRF(1,:,id1));
            Norm_FRf_1(isnan(Norm_FRf_1))=0;

            Norm_FRf_2 = [];

                Norm_FRf_2 = FRF(1,:,id2) ./ max(FRF(1,:,id2));
             Norm_FRf_2(isnan(Norm_FRf_2))=0;


        RDI_L1 = []; RDI_R1 = []; RDI_C1 = [];

                F = load([ROOT.Unit1 '\' ID1 '.mat']); Frdi= F.RDIs_field;
                if size(Frdi,1)<3, Frdi(3,:)=nan; end
                FZ1 = F.thisFieldMap1D.skaggsMap1D{2};
                FB1 = F.thisFieldMap1D.skaggsMap1D{4};
                M = max([FZ1;FB1]); FZ1=FZ1./M; FB1=FB1./M;
                FP1 = F.thisFieldMap1D.skaggsMap1D{3};
                FM1 = F.thisFieldMap1D.skaggsMap1D{5};
                M = max([FP1;FM1]); FP1=FP1./M; FM1=FM1./M;
                FL1 = F.thisFieldMap1D.skaggsMap_left1D;
                FR1 = F.thisFieldMap1D.skaggsMap_right1D;
                M = max([FL1;FR1]); FL1=FL1./M; FR1=FR1./M;
                RDI_L1 = [interp1(1:size(Frdi,2),Frdi(1,:),linspace(1,size(Frdi,2),40))]; RDI_L1(isnan(RDI_L1))=0;
                RDI_R1 = [interp1(1:size(Frdi,2),Frdi(2,:),linspace(1,size(Frdi,2),40))]; RDI_R1(isnan(RDI_R1))=0;
                RDI_C1 = [interp1(1:size(Frdi,2),Frdi(3,:),linspace(1,size(Frdi,2),40))]; RDI_C1(isnan(RDI_C1))=0;


                RDI_L2 = []; RDI_R2 = []; RDI_C2 = [];
                F = load([ROOT.Unit1 '\' ID2 '.mat']); Frdi= F.RDIs_field;
                if size(Frdi,1)<3, Frdi(3,:)=nan; end
                FZ2 = F.thisFieldMap1D.skaggsMap1D{2};
                FB2 = F.thisFieldMap1D.skaggsMap1D{4};
                M = max([FZ2;FB2]); FZ2=FZ2./M; FB2=FB2./M;
                FP2 = F.thisFieldMap1D.skaggsMap1D{3};
                FM2 = F.thisFieldMap1D.skaggsMap1D{5};
                M = max([FP2;FM2]); FP2=FP2./M; FM2=FM2./M;
                FL2 = F.thisFieldMap1D.skaggsMap_left1D;
                FR2 = F.thisFieldMap1D.skaggsMap_right1D;
                M = max([FL2;FR2]); FL2=FL2./M; FR2=FR2./M;
                RDI_L2 = [interp1(1:size(Frdi,2),Frdi(1,:),linspace(1,size(Frdi,2),40))]; RDI_L2(isnan(RDI_L2))=0;
                RDI_R2 = [interp1(1:size(Frdi,2),Frdi(2,:),linspace(1,size(Frdi,2),40))]; RDI_R2(isnan(RDI_R2))=0;
                RDI_C2 = [interp1(1:size(Frdi,2),Frdi(3,:),linspace(1,size(Frdi,2),40))]; RDI_C2(isnan(RDI_C2))=0;



    %%
    figure('position',[135,80,1041,885])


    subplot('position',[.1 .7 .3 .03])
    FL = FRF(1,:,id1)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc(smooth(Norm_FRf_1)'); axis off
%     idf1 = find(strncmp(UT.ID{id1},UTF.ID,12));




    subplot('position',[.5 .7 .3 .03])
    FL = FRF(1,:,id2)'; FL(isnan(FL))=[]; FL(FL==0)=[];
    imagesc(smooth(Norm_FRf_2)'); axis off
    
%     idf2 = find(strncmp(UT.ID{id2},UTF.ID,12));


    colormap jet


    subplot('position',[.85 .72 .1 .05])
    text(0,0,['Sp.Corr = ' jjnum2str((UP.Sp(pid)),3)],'fontsize',15,'FontWeight','b')
    axis off
    %%
    subplot('position',[.1 .8 .3 .141])
    hold on
    mx=max(FRF(1,:,[id1]),[],'all');
    plot(smooth(FRF(1,:,id1)')','linewidth',4,'color','k')
    ylabel('FR(Hz)'); xlabel('Pos(bin)')
       title(ID1)
    set(gca,'FontSize',12,'FontWeight','b')

    subplot('position',[.5 .8 .3 .141])
    hold on
    mx=max(FRF(1,:,[id2]),[],'all');
%     thisField_2 = UTF(find(strncmp(UTF.ID,UT.ID{id2},12)),:);
    plot(smooth(FRF(1,:,id2)')','linewidth',4,'color','k')
        ylabel('FR(Hz)'); xlabel('Pos(bin)')
        title(ID2)
    set(gca,'FontSize',12,'FontWeight','b')
    %%

    ax2 = subplot('position',[.1 .6 .3 .025]);
    imagesc(smooth(RDI_L1)'); axis off
    colormap(ax2,rb); caxis([-1 1])
    title('In-field RDI mean (Left scene)','fontsize',12)

    ax3 = subplot('position',[.5 .6 .3 .025]);
    imagesc(smooth(RDI_L2)'); axis off
    colormap(ax3,rb); caxis([-1 1])
    title('In-field RDI mean (Left scene)','fontsize',12)

    subplot('position',[.1 .55 .3 .025]);
    imagesc(smooth(FZ1)'); axis off; clim([0 1])
    text(-2,1,'Z','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.1 .5 .3 .025]);
    imagesc(smooth(FB1)'); axis off;  clim([0 1])
    text(-2,1,'B','fontsize',12,'fontweight','b','COLOR','b')

        subplot('position',[.5 .55 .3 .025]);
    imagesc(smooth(FZ2)'); axis off;  clim([0 1])
  text(-2,1,'Z','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.5 .5 .3 .025]);
    imagesc(smooth(FB2)'); axis off;  clim([0 1])
    text(-2,1,'B','fontsize',12,'fontweight','b','COLOR','b')


    ax4 = subplot('position',[.1 .4 .3 .025]);
    imagesc(smooth(RDI_R1)'); axis off
    colormap(ax4,rb); caxis([-1 1])
    title('In-field RDI mean (Right scene)','fontsize',12)

    ax5 = subplot('position',[.5 .4 .3 .025]);
    imagesc(smooth(RDI_R2)'); axis off
    colormap(ax5,rb); caxis([-1 1])
    title('In-field RDI mean (Right scene)','fontsize',12)

       subplot('position',[.1 .35 .3 .025]);
    imagesc(smooth(FP1)'); axis off;  clim([0 1])
    text(-2,1,'P','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.1 .3 .3 .025]);
    imagesc(smooth(FM1)'); axis off;  clim([0 1])
    text(-2,1,'M','fontsize',12,'fontweight','b','COLOR','b')

        subplot('position',[.5 .35 .3 .025]);
    imagesc(smooth(FP2)'); axis off;  clim([0 1])
  text(-2,1,'P','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.5 .3 .3 .025]);
    imagesc(smooth(FM2)'); axis off;  clim([0 1])
    text(-2,1,'M','fontsize',12,'fontweight','b','COLOR','b')


    ax6 = subplot('position',[.1 .2 .3 .025]);
    imagesc(smooth(RDI_C1)'); axis off
    colormap(ax6,rb); caxis([-1 1])
    title('In-field RDI mean (Choice)','fontsize',12)

    ax7 = subplot('position',[.5 .2 .3 .025]);
    imagesc(smooth(RDI_C2)'); axis off
    colormap(ax7,rb); caxis([-1 1])
    title('In-field RDI mean (Choice)','fontsize',12)

           subplot('position',[.1 .15 .3 .025]);
    imagesc(smooth(FL1)'); axis off;  clim([0 1])
    text(-2,1,'L','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.1 .1 .3 .025]);
    imagesc(smooth(FR1)'); axis off;  clim([0 1])
    text(-2,1,'R','fontsize',12,'fontweight','b','COLOR','b')

        subplot('position',[.5 .15 .3 .025]);
    imagesc(smooth(FL2)'); axis off;  clim([0 1])
  text(-2,1,'L','fontsize',12,'fontweight','b','COLOR','R')
    subplot('position',[.5 .1 .3 .025]);
    imagesc(smooth(FR2)'); axis off;  clim([0 1])
    text(-2,1,'R','fontsize',12,'fontweight','b','COLOR','b')



    subplot('position',[.83 .53 .1 .051])
    text(0,0,['Nsp.Corr = ' jjnum2str((UP.Nsp_L(pid)),3)],'fontsize',15,'FontWeight','b')
    axis off

       subplot('position',[.83 .33 .1 .051])
    text(0,0,['Nsp.Corr = ' jjnum2str((UP.Nsp_R(pid)),3)],'fontsize',15,'FontWeight','b')
    axis off

      subplot('position',[.83 .13 .1 .051])
    text(0,0,['Nsp.Corr = ' jjnum2str((UP.Nsp_C(pid)),3)],'fontsize',15,'FontWeight','b')
    axis off

%     subplot('position',[.83 .63 .1 .051])
%     text(0,0,['Max. Nsp.Corr = ' jjnum2str((max([UP.Nsp_L(pid) UP.Nsp_R(pid) UP.Nsp_C(pid)])),3)],'fontsize',15,'FontWeight','b')
%     axis off

    %%
    %     subplot('position',[.1 .3 .3 .3]);

    %
    % for f = 1:min(4,size(idf1,1))
    %     subplot('position',[.15 .4-.03*f .2 .015]);
    %     imagesc(smooth(FRF(2,:,idf1(f))')'); axis off
    %     title(['Field ' num2str(f)],'fontsize',12,'color',hex2rgb(Clist_f{f}),'FontWeight','b')
    %     text(-3,1,'P','FontSize',15,'FontWeight','b')
    %     caxis([0 max((FRF(2:3,:,idf1(f))),[],'all')])
    %
    %         subplot('position',[.55 .4-.03*f .2 .015]);
    %     imagesc(smooth(FRF(3,:,idf1(f))')'); axis off
    %     text(-3,1,'M','FontSize',15,'FontWeight','b')
    %     caxis([0 max((FRF(2:3,:,idf1(f))),[],'all')])
    %
    % %     subplot('position',[.41 .69-.1*f .04 .05])
    % %     if abs(UTF.RDI_LScene(idf1(f)))>0.1, c='r'; else, c='k'; end
    % %     text(0,0,jjnum2str(UTF.RDI_LScene(idf1(f)),2),'fontsize',12,'FontWeight','b','color',c)
    % %     axis off
    % end


    %%

    % for f = 1:min(4,size(idf2,1))
    %     subplot('position',[.5 .70-.1*f .3 .025]);
    %     imagesc(smooth(FRF(2,:,idf2(f))')'); axis off
    %     title(['Field ' num2str(f)],'fontsize',12,'color',hex2rgb(Clist_f{f}),'FontWeight','b')
    %     text(-3,1,'Z','FontSize',15,'FontWeight','b')
    %     caxis([0 max((FRF(2:3,:,idf2(f))),[],'all')])
    %
    %         subplot('position',[.5 .66-.1*f .3 .025]);
    %     imagesc(smooth(FRF(3,:,idf2(f))')'); axis off
    %     text(-3,1,'B','FontSize',15,'FontWeight','b')
    %     caxis([0 max((FRF(2:3,:,idf2(f))),[],'all')])
    %
    %     subplot('position',[.81 .69-.1*f .04 .05])
    %     if abs(UTF.RDI_LScene(idf2(f)))>0.1, c='r'; else, c='k'; end
    %     text(0,0,jjnum2str(UTF.RDI_LScene(idf2(f)),2),'fontsize',12,'FontWeight','b','color',c)
    %     axis off
    % end
    %
    % subplot('position',[.85 .56 .1 .05])
    % text(0,0,['Nsp.Corr = ' jjnum2str((UP.Nsp_L(pid)),3)],'fontsize',12)
    % axis off
    %%


    subplot('position',[.83 .9 .1 .051])
    text(0,0,['Co-React. Prob.'],'fontsize',15,'FontWeight','b')
    axis off

    subplot('position',[.83 .85 .1 .05])
    text(0,0,['= ' jjnum2str(UP.pq(pid)/UP.un(pid),3)],'fontsize',15,'FontWeight','b')
    axis off
    %%
    ROOT.Fig = [ROOT.Mother '\Processed Data\units_mat\ProfilingSheet\U8f_sub'];

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