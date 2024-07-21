Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Processed '\ripples_mat\R1'];
ROOT.Rip = [ROOT.Processed '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Processed '\ripples_mat\R4'];
ROOT.Rip5 = [ROOT.Processed '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Processed '\Manuscript\fig2\R0_fig_2'];
ROOT.Units = [ROOT.Processed '\units_mat\U2'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

unit = 5.0000e-04;
EEGpos = [.05 .96 .9 .035];
microSEC = 1e-06;
len = 20000;
thisFRMapSCALE=2;
%%
thisRegion0 = 'SUB';
thisRegion = 'SUB';
thisRegion2 = 'SUB_field';

% thisRegion0 = 'CA1';
% thisRegion = 'CA1';
% thisRegion2 = 'CA1_field';

Sel = 'LScene';
% RipplesTable = readtable([ROOT.Rip '\RipplesTable_Behav_' thisRegion '.xlsx']);
% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);

RipplesTable = readtable([ROOT.Rip0 '\RipplesTable_' thisRegion0 '.xlsx']);
% RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_a = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);

% RipplesTable(RipplesTable.RippleDuration<0.04,:)=[];

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

for sid=1:size(SessionList,1)
% for sid=32:32
    if ~SessionList.include(sid) | ~strcmp(SessionList.experimenter(sid),'LSM'), continue; end
    thisRIDn = SessionList.rat(sid); thisRID=jmnum2str(thisRIDn,3);
    thisSIDn = SessionList.session(sid); thisSID=jmnum2str(thisSIDn,2);
    thisRSID = [thisRID '-' thisSID];
    load([ROOT.Behav '\' thisRSID '.mat']);
    BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);
    Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);

    theseRip = RipplesTable(RipplesTable.rat==thisRIDn & RipplesTable.session==thisSIDn,:);

    if isempty(theseRip), continue; end


    Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
    diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

    diverging_point = diverging_point*0.23;
    stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

    Recording_region_TT = Recording_region(thisRSID,:);
    TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0.03);

    thisTT_table = TT_table(TT_table.rat==thisRIDn & TT_table.session==thisSIDn,:);
    for t=1:size(thisTT_table,1)
        if ~ismember(thisTT_table.TT(t),TargetTT)
            thisTT_table.TT(t)=0;
        end
    end
    thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
    [~,t] = max(thisTT_table.RippleBandMean);
    TargetTT_p = thisTT_table.TT(t);



    EEG = LoadEEGData(ROOT, thisRSID, TargetTT_p,Params,Params_Ripple);


%%

    for bid=1:size(BehavTable,1)-1

        try
            ITIs = BehavTable.xEnd(bid);
            ITIe = BehavTable.start(bid+1);

            thisRip = RipplesTable(RipplesTable.rat==thisRIDn & RipplesTable.session==thisSIDn&...
                RipplesTable.STtime>ITIs & RipplesTable.STtime<ITIe,:);

                        if isempty(thisRip), continue; end

            Tori = theseRip.STtime(1) - theseRip.STindex(1) * unit;
            EPOCHst = floor((BehavTable.start(bid)-Tori) / unit );
            EPOCHrw = ceil((BehavTable.reward(bid)-Tori )/ unit );
            EPOCHist = ceil((BehavTable.xEnd(bid)-Tori )/ unit );
            EPOCHied = ceil((BehavTable.start(bid+1)-Tori )/ unit );
            EPOCHed = ceil((BehavTable.reward(bid+1)-Tori )/ unit );
            thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw;
            thisEEGr = EEG.(['TT' num2str(TargetTT_p)]).Filtered;

            EPOCHw = [EPOCHst:EPOCHrw,EPOCHist:EPOCHed];

            Emax = max(abs(thisEEG(EPOCHw)))*1.1;
            Emaxr = max(abs(thisEEGr(EPOCHw)))*1.5;

            EPOCHedT = (EPOCHed-EPOCHst)*unit;
            xlist = [0:1000:EPOCHrw-EPOCHst]; xlistS={};
            for x=1:length(xlist),xlistS{x}=num2str(xlist(x)*unit); end
            xlist2 = [0:1000:EPOCHied-EPOCHist]; xlistS2={};
            for x=1:length(xlist2),xlistS2{x}=num2str(xlist2(x)*unit); end
            xlist3 = [0:1000:EPOCHed-EPOCHied]; xlistS3={};
            for x=1:length(xlist3),xlistS3{x}=num2str(xlist3(x)*unit); end


             P_base = [EPOCHist:EPOCHied];
             for r=1:size(thisRip,1)
                 P_base(P_base>=thisRip.STindex(r) & P_base<=thisRip.EDindex(r))=[];
             end

             Amp_base = nanmean(abs(thisEEGr(P_base)));
            %%
            clusters = UnitsTable(UnitsTable.rat==thisRIDn & UnitsTable.session==thisSIDn,:);

            cls_all = size(clusters,1);


            FRMaps = LoadFRMap(ROOT,clusters);

            for c=1:size(clusters.ID,1)
                [~,m] = nanmax(FRMaps(1,:,c));
                clusters.MaxBin(c)=m;
            end

            [~,j]=sort(clusters.MaxBin);
            clusters = clusters(j,:);
            %%
            figure('position',[100,100,1907,876]);
            subplot(5,6,1)
            plot(thisEEG(EPOCHst:EPOCHrw),'color','k')


            xticks(xlist)
            xticklabels(xlistS)
            xlim([0 xlist(end)])

            ylabel('Amplitude (uV)')
            ylim([-Emax Emax])

            subplot(5,6,[2 5])
            plot(thisEEG(EPOCHist:EPOCHied),'color','k')
            %         line([1 1]*(BehavTable.xEnd(bid)-Tori) / unit - EPOCHist,[-1 1]*Emax,'color','r')
            %         line([1 1]*(BehavTable.start(bid+1)-Tori) / unit - EPOCHist,[-1 1]*Emax,'color','r')

            %         line([1 1]*(BehavTable.reward(bid)-Tori) / unit - EPOCHst,[-2500 2500],'color','b')
            for rid=1:size(thisRip,1)
                if thisRip.ensemble(rid)>=3,c='y';else,c='k';end

                NormAmp_v3 = nanmean(abs(thisEEGr(thisRip.STindex(rid):thisRip.EDindex(rid))))/Amp_base;

                
                x1 = (thisRip.STtime(rid)-Tori)/ unit - EPOCHist; x2 = (thisRip.EDtime(rid)-Tori)/ unit - EPOCHist;
                patch([x1 x2 x2 x1], [-1 -1 1 1]*Emax,c,'facealpha',0.3,'edgealpha',0)
                text(x1, Emax*round(mod(rid,2)-0.5), thisRip.ID{rid}(end-3:end),'color','r')
%                 text(x1, Emax*round(mod(rid,2)-0.5)*0.8, jjnum2str(RipplesTable.NormFilteredPeak(rid),2),'color','b')
                text(x1, Emax*round(mod(rid,2)-0.5)*0.6, jjnum2str(NormAmp_v3,2),'color',[.5 0 1])

                Rid = find(strcmp(RipplesTable.ID,thisRip.ID(rid)));
                RipplesTable.NormAmp_v3(Rid) = NormAmp_v3;

            end

            xticks(xlist2)
            xticklabels(xlistS2)
            xlim([0 xlist2(end)])
            ylim([-Emax Emax])

            subplot(5,6,6)
            plot(thisEEG(EPOCHied:EPOCHed),'color','k')
            %         line([1 1]*(BehavTable.xEnd(bid)-Tori) / unit - EPOCHist,[-1 1]*Emax,'color','r')
            %         line([1 1]*(BehavTable.start(bid+1)-Tori) / unit - EPOCHist,[-1 1]*Emax,'color','r')

            %         line([1 1]*(BehavTable.reward(bid)-Tori) / unit - EPOCHst,[-2500 2500],'color','b')


            xticks(xlist3)
            xticklabels(xlistS3)
            xlim([0 xlist3(end)])
            ylim([-Emax Emax])
            %%
            subplot(5,6,7)
            plot(thisEEGr(EPOCHst:EPOCHrw),'color','b')
            xticks(xlist)
            xticklabels(xlistS)
            xlim([0 xlist(end)])
            ylabel('Amplitude (uV)')
            ylim([-Emaxr Emaxr])


            subplot(5,6,[8 11])
            plot(thisEEGr(EPOCHist:EPOCHied),'color','b')

            for rid=1:size(thisRip,1)
                if thisRip.ensemble(rid)>=3,c='y';else,c='k';end
                x1 = (thisRip.STtime(rid)-Tori)/ unit - EPOCHist; x2 = (thisRip.EDtime(rid)-Tori)/ unit - EPOCHist;
                patch([x1 x2 x2 x1], [-1 -1 1 1]*Emaxr,c,'facealpha',0.3,'edgealpha',0)
            end

            xticks(xlist2)
            xticklabels(xlistS2)
            xlim([0 xlist2(end)])
            ylim([-Emaxr Emaxr])

            subplot(5,6, 12)
            plot(thisEEGr(EPOCHied:EPOCHed),'color','b')

            xticks(xlist3)
            xticklabels(xlistS3)
            xlim([0 xlist3(end)])
            ylim([-Emaxr Emaxr])
            %%
            subplot(5,6,13)
            thisPos = (Pos.t>=EPOCHst*unit+Tori & Pos.t<=EPOCHrw*unit+Tori);
       
            plot(smooth(max(Pos.y_linearized) - Pos.y_linearized(thisPos),0.2),'k','linewidth',2)
            line([0 sum(thisPos)], [1 1]*max(Pos.y_linearized) -diverging_point,'color','b')
            ylabel('Pos (cm)')
            xlim([0 sum(thisPos)]); ylim([0 max(Pos.y_linearized)-min(Pos.y_linearized)]*1.2)

            subplot(5,6,[14 17])
            thisPos = (Pos.t>=EPOCHist*unit+Tori & Pos.t<=EPOCHied*unit+Tori);
            plot(smooth(max(Pos.y_linearized) - Pos.y_linearized(thisPos),0.2),'k','linewidth',2)
            line([0 sum(thisPos)], [1 1]*max(Pos.y_linearized) -diverging_point,'color','b')
            xlim([0 sum(thisPos)]); ylim([0 max(Pos.y_linearized)-min(Pos.y_linearized)]*1.2)

            subplot(5,6,18)
            thisPos = (Pos.t>=EPOCHied*unit+Tori & Pos.t<=EPOCHed*unit+Tori);
            plot(smooth(max(Pos.y_linearized) - Pos.y_linearized(thisPos)),'k','linewidth',2)
            line([0 sum(thisPos)], [1 1]*max(Pos.y_linearized) -diverging_point,'color','b')
            xlim([0 sum(thisPos)]); ylim([0 max(Pos.y_linearized)-min(Pos.y_linearized)]*1.2)

            %%
            subplot(5,6,[19 25]); hold on
            ylistS={};
            for un = 1:cls_all

                [thisRIDc,thisSIDc,thisTTIDc,thisCLIDc, thisFLIDc] = parsing_clusterID(clusters.ID{un},1);

                Spk = Spike.(['TT' (thisTTIDc)]).(['Unit' thisCLIDc]);

                thisSpks = Spk.t_spk(Spk.t_spk>=EPOCHst*unit+Tori & Spk.t_spk<=EPOCHrw*unit+Tori);

                if clusters.(['Selectivity_' Sel])(un)>=1
                    if (clusters.(['RDI_' Sel])(un))>0
                        c='r';
                    else
                        c='b';
                    end
                else
                    c='k';
                end

                scatter(thisSpks-(EPOCHst*unit+Tori),ones(length(thisSpks),1)*un,c,'|')

                ylistS{un} = [(thisTTIDc) '-' thisCLIDc];
            end
            xlim([0 str2double(xlistS{end})])
            yticks([1:un])
            yticklabels(ylistS)
            ylim([0 un+.5])
            ylabel('Unit')

            subplot(5,6,[20 29]); hold on
            ylistS={};
            for un = 1:cls_all

                [thisRIDc,thisSIDc,thisTTIDc,thisCLIDc, thisFLIDc] = parsing_clusterID(clusters.ID{un},1);

                Spk = Spike.(['TT' (thisTTIDc)]).(['Unit' thisCLIDc]);

                thisSpks = Spk.t_spk(Spk.t_spk>=EPOCHist*unit+Tori & Spk.t_spk<=EPOCHied*unit+Tori);

                if clusters.(['Selectivity_' Sel])(un)>=1
                    if (clusters.(['RDI_' Sel])(un))>0
                        c='r';
                    else
                        c='b';
                    end
                else
                    c='k';
                end

                scatter(thisSpks-(EPOCHist*unit+Tori),ones(length(thisSpks),1)*un,c,'|')

                ylistS{un} = [(thisTTIDc) '-' thisCLIDc];
            end
            xlim([0 str2double(xlistS2{end})])

            for rid=1:size(thisRip,1)
                if thisRip.ensemble(rid)>=3,c='y';else,c='k';end
                x1 = thisRip.STtime(rid)-(EPOCHist*unit+Tori); x2 = thisRip.EDtime(rid)-(EPOCHist*unit+Tori);
                patch([x1 x2 x2 x1], [0 0 1 1]*(un+.5),c,'facealpha',0.3,'edgealpha',0)
            end

            yticks([1:un])
            yticklabels(ylistS)
            ylim([0 un+.5])



            subplot(5,6,[24 30]); hold on
            ylistS={};
            for un = 1:cls_all

                [thisRIDc,thisSIDc,thisTTIDc,thisCLIDc, thisFLIDc] = parsing_clusterID(clusters.ID{un},1);

                Spk = Spike.(['TT' (thisTTIDc)]).(['Unit' thisCLIDc]);

                thisSpks = Spk.t_spk(Spk.t_spk>=EPOCHied*unit+Tori & Spk.t_spk<=EPOCHed*unit+Tori);

                if clusters.(['Selectivity_' Sel])(un)>=1
                    if (clusters.(['RDI_' Sel])(un))>0
                        c='r';
                    else
                        c='b';
                    end
                else
                    c='k';
                end

                scatter(thisSpks-(EPOCHied*unit+Tori),ones(length(thisSpks),1)*un,c,'|')

                ylistS{un} = [(thisTTIDc) '-' thisCLIDc];
            end
            xlim([0 str2double(xlistS3{end})])

            yticks([1:un])
            yticklabels(ylistS)
            ylim([0 un+.5])
            %%
            CxList = {'Zebra','Pebbles','Bamboo','Mountains'};
            CrList={'C','W'};
            cx1 = CxList{BehavTable.context(bid)}; cx2 = CxList{BehavTable.context(bid+1)};
            cr1 = CrList{BehavTable.correctness(bid)}; cr2 = CrList{BehavTable.correctness(bid+1)};
            sgtitle([thisRegion0 ' - Trial ' BehavTable.TrialID{bid} ' (' cx1 '(' cr1 ') -' cx2 '(' cr2 '))'],...
                'interpreter','none')


            ROOT.Fig_en = [ROOT.Fig '\' thisRegion0 ];
            if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
            saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.svg'])
            saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.png'])

            if ~isempty(thisRip)
                ROOT.Fig_en = [ROOT.Fig '\' thisRegion0 ];
                if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
                saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.svg'])
                saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.png'])
            end

            close all
        catch
            close all
        end
    end

end
