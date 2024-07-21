Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Old2 = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Fig1 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R1_sub'];
ROOT.Fig2 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R2_sub'];
ROOT.Fig3 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R3_sub'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U0'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Fig1), mkdir(ROOT.Fig1); end
if ~exist(ROOT.Fig2), mkdir(ROOT.Fig2); end
if ~exist(ROOT.Fig2), mkdir(ROOT.Fig3); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);



thisRegion = 'SUB';

Cluster_List = readtable([ROOT.Units '\ClusterList_SWR_' thisRegion '_filtered.xlsx']);
RipplesTable_ori = readtable([ROOT.Old2 '\RipplesTable_Behav_' thisRegion '.xlsx']);
RipplesTable = readtable([ROOT.Old2 '\RipplesTable_Behav_' thisRegion '_speed_filtered.xlsx'],'ReadRowNames',false);
RipplesTable = RipplesTable(RipplesTable.ensemble>=3,:);
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
%%
for sid=1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(sid,:);
%         if thisRip.rat~=80, continue; end
        Ist = thisRip.STindex-mar;
        Ied = thisRip.EDindex+mar;
        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];

        
        if ~strcmp(thisRSID, thisRSID_old)
            Recording_region_TT = Recording_region(thisRSID,:);
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0.03);
            
            thisTT_table = TT_table(TT_table.rat==thisRip.rat & TT_table.session==thisRip.session,:);
            for t=1:size(thisTT_table,1)
                if ~ismember(thisTT_table.TT(t),TargetTT)
                    thisTT_table.TT(t)=0;
                end
            end
            thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
            [~,t] = max(thisTT_table.RippleBandMean);
            TargetTT_p = thisTT_table.TT(t);
            
          EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
            NumTT = length(fieldnames(EEG));
            Pos = load([ROOT.Raw.Mother '\rat' thisRSID(1:3) '\rat' thisRSID '\ParsedPosition.mat']);
            clusters = Cluster_List(Cluster_List.rat==thisRip.rat & Cluster_List.session==thisRip.session,:);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            load([ROOT.Old '\' thisRSID '.mat'])
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;
        end
        
        figure('position',[317,63,923,915]);
        % title
        subplot(5+fix(NumTT/6),6,1)
        title([thisRip.ID ' ' Exper],'fontsize',15)
        axis off
        
        % trial info
        subplot(4+fix(NumTT/6),6,7)
        if strcmp(Exper,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
        elseif strcmp(Exper, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
        elseif strcmp(Exper, 'JS'), CxtList = {'Forest','City'};
            if thisRip.area==5, thisRip.area=0; end
        end
        cxt = CxtList{thisRip.context};
        if thisRip.correctness, corr='Correct'; else, corr='Wrong'; end
        title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',15)
        axis off
        
        % position
        subplot(5+fix(NumTT/6),6,[5 6 23 24])
        if thisRip.speed>5 || thisRip.area~=0, cl='k'; else, cl='r'; end
        if strcmp(cell2mat(thisRip.experimenter),'JS')
            x=Pos.y; y=Pos.x+300*(~Pos.area(:,1));
        else, x=smooth(Pos.x); y=smooth(Pos.y); end
        scatter(x,y,5,[.7 .7 .7],'filled')
        hold on
        if strcmp(cell2mat(thisRip.experimenter),'JS')
            x=thisRip.PosY; y=thisRip.PosX+300*(~thisRip.area);
            line([0.02 0.04],[300 300], 'color','k');
            ylim([0 330]); xlim([0 .06])
        else, x=thisRip.PosX; y=thisRip.PosY; end
        scatter(x,y,20,cl,'filled')
        title([jjnum2str(thisRip.speed,2) 'cm/s'],'color',cl)
        axis off
        
        % LFP trace
        t2=1; RP1=[]; RP2=[];
        for t1=1:24
            if ~ismember(t1, TargetTT), continue; end
            subplot(5+fix(NumTT/6),6,24+t2)
            axis off
            try
                plot(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied),'k')
                hold on
                if ismember(t1,RippleTT)
                    cl='r';
                    rps = Ripples.(['TT' num2str(t1)]).index;
                    ix = find((rps(:,1)>Ist & rps(:,1)<Ied) | (rps(:,2)>Ist & rps(:,2)<Ied));
                    if ~isempty(ix)
                        
                        for i = 1:size(ix,1)
                            plot([rps(ix(i),1)-Ist:rps(ix(i),2)-Ist],EEG.(['TT' num2str(t1)]).Raw(rps(ix(i),1):rps(ix(i),2)),'r')
                        end
                        
                        RP1(t1) = nanmax(abs(EEG.(['TT' num2str(t1)]).Filtered(thisRip.STindex:thisRip.EDindex)));
                        RP2(t1) = max(rps(ix,2) - rps(ix,1));
                    end
                else, cl='k'; end
                %                 plot(EEG.(['TT' num2str(t1)]).Filtered(Ist:Ied),'b')
                
                x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
                line([x1 x1], [min(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)) max(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied))+50], 'color','k','linestyle','--')
                line([x2 x2], [min(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)) max(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied))+50], 'color','k','linestyle','--')
                title(['TT' num2str(t1)],'color',cl)
                text(x1+1,min(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)),jjnum2str(min(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)),2))
                text(x2+1,max(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)),jjnum2str(max(EEG.(['TT' num2str(t1)]).Raw(Ist:Ied)),2))
                
                axis off
                t2 = t2+1;
                xlim([0 Ied-Ist])
            end
        end
        
        
        % reactivated cells
        %         [~,m] = max(abs(RP1));
        subplot(5+fix(NumTT/6),6,[3,4])
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(Ist:Ied);
        plot(thisEEG,'k')
        %         hold on
        %         plot(EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied),'b')
        
        title(['TT' num2str(TargetTT_p)])
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        xlim([0 dur])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        axis off
        
        %         RP2(m)=0;
        %         [~,m] = max(abs(RP2));
        subplot(5+fix(NumTT/6),6,[9,10])
        x1=0; x2=thisRip.RippleDuration*2e3;
        %         thisEEG = EEG.(['TT' num2str(m)]).Raw(Ist:Ied);
        %         plot(thisEEG,'k')
        %         hold on
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied);
        plot(thisEEG ,'b')
        
        %         title(['TT' num2str(m)])
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        xlim([0 dur])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        axis off
        
        
        subplot(5+fix(NumTT/6),6,[15,16,21,22])
        spks_epoch=[];u=0; Units={};
        cls = size(clusters,1);
        for un = 1:cls
            thisTTID = num2str(clusters.TT(un));
            thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
            
            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            if ~isempty(thisSpks)
                spks_epoch = [spks_epoch;[thisSpks,ones(size(thisSpks,1),1)*un,ones(size(thisSpks,1),1)*u]];
                u=u+1;
            end
            Units = [Units; [thisTTID '-' thisCLID]];
        end
        
        hold on
        for s=1:size(spks_epoch,1)
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+100;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],'k')
        end
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur/2])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        text(x1+1,0,'0')
        text(x2+1,0,[num2str(thisRip.RippleDuration*1e3) 'ms'])
        axis off
        
        fd = [ROOT.Fig1 '\' ['rat' jmnum2str(thisRip.rat,3)] '\rat' thisRSID];
        if ~exist(fd), mkdir(fd); end
        saveas(gca,[fd '\' cell2mat(thisRip.ID) '.png'])
        if thisRip.speed<=5 && thisRip.area==0
            fd = [ROOT.Fig2 '\' ['rat' jmnum2str(thisRip.rat,3)]];
            if ~exist(fd), mkdir(fd); end
            saveas(gca,[fd '\' cell2mat(thisRip.ID) '.png'])
        end
        close all
    catch
        close all
    end
    
end