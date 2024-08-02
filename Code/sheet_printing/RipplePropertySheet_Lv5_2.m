Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig3 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R10_AI'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = 'Choice';
thisRDI = "RDI_LR";
ROOT.Fig = [ROOT.Fig3 '\' dir];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;
%%
for sid=1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(sid,:);
                if thisRip.rat~=561, continue; end
        thisRID = jmnum2str(thisRip.rat,3);
        thisSID = jmnum2str(thisRip.session,2);

        Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
        diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

        diverging_point = diverging_point*0.23;
        stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;


        Ist = thisRip.STindex-mar;
        Ied = thisRip.EDindex+mar;
        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
        %         if ~strcmp(thisRSID(1:3),'561'), continue; end

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

            EEG = LoadEEGData(ROOT, thisRSID, TargetTT_p,Params,Params_Ripple);

            %             Pos = load([ROOT.Raw.Mother '\rat' thisRSID(1:3) '\rat' thisRSID '\ParsedPosition.mat']);
            clusters = UnitsTable(UnitsTable.rat==thisRip.rat & UnitsTable.session==thisRip.session,:);
            clusters_A = UnitsTable_A(UnitsTable_A.rat==thisRip.rat & UnitsTable_A.session==thisRip.session,:);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            load([ROOT.Rip0 '\' thisRSID '.mat'])
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;
        end


        %% Load units
        spks_epoch=[]; spks_epoch_in=[];
        u=0; Units={};UnitsA={};
        Clist=jet(256);

        clusters_A =  clusters(~isnan(clusters.(thisRDI)),:);
        clusters_A = sortrows(clusters_A,thisRDI,'ascend');
        unz = find(clusters_A.(thisRDI)<0,1,'last');

        cls_all = size(clusters_A,1);
        for un = 1:cls_all
            thisTTID = num2str(clusters_A.TT(un));
            thisCLID = num2str(str2double(clusters_A.ID{un}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
            if ~isempty(thisSpks_in)
                s1 = ones(size(thisSpks,1),1); s2 = ones(size(thisSpks_in,1),1);
                spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters_A.SI(un),...
                    s1*clusters_A.RDI_LScene(un), s1*clusters_A.RDI_RScene(un), s1*clusters_A.RDI_LR(un)]];
                spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters_A.SI(un),...
                    s2*clusters_A.RDI_LScene(un), s2*clusters_A.RDI_RScene(un), s2*clusters_A.RDI_LR(un)]];
                u=u+1;
                Units = [Units; [thisTTID '-' thisCLID]];
                UnitsA = [UnitsA; [clusters_A.ID(un)]];
            end

        end

        [~,ia,~] = unique(spks_epoch(:,3),'rows');
        spks_epoch_u = spks_epoch(ia,:);
        [~,ia,~] = unique(spks_epoch_in(:,3),'rows');
        spks_epoch_u_in = spks_epoch_in(ia,:);


        %% Load FRMap

        FRMap = [];
        for cl=1:size(clusters_A.ID,1)
            thisField = table;
            temp=[];

            thisTTID = num2str(str2double(clusters_A.ID{cl}(8:9)));
            thisCLID =  clusters_A.ID{cl}(11:12);

            thisMap.thisFieldMap1D = load([ROOT.Raw.Map '\rat' thisRSID '-' thisTTID '-' thisCLID  '.mat']);

            for i=1:numel(thisMap.thisFieldMap1D.skaggsMap1D)
                FRMap(i,1:length(thisMap.thisFieldMap1D.skaggsMap1D{i}),cl) = thisMap.thisFieldMap1D.skaggsMap1D{i};
            end

            FRMap(6,:,cl) = thisMap.thisFieldMap1D.skaggsMap_left1D;
            FRMap(7,:,cl) = thisMap.thisFieldMap1D.skaggsMap_right1D;

            FRMap(:,:,cl) = FRMap(:,:,cl) / max(max(FRMap(:,:,cl)));


        end
        for cl=1:size(clusters_A.ID,1)
            if sum(FRMap(1:7,:,cl),'all','omitnan')==0, FRMap(1:7,:,cl)=nan; end
        end

        FRMapsm=[];
        for i=1:size(FRMap,3)
            FRMap(:,:,i)  = FRMap(:,:,i)/max(max(max(FRMap(4:5,:,i))));
            for j=1:7
                FRMapsm(j,:,i) = smooth(FRMap(j,:,i));
            end
        end
        FRMapr=FRMapsm;
        for i=1:size(FRMap,3)
            if ~max(i==spks_epoch_u(:,2))
                FRMapr(:,:,i) = nan;
            end
        end

        stem_end_index = stem_end_index*size(FRMap,2)/48;
        %% Set order
        if ~isempty(UnitsA)
            o0=[]; o1=[]; o2=[];
            o0=[1:length(UnitsA)]';
            [~,ord]=sort(spks_epoch_u_in(o0,1),'descend');
            UnitsC = UnitsA(ord);
        else
            ord = [1:size(FRMap,3)];
            continue;
        end
        if isempty(ord), continue; end
        if length(ord)<=1, continue; end
        %%

        for s=1:size(spks_epoch,1)
            spks_epoch(s,8) = find(ord==spks_epoch(s,3)+1)-1;
        end
        [~,ia] = sort(spks_epoch(:,8));
        spks_epoch = spks_epoch(ia,:);
        [~,ia,~] = unique(spks_epoch(:,3),'rows');
        spks_epoch_u = spks_epoch(ia,:);
        Units = Units(ord);
        UnitsA=UnitsA(ord);

        %% load replay
        Replay = load([ROOT.Rip4 '\' RipplesTable.ID{sid} '.mat']);
        %%
        figure('position',[317,63,1200,915],'color','w');
        % title
        subplot(9,6,3)
        title([cell2mat(thisRip.ID) ', ' Exper ', ' dir],'fontsize',15)
        axis off

        % trial info
        subplot(9,6,6)
        if strcmp(Exper,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
        elseif strcmp(Exper, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
        elseif strcmp(Exper, 'JS'), CxtList = {'Forest','','City'};
            if thisRip.area==5, thisRip.area=0; end
        end
        cxt = CxtList{thisRip.context};
        if thisRip.correctness, corr='Correct'; else, corr='Wrong'; end
        title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',15)
        axis off
        %% EEG

        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(Ist:Ied);
        thisEEG_Ripple = EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied);
        hold on

        subplot(9,4,[1 5])
        plot(thisEEG,'k')
        title(['TT' num2str(TargetTT_p)])
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')
        axis off

        subplot(9,4,[9])
        plot(thisEEG_Ripple,'b')

        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')
        axis off

        %%
        subplot(9,4,[13 21])
        hold on
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur2/2])
        ylim([0 max(spks_epoch(:,3))+1])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','r','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','r','linestyle','--')
        axis off
        for s=1:size(spks_epoch,1)
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+mar/2;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,8)+.2 spks_epoch(s,8)+.2 spks_epoch(s,8)+.8 spks_epoch(s,8)+.8],...
                'k','edgecolor','k')
        end
        %%
        marR = mar/Params.Fs/Params.tbinDuration;
        ax1 = subplot(9,4,[25 33]);
        hold on; axis on
        [pbinN, tbinN] = size(Replay.posterior{1});
        imagesc(ax1,Replay.posterior{1});
        colormap(ax1,flipud(gray))

        plot([-marR tbinN+marR]+0.5,[stem_end_index stem_end_index]+0.5, 'k:');
        line([.5 .5], [0 pbinN], 'color','r','linestyle','--')
        line([tbinN tbinN], [0 pbinN], 'color','r','linestyle','--')

        xlim([-marR tbinN+marR]+0.5);
        ylim([0 pbinN]+0.5);
        set(gca, 'xtick', [0 tbinN]+0.5, 'xticklabel', {'0' [jjnum2str(thisRip.RippleDuration*1e3,1) 'ms']}, 'Fontsize', 10);
        set(gca,  'ytick', [0 stem_end_index pbinN]+0.5, 'yticklabel', {'Stbox', 'Dv', 'Fw'}, 'Fontsize', 10);
        xlabel('time(ms)');
        ylabel('Decoded position');
        %             title(['Bayesian decoding']);

        % plot the linear fitting line which has the best goodness of fit
        c_temp = Replay.c{1}([find(Replay.v{1}>0, 1) find(Replay.v{1}<0,1)]);
        v_temp = Replay.v{1}([find(Replay.v{1}>0, 1) find(Replay.v{1}<0,1)]);

        x = [0 tbinN]+0.5;
        for i = 1 : length(c_temp)
            plot(x,v_temp(i) * x + c_temp(i), '-', 'color', 'k', 'linewidth', 1);
            temp = v_temp(i) * x + c_temp(i);

        end
        if temp(1)<0, temp(1)=0; end
        if temp(2)<0, temp(2)=0; end
        slope = abs(temp(1)-temp(2));
        hold off;

        %% Reactivated cells
        mar = mar*2;
        subplot(9,4,[6 34])
        hold off

        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur-mar]/2); ylim([0 cls_all])
        line([x1 x1], [0 cls_all+1], 'color','r','linestyle','--')
        line([x2 x2], [0 cls_all+1], 'color','r','linestyle','--')
        line([0 dur/2], [unz+1 unz+1], 'color','k')

        for s=1:size(spks_epoch_u,1)
            if spks_epoch_u(s,7)>0,cl='r';elseif spks_epoch_u(s,7)<0,cl='b'; else,cl='k'; end
            patch([0 dur/2 dur/2 0], [spks_epoch_u(s,2)+.2 spks_epoch_u(s,2)+.2 spks_epoch_u(s,2)+.8 spks_epoch_u(s,2)+.8],...
                cl,'facealpha',0.2, 'edgealpha',0)
        end

        for s=1:size(spks_epoch,1)
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+mar/2;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,2)+.2 spks_epoch(s,2)+.2 spks_epoch(s,2)+.8 spks_epoch(s,2)+.8],...
                'k','edgecolor','k')
        end

        UnitsR = Units((spks_epoch_u(:,8)+1));

        xticks([x1 x2]); xticklabels({'0',[jjnum2str(thisRip.RippleDuration*1e3,1) 'ms']})
        ylim([1 size(clusters_A,1)+1]); yticks(spks_epoch_u(:,2)+.5); yticklabels(UnitsR)
        axis on
        mar = mar*1/2;
        %%

        ax2 = subplot(9,4,[7 35]);
        fig = popul_FRMap(flip(squeeze(FRMapr(6,:,:))'),Units,[0 stem_end_index 45],{'Stbox', 'Dv','Fw'},0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        clim([0 1])

        ax3 = subplot(9,4,[8 36]);
        fig = popul_FRMap(flip(squeeze(FRMapr(7,:,:))'),Units,[0 stem_end_index 45],{'Stbox', 'Dv','Fw'},0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        clim([0 1])
        colormap(ax2,'jet')
        colormap(ax3,'jet')

        %%
        saveas(gca,[ROOT.Fig '\' cell2mat(thisRip.ID) '_' dir '.svg'])

        close all
    catch
        close all
    end
end