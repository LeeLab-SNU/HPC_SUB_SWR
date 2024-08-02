Initial_SWRFilter_common;
warning off

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1.xlsx']);
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
exper = {'LSM','JS','SEB'};
TargRegion = 'CA1';

TT_table = table;
tid=1;
%%
for sid =1:size(Session_List,1)
    
    if Session_List.include(sid) & ismember(Session_List.experimenter(sid),exper)
        thisRID = jmnum2str(Session_List.rat(sid),3);
        thisSID = jmnum2str(Session_List.session(sid),2);
        thisRSID = [thisRID '-' thisSID];
        r = find(strcmp(Recording_region.SessionID,[thisRID '-' thisSID]));
        [epochST,epochED] = LoadBehavEpoch(ROOT.Raw.Mother,thisRID,thisSID);
        
        EEG=LoadEEGData(ROOT, [thisRID '-' thisSID], [1:24],Params,Params_Ripple);
        load([ROOT.Behav '\' thisRSID '.mat']);
        ideeg = [find(EEG.timestamps*1e6>epochST,1,'first'),find(EEG.timestamps*1e6<epochED,1,'last')];
        ideeg_ITI=[];
        for s = 1:size(Behav.trial_time,1)-1
            temp = [find(EEG.timestamps>Behav.trial_time(s,6),1,'first'),find(EEG.timestamps<Behav.trial_time(s+1,1),1,'last')];
            if ~isempty(temp)
                if length(temp)>=2
                    ideeg_ITI = [ideeg_ITI;[temp(1):temp(2)]'];
                end
            end
        end
        ideeg_Ch=[];
        for s = 1:size(Behav.trial_time,1)
            if Behav.trial_time(s,4)*Behav.trial_time(s,5)==0, continue; end
            temp = [find(EEG.timestamps>Behav.trial_time(s,4),1,'first'),find(EEG.timestamps<Behav.trial_time(s,5),1,'last')];
            if ~isempty(temp)
                if length(temp)>=2
                    ideeg_Ch = [ideeg_Ch;[temp(1):temp(2)]'];
                end
            end
        end
        %%
        for tt = 1:24
            try
                TT_table.rat(tid) = Session_List.rat(sid);
                TT_table.session(tid) = Session_List.session(sid);
                TT_table.experimenter{tid} = Session_List.experimenter(sid);
                TT_table.type{tid} = Session_List.type(sid);
                TT_table.TT(tid) = tt;
                TT_table.region{tid} = Recording_region.(['TT' num2str(tt)]){r};
                
                if isfield(EEG,['TT' num2str(tt)])
                    thisEEG = EEG.(['TT' num2str(tt)]);
                    EEGr = max(abs(thisEEG.Raw));
                    thisEEG.Raw_B = thisEEG.Raw(ideeg(1):ideeg(2));
                    thisEEG.Raw_ITI = thisEEG.Raw(intersect(ideeg_ITI,[ideeg(1):ideeg(2)]));
                    thisEEG.Raw_Ch = thisEEG.Raw(ideeg_Ch);
                    thisEEG.Envelope_smoothed = thisEEG.Envelope_smoothed (ideeg(1):ideeg(2));
                    
                    TT_table.NoiseRatio(tid) = (sum(thisEEG.Raw_B<=-EEGr) + sum(thisEEG.Raw_B>=EEGr))/length(thisEEG.Raw_B);
                    TT_table.NoiseRatio_ITI(tid) = (sum(thisEEG.Raw_ITI<=-EEGr) + sum(thisEEG.Raw_ITI>=EEGr))/length(thisEEG.Raw_ITI);
                    TT_table.NoiseRatio_Ch(tid) = (sum(thisEEG.Raw_Ch<=-EEGr) + sum(thisEEG.Raw_Ch>=EEGr))/length(thisEEG.Raw_Ch);
                    
                    TT_table.RippleBandMean(tid) = nanmean(thisEEG.Envelope_smoothed);
                    TT_table.RippleBandStd(tid) = nanstd(thisEEG.Envelope_smoothed);
                    
                    MyFolderInfo = dir([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' num2str(tt)]);
                    nc=0;
                    for fid = 1: size(MyFolderInfo,1)
                        name = MyFolderInfo(fid).name;
                        if contains(name,'_beh_SS_') && ~contains(name,'old')
                            nc=nc+1;
                        end
                    end
                    TT_table.NumUnits(tid) =  nc;
                    
                    tid = tid+1;
                end
            end
        end
        disp([thisRID '-' thisSID ' is finished!'])
    end
end

 writetable(TT_table,[ROOT.Info '\TT_table.xlsx'],'writemode','replacefile');
%%
TT_table_c=table;
for tid = 1:size(TT_table)
    if ismember( TT_table.experimenter{tid},exper) && strncmp(TT_table.region{tid},'CA1',3)
        TT_table.ratsession{tid} = [jmnum2str(TT_table.rat(tid),3) '-' jmnum2str(TT_table.session(tid),2)];
        TT_table_c = [TT_table_c;TT_table(tid,:)];
        
    end
end
TT_table_c.session_code = TT_table_c.rat*1e2+TT_table_c.session;
boxplot(TT_table_c.NoiseRatio,TT_table_c.session_code)
xticklabels(unique(TT_table_c.ratsession))
xtickangle(90)
ylabel('noise ratio')