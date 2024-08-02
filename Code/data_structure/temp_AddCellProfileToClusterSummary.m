warning off
Initial_SWRFilter_common;


Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);
TargRegion = 'CA1';

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx']);
exper = {'LSM','JS','SEB'};
exper = {'LSM'};


for cid = 1:size(Cluster_List,1)
    if ismember(Cluster_List.experimenter{cid}, exper)
        try
        clusterID = Cluster_List.ID{cid};
 
 [nSPKS, FRRate,withinREFRACPortion, max_width, max_peak, max_amp, peak_ratio, LRATIO, ISODIST,  LogISIPEAKTIME,...
     valley_slope, valley_proportion,SpaInfoScore,onmazeAvgFR,onmazeMaxFR,offmazeAvgFR] = JMGetClusterQuals_lite(clusterID, ROOT,Cluster_List.experimenter{cid});
    
 Cluster_List.nSpks(cid) = nSPKS;
 Cluster_List.withinRef(cid) = withinREFRACPortion;
    Cluster_List.SpkWidth(cid) = max_width;
    Cluster_List.LogISIPeakTime(cid) = LogISIPEAKTIME;
    Cluster_List.AvgPeaktoValley(cid) = max_amp;
    Cluster_List.PeakfromBaseline(cid) = max_peak;
    Cluster_List.PeakRatio(cid) = peak_ratio;
    Cluster_List.AvgFR(cid) = FRRate;
    Cluster_List.onMazeAvgFR(cid) = onmazeAvgFR;
    Cluster_List.onMazeMaxFR(cid) = onmazeMaxFR;
    Cluster_List.offMazeAvgFR(cid) = offmazeAvgFR;
     Cluster_List.SI(cid) = SpaInfoScore;
        catch
            fprintf('\n%s is failed\n', clusterID);
        end
    end
end


%%
Cluster_List(~strcmp(Cluster_List.experimenter,'LSM'),:)=[];
Cluster_List(~strcmp(Cluster_List.region,TargRegion),:)=[];
writetable(Cluster_List,[ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx'],'WriteMode', 'replacefile')
Cluster_List_CA1 = Cluster_List;
% average peak-to-valley amplitude of waveforms ≥ 75uV 
id = Cluster_List_CA1.AvgPeaktoValley>=75;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% proportion of spikes within a 1ms refractory period < 1% (of total spikes)
id = Cluster_List_CA1.withinRef<1;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% average firing rate during the outbound journey on the stem and arms ≥ 1Hz
id = Cluster_List_CA1.AvgFR>=0.5;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% fast-spiking neurons (mean firing rate ≥ 10 Hz; width of the average waveform < 325us) were excluded
id = Cluster_List_CA1.AvgFR>=10 & Cluster_List_CA1.SpkWidth<325;
Cluster_List_CA1 = Cluster_List_CA1(~id,:);


writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_' TargRegion '_filtered.xlsx'],'WriteMode', 'replacefile')

