
Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Processed '\ripples_mat\R1'];
ROOT.Save = [ROOT.Processed '\ripples_mat\R2'];
ROOT.Save3 = [ROOT.Processed '\ripples_mat\R3'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];
if ~exist(ROOT.Save), mkdir(ROOT.Save); end
if ~exist(ROOT.Save3), mkdir(ROOT.Save3); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RegionList = {'CA1','SUB'};

for regID = 2:2
thisRegion0 = RegionList{regID};
thisRegion = RegionList{regID};

RippleList_old = readtable([ROOT.Old '\RipplesTable_' thisRegion '.xlsx'],'ReadRowNames',false);



Experimenter = {'JS','LSM','SEB'};
Experimenter2 = {'LSM'};
RipplesTable = RippleList_old;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';

RipplesTable.RippleDuration=(RipplesTable.EDtime-RipplesTable.STtime);
for sid=1:size(RipplesTable,1)
    % if ~ismember(RipplesTable.experimenter{sid},Experimenter), continue; end

    thisSID = [jmnum2str(RipplesTable.rat(sid),3) '-' jmnum2str(RipplesTable.session(sid),2)];
    if ~strcmp(thisSID,thisSID_old)
        load([ROOT.Behav '\' thisSID '.mat']);
        BehavTable = readtable([ROOT.Behav '\' thisSID '.xlsx']);
        thisSID_old = thisSID;
    end
    %% add to RipplesTable
    Idx = knnsearch(Behav.t,RipplesTable.STtime(sid));
    Idx_t = find(Behav.trial_time(:,1)<RipplesTable.STtime(sid));
    
    if strcmp(RipplesTable.experimenter{sid},'LSM')
       if max(Behav.y) >400
           1
       end
    end
    
    RipplesTable.PosX(sid) = Behav.x(Idx);
    RipplesTable.PosY(sid) = Behav.y(Idx);
    RipplesTable.PosY_linearized(sid) = Behav.y_linearized(Idx);
    RipplesTable.Vx(sid) = Behav.velocity.Vx(Idx);
    RipplesTable.Vy(sid) = Behav.velocity.Vy(Idx);
    RipplesTable.speed(sid) = Behav.velocity.speed(Idx);
    
    if ~isempty(Idx_t)
        Idx_t = max(Idx_t);
        RipplesTable.trial{sid} = [thisSID '-' jmnum2str(Idx_t,3)];
        RipplesTable.context(sid) = Behav.trial_context(Idx_t);
        RipplesTable.correctness(sid) = Behav.trial_correctness(Idx_t);
        RipplesTable.ambiguity(sid) = Behav.trial_ambiguity(Idx_t);
        RipplesTable.area(sid) = Behav.trial_vector(Idx,2);
        if strcmp(RipplesTable.experimenter{sid},'JS') && RipplesTable.area(sid) ==5
            RipplesTable.area(sid) = 0;
        end
        RipplesTable.StartTime_fromTrialEnd(sid) = RipplesTable.STtime(sid) - BehavTable.xEnd(Idx_t);
        RipplesTable.EndTime_fromTrialEnd(sid) = RipplesTable.EDtime(sid) - BehavTable.xEnd(Idx_t);
    else
        RipplesTable.trial{sid} = [thisSID '-000'];
        RipplesTable.StartTime_fromTrialEnd(sid) = nan;
        RipplesTable.EndTime_fromTrialEnd(sid) = nan;
    end
% disp([RipplesTable.ID{sid} ' is finished!'])

end
%   RipplesTable = readtable([ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);
writetable(RipplesTable,[ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx'],'WriteMode', 'replacefile');
RipplesTable =  RipplesTable(RipplesTable.context>0 & RipplesTable.correctness==1,:);
RipplesTable_filtered = RipplesTable(RipplesTable.StartTime_fromTrialEnd>0,:);

writetable(RipplesTable_filtered,[ROOT.Save '\RipplesTable_Behav_' thisRegion '_filtered.xlsx'],'WriteMode', 'replacefile');
RipplesTable_sfiltered = RipplesTable_filtered(RipplesTable_filtered.speed<=5,:);

writetable(RipplesTable_sfiltered,[ROOT.Save '\RipplesTable_Behav_' thisRegion '_speed_filtered.xlsx'],'WriteMode', 'replacefile');

RipplesTable_anal = RipplesTable_filtered(RipplesTable_filtered.ensemble>=3,:);
% RipplesTable_anal = RipplesTable_anal(ismember(RipplesTable_anal.experimenter,Experimenter2),:);
writetable(RipplesTable_anal,[ROOT.Save3 '\RipplesTable_Behav_' thisRegion '.xlsx'],'WriteMode', 'replacefile');

% RipplesTable = readtable([ROOT.Save3 '\RipplesTable_Behav_' thisRegion '.xlsx']);
% RipplesTable(RipplesTable.Overlap==0,:)=[];
% writetable(RipplesTable,[ROOT.Save3 '\RipplesTable_Behav_Overlap' thisRegion '.xlsx'],'WriteMode', 'replacefile');
% RipplesTable_anal(RipplesTable_anal.nFields<5,:)=[];
end