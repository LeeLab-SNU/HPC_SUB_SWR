
Initial_SWRFilter_common

thisRegion = 'SUB';
thisRegion2 = 'SUB';

ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Unit = [ROOT.Save '\units_mat\U2'];
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_field_RDIs_UV.xlsx']);

ROOT.Rip = [ROOT.Save '\ripples_mat\R4'];
ROOT.Fig = [ROOT.Save '\ripples_mat\ProfilingSheet\R4_' thisRegion '(Replay)'];
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


%% decoding parameters
Params.tbinDuration = 0.01;
Params.N_threshold = 3;
Params.fitting_threshold = 1.5;
a_test = 0.05;

%% Ripples' replay table
% RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion '_ReplayP.xlsx']);
pbinN=ones(1,5);

thisFRMapSCALE=2;

for RipNum = 1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(RipNum,:);
        RipID = RipplesTable.ID{RipNum};
        if exist([ROOT.Rip '\' RipID  '.mat'])
            load([ROOT.Rip '\' RipID  '.mat'])
            thisRID = jmnum2str(thisRip.rat,3);
            thisSID = jmnum2str(thisRip.session,2);

            Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
            diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

            diverging_point = diverging_point*0.23;
            stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

            %% Display Title
            figure('position',[317,63,923,600],'color','w');
            % title
            subplot(8,6,1)
            title([cell2mat(thisRip.ID) ', ' thisRip.experimenter],'fontsize',15)
            axis off

            % trial info
            subplot(9,6,6)
            if strcmp(thisRip.experimenter,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
            elseif strcmp(thisRip.experimenter, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
            elseif strcmp(thisRip.experimenter, 'JS'), CxtList = {'Forest','','City',''};
                if thisRip.area==5, thisRip.area=0; end
            end
            cxt = CxtList{thisRip.context};
            if thisRip.correctness, corr='Correct'; else, corr='Wrong'; end
            title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',10)
            axis off
            sceneName=[{'Overall'},CxtList(:)'];
            %% display reconstructed position probability
            pbinN=1;
            for mapRUN = 1 : length(pbinN)

                temp_position = [0.33 0.4 0.12 0.4] + [0.165 0 0 0] * (mapRUN-2);
                ax = subplot('position', temp_position);

                Replay_display_posterior(posterior{mapRUN}, v{mapRUN}, c{mapRUN}, R_actual(mapRUN), R_shuffled{mapRUN},[], ...
                    p_test(mapRUN), 0,stem_end_index, Params.tbinDuration, sceneName{mapRUN}, ax, a_test,[0 .3 0 .2]);

            end
% saveas(gca,[ROOT.Fig '\' cell2mat(thisRip.ID) '.png'])

        close all
        end
    catch
        close all
    end
end