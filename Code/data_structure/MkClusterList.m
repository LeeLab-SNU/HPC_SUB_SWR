Initial_SWRFilter_common;

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_list = table;
exper = {'LSM'};
TargRegion = 'SUB';



for sid = 1:size(Session_List,1)
    if Session_List.include(sid) && ismember(Session_List.experimenter{sid},exper) && Session_List.rat(sid)==295&& Session_List.session(sid)==2
        thisRID = Session_List.rat(sid);
        thisSID = Session_List.session(sid);
        r = find(strcmp(Recording_region.SessionID,[jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)]));
        for thisTTID = 1:24
            loc = [ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\rat' jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '\' 'TT' num2str(thisTTID)];
            %             trans_model_forJS(loc,thisTTID)
            fd = dir(loc);
            for fid = 1: size(fd,1)
                name = fd(fid).name;
                                 TransWinCltoNtt(thisRID,thisSID,thisTTID,name,ROOT)
                if contains(name,'_beh_SS_') && ~contains(name,'old')
                    temp = table;
                    findHYPEN = find(name == '_');
                    thisUID = str2double(name(findHYPEN(end)+1:end-4));
                    temp.ID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' jmnum2str(thisTTID,2) '-' jmnum2str(thisUID,2)];
                    if ~isempty(Cluster_list), if ismember({temp.ID},Cluster_list.ID), continue; end, end
                    temp.rat = thisRID;
                    temp.session = thisSID;
                    temp.TT = thisTTID;
                    temp.region{1} = Recording_region.(['TT' num2str(thisTTID)]){r};
                    temp.experimenter{1} = Session_List.experimenter{sid};
                    temp.date = Session_List.date(sid);
                    temp.session_type{1} = Session_List.type{sid};
                    Cluster_list = [Cluster_list; temp];
                end
                
            end
        end
    end
end

% writetable(Cluster_list,[ROOT.Info '\ClusterList_SWR.xlsx'])
% Cluster_list_CA1 = Cluster_list(strncmp(TargRegion,Cluster_list.region,3),:);
% writetable(Cluster_list_CA1,[ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx'])

function trans_model_forJS(loc,thisTTID)
fd = dir(loc);
            for fid = 1: size(fd,1)
                name = fd(fid).name;
                if (contains(name, ['T' num2str(thisTTID) '.']) && ~contains(name, 'TT')) && str2double(name(end))<100
                    cd(loc)
                d = dir([loc '\' name]);
                movefile(d.name, ['TT' num2str(thisTTID) '_cluster' name(end-1:end)])
                end
            end
end

function TransWinCltoNtt(thisRID,thisSID,thisTTID,name,ROOT)
if contains(name,'_cluster')
    try
        findDOT = find(name == '.');
        clusterID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' num2str(thisTTID) '-' jmnum2str(str2double(name(findDOT+1:end)),2)];
        cl2Ntt(ROOT.Raw.Mother, clusterID, clusterID);
    catch
    end
end
end