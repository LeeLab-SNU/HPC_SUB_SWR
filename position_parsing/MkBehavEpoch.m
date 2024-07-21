% WriteBehaviorEpoch.csv using Parsed Position

function Epoch = MkBehavEpoch(ROOT,Session_List,sid,Epoch)
        clear Behav
        
        thisRID = Session_List.rat(sid);
        thisSID = Session_List.session(sid);
        try
            loc = [ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\rat' jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
            Behav = load([loc '\ParsedPosition.mat']);
%             d = dir([loc '\Behavior\ParsedPosition.mat']);
%             cd([loc '\Behavior'])
%             copyfile(d.name, [loc '\ParsedPosition.mat'])
            Epoch(thisSID,1) = Behav.t(1);
            Epoch(thisSID,2) = Behav.t(end);
        catch
        end
        
        if size(Session_List,1)==sid
                        xlswrite([ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx'],Epoch)
            disp(['behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx is saved!'])
        elseif Session_List.rat(sid+1)~=thisRID
            xlswrite([ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx'],Epoch)
            disp(['behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx is saved!'])
            Epoch=zeros(15,2);
        end
    end

        