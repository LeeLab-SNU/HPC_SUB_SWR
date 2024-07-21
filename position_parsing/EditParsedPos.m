function Behav = EditParsedPos(ROOT,SessionList,sid)
thisSID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];

%% edit Behav structure
Behav = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\Behavior\' 'ParsedPosition.mat']);

if exist([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.csv']])
    Behav.Epoch = csvread([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.csv']]);
else
    Behav.Epoch = table2array(readtable([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.xlsx']]));
end
Behav.SessionST = Behav.Epoch(str2double(thisSID(5:6)),1);
Behav.SessionED = Behav.Epoch(str2double(thisSID(5:6)),end);
Behav.total_trial_number = size(Behav.trial,2);

if strcmp(SessionList.experimenter{sid},'SEB')
    px2cm=0.23;
    if max(Behav.y)>300, Behav.x=Behav.x*px2cm; Behav.y=Behav.y*px2cm; end
    if max(Behav.y)<80, Behav.x=Behav.x/px2cm; Behav.y=Behav.y/px2cm; end
    t=1; cx=360*px2cm; cy =470*px2cm;
    while sum(Behav.trial(t,:))==0
        t=t+1;
    end
    for i=t:size(Behav.trial,1)
        if sum(Behav.trial(i,:))==0
            if Behav.y(i)>370*px2cm
                Behav.trial(i,:) = 0;
                Behav.area(i,1:5)=0;
            else
                Behav.trial(i,:) = Behav.trial(i-1,:);
                Behav.area(i,5)=1;
            end
        else
            if Behav.y(i)<=220*px2cm
                Behav.area(i,4)=1;
            else
                Behav.area(i,1)=1;
            end
        end
    end
elseif strcmp(SessionList.experimenter{sid},'JS')
            px2cm=0.3;
    if max(Behav.y)>900, Behav.x=Behav.x*px2cm; Behav.y=Behav.y*px2cm; end
    cx=1000*px2cm; cy =1;
    if median(Behav.t(2:end) - Behav.t(1:end-1))<0.01
        if median(Behav.t(2:end) - Behav.t(1:end-1))<10^(-7)
            Behav.t = Behav.t * 10^6;
        else
            Behav.t = Behav.t * 10;
        end
    elseif median(Behav.t(2:end) - Behav.t(1:end-1)) >10000
        Behav.t = Behav.t * 10^-6;
    elseif median(Behav.t(2:end) - Behav.t(1:end-1)) >1000
        Behav.t = Behav.t * 10^-5;
    end
elseif strcmp(SessionList.experimenter{sid},'LSM')
        px2cm=0.23;
        if max(Behav.y)>300, Behav.x=Behav.x*px2cm; Behav.y=Behav.y*px2cm; end
    if max(Behav.y)<70, Behav.x=Behav.x/px2cm; Behav.y=Behav.y/px2cm; end
    cx=360*px2cm; cy =470*px2cm;
    Behav.x(Behav.x==0) = cx;
    Behav.y(Behav.y==0) = cy;
end
if ~isfield(Behav,'area')
    Behav.area(:,1) = sum(Behav.trial,2);
    Behav.area(:,5) = sum(Behav.trial,2)==0;
elseif isempty(Behav.area)
    Behav.area(:,1) = sum(Behav.trial,2);
    Behav.area(:,5) = sum(Behav.trial,2)==0;
elseif strcmp(SessionList.experimenter{sid},'JS')
    if sum(Behav.area(:,5))>1000
        Behav.area(:,5) = 0;
        for i=2:size(Behav.area,1)
            if Behav.area(i-1,1)==1 && Behav.area(i,1)==0
                Behav.area(i-1,5)=1;
            end
        end
        Behav.area(Behav.area(:,5)==1,1)=0;
    end
    
end

if ~isfield(Behav,'correctness')
    
    if isfield(Behav,'corr')
      Behav.correctness=  Behav.corr;
    end
end
if ~isfield(Behav,'trial_context')
    Behav.trial_context=[];
    for t=1:size(Behav.trial,2)
        i= min(find(Behav.trial(:,t)));
        if i>0
            Behav.trial_context(t,1) = find(Behav.sc(i,:));
        else
            Behav.trial_context(t,1) = 0;
        end
    end
end
if ~isfield(Behav,'trial_side')
    Behav.trial_side=[];
    for t=1:size(Behav.trial,2)
        i= min(find(Behav.trial(:,t)));
        if i>0
            Behav.trial_side(t,1) = find(Behav.side(i,:));
        else
            Behav.trial_side(t,1) = 0;
        end
    end
end
if ~isfield(Behav,'trial_ambiguity')
    Behav.trial_ambiguity=[];
    for t=1:size(Behav.trial,2)
        i= min(find(Behav.trial(:,t)));
        if i>0
            Behav.trial_ambiguity(t,1) = find(Behav.amb(i,:));
        else
            Behav.trial_ambiguity(t,1) = 0;
        end
    end
end

if ~isfield(Behav,'trial_correctness')
    Behav.trial_correctness=[];
    if isfield(Behav,'corr')
        for t=1:size(Behav.trial,2)
            i= min(find(Behav.trial(:,t)));
            if i>0
                Behav.trial_correctness(t,1) = find(Behav.corr(i,:));
            else
                Behav.trial_correctness(t,1) = 0;
            end
        end
    else
        Behav.trial_correctness(size(Behav.trial,2),:) = 1;
    end
end


for i=1:Behav.total_trial_number
    if isempty(find(Behav.trial(:,i))), Behav.trial_time(i,:) =0; continue; end
    for j=1:5
        if isempty(Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1)))
            Behav.trial_time(i,j) =0;
        else
            Behav.trial_time(i,j) = Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1));
        end
        
    end
    Behav.trial_time(i,6) = Behav.t(min(find(Behav.trial(:,i),1,'last')+1,size(Behav.t,1)));
end
Behav.trial_time(Behav.trial_time(:,1)==0,:) = nan;
for i=1:6
    Behav.trial_time(:,i) = fillmissing(Behav.trial_time(:,i),'linear','SamplePoints',[1:size(Behav.trial_time,1)]);
end

try
    Behav.y_linearized = get_linearized_position(ROOT.Mother,ROOT.Raw.Mother,thisSID);
    Behav.diverging_point = get_divergingPoint(ROOT.Info, thisSID(1:3), thisSID(5:6), 'outbound')*px2cm;
catch
    Behav.y_linearized = Behav.y;
    Behav.diverging_point = floor(min(Behav.y));
end

Behav.trial_vector(:,1) = interp1(Behav.trial_time(:,1),[1:size(Behav.trial_time,1)]',Behav.t, 'previous','extrap');
[~,Behav.trial_vector(:,2)] = max(Behav.area,[],2);
for i=1:size(Behav.trial,1)
    if max(Behav.area(i,:))==0
        Behav.trial_vector(i,2)=0;
    end
end

%% position interpolation

InTrack = Behav.trial_vector(:,2)~=0;

Behav.x = interp_pos(Behav.x,InTrack,cx);
Behav.y = interp_pos(Behav.y,InTrack,cy);

Behav.y_linearized(~InTrack) = cy;
[Behav.y_linearized,~] = fillmissing(Behav.y_linearized,'linear');

%% cal speed
Behav.velocity=table;
bin_size = 5;


v = cal_speed(Behav.t, [Behav.x,Behav.y], bin_size);
if strcmp(SessionList.experimenter{sid},'JS')
    v(v(:,1)<0,1) = 0;
end

Behav.velocity.Vx=v(:,1);
Behav.velocity.Vy=v(:,2);
Behav.velocity.speed = sqrt(v(:,1).^2+v(:,2).^2);


save([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat'],'-struct','Behav');
disp([thisSID ' ParsedPosition editing is finished!'])
end
