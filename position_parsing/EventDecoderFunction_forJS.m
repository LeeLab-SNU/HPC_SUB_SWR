
function EventDecoderFunction_forJS(thisRID,thisSID,ROOT)
session_root = [ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID];  

cd(session_root);
evtfile_listing = dir('Events.csv');

nTable_temp = readtable([session_root '\Events.csv']);
if ~exist('ExtractedEvents','dir')
    mkdir('ExtractedEvents');
end
cd('ExtractedEvents');
% 
% BehStart = find(strcmp("Starting Recording",nTable_temp.Var18));
% BehStop= find(strcmp("Stopping Recording",nTable_temp.Var18));
% 
% n2=[];
% for n=1:size(BehStart,1)
% n2(n) = BehStop(n)-1 - (BehStart(n)+1);
% end
% [m,n] = max(n2);
% nTable = nTable_temp(BehStart(n)+1:BehStop(n)-1,:);


UEtoNlxEvent = ([ROOT.Raw.Mother '\Rat' thisRID '\Rat' thisRID '-' thisSID '\Events.csv']);
filename = [UEtoNlxEvent];
delimiter = ',';
startRow = 2;
formatSpec = '%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Events = table(dataArray{1:end-1}, 'VariableNames', {'Timestamp','Event'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

Timestamp = Events{:, 1};
Event = Events{:, 2};

str0 = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'; %TTL signal from the ARDUINO
str1 = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0010).'; %TTL signal from the ARDUINO
str2 = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0040).'; 
str3 = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0004).';
str4 = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0008).';
str5 = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0050).';

TTL1 = find(strcmp(Event, str1));

TTL1Time = Timestamp(TTL1);

TTL1TimeT = (((TTL1Time-min(TTL1Time)))./10e6);

% Import UE event file
UnrealFile =([ROOT.Raw.Mother '\Rat' thisRID '\Rat' thisRID '-' thisSID '\Rat' thisRID '_' thisSID '.csv']);
delimiter = ',';
formatSpec = '%s%f%*s%f%[^\n\r]';
fileID = fopen(UnrealFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);

X = dataArray{:, 1};
Time = dataArray{:, 2};
Position = dataArray{:, 3};
RewardEvent = 'Reward';
% Position([1:end])=((Position(1:end))./20);

clearvars UnrealFile delimiter formatSpec fileID dataArray ans;

Position(Position<0) = 0;
Filtered = ~isnan(Position);
FilteredPositionIndex = find(Filtered==1);
FilteredTime = Time(FilteredPositionIndex);
FilteredPosition = Position(FilteredPositionIndex);

nTable = table;
nTable.time = TTL1Time./10e6;

if size(nTable.time,1)>size(FilteredPosition,1)
    nTable.pos(1:size(FilteredPosition,1)) = FilteredPosition;
    disp([thisRID '-' thisSID ' unreal vs. neuro mismatch!'])
elseif size(nTable.time,1)<size(FilteredPosition,1)
    nTable.pos = FilteredPosition(1:size(nTable.time,1));
    disp([thisRID '-' thisSID ' unreal vs. neuro mismatch!'])
else
    nTable.pos = FilteredPosition;
end

nTable.cxt(nTable.pos<=10000) = 1;
nTable.cxt(nTable.pos>=10000 & nTable.pos<=20000) = 2;
nTable.cxt(nTable.pos>20000 ) = 0;
posN = size(nTable,1);
trial=[];

t=0;
for i=2:posN
   if nTable.cxt(i)>0 && nTable.cxt(i-1)==0
       t=t+1;
   end
       nTable.trial(i)=t;
end
nTable.trial(nTable.cxt==0)=0;
total_trial_number = max(nTable.trial);
%initialization
reward = zeros(posN, 1);
a = zeros(posN, 1);
t = zeros(posN, 1);
x = zeros(posN, 1);
y = ones(posN, 1);
area=[];
trial = logical(zeros(posN, total_trial_number));

correctness = logical(zeros(posN, 2));
cont = logical(zeros(posN, 6));
side = logical(zeros(posN, 2));
ambiguity = logical(zeros(posN, 3));


% total_trial_number = trialN; % - voidN;
trial_void = logical(zeros(total_trial_number,1)); %it will be not used.
trial_correctness = ones(total_trial_number,1);
trial_context = zeros(total_trial_number,1);
trial_side = zeros(total_trial_number,1);
trial_ambiguity = ones(total_trial_number, 1);

t=nTable.time;
x = nTable.pos-fix(nTable.pos/10000)*10000;
correctness(:,1)=1;
for i=1:2, cont(nTable.cxt==i,i)=1; end
side = cont(:,1:2);
ambiguity(:,1) = 1;
for i=1:total_trial_number
    trial(nTable.trial==i,i)=1;
    trial_context(i) = nTable.cxt(min(find(nTable.trial==i)));
end
trial_side = trial_context;




%% save parsed position results

save([session_root '\ParsedPosition.mat'], 't', 'x', 'y', 'a', 'cont', ...
    'correctness', 'trial', 'side', 'ambiguity', 'area','trial_correctness', 'trial_side', 'trial_context', 'trial_ambiguity', 'total_trial_number'); %void trials are deleted.







