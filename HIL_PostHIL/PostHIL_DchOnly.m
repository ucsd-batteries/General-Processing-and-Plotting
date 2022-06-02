clear all; close all

%% Output Accumulated DCH Only 
% - useful for early terminated tests

test_path = "/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/HIL/PostHIL/Aging15/Aging15_part2";
%   - channelIDs: order of brick IDs (format = brick1, brick2, brick3, brick4)
channelIDs = [1 2 9 10 3 4 11 12 7 6 13 14 5 8 15 16];


%% --------- CHARACTERIZATION PROCESSING -------------
% Process aging data - calculated accumulated discharge amp hours
%
% Output: table of values for DCH2 for each cell
%   - DCH2 = Ah discharged by cell for entire test (characterization + drive cycle)
%
% Assumptions: test file is from Arbin
% 

%% Import data and integrate current
channelStruct = dir(fullfile(test_path, '*.csv'));
channels = {channelStruct.name};
% sort channel files numerically by lowest to highest 
channels = sort_channels(channels);

% for each channel, get Start & End V, DCH1 & DCH2 
DCH2 = zeros(1,length(channels));
days2 = zeros(1,length(channels));
brick_days2 = zeros(1,4);
brick_DCH2 = zeros(1,4);
for c=1:length(channels)
    file_path = fullfile(test_path, channels(c));
    data = readmatrix(file_path{1});
    voltage = data(:,8);
    dch = data(:,11);
    steps = data(:,6);

    timestamp = readmatrix(file_path{1},'range','B:B','output','datetime');
    timestamp = timestamp(2:end);
    timeNum = datenum(timestamp);
    days = timeNum-timeNum(1);

    % calculate cummulative ahdch
    ahdch = zeros(length(dch),1);
    for i=2:length(dch)
       if dch(i) == 0
           ahdch(i) = ahdch(i-1) + dch(i);
       else
           ahdch(i) = ahdch(i-1) + (dch(i)-dch(i-1)); 
       end
    end

    DCH2(c) = ahdch(end); 
    days2(c) = days(end);
    for b=1:4
       ids = channelIDs(b*4-3:b*4);
       brick_days2(b) = max(days2(ids));
       brick_DCH2(b) = max(DCH2(ids));
    end
end


%% Summarize outputs

output = [brick_days2, brick_DCH2, DCH2];
brick_days1=zeros(1,4); brick_DCH1 = zeros(1,4); 
DCH1 = zeros(1,16); Cap = zeros(1,16);
output = [0, brick_days1, brick_days2, brick_DCH1, brick_DCH2, DCH1, DCH2, Cap];

function output_channels = sort_channels(channels)

    % sort channel files numerically by lowest to highest 
    unsorted_order = zeros(1,length(channels));
    for c=1:length(channels)
       str = channels{c};
       expr = 'Channel_\d+_';
       matchStr = regexp(str,expr,'match');
       unsorted_order(c) = str2double(matchStr{1}(9:end-1));
    end
    [~,idx] = sort(unsorted_order);
    for c=1:length(channels)
        channels_sorted{c} = channels{idx(c)};
    end
    output_channels = channels_sorted;
end





