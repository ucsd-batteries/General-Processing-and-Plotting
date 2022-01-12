clear all; close all

%% HIL Data Processing
% Split into two parts: 
%   - CHARACTERIZATION: determines capacity of each cell
%   - HIL Aging: determines Ah discharged during aging for each cell
% 
%  ----------------- USER INPUTS---------------------
%   - charac_path: path to folder for Characterization test data (leave blank if no Charac Data to process)
test_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/HIL/PostHIL/Char4';
% charac_path = '';
%   - channelIDs: order of brick IDs (format = brick1, brick2, brick3, brick4)
channelIDs = [1 2 9 10 3 4 11 12 7 6 13 14 5 8 15 16];

summary_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL/PostHIL_test_summary.csv';
% ------------------------------------------------------
%

%% --------- CHARACTERIZATION PROCESSING -------------
% Process characterization data
%
% Output: table of values for capacity, DCH1, and DCH2 for each cell
%   - capacity = estimated maximum capacity from characterization test
%   - DCH1 = Ah discharged by cell when characterization ends
%   - DCH2 = Ah discharged by cell for entire test (characterization + drive cycle)
%
% Assumptions: Internal Resistance test completed two characterization
% cycles
% 

%% Step 1: Import data
channelStruct = dir(fullfile(test_path, '*.csv'));
channels = {channelStruct.name};
% sort channel files numerically by lowest to highest 
channels = sort_channels(channels);

% for each channel, get Start & End V, DCH1 & DCH2 
StartVs = zeros(1,length(channels));    EndVs = zeros(1,length(channels));
DCH1 = zeros(1,length(channels));  DCH2 = zeros(1,length(channels));
CapDch = zeros(1,length(channels)); 
days1 = zeros(1,length(channels));  days2 = zeros(1,length(channels));
brick_days1 = zeros(1,4); brick_days2 = zeros(1,4);
brick_DCH1 = zeros(1,4); brick_DCH2 = zeros(1,4);
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

    ah0 = 0;
    step1=14; step2=19; 
    flag1 = true; flag2 = true;
    for s=1:length(steps)-1
       if steps(s)==step1 && steps(s+1)>step1 && flag1
           StartVs(c) = voltage(s);
           ah0 = ahdch(s);
           flag1=false;
       elseif steps(s)==step2 && steps(s+1)>step2 && flag2
           EndVs(c) = voltage(s);
           DCH1(c) = ahdch(s);
           CapDch(c) = DCH1(c)-ah0;
           days1(c) = days(s);
           flag2=false;
       end
    end
    DCH2(c) = ahdch(end)-DCH1(c); 
    days2(c) = days(end)-days1(c);
    for b=1:4
       ids = channelIDs(b*4-3:b*4);
       brick_days1(b) = max(days1(ids));
       brick_days2(b) = max(days2(ids));
       brick_DCH1(b) = max(DCH1(ids));
       brick_DCH2(b) = max(DCH2(ids));
    end
end

%% Step 2: Convert Voltages to SOC
SOCcurve = xlsread('SOC_curve.xls');
x=SOCcurve(:,1);
y=SOCcurve(:,2);
SOCstart = zeros(size(StartVs));
for i = 1:length(StartVs)
    V = StartVs(i); %Put Voltage to get SOC
    A = find(y == V);
    if A ~= []
        SOCstart(i) = x(A);
    elseif V >= y(end)
        SOCstart(i) = 100;
    elseif V <= y(1)
        SOCstart(i) = 0;
    else
        A1 = find(y<V);
        A2 = find(y>V);
        SOCstart(i) = ((x(A2(1))-x(A1(end)))/(y(A2(1))-y(A1(end))))*(V-y(A1(end)))+x(A1(end));
    end
end

SOCend = zeros(size(StartVs));
for i = 1:length(EndVs)
    V = EndVs(i); 
    A = find(y == V);
    if A ~= []
        SOCend(i) = x(A);
    elseif V >= y(end)
        SOCend(i) = 100;
    elseif V <= y(1)
        SOCend(i) = 0;
    else
        A1 = find(y<V);
        A2 = find(y>V);
        SOCend(i) = ((x(A2(1))-x(A1(end)))/(y(A2(1))-y(A1(end))))*(V-y(A1(end)))+x(A1(end));
    end
end

%% Step 3: Calculate Capacity for each channel
Cap = zeros(size(SOCstart));
for i =1:length(SOCstart)
   Cap(i) = (100/(SOCstart(i)-SOCend(i)))*CapDch(i); %Capacity
end

%% Step 4: Calculate brick cycles and days
% brick cycles = max(cycles per channel)
Charac_ah = zeros(1,4);
Charac_days = zeros(1,4);
for i=1:4
    ids = channelIDs(4*i-3:4*i);
    Charac_ah(i) = max(DCH2(ids)+DCH1(ids));
    Charac_days(i) = max(days(ids));
end

%% Combine all data into one table

% extract test name
splits = split(test_path,'/');
test_name = splits(end);

%
output = [0, brick_days1, brick_days2, brick_DCH1, brick_DCH2, DCH1, DCH2, Cap];
output = num2cell(output);
output{1} = test_name{1};

input = readcell(summary_path);
output = [input; output];

writecell(output,summary_path,'Delimiter',',')
% writecell(output,fullfile(outpath, 'test.csv'),'Delimiter',',')

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