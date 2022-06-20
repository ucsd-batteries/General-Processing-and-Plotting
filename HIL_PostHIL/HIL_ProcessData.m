%% HIL Data Processing
% Split into two parts: 
%   - CHARACTERIZATION: determines capacity of each cell
%   - HIL Aging: determines Ah discharged during aging for each cell
% 
%  ----------------- USER INPUTS---------------------
%   - charac_path: path to folder for Internal Resistance test data (leave blank if no Charac Data to process)
charac_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/HIL/Internal_Resistance/Internal_Resistance_Test_CC_2021_03_08_172735_RIGHT';
% charac_path = '';
%   - channelIDs: order of brick IDs (format = brick1, brick2, brick3, brick4)
channelIDs = [1 2 9 10 3 4 11 12 7 6 13 14 5 8 15 16];
%   - aging_path: path to HIL Aging tests (leave blank if no Aging Data to process)
% aging_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/HIL/HIL_Aging/Aging14';
aging_path = '';
%   - outpath: path to HIL summary dataset
outpath = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL';
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

if length(charac_path)>0
    %% Step 1: Import data
    channelStruct = dir(fullfile(charac_path, 'Internal_Resistance*'));
    channels = {channelStruct.name};
    % sort channel files numerically by lowest to highest 
    unsorted_order = zeros(1,length(channels));
    for c=1:length(channels)
       str = channels{c};
       expr = '_\d+_';
       matchStr = regexp(str,expr,'match');
       unsorted_order(c) = str2double(matchStr{1}(2:end-1));
    end
    [~,idx] = sort(unsorted_order);
%     channels = channels{idx};
    for c=1:length(channels)
        channels_sorted{c} = channels{idx(c)};
    end
    channels = channels_sorted;

    % for each channel, get Start & End V, DCH1 & DCH2 
    StartVs = zeros(1,length(channels));    EndVs = zeros(1,length(channels));
    DCH1 = zeros(1,length(channels));  DCH2 = zeros(1,length(channels));
    CapDch = zeros(1,length(channels)); 
    days = zeros(1,length(channels));
    for c=1:length(channels)
        channels{c};
        file_path = fullfile(charac_path, channels{c});
        data = readmatrix(file_path);
        voltage = data(:,8);
        dch = data(:,11);
        steps = data(:,6);

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
        for s=1:length(steps)-1
           if steps(s)==step1 && steps(s+1)>step1
               StartVs(c) = voltage(s);
               ah0 = ahdch(s);
           elseif steps(s)==step2 && steps(s+1)>step2
               EndVs(c) = voltage(s);
               DCH1(c) = ahdch(s);
               CapDch(c) = DCH1(c)-ah0;
           end
        end
        DCH2(c) = ahdch(end)-DCH1(c); 
        days(c) = length(dch)/(3600*24);
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
    Charac_cycles = zeros(1,4);
    Charac_days = zeros(1,4);
    CapRating = 56.3;
    for i=1:4
        ids = channelIDs(4*i-1:4*i);
        Charac_cycles(i) = max(DCH2(ids))/CapRating;
        Charac_days(i) = max(days(ids));
    end
    
else
    Cap = zeros(1,16);
    DCH1 = zeros(1,16);
    DCH2 = zeros(1,16);
    Charac_cycles = zeros(1,4);
    Charac_days = zeros(1,4);
end

%% ---------- HIL AGING PROCESSING -------------
% Process characterization data
%
% Output: table of values for Ah discharge for each cell for each test
%
if length(aging_path)>0
    
    %% Step 1: Import Data
    folderStruct = dir(aging_path);
    folderStruct=folderStruct(4:end);
    folders = {folderStruct.name};

    files = [];  all_ah = [];
    hil_cycles = [];    hil_days = [];
    for f=1:length(folders)
       folder_path = fullfile(aging_path, folders{f});
       logStruct = dir(fullfile(folder_path, '*CellParam*'));
       logs = {logStruct.name};
       for l=1:length(logs)
          file = fullfile(folder_path, logs{l});
          data = readmatrix(file);
          timestamp = readmatrix(file,'range','A:A','output','datetime');
          timestamp = timestamp(3:end);
          timeNum = datenum(timestamp);
          
    %% Step 2: Integrate Current to get Ah Discharge 
          ah=zeros(1,16);
          for i=1:16
              current=data(:,(i-1)*6+3);
              c=find(current<=0);
              ah(i) = -trapz(timeNum(c),current(c))*24;
          end
          % sort ah by cell num (currently in brick order)
          ah = ah(channelIDs);
          all_ah = [all_ah; ah];
          files = [files; string(logs{l})];
       end
       
    %% Step 3: Get Brick Cycle Count and Day Count
       brickLogsStruct = dir(fullfile(folder_path, '*BrickParam*'));
       brickLogs = {brickLogsStruct.name};
       for l=1:length(brickLogs)
          file = fullfile(folder_path, brickLogs{l});
          data = readmatrix(file);
          b1 = data(end,2)-data(1,2);
          b2 = data(end,10)-data(1,10);
          b3 = data(end,18)-data(1,18);
          b4 = data(end,26)-data(1,26);
          
          hil_cycles = [hil_cycles; [b1 b2 b3 b4]];
          
          % get number of days elapsed
          times = readmatrix(file, "range","A:A","outputtype", "datetime");
          [~, ~, dys, hrs, mins, secs] = datevec(times(end) - times(3)); 
          hil_days = [hil_days; dys + hrs/24 + mins/(24*60) + secs/(24*3600)];
          
       end
          
    end
    HIL_aging = sum(all_ah);
    HIL_cycles = sum(hil_cycles);
    HIL_days = sum(hil_days);
    
else
    HIL_aging = zeros(1,16);
    HIL_cycles = zeros(1,4);
    HIL_days = 0;
end

%% Combine all data into one table

% extract test name
year = '2021';
path_split = split(charac_path, '/');
test = split(path_split(end),'_');
i = find(contains(test,year));
test_name = strcat(test(i), "_", test(i+1),"_", test(i+2),"_", test(i+3));

%%
output = [0, Charac_days, HIL_days, Charac_cycles, HIL_cycles, Cap, DCH1, DCH2, HIL_aging];
output = num2cell(output);
output{1} = test_name{1};
% output{1} = 'new_aging';

input = readcell(fullfile(outpath, 'HIL_test_summary.csv'));
output = [input; output];

writecell(output,fullfile(outpath, 'HIL_test_summary.csv'),'Delimiter',',')
% writecell(output,fullfile(outpath, 'test.csv'),'Delimiter',',')

% update spreadsheet of channel id SOH values - in channel id order
writematrix(Cap(channelIDs), 'Channel_SOH.csv');


