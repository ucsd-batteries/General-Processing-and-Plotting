%% Nissan HIL Report Plots
clear all; close all;
% Import data and process
data = readmatrix("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL/PostHIL_test_summary.csv");
channelIDs = [1 2 9 10 3 4 11 12 7 6 13 14 5 8 15 16];
ID1=channelIDs(1:4);     ID2=channelIDs(5:8);   
ID3=channelIDs(9:12);    ID4=channelIDs(13:16);
RatedCap = 56.3;
IDs = [ID1; ID2; ID3; ID4];

days1 = data(:,2:5); days2 = data(:,6:9);
brick_DCH1 = data(:,10:13); brick_DCH2 = data(:,14:17);
rows = size(data,1);
brick_days = zeros(rows,4);   brick_days(1,1:4) = days1(1,1:4);
brick_cycles = zeros(rows,4);   brick_cycles(1,1:4) = brick_DCH1(1,1:4)/RatedCap;
% calculate brick-level days and cycles
for i=2:rows
   brick_days(i,1:4) = brick_days(i-1,1:4) + days2(i-1,1:4) + days1(i,1:4);
   brick_cycles(i,1:4) = brick_cycles(i-1,1:4) + (brick_DCH2(i-1,1:4) + brick_DCH1(i,1:4))/RatedCap;
end

cell_DCH1 = data(:,18:33);  cell_DCH2 = data(:,34:49);
Cap = data(:,50:end);      
% calculate cell cycle count
dchAh = zeros(size(Cap));
for c=1:16
    dchAh(1,c)=cell_DCH1(1,c);
    for r=2:rows
        dchAh(r,c)=(dchAh(r-1,c) + cell_DCH2(r-1,c) + cell_DCH1(r,c));
    end
end
cycles=dchAh/RatedCap;

% calculate SOH
SOH = Cap/RatedCap*100;
soh = zeros(rows, 4);
soh_std = zeros(rows, 4);
for b=1:4
    ids = IDs(b,:);
    soh(:,b) = mean(SOH(:,ids),2);
    soh_std(:,b) = std(SOH(:,ids),0,2);
end

%% PLOTS FOR INTERNAL REPORT

save_plot = true;  % Change to True to save plots
outpath = '/Users/quiana/Documents/UCSD/CER/Plots/PostHIL/';  % path that plot should be saved to

% Plot SOH vs cycle count
min_soh = min(min(SOH)); max_soh = max(max(SOH));
range = max_soh-min_soh; low = min_soh - range/10; high = max_soh + range/10;

fig = figure(1); clf; hold on
set(gcf, 'Position',  [100, 100, 1100, 650]);
for j=1:4
    ids = IDs(j,:);
    soh = SOH(:,ids);       % get soh for this brick
    cycs = cycles(:,ids);   % get cycles for this brick
    for i=1:size(soh,2)
        subplot(4,1,j); 
        title(['Brick', ' ', num2str(j), ': Cycle ', num2str(floor(brick_cycles(end,j)))]); 
        hold on; box on;
        plot(cycs(:,i), soh(:,i),'.-','markersize', 18, ...
            'displayname', ['Cell', ' ', num2str(ids(i))],'linewidth',2);
        xlim([0 max(max(cycles))])
        legend('location', 'eastoutside', 'fontsize', 15);
        set(gca, 'fontsize', 15); ylim([low high]);
        if j==4; xlabel("Cycle Count [-]"); end
    end
end
chart=axes(fig,'visible','off');  
chart.YLabel.Visible='on'; chart.Title.Visible='on';
y_pos = ylabel(chart,'State of Health [%]');
y_pos.Position(1) = -.05;
set(gca, 'fontsize', 15);

% ----- Save plot ------
if save_plot
    plot_name = strcat(outpath, strcat('PostHIL Cell Level SOH.png'));
    exportgraphics(gcf, plot_name, 'Resolution', 1000);
end


% SOH Deviation vs brick cycle count and days
figure(2); clf; 
set(gcf, 'Position',  [100, 100, 900, 600]);
for i=1:4
    subplot(1,2,1); hold on; box on; grid on; 
    plot(brick_cycles(:,i), soh_std(:,i), '.-','markersize', 18, 'linewidth', 2, ...
        'displayname', ['Brick ', num2str(i)]);
    legend('location', 'northeast'); ylabel('SOH Deviation [%]'); 
    xlabel('Cycle Count [-]')
    set(gca, 'fontsize', 15);
    subplot(1,2,2); hold on; box on; grid on; 
    plot(brick_days(:,i), soh_std(:,i), '.-','markersize', 18, 'linewidth', 2, ...
        'displayname', ['Brick ', num2str(i)]);
    legend('location', 'northeast'); xlabel('Time [days]');
    set(gca,'fontsize', 15);
end
sgtitle('Nissan Cycling Aging from HIL Testing','fontweight','bold', 'fontsize', 15) 

% ----- Save plot ------
if save_plot
    plot_name = strcat(outpath, strcat('PostHL SOH Std.png'));
    exportgraphics(gcf, plot_name, 'Resolution', 1000);
end

%% Benchmark, HIL, & PostHIL Analysis
clc

save_plot = true;  % Change to True to save plots
outpath = '/Users/quiana/Documents/UCSD/CER/Plots/PostHIL/';  % path that plot should be saved to

HIL_data = readmatrix("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL/HIL_test_summary.csv");
benchmark_data = readmatrix("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/NP6/NP6_test_summary.csv");

% NP6 processing: get cycles, std
DCH1 = benchmark_data(:,2); DCH2 = benchmark_data(:,3);
benchmark_cycles = zeros(1, size(benchmark_data,1));
benchmark_cycles(1) = DCH1(1)/RatedCap;
for i=2:size(benchmark_data,1)
   benchmark_cycles(i) = benchmark_cycles(i-1) + (DCH2(i-1) + DCH1(i))/RatedCap;
end
benchmark_cellCaps = benchmark_data(:,4:end);
cell_i = [1:16]; cell_i = cell_i(randperm(length(cell_i)));
cell_i = [16, 1, 6, 8, 5, 4, 12, 2, 15, 13, 3, 10, 7, 11, 14, 9];
benchmark_cap = zeros(size(benchmark_data,1),4);
benchmark_soh = zeros(size(benchmark_data,1),4);
benchmark_std = zeros(size(benchmark_data,1),4);
benchmark_std_norm = zeros(size(benchmark_data,1),4);
benchmark_std_pD = zeros(size(benchmark_data,1),4);
for b=1:4
    ids = cell_i(b*4-3:b*4);
    benchmark_cap(:,b) = mean(benchmark_cellCaps(:,ids),2);
    benchmark_soh(:,b) = benchmark_cap(:,b)/RatedCap*100;
    benchmark_std(:,b) = std(benchmark_cellCaps(:,ids),0,2);
    benchmark_std_norm(:,b) = std(benchmark_cellCaps(:,ids),0,2)/max(std(benchmark_cellCaps(:,ids),0,2));
    benchmark_std_pD(:,b) = (benchmark_std(:,b)-benchmark_std(1,b))/benchmark_std(1,b);
end


% HIL processing: get cycles, std
charac_cycles = HIL_data(:,7:10); aging_cycles = HIL_data(:,11:14);
HIL_cellCaps = HIL_data(:,15:30);     DCH1 = HIL_data(:,31:46);   
DCH2 = HIL_data(:,47:62);   HIL = HIL_data(:,63:end);
% calculate cumulative brick-levle cycle count and day duration
HIL_cycles = zeros(size(charac_cycles));
HIL_cycles(1,:) = charac_cycles(1,:);
for i=2:size(HIL_cellCaps,1)
    HIL_cycles(i,:) = HIL_cycles(i-1,:) + charac_cycles(i,:) + aging_cycles(i,:);
end
HIL_cap = zeros(size(HIL_data,1),4);
HIL_soh = zeros(size(HIL_data,1),4);
HIL_std = zeros(size(HIL_data,1),4);
HIL_std_norm = zeros(size(HIL_data,1),4);
HIL_std_pD = zeros(size(HIL_data,1),4);
for b=1:4
    ids = channelIDs(b*4-3:b*4)
    HIL_cap(:,b) = mean(HIL_cellCaps(:,ids),2)
    HIL_soh(:,b) = HIL_cap(:,b)/RatedCap*100;
    HIL_std(:,b) = std(HIL_cellCaps(:,ids),0,2);
    HIL_std_norm(:,b) = HIL_std(:,b)/max(HIL_std(:,b));
    HIL_std_pD(:,b) = (HIL_std(:,b)-HIL_std(1,b))/HIL_std(1,b);
end


% postHIL: add HIL cycles to postHIL
postHIL_soh = [HIL_soh(end,:); soh];
postHIL_std = [HIL_std(end,:); soh_std];
postHIL_cycles = [HIL_cycles(end,:); brick_cycles + HIL_cycles(end,:)];
postHIL_std_norm = zeros(rows+1,4);
postHIL_std_pD = zeros(rows+1,4);
for b=1:4
    postHIL_std_norm(:,b) = postHIL_std(:,b)/max(HIL_std(:,b));
    postHIL_std_pD(:,b) = (postHIL_std(:,b)-HIL_std(1,b))/HIL_std(1,b);
end

% plot HIL as one color
HIL_cycles = [HIL_cycles; postHIL_cycles(2:end,:)];
HIL_soh = [HIL_soh; postHIL_soh(2:end,:)];
HIL_std = [HIL_std; postHIL_std(2:end,:)];
HIL_std_norm = [HIL_std_norm; postHIL_std_norm(2:end,:)];
HIL_std_pD = [HIL_std_pD; postHIL_std_pD(2:end,:)];
    
% normalize std
benchmark_std = benchmark_std_pD;
HIL_std = HIL_std_pD;


colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]];
cents = [.495, .485, .548, .553];

% % for weekly update
% for b = 1:4
%     figure(b+2); clf; hold on; box on;
%     P = polyfit(benchmark_cycles,benchmark_soh(:,b),5); x = linspace(min(benchmark_cycles),max(benchmark_cycles),1000);
%     b_soh = polyval(P, x);
% %     b1 = plot(benchmark_cycles, smooth(benchmark_soh(:,b)),'displayname', 'Benchmark SOH','linewidth', 2, 'color', colors(1,:));
%     b1 = plot(x, b_soh,'displayname', 'Benchmark SOH','linewidth', 2, 'color', colors(1,:));
%     P = polyfit(HIL_cycles(:,b),HIL_soh(:,b),4); 
%     x1 = linspace(min(HIL_cycles(:,b)),postHIL_cycles(1,b),1000); x2 = linspace(postHIL_cycles(1,b), max(HIL_cycles(:,b)), 1000);
%     H_soh = polyval(P, x1); pH_soh = polyval(P, x2);
% %     h1 = plot(HIL_cycles(:,b), smooth(HIL_soh(:,b)), 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
% %     h1 = plot(x1, H_soh, 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
% %     ph1 = plot(x2, pH_soh, 'displayname', 'postHIL SOH','linewidth', 2, 'color', colors(3,:));
%     h1 = plot([x1, x2], [H_soh, pH_soh], 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
%     ylabel('State of Health [%]')
%     [low1, high1] = low_high([b_soh; H_soh], [b_soh; H_soh]);
%     ylim([low1 high1])
%     xline(postHIL_cycles(1,b), ':k','linewidth',2)
%     text(postHIL_cycles(1,b),low1+(high1-low1)/15," No Balancing",'FontSize',12,'HorizontalAlignment','left')
%     text(postHIL_cycles(1,b),low1+(high1-low1)/15,"Balancing ",'FontSize',12,'HorizontalAlignment','right')
%     annotation('textarrow',[cents(b),cents(b)+.03],[.41,.41],'linewidth',2)
%     annotation('textarrow',[cents(b)-.01,cents(b)-.04],[.41,.41],'linewidth',2)
%     yyaxis right 
%     ax = gca; ax.YAxis(2).Color = 'black';
%     P = polyfit(benchmark_cycles,benchmark_std(:,b),5);
%     b_std = polyval(P, benchmark_cycles);
%     b2 = plot(benchmark_cycles, b_std,'--', 'displayname', 'Benchmark Std','linewidth', 2, 'color', colors(1,:));
%     P = polyfit(HIL_cycles(:,b), HIL_std(:,b),3);
%     h_std = polyval(P, x1); ph_std = polyval(P,x2);
% %     h2 = plot(x1, h_std, '--', 'displayname', 'HIL Std','linewidth', 2, 'color', colors(2,:));
% %     ph2 = plot(x2, ph_std, '--', 'displayname', 'postHIL Std','linewidth', 2, 'color', colors(3,:));
%     h2 = plot([x1, x2], [h_std, ph_std], '--', 'displayname', 'HIL Std','linewidth', 2, 'color', colors(2,:));
%     ylabel('Imbalance Percent Change [%]')
%     [low, high] = low_high(b_std, h_std);
%     ylim([low, high])
%     xlim([0, 650])
%     title(['Benchmark vs HIL Brick ' num2str(b)])
% %     l = legend([b1, b2, h1, h2, ], 'Benchmark SOH', 'Benchmark Std', 'HIL SOH', 'HIL Std', 'PostHIL SOH', 'PostHIL Std',...
% %         'location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 2)
%     l = legend([b1, b2, h1, h2, ], 'Benchmark SOH', 'Benchmark Std', 'HIL SOH', 'HIL Std',...
%     'location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 2)
% %     l.NumColumns = 1; l.Location = 'eastoutside';
%     xlabel('Cycle Count [-]')
%     set(gca, 'linewidth', 2, 'fontsize', 12)
%     if save_plot
%         plot_name = strcat(outpath, strcat("Brick", num2str(b), " smoothed.png"));
%         exportgraphics(gcf, plot_name, 'Resolution', 1000);
%     end
% end

% raw plotting
for b = 1:4
    figure(b+2); clf; hold on; box on;
    b1 = plot(benchmark_cycles, benchmark_soh(:,b),'displayname', 'Benchmark SOH','linewidth', 2, 'color', colors(1,:));
    h1 = plot(HIL_cycles(:,b), HIL_soh(:,b), 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
    ylabel('State of Health [%]')
    [low1, high1] = low_high([benchmark_soh(:,b); HIL_soh(:,b)], [benchmark_soh(:,b); HIL_soh(:,b)]);
    ylim([low1 high1])
    xline(postHIL_cycles(1,b), ':k','linewidth',2)
    text(postHIL_cycles(1,b),low1+(high1-low1)/15," No Balancing",'FontSize',12,'HorizontalAlignment','left')
    text(postHIL_cycles(1,b),low1+(high1-low1)/15,"Balancing ",'FontSize',12,'HorizontalAlignment','right')
    annotation('textarrow',[cents(b),cents(b)+.04],[.41,.41],'linewidth',2)
    annotation('textarrow',[cents(b)-.01,cents(b)-.05],[.41,.41],'linewidth',2)
    yyaxis right 
    ax = gca; ax.YAxis(2).Color = 'black';
    b2 = plot(benchmark_cycles, benchmark_std(:,b),'--', 'displayname', 'Benchmark Std','linewidth', 2, 'color', colors(1,:));
    h2 = plot(HIL_cycles(:,b), HIL_std(:,b), '--', 'displayname', 'HIL Std','linewidth', 2, 'color', colors(2,:));
    ylabel('Imbalance Percent Change [%]')
    [low, high] = low_high(benchmark_std(:,b), HIL_std(:,b));
    ylim([low, high])
    xlim([0, 650])
    l = legend([b1, b2, h1, h2], 'Benchmark SOH', 'Benchmark Std', 'HIL SOH', 'HIL Std', ...
        'location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 2);
    title(['Benchmark vs HIL Brick ' num2str(b)])
    xlabel('Cycle Count [-]')
    set(gca, 'linewidth', 2, 'fontsize', 12)
    plot_name = strcat("Brick", num2str(b), " raw.png");
    if save_plot
        plot_name = strcat(outpath, strcat("Brick", num2str(b), " raw.png"));
        exportgraphics(gcf, plot_name, 'Resolution', 1000);
    end
end


% % for power point
% for b = 1:1
%     figure(b+2); clf; hold on; box on;
%     P = polyfit(benchmark_cycles,benchmark_soh(:,b),5); x = linspace(min(benchmark_cycles),max(benchmark_cycles),1000);
%     b_soh = polyval(P, x);
% %     b1 = plot(benchmark_cycles, smooth(benchmark_soh(:,b)),'displayname', 'Benchmark SOH','linewidth', 2, 'color', colors(1,:));
%     b1 = plot(x, b_soh,'displayname', 'Benchmark SOH','linewidth', 2, 'color', colors(1,:));
%     P = polyfit(HIL_cycles(:,b),HIL_soh(:,b),4); 
%     x1 = linspace(min(HIL_cycles(:,b)),postHIL_cycles(1,b),1000); x2 = linspace(postHIL_cycles(1,b), max(HIL_cycles(:,b)), 1000);
%     H_soh = polyval(P, x1); pH_soh = polyval(P, x2);
% %     h1 = plot(HIL_cycles(:,b), smooth(HIL_soh(:,b)), 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
% %     h1 = plot(x1, H_soh, 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
%     h1 = plot([x1, x2], [H_soh, pH_soh], 'displayname', 'HIL SOH','linewidth', 2, 'color', colors(2,:));
% %     ph1 = plot(x2, pH_soh, 'displayname', 'postHIL SOH','linewidth', 2, 'color', colors(3,:));
%     [low1, high1] = low_high([b_soh; H_soh], [b_soh; H_soh]);
%     ylim([low1 high1])
%     xline(postHIL_cycles(1,b), ':k','linewidth',2)
%     text(postHIL_cycles(1,b),low1+(high1-low1)/15," No Balancing",'FontSize',12,'HorizontalAlignment','left')
%     text(postHIL_cycles(1,b),low1+(high1-low1)/15,"Balancing ",'FontSize',12,'HorizontalAlignment','right')
%     annotation('textarrow',[.50,.545],[.41,.41],'linewidth',1.5)
%     annotation('textarrow',[.48,.435],[.41,.41],'linewidth',1.5)
%     yyaxis right 
%     ax = gca; ax.YAxis(2).Color = 'black';
%     P = polyfit(benchmark_cycles,benchmark_std(:,b),5);
%     b_std = polyval(P, benchmark_cycles);
%     b2 = plot(benchmark_cycles, b_std,'--', 'displayname', 'Benchmark Std','linewidth', 2, 'color', colors(1,:));
%     P = polyfit(HIL_cycles(:,b), HIL_std(:,b),3);
%     h_std = polyval(P, x1); ph_std = polyval(P,x2);
% %     h2 = plot(x1, h_std, '--', 'displayname', 'HIL Std','linewidth', 2, 'color', colors(2,:));
% %     ph2 = plot(x2, ph_std, '--', 'displayname', 'postHIL Std','linewidth', 2, 'color', colors(3,:));
%     h2 = plot([x1, x2], [h_std, ph_std], '--', 'displayname', 'HIL Std','linewidth', 2, 'color', colors(2,:));
%     [low, high] = low_high(b_std, h_std);
%     ylim([low, high])
%     xlim([0, 650])
% %     title(['Benchmark vs HIL Brick ' num2str(b)])
%     l = legend([b1, b2, h1, h2], 'Control Group SOH', 'Control Group Imbalance', 'Study Group SOH', 'Study Group Imbalance',...
%         'location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 2)
% %     l.NumColumns = 1; l.Location = 'eastoutside';
%     set(gca, 'linewidth', 2, 'fontsize', 14)
%     ylabel('Imbalance Percent Change [%]', 'fontsize', 14)
%     yyaxis left 
%     ylabel('State of Health [%]', 'fontsize', 14)
%     xlabel('Cycle Count [-]', 'fontsize', 14)
%     plot_name = strcat("Brick", num2str(b), " ppt.png");
% end


function [low, high] = low_high(max_vec, min_vec)
    min_std = min(min(min_vec)); max_std = max(max(max_vec));
    range = max_std-min_std; low = min_std - range/4; high = max_std + range/10;
end
