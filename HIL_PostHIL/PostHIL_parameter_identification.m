close all, clear all
%% load data file
data_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/HIL/PostHIL';
summary_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL/postHIL_IR_summary.xlsx';
soc_path = '/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/HIL/SOC_curve.xls';
% testStruct = dir(fullfile(data_path, 'Internal_Resistance_Test*'));
% tests = {testStruct.name};
% tests = ["Char9", "Aging9", "Char10", "Aging10", "Char11", "Aging11", ...
%     "Char12", "Aging12/Aging12_part2", "Char13/Char13_take2", "Aging13/Aging13_part1"];
% tests = ["Char1", "Aging1", "Char2", "Aging2", "Char3", "Aging3", ...
%     "Char4", "Aging4", "Char5", "Aging5", "Char6", "Aging6", "Char7", "Aging7", "Char8", "Aging8"];

tests=["Char16"];

%% initialize arrays for battery parameters
R0_all = zeros(length(tests),16); Rp_all = zeros(length(tests),16);
Cp_all = zeros(length(tests),16); 
RMSE_all = zeros(length(tests),16); MAE_all = zeros(length(tests),16);

for test=1:length(tests)
%     channel_path = '/Users/quiana/Documents/UCSD/CER/ParameterEstimation/Longs/HIL/Internal_Resistance_Test_2021_05_11_143509';
    channel_path = fullfile(data_path, tests(test));
    channelStruct = dir(fullfile(channel_path, '*.CSV'));
    channels = {channelStruct.name};

    % Get battery capacity via integration
    Caps = getCellCap(channels, channel_path, soc_path);


    for l = 1:length(channels)
        %%
        % get channel data
%         l=6;
        channel = channels{l};
        cap = Caps(l);
        data = readmatrix(fullfile(channel_path,channel));

        % get drive cycle dataset
        steps = data(:,6);
        drive_cycle = find(steps==23, 1, 'first');
        t = drive_cycle:length(steps);    % get the data length of drive cycle only
        t = (drive_cycle+round(.1*length(t))):(length(steps)-round(.1*length(t)));    % remove first and last 10% of data

        % get current, voltage, and SOC
        I = -data(t,7);
        V = data(t,8);
        seconds = data(t,3);

        %%
        % integrate current
        ah = zeros(1,length(I));
        for i = 2:length(I)
            ah(i) = ah(i-1) - trapz(seconds(i-1:i),I(i-1:i))/3600;
        end

        %%
        SOC_f = ah/cap;
%         SOC_f = SOC_f + abs(min(SOC_f));
        SOC_f = SOC_f + .95;

        %% 
        % chg = data(:,10);
        % dch = data(:,11);
        % 
        % ahdch = zeros(length(dch),1);
        % for i=2:length(dch)
        %    if dch(i) == 0
        %        ahdch(i) = ahdch(i-1) + dch(i);
        %    else
        %        ahdch(i) = ahdch(i-1) + (dch(i)-dch(i-1)); 
        %    end
        % end
        % 
        % figure(1); clf; plot(chg); title("Ah Charge")
        % figure(2); clf; plot(dch); title("Ah Discharge")
        % figure(3); clf; plot(ahdch); title("Cum Ah Discharge")

        %%
        % t = drive_cycle:length(steps);    % get the data length of drive cycle only
        % t = 1:length(steps);
        % if length(t)<length(V)
        %     V = V(t);
        %     I = I(t);
        %     SOC_f = SOC_f(t);
        % end

        figure(4); plot(I); title("HIL Current")
        figure(5); plot(V); title("HIL Voltage");
        % figure(6); plot(getSOC(V)/100); title("HIL SOC"); 
        figure(7); plot(SOC_f,'linewidth',2); title("HIL SOC Calculated"); ylabel("SOC [x100%]")
        set(gca, 'fontsize', 15)

%         [R0, Rp, Cp, RMSE, MAE] = ParamID(SOC_f, I, V, t, soc_path);
        [R0, Rp, Cp, RMSE, MAE] = ParamID_new(SOC_f, I, V, t,soc_path);

        R0_all(test,l) = R0; Rp_all(test,l) = Rp; Cp_all(test,l) = Cp;
        RMSE_all(test,l) = RMSE; MAE_all(test,l) = MAE;

    end
    disp(['test ', num2str(test), ': ', tests(test), ' completed'])
    
    %% save to ongoing summary 
    % extract test name
    splits = split(channel_path,'/');
    test_name = splits(end);
    
    %
    output = [0, R0_all(test,:), Rp_all(test,:), Cp_all(test,:), RMSE_all(test,:), MAE_all(test,:)];
    output = num2cell(output);
    output{1} = test_name{1};
    
    input = readcell(summary_path);
    output = [input; output];
    
    writecell(output,summary_path)
    % writecell(output,fullfile(outpath, 'test.csv'),'Delimiter',',')

end



% %% Plot R0
% pos_R0=R0_all*1000;
% low = min(min(pos_R0)); high = max(max(pos_R0));
% ohm_range = high - low;
% 
% figure(3); clf; hold on; box on; 
% for j=1:4
%    subplot(4,1,j); hold on; box on;  
%    for i=(j*4-3):(j*4)
%         plot(pos_R0(:,i), '.-', 'markersize', 20,  ...
%             'linewidth', 2,'displayname',['Channel: ', num2str(i)])
%    end
%    legend('location', 'eastoutside')
%    ylabel('Resistance [m\Omega]')
%    set(gca, 'fontsize', 15, 'linewidth', 2)
%    title('HIL Channel Internal Resistance Progression')
%    ylim([low-ohm_range*.1 high+ohm_range*.1])
%    
% end
% xlabel('Test Index [-]')
% 
% %%
% pos_R0=Rp_all*1000;
% low = min(min(pos_R0)); high = max(max(pos_R0));
% ohm_range = high - low;
% figure(4); clf; hold on; box on; 
% for i=1:16
%     if i<9
%         plot(pos_R0(:,i), '.-', 'markersize', 20,  ...
%             'linewidth', 2,'displayname',['Channel: ', num2str(i)])
%     else
%         plot(pos_R0(:,i), '*-', 'markersize', 8,  ...
%             'linewidth', 2,'displayname',['Channel: ', num2str(i)])
%     end
% end
% legend('location', 'eastoutside')
% ylabel('Resistance [m/Omega]')
% set(gca, 'fontsize', 15, 'linewidth', 2)
% title('HIL Channel Internal Resistance')
% ylim([low-ohm_range*.1 high+ohm_range*.1])
% xlabel('Test Index [-]')

%%
% [R0, Rp, Cp, RMSE, MAE] = ParamID(SOC_f, I, V, t)

%%
function [R_0, R_p, C_p, RMSE, MAE] = ParamID(SOC_f, I, V, t)
    %% build OCV-SOC model
    [num1,txt1]=xlsread(soc_path,1);% load open circuit voltage test data
    SOC=num1(1:end,1)/100;% SOC data
    OCV=num1(1:end,2);% OCV data
    p=polyfit(SOC,OCV,9);% fit OCV-SOC model with a 9-degree Power function polynomial. This degree is determined by a comprehensive comparison of the calculation and accuracy of polynomials with different degrees.   
    %% parameter identification initialization
    N=length(t);% the data length
    R0=zeros(1,N); Rp=zeros(1,N); Cp=zeros(1,N);%The Thevenin battery model is applied to depict the battery's dynamic characteristics.Thus, Rp, Cp, R0 represents the polarization resistance, polarization capacitance, the equivalent ohmic resistance, respectively.
    Theta=zeros(N,3);% Theta indicates the parameters to be identified. These are the parameters in the recursive least square (RLS) method instead of the battery parameters.  
    Up=zeros(1,N+1); % Up means the polarization voltage. 
    error=zeros(1,N);% voltage estimation error
    R0_0=0.00143; Rp_0=0.03; Cp_0=1000; % set the initial values of the battery parameters
    Lambda=0.97;% the forgetting factor, and its value usually from 0.95 to 1.The smaller the ¦Ë, the stronger the tracking ability, but the greater the fluctuation.
    dt=1;% sampling interval
    P=1*10^5*eye(3); % error covariance in the RLS
    for i=1:N
        %% Parameter identification with RLS
        if i==1
            R0(i)=R0_0; Rp(i)=Rp_0; Cp(i)=Cp_0;
            % calculate the Theta according to the initial battery parameters
            Theta(i,1)=(2*Rp(i)*Cp(i)-dt)/(2*Rp(i)*Cp(i)+dt);
            Theta(i,2)=(R0(i)*dt+Rp(i)*dt+2*Rp(i)*Cp(i)*R0(i))/(2*Rp(i)*Cp(i)+dt);
            Theta(i,3)=(R0(i)*dt+Rp(i)*dt-2*Rp(i)*Cp(i)*R0(i))/(2*Rp(i)*Cp(i)+dt);
            % observation vector=[OCV - measurement battery voltage]
            z(i)=polyval(p,SOC_f(i))-V(i); % obtain OCV from SOC using OCV-SOC model 
        else
            % observation vector
            z(i)=polyval(p,SOC_f(i))-V(i);
            % RLS
            Phi=[z(i-1) I(i) I(i-1)]';% input and output sequence
            K=P*Phi/(Phi'*P*Phi+Lambda);% calculate algorithm gain
            Theta(i,:)=Theta(i-1,:)+(K*(z(i)-Phi'*Theta(i-1,:)'))'; % calculate the Theta
            P=(P-K*Phi'*P)/Lambda; % update error covariance
            % calculate battery parameters
            R0(i)=(Theta(i,2)-Theta(i,3))/(1+Theta(i,1));
            Cp(i)=(1+Theta(i,1))^2*dt/(4*(Theta(i,1)*Theta(i,2)+Theta(i,3)));
            Rp(i)=2*(Theta(i,1)*Theta(i,2)+Theta(i,3))/((1-Theta(i,1))*(1+Theta(i,1)));
        end
        %% calculate the voltage estimation error
        A=exp(-dt/(Rp(i)*Cp(i))); % calculate factor A
        B=Rp(i)*(1-exp(-dt/(Rp(i)*Cp(i))));% calculate factor B
        Up(i+1)=A*Up(i)+B*I(i);% update Up
        error(i)=polyval(p,SOC_f(i))-Up(i+1)-R0(i)*I(i)-V(i);% calculate voltage estimation error according to the Thevenin battery model and OCV-SOCC model
    end
    %% calculate the voltage estimation error
    R_0 = R0(end);
    R_p = Rp(end);
    C_p = Cp(end);
    RMSE=sqrt(sum(error.^2)/length(error)); % mean absolute error of voltage estimation
    MAE=mean(abs(error)); % root mean square error of voltage estimation
    %% plot
    figure(1)% model parameters 
    clf;
    subplot(3,1,1)
    plot(t/3600,R0,'k','LineWidth',1)
%     ylim([-0.01 0.01]);
    xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    ylabel('R_{0}(\Omega)','Fontname', 'Times new roman','FontSize',8)
    set(gca, 'Fontname', 'Times new roman', 'Fontsize',8,'LineWidth',1.2); 
    hold on
    subplot(3,1,2)
    plot(t/3600,Rp,'k','LineWidth',1)
%     ylim([-50 50])
    xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    ylabel('R_{p}(\Omega)','Fontname', 'Times new roman','FontSize',8)
    set(gca, 'Fontname', 'Times new roman', 'Fontsize',8,'LineWidth',1.2); 
    hold on
    subplot(3,1,3)
    plot(t/3600,Cp,'k','LineWidth',1)
%     ylim([-5000 50000])
    xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    ylabel('C_{p}(F)','Fontname', 'Times new roman','FontSize',8)

    figure(2) % voltage estimation error
    clf;
    plot(t/3600,error,'k','LineWidth',1)
%     ylim([-0.15 0.15]);
    xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    ylabel('voltage estimation error','Fontname', 'Times new roman','FontSize',8)
end

function [R0, Rp, Cp, RMSE, MAE] = ParamID_new(SOC_f,I,V,t,soc_path)
    %% build OCV-SOC model
    [num1,txt1]=xlsread(soc_path,1);% load open circuit voltage test data
    SOC=num1(1:end,1)/100;% SOC data
    OCV=num1(1:end,2);% OCV data
    p=polyfit(SOC,OCV,9);% fit OCV-SOC model with a 9-degree Power function polynomial. This degree is determined by a comprehensive comparison of the calculation and accuracy of polynomials with different degrees.   

    %% parameter identification initialization
    N=length(t);% the data length
    %The Thevenin battery model is applied to depict the battery's dynamic characteristics.Thus, Rp, Cp, R0 represents the polarization resistance, polarization capacitance, the equivalent ohmic resistance, respectively.
    z=zeros(N,1);%observation vector
    Phi=zeros(N-1,3);% input and output sequences
    dt=1;% sampling interval
    Up=zeros(1,N+1); % Up means the polarization voltage. 
    error=zeros(1,N);% voltage estimation error
    Lambda=1;% the forgetting factor, and its value usually from 0.95 to 1.The smaller the ¦Ë, the stronger the tracking ability, but the greater the fluctuation.
    for i=1:N
        if i==1
        % observation vector=[OCV - measurement battery voltage]
        z(i)=(polyval(p,SOC_f(i))-V(i))*Lambda^(N-i); % obtain OCV from SOC using OCV-SOC model
        else
        % observation vector
        z(i)=(polyval(p,SOC_f(i))-V(i))*Lambda^(N-i);
        % input and output sequences
        Phi(i-1,:)=Lambda^(N-i)*[z(i-1) I(i) I(i-1)];
        end
    end
    %% Parameter identification with RLS
    % Theta indicates the parameters to be identified. These are the parameters in the recursive least square (RLS) method instead of the battery parameters.  
    Theta=inv(Phi'*Phi)*Phi'*z(2:N);
    % calculate battery parameters
    R0=(Theta(2)-Theta(3))/(1+Theta(1))
    Rp=2*(Theta(1)*Theta(2)+Theta(3))/((1-Theta(1))*(1+Theta(1)))
    Cp=(1+Theta(1))^2*dt/(4*(Theta(1)*Theta(2)+Theta(3)))
%     R0=-Theta(3)/Theta(1)
%     Rp=(Theta(1)*Theta(2)+Theta(3))/((-1-Theta(1))*Theta(1))
%     Cp=(Theta(1))^2*dt/(Theta(1)*Theta(2)+Theta(3))

    %% calculate the voltage estimation error
    A=exp(-dt/(Rp*Cp)); % calculate factor A
    B=Rp*(1-exp(-dt/(Rp*Cp)));% calculate factor B
    for j=1:N
        Up(j+1)=A*Up(j)+B*I(j);% update Up
        error(j)=polyval(p,SOC_f(j))-Up(j+1)-R0*I(j)-V(j);% calculate voltage estimation error according to the Thevenin battery model and OCV-SOCC model
    end
    % mean absolute error of voltage estimation
    RMSE=sqrt(sum(error.^2)/length(error))
    % root mean square error of voltage estimation
    MAE=mean(abs(error))

    %% plot
    figure(1) % voltage estimation error
    clf;
    plot(t/3600,error,'k','LineWidth',1)
    xlabel('Time(h)','FontSize',15)
    ylabel('voltage estimation error','FontSize',15)
    set(gca, 'fontsize', 15)
end

%% FUNCTIONS
% ------------------------------------------------------------------------------------------
% Function getSOC -- calculates SOC given voltage V from OCV/SOC curve
function SOC = getSOC(V, soc_path)
    SOCcurve = xlsread(soc_path);
    x=SOCcurve(:,1);
    y=SOCcurve(:,2);
    for i = 1:length(V)
        v = V(i); %Put Voltage to get SOC
        A = find(y == v);
        if A ~= []
            SOC(i) = x(A);
        elseif v >= y(end)
            SOC(i) = 100;
        elseif v <= y(1)
            SOC(i) = 0;
        else
            A1 = find(y<v);
            A2 = find(y>v);
            SOC(i) = ((x(A2(1))-x(A1(end)))/(y(A2(1))-y(A1(end))))*(v-y(A1(end)))+x(A1(end));
        end
    end
end

% Function getCellCap -- calculates cell capacities
function Cap = getCellCap(channels, charac_path, soc_path)
    % Inputs:
    % - channels = list of filename for each channel
    % - charac_path = full path to folder holding channel filenames
    
    % Outputs:
    % - Cap = list of capacities for each cell
    
    %% Step 1: Import data
    % for each channel, get Start & End V, DCH1 & DCH2 
    StartVs = zeros(1,length(channels));    EndVs = zeros(1,length(channels));
    CapDch = zeros(1,length(channels));    
    for c=1:length(channels)
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
    end

    %% Step 2: Convert Voltages to SOC
    SOCcurve = xlsread(soc_path);
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
    
end
