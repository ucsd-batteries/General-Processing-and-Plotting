close all, clear all
%% load data file
data_path = 'C:\Users\amirs\OneDrive - UC San Diego\college\research\Dr Tong ESS\Nissan cycle testing\raw data\Tesla';
testStruct = dir(fullfile(data_path, '*.csv'));
tests = {testStruct.name};

%% initialize arrays for battery parameters
R0_all = zeros(length(tests),18); Rp_all = zeros(length(tests),18);
Cp_all = zeros(length(tests),18); 
RMSE_all = zeros(length(tests),18); MAE_all = zeros(length(tests),18);


% test = 1;
%for test=7:length(tests)
for test =1:length(tests)
    %test = 1;
    test_file =fullfile(data_path, tests{test});
%     data = readmatrix(test_file);
        
    
    % Get battery capacity via integration
    [steps, Caps, data] = getCellCap(test_file);
   %%
    % get drive cycle dataset
%     % characterization (steps 24-146)
%     drive_cycle_start = steps(24);
%     drive_cycle_end = steps(146);
    % duty cycle (steps 29 to end)
    drive_cycle_start = steps(26);
    drive_cycle_end = size(data,1);
    
    t = drive_cycle_start:drive_cycle_end;    % get the data length of drive cycle only
    t = (drive_cycle_start+round(.1*length(t))):(drive_cycle_end-round(.1*length(t)));        % remove first and last 10% of data
    
    hours = data(t,1);
    I = data(t,2);
    cell_voltages = data(t,4:end);
    %%
    for l = 1:18
        % get current, voltage, and SOC
        V = cell_voltages(:,l);
        cap = Caps(l);

        
        % integrate current
        ah = zeros(1,length(I));
        for i = 2:length(I)
            ah(i) = ah(i-1) - trapz(hours(i-1:i),I(i-1:i));
        end

        
        SOC_f = ah/cap;
%         SOC_f = SOC_f + abs(min(SOC_f));
        SOC_f = SOC_f + abs(.99-max(SOC_f));
%         SOC_f = SOC_f + .95;


        %figure(4); plot(I); title("HIL Current")
        %figure(5); plot(V); title("HIL Voltage");
        % figure(6); plot(getSOC(V)/100); title("HIL SOC"); 
        %figure(7); plot(SOC_f,'linewidth',2); title("HIL SOC Calculated"); ylabel("SOC [x100%]")
        %set(gca, 'fontsize', 15)
        
%         [R0, Rp, Cp, RMSE, MAE] = ParamID(SOC_f, I, V, t);
        [R0, Rp, Cp, RMSE, MAE] = ParamID_new(SOC_f, I, V, t);

        R0_all(test,l) = R0; Rp_all(test,l) = Rp; Cp_all(test,l) = Cp;
        RMSE_all(test,l) = RMSE; MAE_all(test,l) = MAE;

    end
    disp(['test ', num2str(test), ': ', tests{test}, ' completed'])
end

%% Plot R0
pos_R0=R0_all*1000;
low = min(min(pos_R0)); high = max(max(pos_R0));
ohm_range = high - low;

%figure(3); clf; hold on; box on; 
for j=1:4
   %subplot(4,1,j); hold on; box on;  
   for i=(j*4-3):(j*4)
        %plot(pos_R0(:,i), '.-', 'markersize', 20,  ...
            %'linewidth', 2,'displayname',['Channel: ', num2str(i)])
   end
   %legend('location', 'eastoutside')
   %ylabel('Resistance [m/Omega]')
   %set(gca, 'fontsize', 15, 'linewidth', 2)
   %title('HIL Channel R0 Progression')
   %ylim([low-ohm_range*.1 high+ohm_range*.1])
   
end
%xlabel('Test Index [-]')

%%
low = min(min(Cp_all)); high = max(max(Cp_all));
ohm_range = high - low;
%figure(4); clf; hold on; box on; 
for i=1:18
    if i<9
        %plot(Cp_all(:,i), '.-', 'markersize', 20,  ...
            %'linewidth', 2,'displayname',['Channel: ', num2str(i)])
    else
        %plot(Cp_all(:,i), '*-', 'markersize', 8,  ...
            %'linewidth', 2,'displayname',['Channel: ', num2str(i)])
    end
end
%legend('location', 'eastoutside')
%ylabel('Capacitance [Farad]')
%set(gca, 'fontsize', 15, 'linewidth', 2)
%title('HIL Channel Cp Values')
%ylim([low-ohm_range*.1 high+ohm_range*.1])
%xlabel('Test Index [-]')

%%
Rp_all = R0_all*1000;
low = min(min(Rp_all)); high = max(max(Rp_all));
ohm_range = high - low;
%figure(4); clf; hold on; box on; 
for i=1:18
    if i<9
        %plot(Rp_all(:,i), '.-', 'markersize', 20,  ...
            %'linewidth', 2,'displayname',['Channel: ', num2str(i)])
    else
        %plot(Rp_all(:,i), '*-', 'markersize', 8,  ...
            %'linewidth', 2,'displayname',['Channel: ', num2str(i)])
    end
end
%legend('location', 'eastoutside')
%ylabel('Resistance [m\Omega]')
%set(gca, 'fontsize', 15, 'linewidth', 2)
%title('HIL Channel R0 Values')
%ylim([low-ohm_range*.1 high+ohm_range*.1])
%xlabel('Test Index [-]')

%%
% [R0, Rp, Cp, RMSE, MAE] = ParamID(SOC_f, I, V, t)

%%
function [R_0, R_p, C_p, RMSE, MAE] = ParamID(SOC_f, I, V, t)
    %% build OCV-SOC model
    [num1,txt1]=xlsread('SOC_curve.xls',1);% load open circuit voltage test data
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
    %figure(1)% model parameters 
    %clf;
    %subplot(3,1,1)
    %plot(t/3600,R0,'k','LineWidth',1)
%     ylim([-0.01 0.01]);
    %xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    %ylabel('R_{0}(\Omega)','Fontname', 'Times new roman','FontSize',8)
    %set(gca, 'Fontname', 'Times new roman', 'Fontsize',8,'LineWidth',1.2); 
    %hold on
    %subplot(3,1,2)
    %plot(t/3600,Rp,'k','LineWidth',1)
%     ylim([-50 50])
    %xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    %ylabel('R_{p}(\Omega)','Fontname', 'Times new roman','FontSize',8)
    %set(gca, 'Fontname', 'Times new roman', 'Fontsize',8,'LineWidth',1.2); 
    %hold on
    %subplot(3,1,3)
    %plot(t/3600,Cp,'k','LineWidth',1)
%     ylim([-5000 50000])
    %xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    %ylabel('C_{p}(F)','Fontname', 'Times new roman','FontSize',8)

    %figure(2) % voltage estimation error
    %clf;
    %plot(t/3600,error,'k','LineWidth',1)
%     ylim([-0.15 0.15]);
    %xlabel('Time(h)','Fontname', 'Times new roman','FontSize',8)
    %ylabel('voltage estimation error','Fontname', 'Times new roman','FontSize',8)
end

function [R0, Rp, Cp, RMSE, MAE] = ParamID_new(SOC_f,I,V,t)
    %% build OCV-SOC model
    [num1,txt1]=xlsread('SOC_curve.xls',1);% load open circuit voltage test data
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
    %figure(1) % voltage estimation error
    %clf;
    %plot(t/3600,error,'k','LineWidth',1)
    %xlabel('Time(h)','FontSize',15)
    %ylabel('voltage estimation error','FontSize',15)
    %set(gca, 'fontsize', 15)
end

%% FUNCTIONS
% ------------------------------------------------------------------------------------------
% Function getSOC -- calculates SOC given voltage V from OCV/SOC curve


% Function getCellCap -- calculates cell capacities
function [Step, Cap, data] = getCellCap(fname)
    % Inputs:
    % - fname = complete file path
    
    % Outputs:
    % - Cap = list of capacities for each cell
    % - Steps = vector of step indicies
    
    %% Load BMS File
    %FOR TEST "Nissan_Aging1wk_da2"
    fid = fopen(fname,'r');
    Original = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1,'Delimiter',',');

    % Data Columns
    t = Original{:,1}(1:length(Original{:,63}));
    % SOC0 = Original{:,2};
    Vpack0 = Original{:,3}(1:length(Original{:,63}));
    I0 = Original{:,4}(1:length(Original{:,63}));
    % Climit = Original{:,5};
    % Dlimit = Original(:,6);
    % HighCellID = Original{:,7};
    % HighCellV0 = Original{:,8};
    % LowCellID = Original{:,9};
    % LowCellV0 = Original{:,10};
    % CHAen = Original{:,12};
    % DCHen = Original{:,13};
    % HighT = Original{:,14};
    % LowT = Original{:,15};
    CV01 = Original{:,16}(1:length(Original{:,63}));
    CV02 = Original{:,17}(1:length(Original{:,63}));
    CV03 = Original{:,18}(1:length(Original{:,63}));
    CV04 = Original{:,19}(1:length(Original{:,63}));
    CV05 = Original{:,20}(1:length(Original{:,63}));
    CV06 = Original{:,21}(1:length(Original{:,63}));
    CV07 = Original{:,22}(1:length(Original{:,63}));
    CV08 = Original{:,23}(1:length(Original{:,63}));
    CV09 = Original{:,24}(1:length(Original{:,63}));
    CV10 = Original{:,25}(1:length(Original{:,63}));
    CV11 = Original{:,26}(1:length(Original{:,63}));
    CV12 = Original{:,27}(1:length(Original{:,63}));
    CV13 = Original{:,28}(1:length(Original{:,63}));
    CV14 = Original{:,29}(1:length(Original{:,63}));
    CV15 = Original{:,30}(1:length(Original{:,63}));
    CV16 = Original{:,31}(1:length(Original{:,63}));
    CV17 = Original{:,32}(1:length(Original{:,63}));
    CV18 = Original{:,33}(1:length(Original{:,63}));
    CV1_16 = [CV01 CV02 CV03 CV04 CV05 CV06 CV07 CV08 CV09 CV10 CV11 CV12 CV13 CV14 CV15 CV16 CV17 CV18];

    %% Clear extra measurements
    keep(1) = 1;
    n = 2;
    for i = 2:length(t)
        if (length(t{i}) == length(t{i-1}))
            if prod(t{i} == t{i-1}) == 0
                keep(n) = i;
                n = n+1;
            end
        end
    end
    t0 = t(keep); %Use for analysis
    Vpack = Vpack0(keep);
    I = I0(keep);
    CV = CV1_16(transpose(keep),:);

    %% Convert timestamps to hours (takes a while)
    if  1==1    %str2num(fname(end-4))>1 || ~isempty(str2num(fname(end-5)))||
        for i = 1:length(t0)
            X1 = regexp(t0{i},' ','split');
            X2 = regexp(X1{4},':','split');
            X3 = [X1{3},'-', X1{2}, '-', X1{6}, '-', X1{4}];
            X4 = datenum(X3,'dd-mmm-yyyy HH:MM:SS');
            t1(i) = X4;
        end
        t1 = (transpose(t1)-t1(1))*24;
    else
        t1 = str2double(t0);
    end
    data = [t1 I Vpack CV];
    %[hours amps packV cellV's]
    %writematrix(data,filename) 

    %% Start of Steps
    L = length(data(:,2));
    n = 1;
    for i=1:L-1
        if abs(data(i,2)) > 0
            if data(i+1,2) == 0
                Step(n) = i+1;
                n = n+1;
            end
        else
            if abs(data(i+1,2)) > 0
                Step(n) = i+1;
                n = n+1;
            end
        end
    end
    %Start/End Voltages of Full Discharge
    StartVs = data(Step(17)-1,4:end);
    EndVs = data(Step(23)-1,4:end); %after pause

    %% Convert Voltages to SOH to capcity
    SOCcurve = xlsread('TeslaOCVcurve_matt.xls');
    x=SOCcurve(:,1);
    y=SOCcurve(:,2);
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

    for i = 1:length(EndVs)
    V = EndVs(i); %Put Voltage to get SOC
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
    %% Integrate Discharge AmpHours on BMS
    %get discharge current values

    % get time between steps
    n = 1;
    for i = 2:length(data)
        newt(i) = data(i,1)-data(i-1,1);
    end
    newt = transpose(newt);

    % Use riemann sum to get Ah between steps
    for i = 1:length(Step)-1
        Ah0 = newt.*data(:,2);
           Ah = sum(Ah0(Step(i):Step(i+1)-1));
           StepAh(i) = Ah;
    end

    CapDchAh = 40*(data(Step(18),1)-data(Step(17),1)); %battery DCH capacity from time at constant current

    %% Capacity Calculation
    %CapDchAh = 39.334124690622275;
    for i =1:length(SOCstart)
       Cap(i) = (100/(SOCstart(i)-SOCend(i)))*CapDchAh; %Capacity
    end
    
end
