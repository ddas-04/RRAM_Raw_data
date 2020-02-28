clc
clear all
close all
%%
distribution=["normal", "lognormal", "weibull", "rayleigh", "logistic", "exponential"];
%% data load -----> format  xlsread('filengthame.xlsx','sheet name')
filename='100devices.xlsx';
data_100cycle=xlsread(filename,'100cycles');
data_100devices=xlsread(filename,'devices');

%% Read voltage

target=input('Enter read voltage of RRAM : ');

%% Variable declaration for cycle to cycle
[total_row,total_col]=size(data_100cycle);
n=total_col/2;
LRS_100cycles=zeros(1,n);
HRS_100cycles=zeros(1,n);

%%
icol=1;

i=1;
while icol<=total_col-1
    V_index=find(data_100cycle(:,icol)==target);
    I_val=zeros(1,length(V_index));
    for iv=1:length(V_index)
        icol;
        I_val(iv)=data_100cycle(V_index(iv),icol+1);
    end
    
    if length(V_index)~=2
        print('More/less than two values exists in current, so LRS and HRS are not sufficient')
        break
    else
        I_sorted=sort(I_val);
        HRS_100cycles(i)=abs(target)/I_sorted(1); % as lower current means higher resistance 
        LRS_100cycles(i)=abs(target)/I_sorted(2); % as higher current means lower resistance 
    end
    
    i=i+1;
    icol=icol+2;
end
%         
% figure            
% subplot(3,3,1)
% probplot(distribution(1),LRS_100cycles)
% legend('Fitting',strcat('LRS Resistance at Vread =',num2str(target)),'Location','best')
% grid on
% xlabel('Resistance Cycle to Cycle(LRS)');
% %title('Cycle to Cycle vartiation');
% str=['lognormal'];
% subplot(3,3,2)
% probplot(str,HRS_100cycles)
% %legend('Fitting','HRS Resistance at Vread =-0.1 V','Location','best')
% legend('Fitting',strcat('HRS Resistance at Vread =',num2str(target)),'Location','best')
% grid on
% xlabel('Resistance Cycle to Cycle(HRS)');      
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable declaration for device to device
[total_row,total_col]=size(data_100devices);
n=total_col/2;
LRS_100devices=zeros(1,n);
HRS_100devices=zeros(1,n);
%%
icol=1;

i=1;
while icol<=total_col-1
    V_index=find(data_100devices(:,icol)==target);
    I_val=zeros(1,length(V_index));
    for iv=1:length(V_index)
        I_val(iv)=data_100devices(V_index(iv),icol+1);
    end
    
    if length(V_index)~=2
        print('More/less than two values exists in current, so LRS and HRS are not sufficient')
        break
    else
        I_sorted=sort(I_val);
        HRS_100devices(i)=abs(target)/I_sorted(1); % as lower current means higher resistance 
        LRS_100devices(i)=abs(target)/I_sorted(2); % as higher current means lower resistance 
    end
    i=i+1;
    icol=icol+2;
end

% figure(2)            
% subplot(1,2,1)
% probplot(LRS_100devices)
% legend('Fitting',strcat('LRS Resistance at Vread =',num2str(target)),'Location','best')
% grid on
% xlabel('Resistance Device to Device(LRS)');
% 
% subplot(1,2,2)
% probplot('lognormal',HRS_100devices)
% %legend('Fitting','HRS Resistance at Vread =-0.1 V','Location','best')
% legend('Fitting',strcat('HRS Resistance at Vread =',num2str(target)),'Location','best')
% grid on
% xlabel('Resistance Device to Device(HRS)');


fig_count=1;

%% All distribution LRS plot for cycle to cycle
figure(fig_count)
for i=1:6
    subplot(2,3,i)
    probplot(distribution(i),LRS_100cycles)
    legend('Fitting',strcat('LRS cycle to cycle at Vread =',num2str(target)),'Location','best')
    grid on
    xlabel('Resistance Cycle to Cycle(LRS)');
end
fig_count=fig_count+1;
%% All distribution HRS plot for cycle to cycle
figure(fig_count)
for i=1:6
    subplot(2,3,i)
    probplot(distribution(i),HRS_100cycles)
    legend('Fitting',strcat('HRS cycle to cycle at Vread =',num2str(target)),'Location','best')
    grid on
    xlabel('Resistance Cycle to Cycle(HRS)');
end
fig_count=fig_count+1;
%% All distribution LRS plot for device to device
figure(fig_count)
for i=1:6
    subplot(2,3,i)
    probplot(distribution(i),LRS_100devices)
    legend('Fitting',strcat('LRS device to device at Vread =',num2str(target)),'Location','best')
    grid on
    xlabel('Resistance device to device(LRS)');
end
fig_count=fig_count+1;
%% All distribution HRS plot for device to device
figure(fig_count)
for i=1:6
    subplot(2,3,i)
    probplot(distribution(i),HRS_100devices)
    legend('Fitting',strcat('HRS device to device at Vread =',num2str(target)),'Location','best')
    grid on
    xlabel('Resistance device to device(HRS)');
end
fig_count=fig_count+1;
    