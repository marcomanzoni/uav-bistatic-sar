clear all
close all
clc
%% =============================================================================  LOAD RC MATRIX CUT IN TARGET ZONE
%% COSTANTS
chirp_bw = 56e6*0.9;                                % actual chirp bandwidth (B)
chirp_sr = 56e6;                                % SDR sample rate (fs)
norm_B = chirp_bw / chirp_sr;

experiment_name = "test1_cut1";

radar_folder_name = '..\mat_files\uav_test\20220502\radar\test1\';
drone_track_folder = '..\mat_files\uav_test\20220502\drone track\';

drone_data = load(strcat(drone_track_folder, 'flight1')).drone;
radar_last_mod_date = load(strcat(radar_folder_name,'test1_last_mod_date')).last_mod_date;
targets = load(strcat(drone_track_folder,'target.mat'));

tx_wave = load(strcat('..\tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);
f0 = 1.65e9;
lambda = physconst('LightSpeed')/f0;
PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt);

%% RC Load
addpath(genpath('..\lib' ));       % add path of lib

RC_file = strcat(radar_folder_name,experiment_name);

RC = load_bin(RC_file);
R_margin = size(RC,1)*dR;
R_ax = -R_margin/2:dR:R_margin/2;
tau_ax = 0:PRI:PRI*size(RC,2);
N_PRI_tot = length(tau_ax);

%% taglio sul cross talk
idx_wind_row = 150:300;%180:230;
idx_wind_col = 1:2.4e4;%1:1.4e4;

RC = RC(idx_wind_row,idx_wind_col);
R_ax = R_ax(idx_wind_row); 
tau_ax = tau_ax(idx_wind_col); 

N_PRI = length(tau_ax);
%%
figure
imagesc(abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

figure
imagesc(tau_ax,R_ax,angle(RC))
title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

%%  ======================================================================= SYNCHRONISM ERROR CORRECTION
% Distance tx-rx for each PRI (cross-talk distance)
run('DronePath')

%%
% Synchronism correction algorythm
run('RC_SynchrCorrection.m')

%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on
plot(tau_ax,distance_tx_rx,'r','LineWidth',1.5),
plot(tau_ax,distance_cars(1,:),'g')
legend('crosstalk','car')%,'car 2','human1','human2','human3');

%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on

plot(tau_ax,distance_tx_rx),
plot(tau_ax,distance_cars(1,:),'m'), plot(tau_ax,distance_cars(2,:),'g'), 
plot(tau_ax,distance_humans(1,:),'y'), plot(tau_ax,distance_humans(2,:),'m'), 
plot(tau_ax,distance_humans(3,:),'b')
legend('crosstalk','car','car 2','human1','human2','human3');

%% ===========================================================================  CROSSRTALK REMOVAL - BRUTAL
% cross_talk_size_idxs = 10;
% 
% for n = 1:N_PRI
%     RC(1:cross_talk_idxs_corr(n)+cross_talk_size_idxs,n) = zeros(cross_talk_idxs_corr(n)+cross_talk_size_idxs,1);
% end
% 
% %% Plot
% figure
% imagesc(tau_ax,R_ax,abs(RC))
% title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
% 
% figure
% imagesc(tau_ax,R_ax,angle(RC))
% title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% =========================================================================== FOCUSING

%% FOCUSING BY TIME DOMAIN BACK PROJECTION (TDBP)
% In this case, we choose to focus on the plane x,y  with z = c

%% Set XY grid where to focus
y_min = -30;
y_max = 100;


x_min = -40;
x_max = 200;
pho_az = 2;
dx = pho_az*0.4;
x_ax =  x_min:dx:x_max;
dy = dx;
y_ax =  y_min:dy:y_max;

delta_psi_proc = lambda/pho_az;

[X,Y] = ndgrid(x_ax,y_ax);

z0 = 0;
%% Focusing
wbar = waitbar(0,'Backprojecting');
Ny = length(y_ax);
Nx = length(x_ax);
S = zeros(Nx,Ny);
A = S;
% t = (0:size(RC,1)-1) *dt/OSF ;
t = R_ax./c;

clear Sn    

%Compute path angle
psi_path = atan( (RX_pos(2,end)-RX_pos(2,1)) / (RX_pos(1,end)-RX_pos(1,1)) );
psi_point = psi_path+deg2rad(90);


for n = 1:N_PRI
    waitbar(n/N_PRI,wbar)
   
    % Distance 
    R_tx = sqrt((TX_pos(1,n)-X).^2 + (TX_pos(2,n)-Y).^2  + (TX_pos(3,n)-z0).^2); %  Range distances from the tx antenna [m]
    R_rx = sqrt((RX_pos(1,n)-X).^2 + (RX_pos(2,n)-Y).^2  + (RX_pos(3,n)-z0).^2); %  Range distances from the rx antenna [m]
    distance = R_tx+R_rx; %Total Tx-target-Rx distance [m]
    delay = distance./c;    %Delay
    
    %Compute target angle
    sin_psi = (RX_pos(2,n)-Y)./R_rx;
    psi = asin(sin_psi);
    psi_targ = psi - psi_point;
    
    %Activation function
    Wn = rectpuls(psi_targ/delta_psi_proc);
    cut = find(x_ax>RX_pos(1,n));
    cut = cut(1);
    Wn(1:cut,:) = zeros(size(Wn(1:cut,:)));
%     if not(mod(n,500))
%         imagesc(x_ax,y_ax,Wn.'),hold on
%         plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
%         plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),hold off
% 
%         drawnow
%     end

    % Backprojection of data from a single Radar position 
%     Sn = Wn.*interp1(t,RC(:,n),delay).*exp(1i*2*pi*f0*delay);
    Sn = interp1(t,RC(:,n),delay).*exp(1i*2*pi*f0*delay);
    % Coherent sum over all positions along the trajectory 
    S = S + Sn;
    % Inchoerent sum over all positions along the trajectory (choerent sum)
    A = A + abs(Sn);
end
close(wbar)
%%
Focus = (S./A).';
% Focus = S.';
%Focus = A.';%%

figure,
imagesc(x_ax,y_ax,abs(Focus)), axis xy , 
%title('Final Sensor Position') 
title('Focused image by TDBP ')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start drone track','end drone track','TX','Cars','Humans');
%% Path and pointing vector plot
%     x_vect = 0:dx:x_ax(end)-RX_pos(1,n);
%     pointing = RX_pos(2,n) + (x_vect)*tan(psi_path);
%     pointing = [ones(length(x_ax)-length(pointing),1)', pointing];
%     hold on, plot(x_ax,pointing,'y')
%     
%     pointing = RX_pos(2,n) + (x_vect)*tan(psi_point);
%     pointing = [ones(length(x_ax)-length(pointing),1)', pointing];
%     hold on, plot(x_ax,pointing,'r')
%     
 
%%
figure,
imagesc(x_ax,y_ax,10*log10(abs(Focus))), axis xy , %colormap jet
%title('Initial Sensor Position') 
title('Focused image by TDBP - dB')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start swath','end swath','TX','Cars','Humans');

