function [focus] = focusingTDBP(const,radar,scenario,RX,TX,ang_vec)
%FOCUSINGTDBP compute the focusing on the defined grid with TDBP
%   [focus] = focusingTDBP(const,radar,scenario,RX,TX,ang_vec)

Ny = length(scenario.grid.y_ax);
Nx = length(scenario.grid.x_ax);
t = radar.R_ax./physconst('LightSpeed');

% Processed Antenna aperture
focus.psi_proc = const.lambda / 2 / scenario.grid.pho_az;
focus.R_min = 50;
focus.synt_apert = 2 * tan(focus.psi_proc/2) * focus.R_min;
% Processed wavenumbers
Dk = 2*pi/scenario.grid.pho_az;

%Sqint angle vectors
focus.angle_vec = ang_vec;


wbar = waitbar(0,strcat('Backprojecting n 1/',num2str(length(focus.angle_vec))));

% copy variables for optimizing parfor
idxs = t >=0;
t = t(idxs);
RC = radar.RC(idxs,:);

TX_pos_x = TX.pos(1,:);TX_pos_y = TX.pos(2,:);TX_pos_z = TX.pos(3,:); 
RX_pos_x = RX.pos(1,:);RX_pos_y = RX.pos(2,:);RX_pos_z = RX.pos(3,:); 
RX_speed = RX.speed;
X = scenario.grid.X; Y = scenario.grid.Y; z0 = scenario.grid.z0;
lambda = const.lambda; f0 = const.f0;
x_ax = scenario.grid.x_ax;
median_speed = median(RX_speed);

%Initialize vectors for the result
focus.Focused_vec = zeros(size(X,1),size(X,2),size(focus.angle_vec,1));
focus.not_coh_sum = zeros(size(focus.Focused_vec));
focus.SumCount = zeros(size(focus.Focused_vec));

tic
for ang_idx = 1:length(focus.angle_vec)
    waitbar(ang_idx/length(focus.angle_vec),wbar,strcat('Backprojecting n '...
        ,num2str(ang_idx),"/",num2str(length(focus.angle_vec))));
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/const.lambda); 
 
    S = zeros(Nx,Ny);
    SumCount = zeros(Nx,Ny);
    
    parfor n = 1:size(RC,2)
        [Sn,Wn] = elementFuncTDBP(X,Y,z0,TX_pos_x(n),TX_pos_y(n),TX_pos_z(n),RX_pos_x(n),...
           RX_pos_y(n),RX_pos_z(n),lambda,Dk,RC(:,n),t,f0,k_rx_0,x_ax);
        
        % Give less weight to not moving positions
        speed_norm = RX_speed(n)/median_speed;
        % Count number of summations for each pixel
        SumCount = SumCount + speed_norm.*Wn;
        
        % Coherent sum over all positions along the trajectory 
        S = S + speed_norm .* Sn;
        % Inchoerent sum over all positions along the trajectory
%         A = A + abs(Sn);
    end
    waitbar(ang_idx/length(focus.angle_vec),wbar);
    


% 	Focus = (S./A)';
%     Focus = S;
%     Focus = A';
    focus.Focused_vec(:,:,ang_idx) = S;
    focus.SumCount(:,:,ang_idx) = SumCount;
%     focus.not_coh_sum{ang_idx} = A; 
end

close(wbar)

disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

