function img_corrected = pga(img, params)
%PGA performs one iteration of phase gradient autofocus on SAR image
%   PGA(IMG) for complex-valued SAR image IMG will perform phase gradient
%   autofocus as described in [1]. Images must be rotated such that
%   cross-range is horizontal and down-range is vertical.
%
%   PGA(IMG, PARAMS) where PARAMS is a structure allows tuning of the
%   algorithm. Members needed:
%   - PARAMS.db_down_for_window - dB below peak for automatic windowing to
%     consider (sign doesn't matter). Defaults to 25 dB.
%   - PARAMS.want_plots - 1 if plots of window selected and phase estimated
%     are desired, 0 if not. Defaults to 1.
%
%   FIXME: trivial to extend this for multiple iterations.
%
%   NOTE:  A different, minimum-variance and unbiased, phase estimator is
%   presented in [2] but I can't get it to produce any results. Its
%   implementation is disabled.
%
%   NOTE: In [2], a different automatic windowing system is described: "We
%   have found that a progressively decreasing window width works quite
%   well for low contrast scenes. The initial window width is selected to
%   span the maximum possible blur width (which could be several hundred
%   samples for extremely defocused imagery) and progressively reduced by
%   20% for each iteration." FIXME?
%
%   References:
%
%   [1] C.V. Jakowatz, D.E. Wahl, P.H. Eichel, D.C. Ghiglia, P.A. Thompson,
%   {\em Spotlight-mode Synthetic Aperture Radar: A Signal Processing
%   Approach.} Springer, 1996.
%
%   [2] D.E. Wahl, P.H. Eichel, D.C. Ghiglia, C.V. Jakowatz, "Phase
%   gradient autofocus-a robust tool for high resolution SAR phase
%   correction," IEEE Trans. Aero. & Elec. Sys., vol. 30, no. 3, 1994.
%
% Nominal version 0.01.

if (nargin < 2)
    params.db_down_for_window = 25;
    params.want_plots = 1;
end

center_az_idx = ceil(size(img, 2)/2);

orig_img = img;

%% Center shift largest targets
[tmp maximum_along_az_idx] = max(abs(img), [], 2);
for i = 1:size(img, 1)
    img(i,:) = circshift(transpose(img(i,:)), center_az_idx - maximum_along_az_idx(i));
end

%% Determine window width
noncoh_avg_window = sum(abs(img).^2, 1);

window_cutoff = max(db20(noncoh_avg_window)) - abs(params.db_down_for_window);

if (params.want_plots)
    figure;
    subplot(211)
    plot(db20(noncoh_avg_window)-window_cutoff )
    title('s')
end

leftidx  = find(db20(noncoh_avg_window(1:center_az_idx    )) - window_cutoff<0, 1, 'last' );
rightidx = find(db20(noncoh_avg_window(center_az_idx+1:end)) - window_cutoff<0, 1, 'first');
leftidx = leftidx+1;
rightidx = rightidx + center_az_idx - 1;

noncoh_avg_window = zeros(size(noncoh_avg_window));
noncoh_avg_window(leftidx:rightidx) = 1;

if (params.want_plots)
    hold on
    plot((db20(noncoh_avg_window) - abs(params.db_down_for_window) ).*noncoh_avg_window,'r--');
    title('s')
end

%% Apply window, DFT resulting range lines and, estimate phase error
fft_length_pow_2 = 2^ceil(log2(size(img,2)));

imgF = fft(...
    ifftshift( img .* repmat(noncoh_avg_window, size(img,1), 1), ...
    2),...
    fft_length_pow_2, 2);

if 1
    delta_phase = angle( sum(  conj(imgF(:, 1:end-1  )) .* imgF(:, 2:end) , 1) );
    phase_estim = [0 cumsum(delta_phase)];
else
    % This is the LUMV phase estimator described in [2] above.
    % phase_estim = sum(imag(conj(imgF).*imgF), 1) ./ sum(abs(imgF).^2, 1);
end

%% Remove linear trend in phase error, as suggested in [2], and remove
linear_coefs = polyfit(1:length(phase_estim), phase_estim, 1);
%phase_estim = unwrap(phase_estim - polyval(linear_coefs, 1:length(phase_estim)));

if (params.want_plots)
    subplot(212)
    plot(phase_estim);
    title('\phi hat')
end

% Save memory
clear img

imgF = fft(...
    ifftshift( orig_img, 2),...
    fft_length_pow_2, 2);

%imgF = imgF .* repmat(exp(-1j*phase_estim), fft_length_pow_2, 1);
imgF = imgF .* repmat(exp(-1j*phase_estim), size(imgF,1), 1);
img_corrected = fftshift(ifft(imgF, [], 2), 2);

end


function foo = db20(bar)
%DB20 returns 20*log10(abs()) of input. This can be slow---so profile!
foo=20*log10(abs(bar));
end