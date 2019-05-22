%STA_LTA: Short-Time-Average/Long-Time-Average event detector.
%
%USAGE: events = sta_lta(wave,prop_name_1,prop_val_1,...)
%
%DIAGRAM:
%                       /\      /\/\        /\
%              /\  /\/\/  \/\  /    \/\    /  \/\
%                \/          \/        \/\/      \
%                                           |-STA-| -->
%         --> |-----------------LTA---------------| -->
%
%INPUTS: wave     - A waveform object containing events (maybe)
%        varargin - User-defined parameter name/value pairs (below)
%
%VALID PROP_NAME: 
%
%  'edp'      - Event Detection Parameters (Detailed in 'VALID PROP_VAL')
%  'lta_mode' - LTA window behavior following a 'trigger on' flag
%  'skip'     - Times for STA/LTA to skip over (if known), examples are 
%               data gaps, calibration pulses, periods of excessive noise
%  'eot'      - Event output type (Detailed in 'VALID PROP_VAL')
%  'fix'      - Output event length fixed or variable?
%  'pad'      - Pad output events by a fixed amount of time before start 
%               time (trigger on) and after stop time (trigger off).
%
%VALID PROP_VAL: 
%
%  'edp'(1x6 numeric) ... default: [1 8 2 1.6 0 3]
%     --> [l_sta l_lta th_on th_off min_sep min_dur]
%        --> l_sta    - STA window length (s)
%        --> l_lta    - LTA window length (s)
%        --> th_on    - STA/LTA trigger on threshold
%        --> th_off   - STA/LTA trigger off threshold
%        --> min_sep  - Minimum seperation between trigger off time of an 
%                       event and trigger on time of the next event (s)
%        --> min_dur  - Minimum event duration to be recorded(s)
%
%	'lta_mode' (string) ... default: 'continuous'
%     --> 'frozen' - LTA window fixed in place after trigger is turned on 
%                    while STA window continues forward.
%     --> 'continuous' - LTA window continues w/ STA window after trigger 
%                        is turned on (Same behavior as before trigger)
%     --> 'grow' - LTA left edge fixed, right edges continues w/ STA right 
%                  edge after trigger on, thus size of LTA window is 
%                  growing until trigger off, at which point LTA window 
%                  returns to original size.
%
%	'skip' (nx2 numeric) ... default: [] (don't skip anything)
%     --> skip_val - List of start/stop times to skip
%
%	'eot' (string) ... default: 'sst'
%     --> 'st'  - Return event start times (nx1 numeric)
%                 (a.k.a. trigger on times)
%     --> 'sst' - Return event start/stop times (nx2 numeric)
%                 (a.k.a. trigger on/off times) 
%     --> 'ssd' - Return event start/stop datapoints (nx2 integer)
%                 (referenced from first datapoint in 'wave')
%     --> 'wfa' - Return event waveform array (1xn waveform)
%
%   'fix' (1x1 numeric) ... default: 0
%     --> 'fix_val' - If 0, output events will be variable length 
%                     (trigger on time to trigger off time)
%                     If a positive numeric value N, output events will be
%                     fixed length (trigger on time to fixed length N)
%                     Units of N is seconds if non-zero
%
%   'pad' (1x2 numeric) ... default: [0 0] (i.e. no padding)
%     --> 'pad_val' - If [0 0] events will have no time padding
%                     [1 2] - pad by 1s before, 2s after event
%                     this will pad events regarless of the value of fix:
%                     i.e. fix = 7, pad = [1 2], returned events are 10s
%                         
%OUTPUTS: events - events in format specified by 'return'

% Author: Dane Ketner, Alaska Volcano Observatory
% $Date$
% $Revision$

%% Initialize Waveform Variables
clear;clc;close;
%%
% load mic_syn.mat
load Data/data0_syn2.mat

dt = 0.001;
T = 3.1;
nt = size(d, 1);
nx = size(d, 2);
dx = 7.5;
x = 0 : dx : (nx-1)*dx;

%% interpolation
dt_fine = 0.0005;
t = 0 : dt : T;
t_fine = 0 : dt_fine : T;
nt_fine = length(t_fine);
d_fine = zeros(nt_fine, nx);
for ix = 1 : nx
    d_fine(:, ix) = interp1(t, d(:, ix), t_fine, 'spline');
end

dt = dt_fine;
d = d_fine;
t = t_fine;
nt = nt_fine;
    
%% add random noise
d_noise = awgn(d, 0.1,'measured','linear');

figure;
imagesc(d_noise);

%%

wave.data = d_noise(:);
nt = size(wave.data,1);
t = 0 : dt : (nt-1)*dt;
wave.freq = 1 / dt;
wave.data_length = nt;
wave.timevector = t;

%% Check waveform variable
Fs = wave.freq;         % Sampling frequency
l_v = wave.data_length; % Length of time series
tv = wave.timevector;  % Time vector of waveform

%% Set all default parameters
l_sta = 58*2;     % STA window length
l_lta = 8*l_sta;     % LTA window length
th_on = 1.2;        % Trigger on when sta_to_lta exceeds this theshold
th_off = 0.5;     % Trigger off when sta_to_lta drops below threshold
min_sep = 0*Fs;   % Skip ahead after end of event
min_dur = 58*2;   % Any triggers shorter than min_dur are discarded

lta_mode = 'continuous'; % Post trigger-on LTA behavior
skip = [];               % List of start/stop times to skip
eot = 'wfa';             % Event output type 
fix = 0*Fs;              % Fix event length
pad = [0 0]*Fs;          % Pad event start/stop times



%% Initialize waveform data
% wave = set2val(wave,skip,NaN);   % NaN skip times
wave.data = zero2nan(wave.data,10);        % NaN any data gaps
v = wave.data;            % Waveform data
abs_v = abs(v);                  % Absolute value of time series

%% Initialize flags and other variables
lta_calc_flag = 0;       % has the full LTA window been calculated?
ntrig = 0;               % number of triggers
trig_array = zeros(1,2); % array of trigger times: [on,off;on,off;...]

%% Loops over data
% i is the primary reference point (right end of STA/LTA window)
sta_to_lta_array = zeros(nt, 1);
i = l_lta+1;
while i <= l_v % START STA_LTA MAIN LOOP

%% Skip data gaps (NaN values in LTA window)?
   if any(isnan(abs_v(i-l_lta:i)))
      gap = 1;
      lta_calc_flag = 0; % Force full claculations after gap
      while (gap == 1) && (i < l_v)
         i = i+1;
         if ~any(isnan(abs_v(i-l_lta:i)))
            gap = 0;
         end
      end
   end

%% Calculate STA & LTA Sum (Do Full Calculation?)
   if (lta_calc_flag == 0)
      lta_sum = 0;
      sta_sum = 0;
      for j = i-l_lta:i-1              % Loop to compute LTA & STA
         lta_sum = lta_sum + abs_v(j); % Sum LTA window
         if (i - j) <= l_sta           % Sum STA window (right side of LTA)
            sta_sum = sta_sum + abs_v(j);
         end
      end
      lta_calc_flag = 1;
   else
      
%% Calculate STA & LTA Sum (Single new data point if not Full) 
      lta_sum = lta_sum - abs_v(i-l_lta-1) + abs_v(i-1);
      sta_sum = sta_sum - abs_v(i-l_sta-1) + abs_v(i-1);
   end

%% Calculate STA & LTA
   lta = lta_sum/l_lta;
   sta = sta_sum/l_sta;

%% Calculate STA/LTA Ratio
   sta_to_lta = sta/lta;
   sta_to_lta_array(i) = sta_to_lta;
% %% Trigger on? (Y/N)
%    if (sta_to_lta > th_on)
%       j = i;   % Set secondary reference point = primary
%       g = 0;   % l_lta growth, only used if LTA growing
%       while (sta_to_lta > th_off)
%          j = j+1;
%          if j < l_v
%             sta_sum = sta_sum - abs_v(j-l_sta-1) + abs_v(j-1);
%             switch lta_mode
%                case 'frozen'
%                   % LTA is good just the way it is
%                case 'continuous'
%                   % Add new data point, remove oldest data point
%                   lta_sum = lta_sum - abs_v(j-l_lta-1) + abs_v(j-1);
%                case 'grow'
%                   % Add new data point, increase
%                   lta_sum = lta_sum + abs_v(j-1);
%                   l_lta = l_lta + 1;
%                   g = g+1;
%             end
%             sta = sta_sum/l_sta;
%             lta = lta_sum/l_lta;
%             sta_to_lta = sta/lta;
%             if any(isnan(abs_v(j-l_sta:j))) % NaN gaps to skip?
%                sta_to_lta = 0; % Force trigger off (data gap in STA window)
%             end
%          else
%             sta_to_lta = 0; % Force trigger off (end of data)
%          end
%       end
%       duration = (j-i); % span from trigger on to trigger off
%       l_lta = l_lta-g;
% 
% %% Triggered period long enough? (Y/N)
%       if duration > min_dur % If duration < min_dur then skip it
%          trig_t = i-l_sta;  % Beginning of STA window during trigger on
%          end_t  = j;        % End of STA window during trigger off
%          ntrig = ntrig + 1; % Event counter
%          trig_array(ntrig,:) = [trig_t, end_t];
%       end
%       i = j + min_sep;   % Skip ahead by minimum event seperation
%       lta_calc_flag = 0; % Reset LTA calc flag to force new computation
%    end
   i = i + 1;
end % END STA_LTA MAIN LOOP

% test = zeros(nt, 1);
% for i = 1 : size(trig_array, 1)
%     test(trig_array(i, 1) : trig_array(i, 2)) = wave.data(trig_array(i, 1) : trig_array(i, 2));
% end
sta_to_lta_toshow = reshape(sta_to_lta_array, nt_fine, nx);
figure;imagesc(sta_to_lta_toshow)
figure(1);imagesc((sta_to_lta_toshow>1.05) .* d_noise );
% (sta_to_lta_toshow>1.01) .* d_noise


load cgray.mat
fig=figure(1)
clf
%ha = tight_subplot(Nh, Nw, gap, marg_ver, marg_hor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha = tight_subplot(1,1,[0.01 0.01],[0.15 0.05],[0.2 0.1]);
axes(ha(1)); imagesc(x(6:end-5),t,d_noise(:, 6:end-5).*(sta_to_lta_toshow(:, 6:end-5)>1.07), [-0.15 0.15]);colorbar('eastoutside');colormap(cgray)

axes(ha(1));
ylabel('Time [s]')
xlabel('x [m]')
set(gca,'fontsize',30); 

set(fig,'paperposition',[0 0.1 10 10])


print -depsc Fig/d2_ltasta_noise_n01.eps
print -djpeg Fig/d2_ltasta_noise_n01.jpg
