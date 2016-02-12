function [wavepackets] = findwavepackets(x,params)
%% Wave Packet Finder
% DE Frederick
% Kay Lab (kaylab.uchicago.edu/)
% 2014
% 
% Designed for Frederick et al., 2016
% 
% Set of greedy algorithms backed by wavelet analysis for finding wave
% packets. 
% 
% Builds upon the concept that a neural oscillation can be defined as a
% wave packet, something that is discrete. 
% 
% These packets reflect discrete coordinated (or synchronus, or
% co-temporaneous) events above background noise. 
% 
% Freeman WJ (2003) The wave packet: an action potential for the 21st century. 
% Journal of Integrative Neuroscience 2:3–30.
% 
% 
% Dependencies
% MATLAB v2013b or higher.
% |--- table data type
% 
% MATLAB Wavelet Toolbox
% MATLAB Curve Fitting Toolbox
% 
% GNU Statement
% 
% $Verison History$
% $v.0.0.1 - Feb, 2014

% Input
% data <- a row vector of (filtered) data

%% Params package
% .drawwavepackets = true,false
% .exitcondition = 1,2
% .thres = int, default 3
% .maxwidth
% .minwidth
% .lag
% .sf
% .freqrange

%% Make sure data is a vector, not a matrix
if size(x,1)>1 && size(x,2)>1
	error('Data must be a vector')
end

% Change to row if column
if iscolumn(x)
	x=x';
end

%% Params
if isfield(params,'drawwavepackets')
    DRAWWAVEPACKETS = params.drawwavepackets;
else
    DRAWWAVEPACKETS= 0; %default
end


%set type of exit condition;
if isfield(params,'exitcondition')
    EXITCONDITION = params.exitcondition;
else
    EXITCONDITION = 1; %1,2
end

%set threshold
if isfield(params,'thres')
    THRES = params.thres;
else
    THRES = 3;
end

%set max width in ms
if isfield(params,'maxwidth')
    MAXWIDTH=params.maxwidth;
else
    MAXWIDTH = 0.080; % max width in ms; % @2020.2 fs, 10pts~5ms, 50pts~25ms  100pts~50ms
end

%set min width
if isfield(params,'minwidth')
    MINWIDTH = params.minwidth;
else
    MINWIDTH = 0.030; % min width in ms
end

%set step size
if isfield(params,'stepsize')
    STEPSIZE = params.stepsize;
else
    STEPSIZE = 2;     % initial step (in pts) to take left (or right)
end

%set lag step size
if isfield(params,'lag')
    LAG = params.lag;
else
    LAG = 5;          % number of points for edge detection
end

%set sampling frequency
if isfield(params,'sf')
    SF = params.sf;
else
    SF = 2020.2;      % sampling frequency
end


%% Wavelet parameters
wname='morl'; %select wavelet type; morlet works well
fc=centfrq(wname);

%% Setup wavepackets table to store waves and wave stats -----------------%
wavepackets = table(); 
wavepackets.record = 0; % surrogate key
wavepackets.start  = 0; % start of wave; relative to total vector, x
wavepackets.peak   = 0; % peak of wave;  relative to total vector, x
wavepackets.stop   = 0; % end of wave;   relative to total vector, x
wavepackets.length = 0;
wavepackets.ifreq  = 0;
wavepackets.iamp   = 0;
wavepackets.wvfreq = 0; % peak wavelet frequency
wavepackets.wvpow  = 0; % average wavelet power in segment 
wavepackets.snr    = 0; % signal-to-noise 
wavepackets.wave   = {0};
wavepackets = repmat(wavepackets,1000,1);
null_wavepackets = wavepackets; %save to grow table, if needed

% Convert max width from ms to pts based upon sampling frequency
MAXWIDTH=round(MAXWIDTH*SF);
if mod(MAXWIDTH,2)==1 %force even
    MAXWIDTH=MAXWIDTH+1;
end
% Convert min width from ms to pts based upon sampling frequency
MINWIDTH = round(MINWIDTH*SF);
if mod(MINWIDTH,2)==1 %force even
    MINWIDTH=MINWIDTH+1;
end


% Set wavelet parameters
if ~isfield(params, 'freqrange')
	error('Need frequency range.')
end
freqrange  = params.freqrange;
scalerange = fc./(freqrange*(1/SF));
scales = scalerange(end):0.2:scalerange(1);


% Do Wavelet (CWT)
coefs     = cwt(x,scales,wname); %
coefs     = abs(coefs.*coefs);    
coefs     = flipud(coefs);    % for some reason, need to flip the coefs mat
energyProfile = sum(coefs,1); % sum down the rows (@time pt)


% Start constructing the plots
if DRAWWAVEPACKETS
	% Find biggest peak and then get it's matrix location
	%[m1,m2]=max(coefs(:));
	%[ij,ji]=ind2sub(size(coefs),m2);

	figure;hold on;
	subplot(4,1,1:2);hold on;
%     xt=1:size(coefs,2);
%     xt=xt./SF;
    psfrq = fliplr(scal2frq(scales,'morl',1/SF));
	pcolor(1:size(coefs,2),psfrq, coefs); shading interp;
    
    xlim([1 size(coefs,2)])
	ylim([min(psfrq) max(psfrq)])

	subplot(4,1,4);hold on;
	plot(x);
	xlim([0 numel(x)])

	subplot(4,1,3);hold on;
	%             plot( s )
end

% Local interpolation of local maxima using cubic splines
% need this to create a smoother function. 
s3  = energyProfile;   % duplicate energyProfile in order to alter this vector
s3m = zeros(size(s3)); % logical for maxima
for is=2:(numel(s3)-1)
	if (s3(is-1)<s3(is)) && (s3(is)>s3(is+1)) %local maxima
		s3m(is)=is;
	end
end
s3mpts = s3m(s3m>0);
xy = s3mpts;
y  = s3(logical(s3m));
xx = 1:numel(x);
splineEnergyProfile = spline(xy,y,xx);

% Continue drawing if marked
if DRAWWAVEPACKETS
	area(splineEnergyProfile) %draw engery profile as shaded area. perty.
end

% Now apply the greedy algorithm on the splined data    
pk = zeros(size(splineEnergyProfile));
tmpSplineEnergyProfile=splineEnergyProfile; %duplicate energy profile to use within inner loop        
chckdpts = ones(size(x));%vector for checked points

% MAIN ALGORITHM LOOP --------------------------------------------------- %
ipeaks=0;  %set counter to label peaks
while true %run until no more significant peaks

	% make sure there is still data to be checked
	if all(chckdpts==0)
		break;
	end

	% iterate peak counter; used to label (number) found peaks    
	ipeaks=ipeaks+1;

	% kills peaks that have already been found. on first pass, this doesn't do anything.      
	tmpSplineEnergyProfile(pk~=0)=0;

	% find biggest (current) peak
	[maxPeakVal,maxPeakLoc]= max(tmpSplineEnergyProfile);
	%chckdpts(maxPeakLoc)=0;

	% Exit condition. Once the current peak falls below a threshold, which 
	% you can think of as a background or noise level, then the algorithm ends.
	% Sensitivity of algorithm depends on the rule for the exit condition.
	%
	% TODO: should we use the splineEnergyProfile or the energyprofile for
	% computing RMS? 
	switch EXITCONDITION
		case 1
			snr = THRES*rms(splineEnergyProfile(pk<1)); %pk==0
			%             snr = THRES*rms(energyprofile(pk==0));
			if maxPeakVal< snr
				ipeaks=ipeaks-1;
				break;
			end

		case 2 %this keeps the previous packet in for rms calculation
			snr = THRES*rms(splineEnergyProfile(pk==0 | pk==(ipeaks-1) ));
			%             snr = THRES*rms(energyprofile(pk==0 | pk==(ipeaks-1) ));
			if maxPeakVal< snr
				ipeaks=ipeaks-1;
				break;
			end
		otherwise
			error('Need to specify an exit condition.')
	end

	%march border left
	left = maxPeakLoc-STEPSIZE;
	while left>2 && left<(numel(x)-LAG) 
		if all( splineEnergyProfile(left-1) > splineEnergyProfile(left:(left+LAG)) )

			if splineEnergyProfile(left-1) > splineEnergyProfile(maxPeakLoc)/2 		
			%still greater than half-energy. likely still part of same oscillation
			else
			    break;
			end
		end        
		%collision detection
		if pk(left)>0
			left=left+2;            
			break
		end        
		%iterate
		left=left-1;        
	end
	%check width condition
	if (maxPeakLoc-left)> MAXWIDTH/2
		left = maxPeakLoc - MAXWIDTH/2;
	end
	%check boundary condition
	if left<=1
		left=1;
	end


	%march border right
	right=maxPeakLoc+STEPSIZE;
	while (right-LAG>0) && right<(numel(x)-LAG) 
		if all( splineEnergyProfile(right+1) > splineEnergyProfile( (right-LAG):right ) )

		%still greater than half-energy. likely still part of same oscillation
			if splineEnergyProfile(right+1) < splineEnergyProfile(maxPeakLoc)/2 
			    break;
			end
		end
		%collision detection
		if pk(right)>0
			right=right-2;            
			break
		end        
		%iterate
		right=right+1;
	end
	%check width condition
	if (right-maxPeakLoc) > MAXWIDTH/2
		right = maxPeakLoc+ MAXWIDTH/2;
	end
	%check boundary condition
	if right>=numel(x)
		right=numel(x);
	end


	%mark found peak with a peak number. 
	%can use logical(pk) to get back peaks
	pk(left:right)=ipeaks;

	%check minimal width condition
	if (right-left) < MINWIDTH
		ipeaks=ipeaks-1;
		pk(left:right)=-1;
		%break;
	else

		%make sure we don't need to increase wavepackets table size
		if ipeaks>size(wavepackets,1)
			wavepackets = [wavepackets; null_wavepackets]; %#ok<AGROW>
		end

		%compute wave stats, save to wavepacket table
		wavepackets.peak(ipeaks)   = maxPeakLoc;
		wavepackets.start(ipeaks)  = left;
		wavepackets.stop(ipeaks)   = right;
		wavepackets.length(ipeaks) = right-left+1;
		wavepackets.snr(ipeaks)    = maxPeakVal/snr;

		%compute hilbert frequency, etc 
		hy     = x(left:right);
		tvec   = (1:numel(hy))./SF;
		hy     = hilbert(hy);
		uphase = unwrap(angle(hy));
		ifreq  = diff(uphase)/diff(tvec)/(2*pi);
		iamp   = abs(hy);

		[~,mf]=max(coefs(:,maxPeakLoc));

		wavepackets.wvfreq(ipeaks) = scal2frq(scales(mf),'morl',1/SF);
		wavepackets.wvpow(ipeaks)  = mean(mean(coefs(:,left:right)));

		wavepackets.ifreq(ipeaks)  = ifreq;
		wavepackets.iamp(ipeaks)   = mean(iamp);

		%store wave
		wavepackets.wave{ipeaks}   = x(left:right);

		wavepackets.record(ipeaks) = ipeaks;


	end

end %end 

if DRAWWAVEPACKETS
	stairs(max(splineEnergyProfile)*logical(pk>0),'k')
	xlim([0 numel(x)])
end

if (ipeaks+1)< size(wavepackets,1)
    wavepackets((ipeaks+1):end,:)=[];
end



