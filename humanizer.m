function [ notes_humanized, seed, sort_ind ] = humanizer( input_filename, sigma, alpha, type, seed, verbose, sigma2, alpha2, W)
	% Couples two audio tracks ("as if two musicians played together") or
	% Humanizes an audio track ("as if played by a musician") in a MIDI file.
	% Written by Holger Hennig (2014)
	% Contact/feedback: holger.hen (at) gmail.com
    % Based on the article: H. Hennig, PNAS 111,12974 (2014)
    % Free article download: www.nld.ds.mpg.de/~holgerh/download
    % License: This work is licensed under a Creative Commons 
    %   Attribution-NonCommercial-ShareAlike 4.0 International License.
    
    % notes columns:
	% 1: track number1
	% 2: channel number
	% 3: pitch
	% 4: velocity
	% 5: start time (seconds)
	% 6: end time (seconds)
	% 7: message number of note_on, 8: message number of note_off
	
	% channel # and pitch define the "sound". Examples:
	% drum: channel 9
	% pitch 36: bass drum, 40: snare, 42: hi-hat
	% bass: channel 1
	% pitch 42: F, 37: C#, 40: E, 35: F
    
	%% check input
	
	% parameters for coupling / humanizing
	% alpha: strength of correlations between errors ("How strongly an error depends on previous errors")
	% sigma: standard deviation of introduced errors in milliseconds ("How large is an error in average. Deviations will be cut off at 2*sigma (to exclude rare extreme events)")
	if ~exist('sigma','var') || isempty(sigma), sigma=10; end
	if ~exist('alpha','var') || isempty(alpha), alpha=1; end
	if ~exist('type','var') || isempty(type), type=1; end
	if ~exist('sigma2','var') || isempty(sigma2), sigma2=sigma; end
	if ~exist('alpha2','var') || isempty(alpha2), alpha2=alpha; end
	if ~exist('W','var') || isempty(W), W=0.5; end %maximum coupling
	if ~exist('seed','var') || isempty(seed), seed = floor(1e6*rand); end
	if ~exist('verbose','var') || isempty(verbose), verbose=1; end %be talkative by default and write misc files and figures to disk
	
	% song prefix for output filenames, e.g. 'mysong_'	
	prefix = 'mysong_';
	
    % extract extension from input_filename
    [input_path,input_name,extension] = fileparts(input_filename);
    
    % convert mid to txt
    switch extension
        case '.mid'
        
            % convert midi to textfile: mf2t song.mid song.txt
            output_filename = [input_path '/' input_name '.txt'];
            midi2txt = ['mf2t ' input_filename ' ' output_filename];
            unix(midi2txt);
            fprintf('Converting MIDI to ascii txt file. Wrote file %s to disk\n', output_filename);
            error('Humanizer:FileFormat', 'Currently only txt files are supported. Your MIDI file was converted to a *.txt file. See documentation how to proceed with txt file.')
    end
    
    % only work with txt file from here on
    if ~strcmpi(extension, '.txt')
        error('Humanizer:InputCheck', 'Input file extension needs to be txt. Example: mysong.txt')
    end
    
    %notes_start is an index pointing to starting notes in C, same length as beats_start (which are times in ms)
	[notes, C, notes_start, notes_end] = humanizer_midiread(input_filename);
    
	%% find tatum ("smallest" unit, smallest periodic interbeat interval)
	offset = notes(1,5);
    %offset = notes(2,5)
	res = 10; % increase temporal resolution by res (e.g., res = 10). Useful to imprint small (submillisecond) deviations
	
	% Billie Jean has 117 BPM, the exact tatum is
	tatum = 60/117/2 * 10/res;
	
	% Prepare for modifying the tatums. Modify start (col 5) and end times (col 6)
	beats_start = notes(:,5);
	beats_end = notes(:,6);
	channels = notes(:,2);
	
    
	%% identify instruments
	
	bass = channels==2;
	bass(1)=0; %don't humanize the first inaudible bass sound (simply marks the beginning of the song)
	drum = channels==10;
	keyboard = channels==3;
	percussion = channels==11;
	
	% merge drum with other instruments (keyboard, percussion) into one track for coupling
	drum = drum | keyboard;
	drum = drum | percussion;
	
    
	%% find beats on grid
	
	%beats_start in tatum grid units, e.g. [0,1,1,2,3,3.1,4,4,4,5.8,6,6,8,...]
	%beats_start_grid = (beats_start - offset) / tatum;
	beats_start_grid = round((beats_start - offset) / tatum);
	% for humanizing: use only the beats_on_grid
	beats_on_grid = abs(beats_start_grid - round(beats_start_grid)) < 3e-3;
	
	beats_start_humanized = beats_start;
	beats_end_humanized = beats_end;
		
	len_humanize = round(beats_start_grid(end))+1; %length(beats_start(beats_on_grid));	
	%len_humanize = round(beats_start_grid(end))+10; 
    
	%% get ready: Humanize!
	switch type
		
        % exact version
        case 0 
			
            fprintf('Generating exact version.\n');
			typename = 'exact';
			x = zeros(len_humanize,1);
			y = x;
		
        % "Group humanize": couple instruments using two-component ARFIMA process    
        case 1 
			
            fprintf('Group Humanizer: Coupling two instruments (2-component ARFIMA).\n');
			typename = 'couple';
			
            % Coupling: generate two-component ARFIMA process for interbeat intervals
        	[x, y] = arfima_model(alpha-0.5, alpha2-0.5, W, len_humanize, seed);
            x = x'; 
            y = y';
	
            % renormalize. ARFIMA outputs x_couple where standard deviation of x_couple is only 3.5 for sigma=5
            x = x/std(x) * sigma/1000;
        	y = y/std(y) * sigma2/1000;
            x = x - mean(x);
            y = y - mean(y);
		
        % humanize, SAME error on each tantum for all instruments
		case 2 
			
            fprintf('Humanizing two instruments with the same deviations.\n');
			typename = 'same';
		
            % Humanizing: generate 1/f noise for interbeat intervals
            x = noise(alpha, len_humanize, seed)' * sigma/1000;
			y = x * sigma2/1000; %because correlations in x and y shall be identical here (but stddev can differ)!
			
        % humanize instruments separately
		case 3 
			
            fprintf('Humanizing two instruments separately.\n');
			typename = 'sep';
            x = noise(alpha, len_humanize, seed)' * sigma/1000;
            y = noise(alpha, len_humanize, seed+1)' * sigma2/1000;

		
        % "Group humanize": couple instruments using MICS model		
		case 4 
            
            fprintf('Group Humanizer: Coupling two instruments based on MICS model\n');
			typename = 'mics';	
            [x,y] = mics_model(alpha, sigma, len_humanize, W, [], 0, seed, alpha2, sigma2);
            x = x'; 
            y = y';
            
            % renormalize to get the standard deviation right
            x = x/std(x) * sigma/1000;
        	y = y/std(y) * sigma2/1000;
            x = x - mean(x);
            y = y - mean(y);            
    end

    % renormalize deviations (except for exact version)
    if type ~= 0 
        
        % don't modify the first beat!
        x(1) = 0;
        y(1) = 0;
	
        % reduce rare large changes in the deviations (outside 2*sigma) by factor of 2
        xcounter = 0;
        ycounter = 0;
        for kk = 2:length(x)
            if abs(x(kk) - x(kk-1)) > 2*sigma/1000
                x(kk) = ( x(kk) + x(kk-1) ) /2;
                xcounter = xcounter + 1;
            end
        end
	
        for kk = 2:length(y)
            if abs(y(kk) - y(kk-1)) > 2*sigma2/1000
                y(kk) = ( y(kk) + y(kk-1) ) /2;
                    ycounter = ycounter + 1;
            end
        end
    
        x = x / std(x) * sigma/1000;
        y = y / std(y) * sigma2/1000;
    end
  
    
	%% calculate deviations from interbeat intervals
	% interbeat intervals: I = x + tatum;
	% deviations d_n = I_1 + I_2 + ... + I_n - n * tatum
	% therefore: d_n = cumsum(I - tatum) = cumsum(x);
	% BUT: need to work on a grid, otherwise drum and bass are not playing together at all!
	% therefore consider x,y directly as the deviations (as if a metronome was present). 
	
	dev_x = x; %cumsum(x);
	dev_y = y; %cumsum(y);

    
	%% Couple / Humanize!
	beats_start_humanized(beats_on_grid & drum) = beats_start(beats_on_grid & drum) + dev_x(1+round(beats_start_grid(beats_on_grid & drum)));
	beats_end_humanized(beats_on_grid & drum)   = beats_end(beats_on_grid & drum)   + dev_x(1+round(beats_start_grid(beats_on_grid & drum)));
	beats_start_humanized(beats_on_grid & bass) = beats_start(beats_on_grid & bass) + dev_y(1+round(beats_start_grid(beats_on_grid & bass)));
	beats_end_humanized(beats_on_grid & bass)   = beats_end(beats_on_grid & bass)   + dev_y(1+round(beats_start_grid(beats_on_grid & bass)));
	
	dev_drum_bass = zeros(len_humanize,2);
	dev_drum_bass(:,1) = x;
	dev_drum_bass(:,2) = y;
	
    
	%% consistency check
	% make sure all beats are at times >=0. Too large devs can result in negative times for the very first beat.
	
	if min(beats_start_humanized)<0
		beats_start_humanized(beats_start_humanized<0)=0;
		beats_end_humanized(beats_end_humanized<0)=0;
		fprintf('Warning, there are beats at negative times. Shifted the beats to time zero.\n')
	end
	
	notes_humanized = notes;
	notes_humanized(:,5) = beats_start_humanized;
	notes_humanized(:,6) = beats_end_humanized;
	
    
	%% Graphical output
	
	filename_plot = [prefix typename '_s' num2str(sigma) '_a' num2str(alpha) '_dev'];
		
	drum_ind = 1 + unique(round(beats_start_grid(beats_on_grid & drum)));
	drumdev_used = zeros(length(x),1);
	drumdev_used(drum_ind) = 1;
	drumdev_used = logical(drumdev_used);
	%
	bass_ind = 1 + unique(round(beats_start_grid(beats_on_grid & bass)));
	bassdev_used = zeros(length(y),1);
	bassdev_used(bass_ind) = 1;
	bassdev_used = logical(bassdev_used);
	
    % plot time series of deviations
	figure(1)
	xaxis_d = 1:length(x);
	xaxis_b = 1:length(y);	
	h = plot(x*1000,'k-','linewidth',1.5); %drum
	hold on
	plot(y*1000 + 5*sigma,'b-','linewidth',1.5) %bass
	legend('drum dev.','bass dev.')	
	xlim([1, length(x)]); ylim([min(x*1000), max(y*1000) + 7.5*sigma])
	xlabel('Beat no.','Fontsize',24); ylabel('Deviations /ms ','Fontsize',24)
	set(gca,'FontSize',20,'Linewidth',1)
	
    % optional: write time series plot to disk
	if verbose>0
		saveas(h, [filename_plot '.fig'])
		saveas(h, [filename_plot '.eps'],'epsc2')
    end
	
    % mark used deviations in time series plot
	plot(xaxis_d(drumdev_used), x(drumdev_used)*1000,'ro')
	plot(xaxis_b(bassdev_used), y(bassdev_used)*1000 + 5*sigma,'ro')
	ylim([min(x*1000), max(y*1000) + 9*sigma])
	legend('drum dev.','bass dev.','dev. used')
	set(h,'MarkerSize',5)
	
    % optional: write time series plot (including used deviations) to disk
	if verbose>0
		saveas(h, [filename_plot '_used.fig'])
		saveas(h, [filename_plot '_used.eps'],'epsc2')
	end
	hold off
		
    
	%% Export data: write time series of deviations to disk	
	
	if verbose>0
        filename_dev = [prefix typename '_s' num2str(sigma) '_a' num2str(alpha) '_dev_seed' num2str(seed) '.txt'];
        save(filename_dev,'dev_drum_bass','-ascii','-tabs')
    end
	
	% add deviations to cell C
	% extract times
	time_cell = C(:,1); % time in midiunits
	time = cellfun(@str2num,time_cell); %double array
	time_h = time; 
	
	% MIDI File header: Mfile <format> <ntrks> <division>
	% <division> = TPB = ticks per beat. TPB is the time resolution and is identical to pulses per quarter (PPQ) note
	%192 is typical value, but higher values are possible (e.g., 1920)
	
	TPB = 192; %ticks per beat
	tempo = 512820; % in 10^-6 seconds (per beat)
	BPM_u = 60 / (tempo*1e-6); % which is 117 BPM, yay!
	% realtime = 60 / (TPB * BPM_u)*time; %realtime in seconds
	
	%convert realtime (s) to miditime. time(midi) = time(s) * (TPB * BPM_u) / 60 ;
	sec2miditime = (TPB * BPM_u) / 60; %374
	time_h(notes_start) = notes_humanized(:,5) * sec2miditime;
	time_h(notes_end) = notes_humanized(:,6) * sec2miditime;
		
	%increase TPB (time resolution)... TPB' = TPB * res in header txt file.
	%res = 10;  
	time_h_str = num2str(round(time_h*res));	
	
	% write modified times in cell C_h(:,1)
    C_tmp = C;
    
    % write modified onsets in cell
	if abs(type)<1e-6
        % exact version for type=0
        C_tmp(:,1) = cellstr(num2str(round(time*res))); % use this line for not adding temporal modifications!
    else
        C_tmp(:,1) = cellstr(time_h_str); 
    end
	
    %sort the cell based on times. MIDI requires time-ordered entries (ascending order)
	[C_h, sort_ind] = sortrows(C_tmp,1); 
    C_h(:,1) = strtrim( C_h(:,1) ); %remove whitespace (warning: remove whitespace only after sorting, otherwise sorting won't work)
	
	
	%% select instruments for export
		
	% convert channel to number (e.g. 'ch=10' to 10)
	ch = C_h(:,3);
	ch =  regexp(ch,'\d*','match');
	ch = cat(2,ch{:});
	ch = cellfun(@str2num,ch)';
	
	ch_bass = 2; % bass channel is 2, but depending on MIDI sequence may need to be adjusted (e.g. ch_bass=1)
	ch_drum = 10;
	ch_keyboard = 3;
	ch_extra = 11; % misc sounds

	% select channels
	instr = (ch==ch_drum | ch==ch_bass | ch==ch_keyboard | ch==ch_extra);
	C_h = C_h(instr,:);
	
	%% optional: choose length of silence at the beginning of the song
	stime_cell = C_h(:,1);
	stime = cellfun(@str2num,stime_cell); %double array
	silence = 7200-750; %7200 is current first note on. 750 is wanted first note on (750*midiunit=200 ms)
	nosilence = stime - silence;
	
	% write time into cell 
	time_h_str = num2str(nosilence);
	C_h(:,1) = cellstr(time_h_str); % adjust silence
	C_h(:,1) = strtrim( C_h(:,1) ); % remove whitespace
	    
	%% Export MIDI using t2mf (UNIX, see http://code.google.com/p/midi2text/)
	
	if abs(type) < 1e-6
		filename = [prefix 'exact'];
		f_core = [prefix 'exact_core.txt'];	
	else
		filename = [prefix typename '_s' num2str(sigma) '_a' num2str(alpha)];
		f_core = [prefix typename '_s' num2str(sigma) '_a' num2str(alpha) '_core.txt'];
    end
	
    % write cell C_h as ascii textfile f_core to disk
	fprintf('Writing midi file %s.mid\n',filename);
	dlmcell(f_core, C_h, ' ');
        
	% header + core + footer = txtfile
	core2txt = ['cat org/org_bj_header.txt ' f_core ' org/org_bj_footer.txt > ' filename '.txt'];
	unix(core2txt);
	
	% delete core file
	deletecore = ['rm ' f_core];
	unix(deletecore);
	
	% txt -> midi
    txt2midi = ['t2mf ' filename '.txt ' filename '.mid'];
	%txt2midi = ['/Users/holgerh/bin/t2mf ' filename '.txt ' filename '.mid'];
	unix(txt2midi);
	
	% clean up
	if verbose<=0
		unix(['rm ' filename '.txt ']);
    end	
	
	fprintf('Done\n\n')
end

