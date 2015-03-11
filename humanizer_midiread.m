function [notes, C0, notes_start, notes_end] = humanizer_midiread( textfile, row_start )
	% Read midi file and generate matrix "notes" in a readable format with beat start and end times.
	% Written by Holger Hennig (2013)
	% The output matrix "notes" is compatible with Schutte's readmidi (but the latter isn't exact)
	
	% notes columns:
	% 1: track number % set to zero
	% 2: channel number
	% 3: pitch
	% 4: velocity % set to zero
	% 5: start time (seconds)
	% 6: end time (seconds)
	
    if ~exist('row_start','var') || isempty(row_start), row_start=1; end
	
	%% read data in textfile into cell C
    % Open the file
	fid = fopen(textfile,'r');                        
	
    % read data as strings
    data = textscan(fid,'%s %s %s %s %s','CollectOutput',true);
	fclose(fid);
    
    % write data into cell C
	C0 = data{1};
	C = C0(row_start:end,:);
    
	%% prepare translating cell to notes matrix
	time_cell = C(:,1);
	onoff_cell = C(:,2);
	channel_cell = C(:,3);
	pitch_cell = C(:,4);
	
	%convert cell string to number array
	time = cellfun(@str2num,time_cell);
	
	%convert cell string to string array
	onoff = char(onoff_cell);
	
	%convert channel to number (e.g. 'ch=5' to 5)
	ch_tmp =  regexp(channel_cell,'\d*','match');
	ch_tmp = cat(2,ch_tmp{:});
	ch = cellfun(@str2num,ch_tmp)';
	
	%convert pitch to number (e.g. 'n=11' to 11)
	p_tmp =  regexp(pitch_cell,'\d*','match');
	p_tmp = cat(2,p_tmp{:});
	pitch_unix = cellfun(@str2num,p_tmp)';
	
	len_time = length(time);
	notes_start = zeros(round(len_time/2),1);
	notes_end = zeros(round(len_time/2),1);
	
	%% Find beats_start and beats_end
	for k=1:len_time
		str = onoff(k,:);

        if strcmpi(str, 'On  ')
			%note_start found
			notes_start(k)=k;
			k2 = k;
			cand = 0;
			
			%check candidate for note_end
			while cand==0
				k2 = k2+1;
				if strcmpi(onoff(k2,:), 'Off ')
					
					%does note_end have the same pitch and channel?
					if pitch_unix(k2)==pitch_unix(k) && ch(k2)==ch(k)
						cand = 1; % note_end found!
						notes_end(k) = k2;
					end
				end
			end
		end
	end
	
	% remove zeros from notes_start and notes_end
	notes_start(notes_start==0) = [];
	notes_end(notes_end==0) = [];

	%% convert time to realtime
	% File header: Mfile <format> <ntrks> <division>
	% TPB = <division> is the time resolution in clicks per quarter note
	TPB = 192; %ticks per beat
	tempo = 512820; % in 10^-6 seconds. tempo is 0.51s per beat
	BPM_u = 60 / (tempo*1e-6); % which is 117 BPM, yay!
	realtime = 60 / (TPB * BPM_u)*time; %realtime in seconds
	%realtime = tempo*1e-6 / TPB; % one unit in midi is 0.0027 s = 2.7ms

	notes_start_tmp = realtime(notes_start);
	notes_end_tmp = realtime(notes_end);
	num_beats = length(notes_start_tmp);
	
	%% generate matrix notes
	notes = zeros(num_beats,6);
	notes(:,2) = ch(notes_start);
	notes(:,3) = round(pitch_unix(notes_start));
	notes(:,5) = notes_start_tmp;
	notes(:,6) = notes_end_tmp;
	
	
end

