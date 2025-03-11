%% Wrapper function to loop over datasets and compute timescales

clear

% Determine the computer name to set the correct file paths
[~, comp_name] = system('hostname');

% Define paths based on workstation name
if strcmp(comp_name(1:5), 'thebr') || strcmp(comp_name(1:5), 'pinky') % Workstations: Pinky + Brain
    path2go = '/home/fred/TIMESCALES/data/';
    ExistingPool = gcp('NoCreate'); % Check if a parallel pool exists
    
    if isempty(ExistingPool)
        parpool('local', 36) % Start a parallel pool with 36 workers
    end
elseif strcmp(comp_name(1:5), 'DESKT') % Desktop machine
    path2go = 'S:\TIMESCALES\data\';
else
    error('Define path where the datasets/processed folders are') % Error if path is not specified
end

% Parameters for autocorrelation fitting
nTau = 1; % Set to 1 for all timestamps, otherwise use 50
nSpk = 1000; % Number of spikes used if nTau > 1
fol = 'revision'; % Folder name, can be 'draft', 'final', 'revision' 

% Parameters for spike count analysis
param_ac.tot_time = 1; % Time window before trial start (1 sec)
param_ac.bins = 0.05; % Bin size in seconds
param_ac.bins_start = 0:param_ac.bins:param_ac.tot_time - param_ac.bins;
param_ac.bins_end = param_ac.bins:param_ac.bins:param_ac.tot_time;
param_ac.evt = 't_on'; % Event timestamp

% Load dataset file list
list = dir([path2go 'datasets_' fol '/*.mat']);

% Loop through datasets
for d = 1:length(list)
    % Check if the processed file already exists
    if exist([path2go 'processed_' fol '/ISI_Timescales_' list(d).name(1:end-4) '.mat'], 'file') == 0
        
        % Load dataset
        load([path2go 'datasets_' fol '/' list(d).name])
        disp(list(d).name)
        unt = find(list(d).name == '_');

        % Determine which analyses can be run based on available data
        if exist('spikes_rest', 'var')
            run_rest = true;
        else
            run_rest = false;
            spikes_rest = cell(size(spikes));
        end
        if exist('spikes_task', 'var')
            run_task = true;
        else
            run_task = false;
            spikes_task = cell(size(spikes));
        end
        if isfield(cell_info{1}, param_ac.evt)
            run_ac = true;
        else
            run_ac = false;
        end

        % Initialize output structures
        out = {}; % ISI method with all spikes
        out_rest = {}; % ISI method only rest periods
        out_task = {}; % ISI method only task periods
        out_ac = {}; % Autocorrelation method (a la Murray et al.)

        % Loop through neurons
        for p = 1:nTau
            parfor n = 1:length(spikes)
                disp(['Processing area ' num2str(d) '/' num2str(length(list)) ' - neuron ' num2str(n) '/' num2str(length(spikes))])
                
                if nTau > 1 % If subsampling is enabled
                    if length(spikes{n}) > nSpk
                        % Select a random subset of spike times
                        rdperm = randperm(length(spikes{n}));
                        [out{n, p}] = ISI_timescales(sort(spikes{n}(rdperm(1:nSpk))));
                    else
                        [out{n, p}] = ISI_timescales([]);
                    end
                else
                    % Compute timescales for full spike data
                    [out{n, p}] = ISI_timescales(spikes{n});
                    if run_rest
                        [out_rest{n, p}] = ISI_timescales(spikes_rest{n});
                    end
                    if run_task
                        [out_task{n, p}] = ISI_timescales(spikes_task{n});
                    end
                    if run_ac
                        tr_start = cell_info{n}.(param_ac.evt); % Extract event times
                        
                        % Construct spike count matrix based on dataset type
                        if ~ismember({list(d).name(1:unt-1)}, 'minxha')
                            spkcount = [];
                            for t = 1:length(tr_start)
                                for b = 1:length(param_ac.bins_start)
                                    spkcount(t, b) = sum(spikes{n} >= (tr_start(t) - param_ac.tot_time + param_ac.bins_start(b)) & spikes{n} <= (tr_start(t) - param_ac.tot_time + param_ac.bins_end(b)));
                                end
                            end
                        else
                            % For human dataset, use post-stimulus period (1 sec after stimulus)
                            spkcount = [];
                            for t = 1:length(tr_start)
                                for b = 1:length(param_ac.bins_start)
                                    spkcount(t, b) = sum(spikes{n} >= (tr_start(t) + param_ac.bins_start(b)) & spikes{n} <= (tr_start(t) + param_ac.bins_end(b)));
                                end
                            end
                        end

                        % Compute autocorrelation timescales
                        [out_ac{n, p}] = AC_timescales(spkcount, param_ac.bins * 1000);
                    end
                end
            end
        end
        
        % Ensure proper output formatting
        if nTau == 1
            out = out';
            out_rest = out_rest';
            out_task = out_task';
            out_ac = out_ac';
        end

        % Save processed output
        save([path2go 'processed_' fol '/ISI_Timescales_' list(d).name(1:end-4) '.mat'], 'out', 'out_rest', 'out_task', 'out_ac', 'spikes', 'spikes_rest', 'spikes_task')
        clear out out_rest out_task out_ac spikes spikes_rest spikes_task
    end
end

%% Post processing - create a csv file with all the required infos!

clearvars -except path2go fol

savenames = {'full' 'rest' 'task' 'ac'};
mat2take = {'out' 'out_rest' 'out_task' 'out_ac'};
spk2take = {'spikes' 'spikes_rest' 'spikes_task' 'spikes'};

list = dir([path2go 'processed_' fol '/ISI_Ti*.mat']);

for mm = 1 : length(mat2take)
    data = [];
    for d = 1 : length(list)
        load([path2go 'processed_' fol '/' list(d).name])
    
        %- find dataset name and area
        unde = find(list(d).name=='_');
        PIname = {list(d).name(unde(2)+1:unde(3)-1)};
        Area = {list(d).name(unde(3)+1:end-4)};
        
        disp(['Processing: ' PIname{1} ' ' Area{1} ' - ' savenames{mm} ' ....'])
        %- if FROOT, load the FR from the initial matrix (cell_info)
        if strcmp(PIname,'minxha') | strcmp(PIname,'stoll') 
            load([path2go 'datasets_' fol '/' list(d).name(unde(2)+1:end)],'cell_info');
        end
            
        eval(['curr_mat = ' mat2take{mm} ';']); % pick the right mat
        eval(['curr_spk = ' spk2take{mm} ';']); % pick the right mat
        
        if isempty(curr_mat)
            curr_mat{1}.f0=[];
            curr_mat = repmat(curr_mat,1,length(spikes));
        end
        timescale_data=[];
        for i = 1 : length(spikes)
            if ~isempty(curr_mat{i}.f0)
                if strcmp(PIname,'minxha') 
                    FR = cell_info{i}.fr;
                else
                    FR = length(curr_spk{i})/(curr_spk{i}(end)-curr_spk{i}(1));
                end
                %- FR calculated only for main analysis
                if mm ~= 1
                    FR = 0;
                end
                if strcmp(PIname,'stoll') 
                    timescale_data = [timescale_data ; length(curr_spk{i}) curr_mat{i}.lat curr_mat{i}.f0.time_constant curr_mat{i}.f0.A curr_mat{i}.f0.B curr_mat{i}.sse curr_mat{i}.r2 FR i cell_info{i}.loc strcmp(cell_info{i}.session(1),'X')+1];
                else
                    timescale_data = [timescale_data ; length(curr_spk{i}) curr_mat{i}.lat curr_mat{i}.f0.time_constant curr_mat{i}.f0.A curr_mat{i}.f0.B curr_mat{i}.sse curr_mat{i}.r2 FR i 0 0 0 0];
                end
                    
                else
                timescale_data = [timescale_data ; 0 0 0 0 0 0 0 0 0 0 0 0 0];
            end
        end
        %- criteria to keep units:
        keep = timescale_data(:,1)>250 & (timescale_data(:,3)>0 & timescale_data(:,3)<1000) & timescale_data(:,4)>0 & timescale_data(:,7)>=0.5;
        timescale_data(:,end+1) = keep;
        % timescale_data(~keep,:)=[];
    
        PI = repmat(PIname,length(timescale_data(:,1)),1);
        Ar = repmat(Area,length(timescale_data(:,1)),1);
        
        if strcmp(PIname,'stein')
            species = {'mouse'};
        elseif strcmp(PIname,'meg')  | strcmp(PIname,'froot')  | strcmp(PIname,'wirth')  | strcmp(PIname,'stoll')  | strcmp(PIname,'fontanier') 
            species = {'monkey'};
        elseif strcmp(PIname,'kepecs') | strcmp(PIname,'buzsaki') | strcmp(PIname,'feierstein') | strcmp(PIname,'peyrache') | strcmp(PIname,'lemerre')
            species = {'rat'};
        elseif strcmp(PIname,'minxha') | strcmp(PIname,'chandravadia')  | strcmp(PIname,'faraut') 
            species = {'human'};
        end
        Sp = repmat(species,length(timescale_data(:,1)),1);
    
        clear T
    
        %- put that in a table
        T = array2table(timescale_data(:,[9 3 2 4:7 1 8 10 11 12 13 14])); % different order but match the variable names below
        T = [T PI Ar Sp];
        T.Properties.VariableNames = {'unitID' 'tau' 'lat' 'A' 'B' 'sse' 'r2' 'nbSpk' 'FR' 'ML' 'VD' 'AP' 'animID' 'keep' 'name' 'area' 'species' };
        
        data = [data ; T];
        % [sum(data.keep)/height(data)]
        % [mean(data.tau(data.keep==1)) std(data.tau(data.keep==1)) ]
    end
    disp([num2str((sum(data.keep)/height(data))*100) '% of neurons processed'])
    save([path2go 'ISI_timescales_' savenames{mm} '.mat'],'data')

end




