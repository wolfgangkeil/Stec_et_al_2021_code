function worm = read_single_worm_molting_data(wormfile)
%
% function worm = read_single_worm_molting_data(wormfile)
%
% This function reads a worm lineage textfile and returns a structure
% called worm, with molting times etc.,
% in the text file, NC means 'not captured',
% happened
% 
% a typical worm structure might look like this
%
% worm = 
% 
%                   wormfile: '09_Jun_2015_2.txt'
%                   genotype: 'N2'
%             imaging_method: 'automated'
%                temperature: '2015-06-09'
%                 synch_date: '2015-06-09'
%                 synch_time: '15:15:00'
%           growth_condition: 'liquid'
%                 mount_date: '2015-06-10'
%                 mount_time: '16:40:00'
%     end_of_experiment_date: '2015-06-11'
%     end_of_experiment_time: '12:00:00'
%          l1_lethargus_date: 'NC'
%          l1_lethargus_time: 'NC'
%               l1_molt_date: 'NC'
%               l1_molt_time: '00:00:00'
%          l2_lethargus_date: 'NC'
%          l2_lethargus_time: 'NC'
%               l2_molt_date: 'NC'
%               l2_molt_time: '01:20:00'
%          l3_lethargus_time: 'NC'
%               l3_molt_date: 'NC'
%               l3_molt_time: 'NC'
%               l4_molt_date: 'NC'
%               l4_molt_time: 'NC'
%
%
% wolfgang.keil@curie.fr, 2020

    worm = [];
    
    fid = fopen(wormfile);
    
    if fid > 0
        closeTheFile = onCleanup(@() fclose(fid));

        % Skip the header 
        tmp = textscan(fid, '%s %*[^\n]',1);
        while ~strcmp(tmp{1}, '##')
            tmp = textscan(fid, '%s %*[^\n]',1);
        end

        worm.wormfile = wormfile;
        %%%%%%%%%  THIS READS EXPERIMENT INFORMATIONS        
        tmp  = {''};
        while ~strcmp(tmp{1}, '##');            
            tmp = textscan(fid, '%s',1);
            if strcmpi(tmp{1}, 'genotype')
                t  = textscan(fid, '%s %*[^\n]',1);
                worm.genotype =char(t{1});
            elseif strcmpi(tmp{1}, 'temperature')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.temperature = char(t{1});
            elseif strcmpi(tmp{1}, 'synch_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.synch_date = char(t{1});
                worm.temperature = char(t{1});
            elseif strcmpi(tmp{1}, 'synch_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.synch_time =char(t{1});
            elseif strcmpi(tmp{1}, 'imaging_method')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.imaging_method = char(t{1});        
            elseif strcmpi(tmp{1}, 'growth_condition')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.growth_condition = char(t{1});            
            elseif strcmpi(tmp{1}, 'mount_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.mount_date =char(t{1});
            elseif strcmpi(tmp{1}, 'mount_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.mount_time = char(t{1});
            elseif strcmpi(tmp{1}, 'end_of_experiment_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.end_of_experiment_date = char(t{1});
            elseif strcmpi(tmp{1}, 'end_of_experiment_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.end_of_experiment_time = char(t{1});
            elseif strcmpi(tmp{1}, 'experiment_folder')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.experiment_folder = char(t{1});

                if ~ismac
                    % Code to run on Linux plaform
                    worm.experiment_folder = strrep(worm.experiment_folder, '/Volumes', '/media/wolfgang');
                end
            elseif strcmpi(tmp{1}, 'fluorescence_peaks')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.fluorescence_peaks = char(t{1});
            elseif strcmpi(tmp{1}, 'experiment_type')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.experiment_type = char(t{1});
            elseif strcmpi(tmp{1}, 'worm_development_arrested')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.worm_development_arrested = char(t{1});
            end

        end

        %%%%%%%%%  THIS SECTION THEN READS LETHARGUS INFORMATIONS
        tmp  = {''};        
        while ~strcmp(tmp{1}, '##');            
            tmp = textscan(fid, '%s',1);

            %%% L1/L2 molts        
            if strcmpi(tmp{1}, 'l1_lethargus_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_lethargus_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l1_lethargus_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_lethargus_time = char(t{1});
            elseif strcmpi(tmp{1}, 'l1_mouth_cap_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_mouth_cap_date = char(t{1});            
            elseif strcmpi(tmp{1}, 'l1_mouth_cap_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_mouth_cap_time = char(t{1});            
            elseif strcmpi(tmp{1}, 'l1_molt_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_molt_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l1_molt_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l1_molt_time = char(t{1});

                %%% L2/L3 molts
            elseif strcmpi(tmp{1}, 'l2_lethargus_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_lethargus_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l2_lethargus_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_lethargus_time = char(t{1});
            elseif strcmpi(tmp{1}, 'l2_mouth_cap_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_mouth_cap_date = char(t{1});            
            elseif strcmpi(tmp{1}, 'l2_mouth_cap_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_mouth_cap_time = char(t{1});                        
            elseif strcmpi(tmp{1}, 'l2_molt_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_molt_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l2_molt_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l2_molt_time = char(t{1});

                %%% L3/L4 molts
            elseif strcmpi(tmp{1}, 'l3_lethargus_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_lethargus_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l3_lethargus_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_lethargus_time = char(t{1});
            elseif strcmpi(tmp{1}, 'l3_mouth_cap_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_mouth_cap_date = char(t{1});            
            elseif strcmpi(tmp{1}, 'l3_mouth_cap_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_mouth_cap_time = char(t{1});                        
            elseif strcmpi(tmp{1}, 'l3_molt_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_molt_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l3_molt_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l3_molt_time = char(t{1});

                %%% L4/adult molts
            elseif strcmpi(tmp{1}, 'l4_lethargus_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_lethargus_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l4_lethargus_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_lethargus_time = char(t{1});
            elseif strcmpi(tmp{1}, 'l4_mouth_cap_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_mouth_cap_date = char(t{1});            
            elseif strcmpi(tmp{1}, 'l4_mouth_cap_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_mouth_cap_time = char(t{1});            
            elseif strcmpi(tmp{1}, 'l4_molt_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_molt_date = char(t{1});
            elseif strcmpi(tmp{1}, 'l4_molt_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.l4_molt_time = char(t{1});            

               %%%% Time of first egg-laying
            elseif strcmpi(tmp{1}, 'first_egg_date')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.first_egg_date = char(t{1});
            elseif strcmpi(tmp{1}, 'first_egg_time')
                t = textscan(fid, '%s %*[^\n]',1);
                worm.first_egg_time = char(t{1});            
            end

        end

        delete(closeTheFile)
        
    else
        disp(['Could not load data from ' wormfile '. Skipping...']);
    end
   
end
