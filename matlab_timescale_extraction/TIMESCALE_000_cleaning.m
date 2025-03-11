
%% cleanup of some matrices

% clear
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% load([path2go 'meg_amygdala.mat'])
% for u = 1 : length(spikes)
%      spikes{u} = spikes{u}/1000;
%      spikes_rest{u} = spikes_rest{u}/1000;
%      spikes_task{u} = spikes_task{u}/1000;
%      cell_info{u}.t_on = cell_info{u}.t_on /1000;
%      cell_info{u}.t_off = cell_info{u}.t_off /1000;
% end
% save([path2go 'meg_amygdala.mat'],'spikes','spikes_rest',"spikes_task","cell_info")
% 
% clear
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% load([path2go 'meg_scACC.mat'])
% for u = 1 : length(spikes)
%      spikes{u} = spikes{u}/1000;
%      spikes_rest{u} = spikes_rest{u}/1000;
%      spikes_task{u} = spikes_task{u}/1000;
%      cell_info{u}.t_on = cell_info{u}.t_on /1000;
%      cell_info{u}.t_off = cell_info{u}.t_off /1000;
% end
% save([path2go 'meg_scACC.mat'],'spikes','spikes_rest',"spikes_task","cell_info")
% 

% clear
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% w1 = load([path2go 'wirth1_hc.mat']);
% w2 = load([path2go 'wirth2_hc.mat']);
% spikes = [w1.spikes w2.spikes];
% spikes_rest = [w1.spikes_rest w2.spikes_rest];
% spikes_task = [w1.spikes_task w2.spikes_task];
% cell_info = [w1.cell_info w2.cell_info];
% 
% remove = [];
% for u = 1 : length(spikes)
%      remove(u) = isempty(cell_info{u});
% end
% spikes(find(remove==1)) = [];
% spikes_rest(find(remove==1)) = [];
% spikes_task(find(remove==1)) = [];
% cell_info(find(remove==1)) = [];
% save([path2go 'wirth_hc.mat'],'spikes','spikes_rest',"spikes_task","cell_info")
% % 
% clear
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% load([path2go 'minxha_ha.mat'])
% x=[0 0];
% for u = 1 : length(spikes)
%     if ismember(cell_info{u}.BrainArea,'amygdala')
%         x(1) = x(1) + 1;
%         spikes_amg{x(1)} = spikes{u};
%         cell_info_amg{x(1)} = cell_info{u};
%     elseif ismember(cell_info{u}.BrainArea,'hippocampus')
%         x(2) = x(2) + 1;
%         spikes_hp{x(2)} = spikes{u};
%         cell_info_hp{x(2)} = cell_info{u};
%     end
% end
% spikes = spikes_hp;
% cell_info = cell_info_hp;
% save([path2go 'minxha_hippocampus.mat'],'spikes',"cell_info")
% 
% spikes = spikes_amg;
% cell_info = cell_info_amg;
% save([path2go 'minxha_amygdala.mat'],'spikes',"cell_info")
% % 
% 
% clear
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% load([path2go 'minxha_mfc.mat'])
% x=0;
% for u = 1 : length(spikes)
%     if ismember(cell_info{u}.BrainArea,'dorsal ACC')
%         x(1) = x(1) + 1;
%         spikes_acc{x(1)} = spikes{u};
%         cell_info_acc{x(1)} = cell_info{u};    
%     end
% end
% spikes = spikes_acc;
% cell_info = cell_info_acc;
% save([path2go 'minxha_dACC.mat'],'spikes',"cell_info")

% create Fontanier dataset:
% clear
% path2go = 'F:\PROCYK\DECCA\UNITS COMBINED\ACC\';
% cd(path2go)
% list = dir('*.mat');
% 
% for u = 1 : length(list)
%     clear SPK keep_evt sub_evt sub_evtcode evts
%     load([path2go list(u).name])
% 
%      if isfield(SPK,'originalevent')
%           evts = SPK.originalevent;
%      else
%           evts = SPK.event;
%      end
%      keep_evt = evts.timestamp>SPK.clip(1) & evts.timestamp<SPK.clip(2) ;

%      sub_evt = evts.timestamp(keep_evt);
%      sub_evtcode = evts.data(keep_evt);
% 
%      t_on = sub_evt(sub_evtcode == 100);
%      t_off = sub_evt(sub_evtcode == 101)-0.5; % remove 500ms to account for some of the ITI occuring before this code
%      t_off(t_off<=0)=[];
%      [~, loc] = min([t_on(1) t_off(1)]);
%      if loc == 2
%           t_off = t_off(2:end);
%      end
%      fake = find(diff(t_on)<1);
%      if ~isempty(fake)
%           t_on(fake) = [];
%      end
% 
%      nbtr = min([length(t_on) length(t_off)]);
%      t_on = t_on(1:nbtr);
%      t_off = t_off(1:nbtr);
%     temp =  t_on(2:end)-t_off(1:end-1);
%     disp([num2str(u) ' --- ' num2str(temp(end))])
% 
%      behav.task = [t_on t_off];
%      behav.rest = [t_off(1:end-1) t_on(2:end)];
% 
%      cell_info{u}.t_on = t_on;
%      cell_info{u}.t_off = t_off;
%      cell_info{u}.animal_id =  list(u).name(1);
%      cell_info{u}.unit =  list(u).name;
%      cell_info{u}.brain_area =  'MCC';
%      
%      spikes{u} = SPK.clipped.timestamp';
% 
%                spikes_rest{u}= [] ;
%                for tr = 1 : length(behav.rest(:,1))
%                     spikes_rest{u} = [spikes_rest{u} , single(SPK.clipped.timestamp(SPK.clipped.timestamp>=behav.rest(tr,1) & SPK.clipped.timestamp<=behav.rest(tr,2)))' ] ;
%                end
%                
%                spikes_task{u}= [] ;
%                for tr = 1 : length(behav.task(:,1))
%                     spikes_task{u} = [spikes_task{u} , single(SPK.clipped.timestamp(SPK.clipped.timestamp>=behav.task(tr,1) & SPK.clipped.timestamp<=behav.task(tr,2)))' ] ;
%                end
%     
% 
% end
% 
% path2go = 'S:\TIMESCALES\data\datasets_revision\';
% save([path2go 'fontanier_mcc.mat'],'spikes','spikes_rest',"spikes_task","cell_info") 
% 

