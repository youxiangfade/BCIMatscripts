% �����2024.09.21
% ���ߣ� ��ع��

EEG.etc.eeglabvers = '2024.2'; % this tracks which version of EEGLAB is being used, you may ignore it

rootPath = 'D:\WorkSpace\MyDataSets\BCICIV_2a_216\LRP'
savefolderPath = 'D:\WorkSpace\SelfDataSet\BCICIV_2a_mat_216\'
Daypath = 'D2\'
fileInfo = dir(fullfile(rootPath, Daypath,'run*_Raw_Data.mat'));
disp(fileInfo);

% ��ȡ�ļ���
fileNames = {fileInfo.name};

filesPath = {};
for sub = 1:length(fileNames)
    filesPath{sub} = fullfile(rootPath, Daypath, fileNames{sub});
end

function [eegdata_sub, event_Raw] = Read_mat (sub, filesPath, isTrain)
    if isTrain==true
        EEG = pop_loadbva(filesPath{sub});
    else
    end
    EEG=pop_chanedit(EEG, 'changefield',{64,'labels','FCz'}); % ��lzͨ������ΪFCz
    EEG=pop_chanedit(EEG, 'changefield',{64,'theta','0'},'changefield',{64,'radius','0.17222'},'changefield',{64,'X','0.37049'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.85717'},'changefield',{64,'sph_phi','59'},'changefield',{64,'sph_theta','0'},'changefield',{64,'radius','0.125'}); % ����FCz����
    % EEG = pop_select( EEG,'channel',{'FC1','Fz','C3','CP1','Pz','CP2','Cz','C4','FC2','FC3','C1','C5','CP3','P1','POz','P2','CPz','CP4','C6','C2','FC4','Iz'});
    %  ���Ĳ��ܸı�˳��
    EEG = pop_select( EEG, 'channel',{'Fz','FC1','C3','CP1','Pz','CP2','Cz','C4','FC2','FC3','C1','C5','CP3','P1','POz','P2','CPz','CP4','C6','C2','FC4','FCz'}); % ֻ������Ҫ��2a���ݼ���22��EEGͨ��
    EEG = pop_resample( EEG, 250); % ��������250/sec
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',32,'plotfreqz',0); % ��ͨ�˲���Ƶ��Ϊ(0.5~32)�������˲�Ƶ��ͼ
    % ��EEG�ṹ������ȡ 'data', 'event', 'chanlocs'���ֶ�
    eegdata_sub = EEG.data;
    event_Raw = EEG.event;
    % originalchs = EEG.chanlocs;
end

function [event] = Event_Process (event_Raw)
    % ����Ҫɸѡ�ı�ǩ�����Ӧ������
    type_map = containers.Map({'S 10', 'S 20', 'S 30', 'S 40'}, [1, 2, 3, 4]);
    
    % ɸѡ�������ض���ǩ�Ľṹ��Ԫ��
    selected_tags = {'S 10', 'S 20', 'S 30', 'S 40'};
    filtered_events = event_Raw(ismember({event_Raw.type}, selected_tags));
    
    % Ϊɸѡ��Ľṹ������������ֶ� 'label' ����ֵ
    for i = 1:length(filtered_events)
        if type_map.isKey(filtered_events(i).type)
            filtered_events(i).label = type_map(filtered_events(i).type);
        else
            filtered_events(i).label = NaN; % ���û��ƥ������ͣ����Ը�ֵΪ NaN �������ʵ���ֵ
        end
        % �� 'latency' �ֶ����ݽ���ȡ��
        filtered_events(i).latency = round(filtered_events(i).latency);
    end
    
    % ����һ���µĽṹ������ event_tmp��ֻ���� label �� latency �ֶ�
    event = struct('label', {}, 'latency', {});
    
    % ���� filtered_events���� label �� latency �ֶε�ֵ���� event_tmp
    for i = 1:length(filtered_events)
        event(i).label = filtered_events(i).label;
        event(i).latency = filtered_events(i).latency;
    end
    
    % ��ʾɸѡ��������ֶβ�ȡ����Ľṹ������
    disp(event);
end

function [data_tmp, eeglabel_tmp] =Data_Process(eegdata_sub, event, event_duration) 
    % ��ʼ�� data �� eeglabel_tmp
    j = 1;
    num_events = length(event); % ���� event ��һ���ṹ������
    max_duration = size(eegdata_sub, 2); % ���� eegdata_Raw ��һ�� 2D ����
    data_tmp = zeros(num_events, size(eegdata_sub, 1), event_duration); % Ԥ����ռ�
    eeglabel_tmp = zeros(num_events, 1); % Ԥ����ռ�
    
    for i = 1:num_events
        if any(event(i).label == [1, 2, 3, 4])
            % ȷ�����ᳬ�� eegdata_Raw �ı߽�
            start_idx = event(i).latency;
            end_idx = min(start_idx + event_duration - 1, max_duration);
            
            % ��ȡ����Ƭ��
            data_tmp(j, :, 1:end_idx - start_idx + 1) = eegdata_sub(:, start_idx:end_idx);
            
            % �洢��ǩ
            eeglabel_tmp(j) = event(i).label;
            
            % ��������
            j = j + 1;
        end
    end
    
    % �Ƴ�δʹ�õĿռ䣨����б�Ҫ��
    data_tmp = data_tmp(1:j-1, :, :);
    eeglabel_tmp = eeglabel_tmp(1:j-1);
    
    % % ��ʾ���
    % disp(data_tmp);
    % disp(eeglabel_tmp);
end

% data = zeros(288, 22, 1000);
% label = zeros(288, 1);
data = zeros(144, 22, 1000);
label = zeros(144, 1);
% for sub = 1:length(filesPath)
% for sub = 1:3
for sub = 4:6
    event_duration = 1000; % ������,�����빲����1000/250=4s������
    isTrain = true;
    disp('Read mat ...')
    [eegdata_sub, event_Raw] = Read_mat (sub, filesPath, isTrain);
    disp('Process Event ...')
    [event] = Event_Process (event_Raw);
    disp('Process Data ...')
    [data_tmp, eeglabel_tmp] = Data_Process(eegdata_sub, event, event_duration);
    % data((sub-1)*48+1:(sub-1)*48+48, :, :) = data_tmp(:,:,:);
    % label((sub-1)*48+1:(sub-1)*48+48,1) = eeglabel_tmp(:);
    data((sub-4)*48+1:(sub-4)*48+48, :, :) = data_tmp(:,:,:);
    label((sub-4)*48+1:(sub-4)*48+48,4) = eeglabel_tmp(:);
end

save_file_name = "A02E.mat";
save_path = fullfile(savefolderPath, save_file_name);
save(save_path, 'data', 'label') 
