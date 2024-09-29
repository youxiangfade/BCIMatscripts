% 完成于2024.09.27
% 相较于ReadF2162a2mat，增加了对电极排序和利用FC5，FC6插值得到FCz数据
% 作者： 王毓晨

EEG.etc.eeglabvers = '2024.2'; % this tracks which version of EEGLAB is being used, you may ignore it

rootPath = 'D:\WorkSpace\MyDataSets\BCICIV_2a_216\CZC';
savefolderPath = 'D:\WorkSpace\SelfDataSet\BCICIV_2a_mat_216\'
Daypath = 'D2\'
fileInfo = dir(fullfile(rootPath, Daypath,'run*_Raw_Data.mat'));
disp(fileInfo);

% 提取文件名
fileNames = {fileInfo.name};

filesPath = {};
for sub = 1:length(fileNames)
    filesPath{sub} = fullfile(rootPath, Daypath, fileNames{sub});
end

function [eegdata_sub, event_Raw, originalchs] = Read_mat (sub, filesPath, isTrain)
    if isTrain==true
        EEG = pop_loadbva(filesPath{sub});
    else
    end
    EEG=pop_chanedit(EEG, 'changefield',{64,'labels','FCz'}); % 将lz通道改名为FCz
    EEG=pop_chanedit(EEG, 'changefield',{64,'theta','0'},'changefield',{64,'radius','0.17222'},'changefield',{64,'X','0.37049'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.85717'},'changefield',{64,'sph_phi','59'},'changefield',{64,'sph_theta','0'},'changefield',{64,'radius','0.125'}); % 更改FCz坐标
    % EEG = pop_select( EEG, 'channel',{'Fz','FC1','C3','CP1','Pz','CP2','Cz','C4','FC2','FC3','C1','C5','CP3','P1','POz','P2','CPz','CP4','C6','C2','FC4','FCz'}); % 只保留需要的2a数据集的22个EEG通道
    EEG = pop_select( EEG, 'channel',{'Fz','FC1','C3','CP1','Pz','CP2','Cz','C4','FC2','FC3','C1','C5','CP3','P1','POz','P2','CPz','CP4','C6','C2','FC4','FCz','FC5','FC6','FCz'}); % 只保留需要的2a数据集的21个EEG通道没有FCz使用FC5和FC6插值并存入原Oz中,一共24个通道
    if EEG.srate ~= 250
        EEG = pop_resample( EEG, 250); % 降采样到250/sec
        disp('采样率更改为250/sec')
    else
        disp('采样率为250/sec')
    end
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',32,'plotfreqz',0); % 带通滤波的频带为(0.5~32)，不画滤波频带图
    % 从EEG结构体中提取 'data', 'event', 'chanlocs'等字段
    eegdata_sub = EEG.data;
    event_Raw = EEG.event;
    originalchs = EEG.chanlocs;
end

function [rearrangedData] = rearrangeData(originalchs, target_order, eegdata_Raw)
    Raw_label_order = {originalchs.labels};
    % % FC5,FC6用于插值FCz
    % target_order = {'Fz','FC3','FC1','FCz','FC2','FC4','C5','C3','C1','Cz','C2','C4','C6','CP3','CP1','CPz','CP2','CP4','P1','Pz','P2','POz','FC5','FC6'}; 
    % 获取对应的序号
    target_label_indices = zeros(length(target_order), 1); % 注意变量名的一致性
    for i = 1:length(target_order)
        index = find(strcmp(Raw_label_order, target_order{i}));
        if ~isempty(index) % 确保找到了索引
            target_label_indices(i) = index(1); % 使用找到的第一个索引
        else
            target_label_indices(i) = NaN; % 如果找不到，赋值为NaN或其他标记值
        end
    end
    
    % 获取特定电极的索引
    fc5_index = find(strcmp(Raw_label_order, 'FC5'));
    fc6_index = find(strcmp(Raw_label_order, 'FC6'));
    fcz_index = find(strcmp(Raw_label_order, 'FCz'));
    
    % 计算 FCz 的数据
    if ~isnan(fcz_index) && ~isnan(fc5_index) && ~isnan(fc6_index)
        eegdata_Raw(fcz_index, :) = (eegdata_Raw(fc5_index, :) + eegdata_Raw(fc6_index, :)) / 2;
    end
    target_label_indices = target_label_indices(1:end-2);
    % 重新排列数组
    rearrangedData = eegdata_Raw(target_label_indices, :);
    % rearrangedData = eegdata_Raw(target_label_indices, :);
    
    % 显示重新排列后的数组大小，确认操作成功
    disp(size(rearrangedData));
end

function [event] = Event_Process (event_Raw)
    % 定义要筛选的标签及其对应的数字
    type_map = containers.Map({'S 10', 'S 20', 'S 30', 'S 40'}, [1, 2, 3, 4]);
    
    % 筛选出包含特定标签的结构体元素
    selected_tags = {'S 10', 'S 20', 'S 30', 'S 40'};
    filtered_events = event_Raw(ismember({event_Raw.type}, selected_tags));
    
    % 为筛选后的结构体数组添加新字段 'label' 并赋值
    for i = 1:length(filtered_events)
        if type_map.isKey(filtered_events(i).type)
            filtered_events(i).label = type_map(filtered_events(i).type);
        else
            filtered_events(i).label = NaN; % 如果没有匹配的类型，可以赋值为 NaN 或其他适当的值
        end
        % 对 'latency' 字段数据进行取整
        filtered_events(i).latency = round(filtered_events(i).latency);
    end
    
    % 创建一个新的结构体数组 event_tmp，只包含 label 和 latency 字段
    event = struct('label', {}, 'latency', {});
    
    % 遍历 filtered_events，将 label 和 latency 字段的值赋给 event_tmp
    for i = 1:length(filtered_events)
        event(i).label = filtered_events(i).label;
        event(i).latency = filtered_events(i).latency;
    end
    
    % 显示筛选、添加新字段并取整后的结构体数组
    disp(event);
end

function [data_tmp, eeglabel_tmp] =Data_Process(eegdata_sub, event, event_duration) 
    % 初始化 data 和 eeglabel_tmp
    j = 1;
    num_events = length(event); % 假设 event 是一个结构体数组
    max_duration = size(eegdata_sub, 2); % 假设 eegdata_Raw 是一个 2D 矩阵
    data_tmp = zeros(num_events, size(eegdata_sub, 1), event_duration); % 预分配空间
    eeglabel_tmp = zeros(num_events, 1); % 预分配空间
    
    for i = 1:num_events
        if any(event(i).label == [1, 2, 3, 4])
            % 确保不会超出 eegdata_Raw 的边界
            start_idx = event(i).latency;
            end_idx = min(start_idx + event_duration - 1, max_duration);
            
            % 提取数据片段
            data_tmp(j, :, 1:end_idx - start_idx + 1) = eegdata_sub(:, start_idx:end_idx);
            
            % 存储标签
            eeglabel_tmp(j) = event(i).label;
            
            % 更新索引
            j = j + 1;
        end
    end
    
    % 移除未使用的空间（如果有必要）
    data_tmp = data_tmp(1:j-1, :, :);
    eeglabel_tmp = eeglabel_tmp(1:j-1);
    
    % % 显示结果
    % disp(data_tmp);
    % disp(eeglabel_tmp);
end

% data = zeros(288, 22, 1000);
% label = zeros(288, 1);
data = zeros(144, 22, 1000);
label = zeros(144, 1);
target_order = {'Fz','FC3','FC1','FCz','FC2','FC4','C5','C3','C1','Cz','C2','C4','C6','CP3','CP1','CPz','CP2','CP4','P1','Pz','P2','POz','FC5','FC6'}; 
% for sub = 1:length(filesPath)
% for sub = 1:3
for sub = 1:6
    event_duration = 1000; % 样本数,本代码共保存1000/250=4s的数据
    isTrain = true;
    disp('Read mat ...')
    [eegdata_sub, event_Raw, originalchs] = Read_mat (sub, filesPath, isTrain);
    [eegdata_sub] = rearrangeData(originalchs, target_order, eegdata_sub);
    disp('Process Event ...')
    [event] = Event_Process (event_Raw);
    disp('Process Data ...')
    [data_tmp, eeglabel_tmp] = Data_Process(eegdata_sub, event, event_duration);
    data((sub-1)*48+1:(sub-1)*48+48, :, :) = data_tmp(:,:,:);
    label((sub-1)*48+1:(sub-1)*48+48,1) = eeglabel_tmp(:);
    % data((sub-4)*48+1:(sub-4)*48+48, :, :) = data_tmp(:,:,:);
    % label((sub-4)*48+1:(sub-4)*48+48,4) = eeglabel_tmp(:);
end

save_file_name = "A03E.mat";
save_path = fullfile(savefolderPath, save_file_name);
save(save_path, 'data', 'label') 
