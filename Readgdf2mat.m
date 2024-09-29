% 需要在eeglab中安装BIOSIG的插件，通过插件将gdf的数据读进工作区
% 完成于2024.07.03
% 作者： 王毓晨


% 定义文件夹路径
folderPath = 'D:\WorkSpace\MyDataSets\BCICIV_2a_gdf\'
savefolderPath = 'D:\WorkSpace\SelfDataSet\BCICIV_2a_mat\'
% savefolderPath = 'D:\WorkSpace\SelfDataSet\BCICIV_2a_mat_400ms\'
if ~exist(savefolderPath, 'dir')
    mkdir(savefolderPath);
end
% 读取训练集true，测试集false
isTrain = false

channelNum = 22

% % 使用dir函数读取文件夹内容，设置属性为'file'以获取文件信息
% fileInfo = dir(fullfile(folderPath, '*.*'));  % 读取所有文件
% 
% % 过滤掉子文件夹和隐藏文件（如果需要）
% fileInfo(subs(fileInfo.isdir)) = [];
% fileInfo([fileInfo.isdir] | [fileInfo.isHidden]) = [];

% 只读取训练集 .gdf 文件
% fileInfo = dir(fullfile(folderPath, '*T.gdf'));
% 只读取测试集 .gdf 文件
% fileInfo = dir(fullfile(folderPath, '*E.gdf'));
if isTrain == true
    fileInfo = dir(fullfile(folderPath, '*T.gdf'));
else
    fileInfo = dir(fullfile(folderPath, '*E.gdf'));
end

% 提取文件名
fileNames = {fileInfo.name};

filesPath = {};

% 为提升处理速度,假定持续时间全为313，根据之前处理经验提前定好数据数组尺寸
% event_duration = 400;
event_duration = 1000;
data = {};
eeglabel_tmp = {};
indicesToKeep = {};
eeglabel = {};
chs = {};

eegdata = {};
event = {};

function [eegdata, event, chs] = read_gdf (sub, filesPath, isTrain, chs)
    EEG.etc.eeglabvers = 'dev'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_biosig(filesPath);
    if isTrain == true
        EEG.setname = sprintf('A0%dT', sub); % sub<=9
    else
        EEG.setname = sprintf('A0%dE', sub); % sub<=9
    end
    eegdata = EEG.data;
    event = EEG.event;
    originalchs = EEG.chanlocs;
    for i = 1:size(originalchs)
        chs{end+1} = originalchs(i).labels;
    end
    disp([EEG.setname, 'The file was successfully read'])
end

function [data, eeglabel_tmp, indicesToKeep] = init_para_size (isTrain, event_duration, event)
    data = zeros(288, 25, event_duration);
    eeglabel_tmp = zeros(size(event));
    % 找出 edftype 字段等于 769、770、771、772 的结构体元素的索引
    if isTrain == true
    end
    indicesToKeep = false(size(event));
    disp('参数初始化...')
end

function [data, label] =TrainDataProcess(eegdata, event, isTrain, indicesToKeep, event_duration, data, eeglabel)
    j = 1;
    disp('data读取中...')
    if isTrain == true
        %%% 训练集
        for i = 1:length(event)
            if any(ismember([769, 770, 771, 772], event(i).edftype))
                indicesToKeep(i) = true;
                data(j,:,:) = eegdata(:, event(i).latency:event(i).latency+event_duration-1);
                eeglabel_tmp(j) = event(i).edftype;
                j = j + 1;
            end
        end

        for k = 1:length(eeglabel_tmp) % 确保使用正确的数组长度
            if eeglabel_tmp(k) ~= 0 % 检查元素是否非零
                % eeglabel{end+1} = eeglabel_tmp(k) - 769; % 存储差值
                eeglabel{end+1} = eeglabel_tmp(k) - 768; % 存储差值
            end
        end
    else
        %%% 测试集
        for i = 1:length(event)
            if any(ismember(783, event(i).edftype))
                indicesToKeep(i) = true;
                data(j,:,:) = eegdata(:, event(i).latency:event(i).latency+event_duration-1);
                j = j + 1;
                eeglabel_tmp(j) = event(i).edftype;
            end
        end
        % for k = 1:length(eeglabel_tmp) % 确保使用正确的数组长度
        %     if eeglabel_tmp(k) ~= 0 % 检查元素是否非零
        %         eeglabel{end+1} = eeglabel_tmp(k) - 783; % 存储差值
        %     end
        % end
    end
    disp(['data size: ', size(data), ''])
    % 将 cell 数组转换为数值数组（如果需要）
    label = cell2mat(eeglabel);
    disp(['label size: ', size(label), ''])
end

%%% main function

    % data size: 288x25x313 double
    % label size: 1x288 double
    % chs size: 1x25 cell

% for sub = 1 % 调试用
for sub = 1:length(fileNames)
    filesPath = fullfile(folderPath, fileNames{sub});
    [eegdata, event, chs] = read_gdf (sub, filesPath, isTrain, chs);
    [data, eeglabel_tmp, indicesToKeep] = init_para_size(isTrain, event_duration, event);
    [data, label] =TrainDataProcess(eegdata, event, isTrain, indicesToKeep, event_duration, data, eeglabel);
    if isTrain == true
        save_file_name = sprintf('A0%dT.mat', sub)
        label = int8(label');
    else
        save_file_name = sprintf('A0%dE.mat', sub)
        loaf_file_path = fullfile('D:\WorkSpace\SelfDataSet\BCICIV_2a_mat\testlabels\', save_file_name)
        load(loaf_file_path, 'classlabel');
        i = 1;
        while i <= length(classlabel)
            % label(i) = classlabel(i) - 1;
            label(i) = classlabel(i);
            i = i + 1;
        end
        label = int8(label'); % label [1, 2, 3, 4]
    end
    data = data(:, 1:channelNum, :);

    save_path = fullfile(savefolderPath, save_file_name);
    save(save_path, 'data', 'label','chs') 
end


