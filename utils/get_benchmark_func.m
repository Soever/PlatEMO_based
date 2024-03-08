function func_handles = get_benchmark_func(name)
func_handles = {};
if name == "CEC2008"
    for i = 1:7
        func_name = "CEC2008_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end 
elseif name == "CEC2010"
    for i = 1:18
        func_name = "CEC2010_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end 
elseif name == "CEC2013"
    for i = 1:15
        func_name = "CEC2013_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
elseif name == "CEC2017"
    for i = 1:28
        func_name = "CEC2017_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
elseif name == "CEC2020"
    for i = 1:10
        func_name = "CEC2020_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
else
    for i = 1:7
        func_name = "CEC2008_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
    for i = 1:18
        func_name = "CEC2010_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
    for i = 1:15
        func_name = "CEC2013_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
    for i = 1:28
        func_name = "CEC2017_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end
    for i = 1:10
        func_name = "CEC2020_F"+ num2str(i); % 创建函数名字符串
        func_handles{end+1} = str2func(func_name); % 创建函数句柄
    end


end






