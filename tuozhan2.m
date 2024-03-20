
%need m_map
% 加载 MAT 文件
argo_data = load('sst_Argo_201511.mat');

% 提取数据
lat = argo_data.lat(:);  % 纠正形状不一致的问题
lon = argo_data.lon(:);
temp = argo_data.temp(:);
% 筛选 10-50N，120-170E 区域的数据
east_sea_mask = (lat >= 10) & (lat <= 50) & (lon >= 120) & (lon <= 170);
lat_east_sea = lat(east_sea_mask);
lon_east_sea = lon(east_sea_mask);
temp_east_sea = temp(east_sea_mask);
% 创建网格
[lon_grid, lat_grid] = meshgrid(120:1:170, 10:1:50);

% 找出 temp_east_sea 中不是 NaN 的索引
validIdx = ~isnan(temp_east_sea);

% 使用这些索引来过滤数据，只保留有效的数据点
lat_valid = lat_east_sea(validIdx);
lon_valid = lon_east_sea(validIdx);
temp_valid = temp_east_sea(validIdx);




% 执行 IDW 插值
temp_grid = idw_interpolation(lon_valid, lat_valid, temp_valid, lon_grid, lat_grid);
% 绘制网格化数据的空间分布图
figure;
contourf(lon_grid, lat_grid, temp_grid, 100, 'LineColor', 'none');
colorbar;
hold on;
scatter(lon_east_sea, lat_east_sea, 36, 'w', 'filled');
title('IDW Interpolation of Sea Surface Temperature in the Expanded Region (2015-11)');
xlabel('Longitude');
ylabel('Latitude');
hold off;

% 加载 NETCDF 文件
ersst_data = ncread('ersst.v5.201511.nc', 'sst');
lon=ncread('ersst.v5.201511.nc','lon');
lat=ncread('ersst.v5.201511.nc','lat');


lon_range=[120 170];
lat_range=[10 50];
i_inx=find(lon<=lon_range(end) & lon>=lon_range(1));
j_inx=find(lat<=lat_range(end) & lat>=lat_range(1));


ersst=ersst_data(i_inx,j_inx);
% 计算平均海表面温度（假设已经在正确的维度上）
 % 调整这里以匹配数据的时间和深度维度






 figure
    m_proj('Equidistant Cylindrical','lon',lon_range,'lat',lat_range);
    s=m_pcolor(lon,lat,(ersst(:,:,i))');
    s.FaceColor = 'interp';
    hold on
    m_coast('line','color','k','linewidth',1,'linestyle','-');
    m_grid('linestyle','none','box','off','tickdir','out','LineWidth',0.5);
    hold on
    caxis([-2,2]);
    colormap(m_colmap('diversing'));
    %shading flat
    colorbar
    title(datestr(sstmeantime(1489+(i-1))))
    print(gcf,datestr(sstmeantime(1489+(i-1))),'-dpng')




















function Z = idw_interpolation(x, y, z, XI, YI, p)
    % 检查是否提供了 p 参数，如果没有，则默认为 2
    if (nargin < 6)
        p = 2;
    end

    % 初始化输出网格的大小
    Z = zeros(size(XI));
    
    % 遍历输出网格中的每一个点
    for i = 1:numel(XI)
        % 计算网格点与每一个数据点的距离
        d = sqrt((XI(i) - x).^2 + (YI(i) - y).^2);
        
        % 防止除以0的情况，将距离为0的位置设置为一个非常小的数
        d(d == 0) = eps;
        
        % 计算权重
        w = 1 ./ d.^p;
        
        % 计算插值结果
        Z(i) = sum(w .* z) / sum(w);
    end
end
