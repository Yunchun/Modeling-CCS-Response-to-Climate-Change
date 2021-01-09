long = ncread('CCS_7k_0-360_fred_grd.nc','lon_rho');
lat = ncread('CCS_7k_0-360_fred_grd.nc','lat_rho');
lon = find(235 <= long & long <= 238);
la = find(32 <= lat & lat <= 36);
index = zeros(1,1);
for l = 1:length(lon)
    for p = 1:length(la)
        if lon(l,1)==la(p,1)
            index1 = la(p,1);
            index = cat(1,index,index1);
        end
    end
end
index = index(2:length(index),1);

%get the row and column indices of long and lat elements that are within
%the range we want
grid_table = zeros(length(index),2);
for i = 1:length(index)
    grid_table(i,2)= fix(index(i)/182)+1;
    grid_table(i,1) = rem(index(i),182);
end

N = 50;
s_rho = ncread('CCS_7k_0-360_fred_grd.nc','s_rho');
hc = ncread('CCS_7k_0-360_fred_grd.nc','hc');
h = ncread('CCS_7k_0-360_fred_grd.nc','h');
Cs_r = ncread('CCS_7k_0-360_fred_grd.nc','Cs_r');
 
% need to loop through this f%
zeta = ncread('2007-01-12.nc','zeta'); % free surface elevation (zet ; 
%You will have to replace 'CCS1-RD.HCof03R_avg_2006-12-31T00:00:00.nc' for
%the name of your file.
zeta=squeeze(zeta(:,:,10));

for  k = 1:N;
    z0 = (hc.*s_rho(k) + h * Cs_r(k))./(hc + h);
    z_r(k,:,:) = zeta + (zeta + h).*z0; % this is the 3D matrix of depth at every grid point
end

%get the table where each element is the actual depth of lon and lat that
%are within the range we want
depth_table = zeros(50,length(index));
for i = 1:length(index)
    depth_table(:,i) = squeeze(z_r(:,grid_table(i,1),grid_table(i,2)));
end

avg_depth = zeros(1,1);
projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

%loop through all files
for K = 1 : num_files
    this_file = filenames{K};
    %read mesozoo value from each file
    meso = ncread(this_file, 'mesozooplankton');
    avgs = zeros(1,1);
  
    for i = 39:50
        depth = squeeze(meso(:,:,i,:));
        %get monthly average mesozoo for a horizontal profile
        for z = 1: 12
            subset = depth(:,:,z);
            meso_avgs = zeros(1,1);
            count = 0;
            for m = 1:length(index)
                meso_avgs = cat(1,meso_avgs, subset(index(m,1)));
                if subset(index(m,1)) ~= NaN
                    count = count + 1;
                end
            end
            meso_avgs = nansum(meso_avgs)/count;
            avgs = cat(1,avgs,meso_avgs);
        end
    end
    avgs = avgs(2:length(avgs),1);
    
    %sort monthly average mesozoo into a table where each column is the
    %monthly average mesozoo at one depth
    meso_avgs = zeros(12,12);
    for i = 0:11
        meso_avgs(:,i+1) = avgs((12*i+1):12*(i+1),1);
    end
    
    %get actual depths from row 37 to row 50(where actual depths are above
    %100 meters)
    depth = depth_table(37:50,1);
    interval = zeros(12,1);
    avgs = zeros(12,1);
    
    %i is the row index indicating month
    for i = 1:12
        sum = 0;
        %j is the column index indicating the depth level(39-50)
        for j = 1:12
            %get the depth interval value 
            interval(j,1) = (depth(j+2,1)-depth(j,1))/2 ;
            %integrate over all depths
            sum = sum + interval(j,1)*meso_avgs(i,j);
        end
        %divide the sum of mesozoo*depth interval by approx 100m
        avgs(i,1) = sum/nansum(interval);
    end
    %concatenate each year's monthly mesozoos averaged over 100m depths to
    %avg_depth
    avg_depth = cat(1,avg_depth,avgs);
end

avg_depth = avg_depth(2:length(avg_depth),1);

year = 1959*ones(12,1);
for i = 1:48
    year = cat(1,year,(1959+i)*ones(12,1));
end

month = [1:12]';
%Make a column vector consisted of 49 copies of month
month = repmat(month, 49,1);

meso_table = [year,month,avg_depth];

years = 1959:2007;
months = 1:12;
%sort the data into table where row index indicates month and column index
%indicates year
avg_depth = zeros(length(months),length(years));
for y = 1:length(years)
    row_index1 = find(meso_table(:,1) == years(y));
    subset = meso_table(row_index1(1):(row_index1(1)+length(row_index1)-1),:);
    for m = 1:length(months)
        row_index2 = find(subset(:,2) == months(m));
        avg_depth(m,y) = subset(row_index2,3);
    end
end

%add mesos of one month of all years together and divide by the
%number of years
avg_mean = zeros(1,12);
for m = 1:12
    avg_mean(1,m) = nanmean(avg_depth(m,:));
end

avg_anom = zeros(12,49);
for i = 1:12
    for j = 1:49
        avg_anom(i,j) = avg_depth(i,j)-avg_mean(1,i);
    end
end

figure(1)
subplot(2,1,1);
z = linspace(1959,2008,588);
plot(z,avg_anom(:),'-');
xticks(1959:2:2008);
refline([0 0])
xlabel('Year'),ylabel('millimole nitrogen meter^-3');
title('0-100m Mesozooplankton Monthly-mean Anomalies over Southern CCS')
hold on




lon = find(232 <= long & long <= 235);
la = find(38 <= lat & lat <= 42);
index = zeros(1,1);
for l = 1:length(lon)
    for p = 1:length(la)
        if lon(l,1)==la(p,1)
            index1 = la(p,1);
            index = cat(1,index,index1);
        end
    end
end
index = index(2:length(index),1);

grid_table = zeros(length(index),2);
for i = 1:length(index)
    grid_table(i,2)= fix(index(i)/182)+1;
    grid_table(i,1) = rem(index(i),182);
end

N = 50;
s_rho = ncread('CCS_7k_0-360_fred_grd.nc','s_rho');
hc = ncread('CCS_7k_0-360_fred_grd.nc','hc');
h = ncread('CCS_7k_0-360_fred_grd.nc','h');
Cs_r = ncread('CCS_7k_0-360_fred_grd.nc','Cs_r');
 
% need to loop through this f%
zeta = ncread('2007-01-12.nc','zeta'); % free surface elevation (zet ; 
%You will have to replace 'CCS1-RD.HCof03R_avg_2006-12-31T00:00:00.nc' for
%the name of your file.
zeta=squeeze(zeta(:,:,1));

for  k = 1:N;
    z0 = (hc.*s_rho(k) + h * Cs_r(k))./(hc + h);
    z_r(k,:,:) = zeta + (zeta + h).*z0; % this is the 3D matrix of depth at every grid point
end

depth_table = zeros(50,length(index));
for i = 1:length(index)
    depth_table(:,i) = squeeze(z_r(:,grid_table(i,1),grid_table(i,2)));
end

avg_depth = zeros(1,1);
projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

for K = 1 : num_files
    this_file = filenames{K};
    meso = ncread(this_file, 'mesozooplankton');
    avgs = zeros(1,1);
    for i = 39:50
        depth = squeeze(meso(:,:,i,:));
        for z = 1: 12
            subset = depth(:,:,z);
            meso_avgs = zeros(1,1);
            count = 0;
            for m = 1:length(index)
                meso_avgs = cat(1,meso_avgs, subset(index(m,1)));
                if subset(index(m,1)) ~= NaN
                    count = count + 1;
                end
            end
            meso_avgs = nansum(meso_avgs)/count;
            avgs = cat(1,avgs,meso_avgs);
        end
    end
    avgs = avgs(2:length(avgs),1);
    
    
    meso_avgs = zeros(12,12);
    for i = 0:11
        meso_avgs(:,i+1) = avgs((12*i+1):12*(i+1),1);
    end
    
    depth = depth_table(37:50,1);
    interval = zeros(12,1);
    avgs = zeros(12,1);
    
    for i = 1:12
        sum = 0;
        for j = 1:12
            interval(j,1) = (depth(j+2,1)-depth(j,1))/2 ;
            sum = sum + interval(j,1)*meso_avgs(i,j);
        end
        avgs(i,1) = sum/nansum(interval);
    end
    avg_depth = cat(1,avg_depth,avgs);
end

avg_depth = avg_depth(2:length(avg_depth),1);

year = 1959*ones(12,1);
for i = 1:48
    year = cat(1,year,(1959+i)*ones(12,1));
end

month = [1:12]';
%Make a column vector consisted of 49 copies of month
month = repmat(month, 49,1);

meso_table = [year,month,avg_depth];

years = 1959:2007;
months = 1:12;
avg_depth = zeros(length(months),length(years));
for y = 1:length(years)
    row_index1 = find(meso_table(:,1) == years(y));
    subset = meso_table(row_index1(1):(row_index1(1)+length(row_index1)-1),:);
    for m = 1:length(months)
        row_index2 = find(subset(:,2) == months(m));
        avg_depth(m,y) = subset(row_index2,3);
    end
end

avg_mean = zeros(1,12);
for m = 1:12
    avg_mean(1,m) = nanmean(avg_depth(m,:));
end

avg_anom = zeros(12,49);
for i = 1:12
    for j = 1:49
        avg_anom(i,j) = avg_depth(i,j)-avg_mean(1,i);
    end
end

subplot(2,1,2);
z = linspace(1959,2008,588);
plot(z,avg_anom(:),'-');
% ylim([-0.05 0.07]);
xticks(1959:2:2008);
refline([0 0]);
xlabel('Year'),ylabel('millimole nitrogen meter^-3');
title('0-100m Mesozooplankton Monthly-mean Anomalies over Northern CCS');
hold off