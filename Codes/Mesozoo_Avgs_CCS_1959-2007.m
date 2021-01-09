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

avgs = zeros(1,1);

projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

for K = 1 : num_files
  this_file = filenames{K};
  mesoz = ncread(this_file, 'mesozooplankton');
  mesoz_depth50 = squeeze(mesoz(:,:,50,:));
  
  for z = 1: 12
      mesoz_subset = mesoz_depth50(:,:,z);
      
      mesoz_avgs = zeros(1,1);
      count = 0;
      for m = 1:length(index)
          mesoz_avgs = cat(1,mesoz_avgs, mesoz_subset(index(m,1)));
          if mesoz_subset(index(m,1)) ~= NaN
              count = count + 1;
          end
      end
      mesoz_avgs = nansum(mesoz_avgs)/count;
      
      avgs = cat(1,avgs,mesoz_avgs);    
  end
  
end
avgs= avgs(2:589,1);

year = 1959*ones(12,1);
for i = 1:48
    year = cat(1,year,(1959+i)*ones(12,1));
end

month = [1:12]';
%Make a column vector consisted of 49 copies of month
month = repmat(month, 49,1);

mesoz_table = [year,month,avgs];

%Calculate Mesozooplankton Monthly Averages from 1959 to 2007 
%over CCS at depth 50 
years = 1959:2007;
months = 1:12;
mesoz = zeros(length(months),length(years));
for y = 1:length(years)
    row_index1 = find(mesoz_table(:,1) == years(y));
    mesoz_subset = mesoz_table(row_index1(1):(row_index1(1)+length(row_index1)-1),:);
    for m = 1:length(months)
        row_index2 = find(mesoz_subset(:,2) == months(m));
        mesoz(m,y) = mesoz_subset(row_index2,3);
    end
end

figure(1)
%Plot Mesozooplankton Monthly Averages from 1959 to 2007 over CCS at depth 50
subplot(4,1,1);
z = linspace(1959,2007,588);
plot(z,mesoz(:),'-');
xticks(1959:3:2007);
xlabel('Year');
title('Mesozooplankton Monthly Averages over CCS at Depth 50')
hold on

%Calculate Mesozooplankton Monthly-Mean over CCS at Depth 50
mesoz_Mean = zeros(1,12);
for m = 1:length(months)
    mesoz_Mean(1,m) = nanmean(mesoz(m,:));
end

%Calculate Mesozooplankton Monthly-Anomalies over CCS at Depth 50
mesoz_Mean_Anom = zeros(12,49);
for i = 1:12
    for j = 1:49
        mesoz_Mean_Anom(i,j) = mesoz(i,j)-mesoz_Mean(1,i);
    end
end

%Plot Mesozooplankton Monthly-Anomalies over CCS at Depth 50
subplot(4,1,2);
z = linspace(1959,2007,588);
plot(z,mesoz_Mean_Anom(:),'-');
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
title('Mesozooplankton Monthly-mean Anomalies over CCS')
hold on

%Plot Mesozooplankton Mean Anomaly in January
subplot(4,1,3);
plot(1959:2007,mesoz_Mean_Anom(1,:),'-*');
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
title('Mesozooplankton Mean Anomaly in January over CCS')
hold on

%Calculate Mesozooplankton Monthly-Average Standard Deviation 
mesoz_std = zeros(1,12);
for i = 1:12
    mesoz_std(1,i) = nanstd(mesoz(i,:));
end

%Calculate Upper Bound of Mesozooplankton Monthly Averages
mesoz_aveUpper = zeros(1,12);
for i = 1:12
    mesoz_aveUpper(1,i) = mesoz_Mean(1,i)+mesoz_std(1,i);
end

%Calculate Lower Bound of Mesozooplankton Monthly Averages
mesoz_aveLower = zeros(1,12);
for i = 1:12
    mesoz_aveLower(1,i) = mesoz_Mean(1,i)-mesoz_std(1,i);
end

%Plot Range of Mesozooplankton Monthly Averages
subplot(4,1,4)
plot(1:12,mesoz_Mean(1,:),'*');
hold on
plot(1:12,mesoz_aveUpper(1,:),'*');
hold on
plot(1:12,mesoz_aveLower(1,:),'*');
hold on

fill([1:12 fliplr(1:12)],[mesoz_aveUpper fliplr(mesoz_Mean)],'g');
fill([1:12 fliplr(1:12)],[mesoz_aveLower fliplr(mesoz_Mean)],'g');
xlabel('Month');
set(gca,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
title('Range of Mesozooplankton Monthly Averages over CCS')
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

avgs = zeros(1,1);

projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

for K = 1 : num_files
  this_file = filenames{K};
  mesoz = ncread(this_file, 'mesozooplankton');
  mesoz_depth50 = squeeze(mesoz(:,:,50,:));
  
  for z = 1: 12
      mesoz_subset = mesoz_depth50(:,:,z);
      
      mesoz_avgs = zeros(1,1);
      count = 0;
      for m = 1:length(index)
          mesoz_avgs = cat(1,mesoz_avgs, mesoz_subset(index(m,1)));
          if mesoz_subset(index(m,1)) ~= NaN
              count = count + 1;
          end
      end
      mesoz_avgs = nansum(mesoz_avgs)/count;
      
      avgs = cat(1,avgs,mesoz_avgs);    
  end
  
end
avgs= avgs(2:589,1);

year = 1959*ones(12,1);
for i = 1:48
    year = cat(1,year,(1959+i)*ones(12,1));
end

month = [1:12]';
%Make a column vector consisted of 49 copies of month
month = repmat(month, 49,1);

mesoz_table = [year,month,avgs];

%Calculate Meszooplankton Monthly Averages from 1959 to 2007 
%over NCCS at depth 50 
years = 1959:2007;
months = 1:12;
mesoz = zeros(length(months),length(years));
for y = 1:length(years)
    row_index1 = find(mesoz_table(:,1) == years(y));
    mesoz_subset = mesoz_table(row_index1(1):(row_index1(1)+length(row_index1)-1),:);
    for m = 1:length(months)
        row_index2 = find(mesoz_subset(:,2) == months(m));
        mesoz(m,y) = mesoz_subset(row_index2,3);
    end
end

figure(1)
%Plot Mesozooplankton Monthly Averages from 1959 to 2007 over NCCS at depth 50
subplot(4,1,1);
z = linspace(1959,2007,588);
plot(z,mesoz(:),'-');
xticks(1959:3:2007);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8);
legend('Mesoz Avgs over SCCS','Mesoz Avgs over NCCS');
title('Mesozooplankton Monthly Averages over CCS at Depth 50')
hold on

%Calculate Mesozooplankton Monthly-Mean over NCCS at Depth 50
mesoz_Mean = zeros(1,12);
for m = 1:length(months)
    mesoz_Mean(1,m) = nanmean(mesoz(m,:));
end

%Calculate Mesozooplankton Monthly-Anomalies over NCCS at Depth 50
mesoz_Mean_Anom = zeros(12,49);
for i = 1:12
    for j = 1:49
        mesoz_Mean_Anom(i,j) = mesoz(i,j)-mesoz_Mean(1,i);
    end
end

%Plot Mesozooplankton Monthly-Anomalies over NCCS at Depth 50
subplot(4,1,2);
z = linspace(1959,2007,588);
plot(z,mesoz_Mean_Anom(:),'-');
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8);
legend('Mesoz Anoms over SCCS','refline 0','Mesoz Anoms over NCCS');
title('Mesozooplankton Monthly-mean Anomalies over CCS')
hold on

%Plot Mesozooplankton Mean Anomaly in January
subplot(4,1,3);
plot(1959:2007,mesoz_Mean_Anom(1,:),'-*');
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8);
legend('Mesoz Anoms over SCCS','refline 0','Mesoz Anoms over NCCS');
title('Mesozooplankton Mean Anomaly in January over CCS')
hold on

%Calculate Mesozooplankton Monthly-Average Standard Deviation 
mesoz_std = zeros(1,12);
for i = 1:12
    mesoz_std(1,i) = nanstd(mesoz(i,:));
end

%Calculate Upper Bound of Mesozooplankton Monthly Averages
mesoz_aveUpper = zeros(1,12);
for i = 1:12
    mesoz_aveUpper(1,i) = mesoz_Mean(1,i)+mesoz_std(1,i);
end

%Calculate Lower Bound of Mesozooplankton Monthly Averages
mesoz_aveLower = zeros(1,12);
for i = 1:12
    mesoz_aveLower(1,i) = mesoz_Mean(1,i)-mesoz_std(1,i);
end

%Plot Range of Mesozooplankton Monthly Averages
subplot(4,1,4)
plot(1:12,mesoz_Mean(1,:),'*');
hold on
plot(1:12,mesoz_aveUpper(1,:),'*');
hold on
plot(1:12,mesoz_aveLower(1,:),'*');
hold on

fill([1:12 fliplr(1:12)],[mesoz_aveUpper fliplr(mesoz_Mean)],'y');
fill([1:12 fliplr(1:12)],[mesoz_aveLower fliplr(mesoz_Mean)],'y');
xlabel('Month');
ylabel('millimole nitrogen meter^-3','Fontsize',8);
legend('Mesoz Avgs Range over NCCS','Mesoz Avgs Range over SCCS');
set(gca,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
title('Range of Mesozooplankton Monthly Averages over CCS')
hold off
