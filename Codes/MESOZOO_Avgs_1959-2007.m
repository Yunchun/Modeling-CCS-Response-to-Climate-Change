projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );
mesozoo = zeros(1,1);
for K = 1 : num_files
  this_file = filenames{K};
  mesozoo_total = ncread(this_file, 'mesozooplankton');
  mesozoo_grid_depth50 = squeeze(mesozoo_total(50,100,50,:));
  mesozoo = cat(1,mesozoo,mesozoo_grid_depth50);
end
mesozoo = mesozoo(2:589,1);

%Testing
%MESOZOO = ncread('2007-01-12.nc','mesozooplankton');
%MESOZOO = squeeze(MESOZOO(50,100,50,:));

year = 1959*ones(12,1);
for i = 1:48
    year = cat(1,year,(1959+i)*ones(12,1));
end

month = [1:12]';
%Make a column vector consisted of 49 copies of month
month = repmat(month, 49,1);

mesozoo_table = [year,month,mesozoo];

%Calculate mesozooplankton Monthly Averages from 1959 to 2007 
%over a grid cell at depth 50 in CCS
years = 1959:2007;
months = 1:12;
mesozoo = zeros(length(months),length(years));
for y = 1:length(years)
    row_index1 = find(mesozoo_table(:,1) == years(y));
    mesozoo_subset = mesozoo_table(row_index1(1):(row_index1(1)+length(row_index1)-1),:);
    for m = 1:length(months)
        row_index2 = find(mesozoo_subset(:,2) == months(m));
        mesozoo(m,y) = mesozoo_subset(row_index2,3);
    end
end

figure(1)
%Plot mesozooplankton Monthly Averages from 1959 to 2007 over a grid cell at depth 50
subplot(4,1,1);
z = linspace(1959,2007,588);
plot(z,mesozoo(:),'-');
ylim([0 0.1]);
xticks(1959:3:2007);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8),
title('Mesoozooplankton Monthly Averages Over a Grid Cell at Depth 50')
hold on

%Calculate Mesozooplankton Monthly-Mean over a grid cell at Depth 50
mesozoo_Mean = zeros(1,12);
for m = 1:length(months)
    mesozoo_Mean(1,m) = nanmean(mesozoo(m,:));
end

%Calculate Mesozooplankton Monthly-Anomalies over a grid cell at Depth 50
mesozoo_Mean_Anom = zeros(12,49);
for i = 1:12
    for j = 1:49
        mesozoo_Mean_Anom(i,j) = mesozoo(i,j)-mesozoo_Mean(1,i);
    end
end

%Plot Mesozooplankton Monthly-Anomalies over a grid cell at Depth 50
subplot(4,1,2);
z = linspace(1959,2007,588);
plot(z,mesozoo_Mean_Anom(:),'-');
ylim([-0.02 0.05]);
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8),
title('Mesozooplankton Monthly-mean Anomalies over a Grid Cell at Depth 50')
hold on

%Plot Mesozooplankton Mean Anomaly in January
subplot(4,1,3);
plot(1959:2007,mesozoo_Mean_Anom(1,:),'-*');
ylim([-0.02 0.05]);
xticks(1959:3:2007);
refline([0 0]);
xlabel('Year');
ylabel('millimole nitrogen meter^-3','Fontsize',8),
title('Mesozooplankton Mean Anomaly in January over a Grid Cell at Depth 50')
hold on

%Calculate Mesozooplankton Monthly-Average Standard Deviation 
mesozoo_std = zeros(1,12);
for i = 1:12
    mesozoo_std(1,i) = nanstd(mesozoo(i,:));
end

%Calculate Upper Bound of Mesozooplankton Monthly Averages
mesozoo_aveUpper = zeros(1,12);
for i = 1:12
    mesozoo_aveUpper(1,i) = mesozoo_Mean(1,i)+mesozoo_std(1,i);
end

%Calculate Lower Bound of Mesozooplankton Monthly Averages
mesozoo_aveLower = zeros(1,12);
for i = 1:12
    mesozoo_aveLower(1,i) = mesozoo_Mean(1,i)-mesozoo_std(1,i);
end

%Plot Range of Mesozooplankton Monthly Averages
subplot(4,1,4)
plot(1:12,mesozoo_Mean(1,:),'o');
hold on
plot(1:12,mesozoo_aveUpper(1,:),'o');
hold on
plot(1:12,mesozoo_aveLower(1,:),'o');
hold on

fill([1:12 fliplr(1:12)],[mesozoo_aveUpper fliplr(mesozoo_Mean)],'y');
fill([1:12 fliplr(1:12)],[mesozoo_aveLower fliplr(mesozoo_Mean)],'y');
xlabel('Month');
ylabel('millimole nitrogen meter^-3','Fontsize',8),
set(gca,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
title('Range of Mesozooplankton Monthly Averages over at Depth 50')
hold off