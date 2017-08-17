clear;

% list all files with .nc (netcdf extension) and store in varible- 'files'
files = dir('*.nc');

% determine the dimensions of the files structure which indicates the
% number of files i.e. 432 files corresponding to 36 years of monthly mean
% data (36*12=432), for the period from January 1980 to December 2015
dim_files = size(files);

% each monthly mean temperature map has 361 latitudes and 576 longitudes
lat1=361;
lon1=576;

% initialize a 3D matrix (temp_data) with the 3 dimensions namely lat1,
% lon1, and number of total months (432)
temp_data = zeros(lat1,lon1,dim_files(1));
temp_data = temp_data*NaN;

for i = 1:dim_files(1)
    
    % data is a temporary variable to store intermediate temperature data.
    % 'M2TMNXSLV_5_12_4_T2M' is the temperature variable within the netcdf
    % files
    data = ncread(files(i).name,'M2TMNXSLV_5_12_4_T2M');
    data = single(data);
    
    % the data needs to be rotated once and the resulting matrix is stored
    % in temp_data. After the for loop is completed running, temp_data contains the
    % temperature data for each month
    temp_data(:,:,i) = rot90(data);

end

% create figure/image of the spatial distribution of temperature (January
% 1980, index of January 1980 is 1). Units are in degree Kelvin.
imagesc(temp_data(:,:,1));colorbar

% algorithm Implementation
[row,column,month1]=size(temp_data);
Month1=zeros(432,1);
for i=1:month1
    Month1(i)=i;
end
MonthsMean=mean(Month1);
MonthsNew=Month1-MonthsMean;

slopeOffset=zeros(row,column,2);
CoefficentOfDetermination=zeros(row,column,1);
CoefficentOfDeterminationInbuilt=zeros(row,column,1);
slopeOffsetInbuilt=zeros(row,column,2);
for i=1:row
    for j=1:column
        TempData=squeeze(temp_data(i,j,:));
        TempDataMean=mean(TempData);
       % TempdataNew represents (y-y')
        TempdataNew=TempData-TempDataMean;
        slopeOffset(i,j,1)=sum(TempdataNew.*MonthsNew)/sum(MonthsNew.^2);
        slopeOffset(i,j,2)=TempDataMean-slopeOffset(i,j,1)*MonthsMean;
        CoefficentOfDetermination(i,j)=((slopeOffset(i,j,1)*sqrt(sum(MonthsNew.^2)))/sqrt(sum(TempdataNew.^2))).^2;
        % algorithm Implememted
        %
        [p,S] = polyfit(Month1,TempData,1);
        slopeOffsetInbuilt(i,j,:)=p;
        TempDataPredacted=polyval(p,Month1);
        TempSSResid = sum((TempData - TempDataPredacted).^2);
        TempSStotal = (length(TempData)-1) * var(TempData);
        CoefficentOfDeterminationInbuilt(i,j)= 1-TempSSResid/TempSStotal;
        %
        %Inbuilt Implementation This implementation gives all parameters
        %but this is sluggis.
        %press cntrl+t  to uncomment and test
        
        
%         Model = LinearModel.fit(Month1,TempData);
%         slopeOffsetInbuilt(i,j,:)=Model.Coefficients.Estimate(:);
%         CoefficentOfDeterminationInbuilt(i,j)=Model.Rsquared.Ordinary;
        

% Inbuilt Finish
    end
end
%Root mean square error in slope
RmseErrorSlope=sqrt(mean(mean((slopeOffsetInbuilt(:,:,1)-slopeOffset(:,:,1)).^2)));
RmseErrorRsquare=sqrt(mean(mean((CoefficentOfDeterminationInbuilt-CoefficentOfDetermination).^2)));

%Plotting algorithm outputs
figure;
subplot(1,2,1)
imagesc(slopeOffset(:,:,1));colorbar
title('Slope Calculated from Algorithm');
subplot(1,2,2)
imagesc(CoefficentOfDetermination);colorbar
title('coefficient of determination from Algorithm');
subplot(1,2,1)
figure
%Plotting Inbuilt Parameters
subplot(1,2,1)
imagesc(slopeOffsetInbuilt(:,:,1));colorbar
title('Slope Calculated from Inbuilt');
subplot(1,2,2)
imagesc(CoefficentOfDeterminationInbuilt);colorbar
title('coefficient of determination from Inbuilt');
subplot(1,2,1)

% Plotting Difference Images
figure
imagesc(slopeOffsetInbuilt(:,:,1)-slopeOffset(:,:,1));colorbar;
title('slope inbuilt -slope algorithm');
figure
imagesc(CoefficentOfDeterminationInbuilt-CoefficentOfDetermination);colorbar;
title('r^2 inbuilt - r^2 algorithm')
RmseErrorSlope=sqrt(mean(mean((slopeOffsetInbuilt(:,:,1)-slopeOffset(:,:,1)).^2)))
RmseErrorRsquare=sqrt(mean(mean((CoefficentOfDeterminationInbuilt-CoefficentOfDetermination).^2)))


% Q4 extracting Particular Pixel and Plotting
x = input('Do you Want to extract a particular Pixel Yes then enter 1 x = ');
if x==1
    RowNo=input('Row No. = ');
    ColumnNo = input('Column No =');
    SelectedPixel=squeeze(temp_data(RowNo,ColumnNo,:));
    p = polyfit(Month1,SelectedPixel,1);
   TempDataPredacted=polyval(p,Month1);
   TempSSResid = sum((SelectedPixel - TempDataPredacted).^2);
   TempSStotal = (length(SelectedPixel)-1) * var(SelectedPixel);
   CoefficentOfDeterminationforPixel= 1-TempSSResid/TempSStotal
   figure
    plot(Month1,SelectedPixel,Month1,TempDataPredacted);    
    xlabel('Time');
    ylabel('Temperature');
    hold on;
end
