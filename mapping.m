
temp_data_jan1980 = temp_data(:,:,1);

latmax = 90;
latmin = -90;
lonmin = -180;
lonmax = 180;

lonrange = lonmax - lonmin;
latrange = latmax -latmin;

extent  = size(temp_data_jan1980);
refvecrange = extent(2)/lonrange;
refvec = [refvecrange latmax lonmin];
squeeze 
% lat_new=180;
% long_new=360;
%  
% lat11=zeros(lat_new,long_new);
% x=latmax-0.5;
% for i=1:lat_new
%     for j=1:long_new
%         lat11(i,j)=x; 
%     end
%     x=x-1;
% end
    
% long11=zeros(lat_new,long_new);
% y=lonmin+0.5;
% for i=1:long_new 
%     for j=1:lat_new
%         long11(j,i)=y;
%     end
%     y=y+1;  
% end

figure
% subplot(2,2,4)
axesm miller
% worldmap([-0.5 40.5],[59.5 100.5])
setm(gca, 'MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...1
'Frame','on','FontSize',12)
setm(gca, 'MlabelLocation', 30, 'PlabelLocation',30,...
'MLabelParallel','north', 'MeridianLabel','on',...
'ParallelLabel','on')
% setm(gca,'MeridianLabel','on');
framem on

geoshow(gca, flipud(temp_data_jan1980), refvec, 'DisplayType', 'texturemap')
S = shaperead('cntry02.shp','UseGeoCoords',true);
geoshow([S.Lat], [S.Lon],'color','black','linewidth',0.2);
colorbar('horizon');
h_bar=findobj(gcf,'Tag','Colorbar');
initpos=get(h_bar,'Position');

set(h_bar,'Position',[initpos(1)+initpos(3)*0.25 initpos(2)-0.04+initpos(4)*0.25 initpos(3)*0.5 initpos(4)*0.5]);
set(gca,'linewidth',2,'fontsize',14);
% colormap(color_map);

% caxis([0 0.8])
set(gcf,'Color',[1,1,1])
axis off;
box off;