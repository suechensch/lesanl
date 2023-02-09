clf;
clear all; close all;
mytitle={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)', ...
        '(k)','(l)', '(m)','(n)','(o)','(p)','(q)'};
varname={'tke','tke buoyancy production','tke shear production','tke dissipation'};
oname={'tke_bt','tke_bp','tke_sh','tke_dp'};
cmin=[-2*10-5,-4*10^-6,-4*10^-5,-2*10^-4];
cmax=[0.5*10^-5,1.5*10^-4, 5*10^-5, 5*10^-5];

subn='fig8'; %output figure number
lenp=4;
lenv=length(varname);
figure;
hold all
TT=tiledlayout(lenv,lenp);
TT.Padding = 'compact';
TT.TileSpacing = 'compact';
data_path='../data/';

for npp=1:lenp
expn={'L1','L1sst','L2','L2sst'};
exp=expn{npp}
dataDIR=[data_path exp '/'];
matDataDIR=['../data/' exp '/'];
runname=expn{npp};

if strcmp(exp,'L2') == 1
    dateStamp='2723';
    dtg='2019092723';
    lt=23;
    tlen=1620;
    fcst=24;
elseif strcmp(exp,'L1') == 1
    dateStamp='2721';
    dtg='2019092721';
    lt=21;
    tlen=1620;
    fcst=24;
elseif strcmp(exp,'L1sst') == 1
    dateStamp='2721';
    dtg='2019092721';
    lt=21;
    tlen=1260;
    fcst=21;
elseif strcmp(exp,'L2sst') == 1
    dateStamp='2723';
    dtg='2019092723';
    lt=23;
    tlen=1258;
    fcst=20;
end
    tlen=1260;
    fcst=21;

filenameSave=[matDataDIR runname '_' dateStamp '.mat']; 
outPath=matDataDIR


if exist(outPath)~=7 % check if outPut folder exists, else create
   mkdir(outPath);
end
load (filenameSave);  
%%%

GMT2LT=8; % conversion from GMT to local time
sunrise_lt=6; % model forecast hour + initial hour + GMT to LT 
sunset_lt=18; % model forecast hour + initial hour + GMT to LT
t1=4;
t2=27;
add_lt=8;
len=length(lh);
tint=1/60; % in hour

kk=length(zm);
k2=kk-80;
datevec(datenum(2019,09,27,lt+4+add_lt,0,0));
for k=1:k2
  ic=0;
  for t=t1:t2
    ic=ic+1;
    zm2d(k,ic)=zm(k)*0.001;
    datev(npp,k,ic)=datenum(datevec(datenum(2019,09,27,lt+t+add_lt,0,0)));
  end
end

for nv=1:length(varname)
  ;
  colormap(jet)

% average to hourly

  for k = 1:size(radht_lw,1)
     if nv == 1
       fld(nv,k,:)=time_smooth(tke(k,1:tlen),60);
     elseif nv == 2
       fld(nv,k,:)=time_smooth(tke_sh(k,1:tlen),60);
     elseif nv == 3
       fld(nv,k,:)=time_smooth(tke_tp(k,1:tlen),60);
     elseif nv == 4
       fld(nv,k,:)=time_smooth(tke_dp(k,1:tlen),60);
     end
  end

end % nv
cmin=[0,-9*10^-5,-4*10^-5,-2*10^-4];
cmax=[0.1,6.4*10^-5, 5*10^-5, 1*10^-4];
  for nv=1:4
  ax = nexttile;
  fld2d=squeeze(fld(nv,1:k2,60:60:tlen));
  pcolor(squeeze(datev(npp,:,1:18)),zm2d(:,1:18),fld2d(:,4:fcst));
  shading interp;
  datetick('x','hh','keeplimits','keepticks')
  axis tight
  fldmx=max(max(fld2d(:,4:fcst)))
  fldmin=min(min(fld2d(:,4:fcst)))
  caxis([cmin(nv) cmax(nv)])
  if npp == 4
    colorbar('location','southoutside')
    xlabel('Local Time (h)');
  end

  if nv ==1 
    ylabel ('Height (km)');
  end
  jj=4*(npp-1)+nv;
  title(mytitle{jj});      
  ax.TitleHorizontalAlignment = 'left';
  clear ax
  box on

end % nv
end % number of exp

  grid on;
  fname=['../figures/fig8_all_tkeTend_hgt.png'];
  print('-dpng',fname);

