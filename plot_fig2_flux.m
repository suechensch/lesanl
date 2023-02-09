clf;
clear all; close all;
exp={'L1','L1sst','L2','L2sst'};
outPath=['../figures/'];
nexp=length(exp);

for np=1:nexp
  expn=exp{np};
  runname=expn;

dataDIR=['../data/' expn '/'];
matDataDIR=['../data/' expn '/'];

if strcmp(expn,'L2') == 1
    dateStamp='2723';
    dtg='2019092723';
elseif strcmp(expn,'L1') == 1
    dateStamp='2721';
    dtg='2019092721';
elseif strcmp(expn,'L1sst') == 1
    dateStamp='2721';
    dtg='2019092721';
elseif strcmp(expn,'L2sst') == 1
    dateStamp='2723';
    dtg='2019092723';
end
len2=1260; % size of L1sst experiment
filenameSave=[matDataDIR runname '_' dateStamp '.mat']; 
createMovie=false;


if exist(outPath)~=7 % check if outPut folder exists, else create
   mkdir(outPath);
end
load (filenameSave);  

% generate figures
%

varname={'surface sensible heat flux','surface latent heat flux','total surface heat flux'};
oname={'sh','lh','tflx'};
GMT2LT=8; % conversion from GMT to local time
sunrise_lt=6; % model forecast hour + initial hour + GMT to LT 
sunset_lt=18; % model forecast hour + initial hour + GMT to LT
t1=0;
t2=20;
  add_lt=8;
len=length(lh);
tint=1/60; % in hour
if strcmp(expn,'sstSep28') == 1 | strcmp(expn,'camp2xSep28-t2') == 1
  for t=1:len
    datev(t,np)=datenum(datevec(datenum(2019,08,27,23+(t-1)*tint+add_lt,0,0)));
  end
else
  for t=1:len
    datev(t,np)=datenum(datevec(datenum(2019,08,27,21+(t-1)*tint+add_lt,0,0)));
  end
end
  lflx(1:len2,np)=lh(1:len2);
  sflx(1:len2,np)=sh(1:len2);
  tflx(1:len2,np)=lh(1:len2)+sh(1:len2);
  tco(1:len2,np)=cld_cover(1:len2);
  asst(1:len2,np)=sst(1:len2);

end % np

%%
mylabel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

  figure;
  set(gca,'fontsize',15)
  colormap(jet)
  subplot(2,1,1)
  hold all
% average to hourly
  colororder({'k','r'})
  yyaxis left
  
  plot(datev(240:60:len2,1),squeeze(tflx(240:60:len2,1)),'b','linewidth',2);
  plot(datev(240:60:len2,2),squeeze(tflx(240:60:len2,2)),'--b','linewidth',2);
  plot(datev(240:60:len2,3),squeeze(tflx(240:60:len2,3)),'k','linewidth',2);
  plot(datev(240:60:len2,4),squeeze(tflx(240:60:len2,4)),'--k','linewidth',2);
  ylim([0 65])
  ylabel ('Surface Flux (W/m^2)');
  
  yyaxis right
  plot(datev(240:60:len2,1),squeeze(asst(240:60:len2,1)-273.15),'r','linewidth',2);
  plot(datev(240:60:len2,2),squeeze(asst(240:60:len2,2)-273.15),'--r','linewidth',2);
  ylim([27.5 31.5])
  ylabel ('SST (^oC)');
  legend('L1','L1sst','L2','L2sst');
  datetick('x','hh')
  title(mylabel{1})
  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
  xlabel('Local Time (h)');

  subplot(2,1,2)
  hold all
  plot(datev(240:60:len2,1),squeeze(tco(240:60:len2,1)),'b','linewidth',2);
  plot(datev(240:60:len2,2),squeeze(tco(240:60:len2,2)),'--b','linewidth',2);
  plot(datev(240:60:len2,3),squeeze(tco(240:60:len2,3)),'k','linewidth',2);
  plot(datev(240:60:len2,4),squeeze(tco(240:60:len2,4)),'--k','linewidth',2);
  datetick('x','hh')
  
  title(mylabel{2})
  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
  ylabel ('Cloud Fraction','fontsize',15);
  xlabel('Local Time (h)');
  legend('L1','L1sst','L2','L2sst');
  fname=[outPath 'fig2_flx_cld_fraction.png'];
  print('-dpng',fname);

