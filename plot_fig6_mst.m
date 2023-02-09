% matlab script to:
% 1. read COAMPS atmos sig filed
% 2. interpolate to the obs location
% 3. plot time-height cross section of the fields
clear
close all

mylabel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
exps={'L1','L2','L1sst','L2sst'};
nexp=length(exps);

tint=5/60;

figure
icc=0;
for np1=1:nexp
    icc=icc+1;
    subplot(2,2,icc)
    hold all;
    if strcmp(exps{np1},'L2') == 1
      dtg='2019092723';
    elseif strcmp(exps{np1},'L2sst') == 1
      dtg='2019092723';
    elseif strcmp(exps{np1},'L1') == 1
      dtg='2019092721';
    elseif strcmp(exps{np1},'L1sst') == 1
      dtg='2019092721';
    end

   kpath=['../data/' exps{np1} '/'];
   sfile=[kpath exps{np1} '_mse_cld_entr.mat'];
   load(sfile)

   yy=str2num(dtg(1:4));
   mm=str2num(dtg(5:6));
   dd=str2num(dtg(7:8));
   hh=str2num(dtg(9:10))

   ntime=size(mthilq_u,2);
   nt_hr=int8(ntime/12);
   hr_mthilq_u(1:kka,1:nt_hr)=0;
   ic2=0;
   t1=1;
   t2=ntime;
   for i=t1:12:t2
     ic2=ic2+1;
     if ic2 > ntime
      break;
     end
     ii=i+11;
     if ii > ntime
       break
     end
     hr_mthilq_u(:,ic2)=nanmean(mthilq_u(:,i:ii),2);  
     for k=1:kka
       z2d(k,ic2)=Z(k)/1000.;
       datev2d(k,ic2)=datenum(yy,mm,dd,hh+GMT2LT+ic2-1+4,0,0);
     end
   end % end i

   k1=1;
   k2=125;
   colormap(jet)
   pcolor(datev2d(k1:k2,:),z2d(k1:k2,:),hr_mthilq_u(k1:k2,:))
   shading interp
   caxis([338 356]);
   colorbar
   datetick('x','hh')
   axis tight
   box on
   titlestring=[mylabel{np1} ' ' exps{np1} ' updraft moist static energy'];
   title(titlestring);
   ax = gca;
   ax.TitleHorizontalAlignment = 'left';
   xlabel('Local Time (hr)');
   ylabel('Height (km)');

end % end np1

fpath=['../figures/'];
fname=[fpath 'fig6_mthilq_u.png']
print('-dpng',fname)


    


