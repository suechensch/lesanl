% function to calculate shallow cumulus cloud entrainment/detrainment rates
clear all
close all

for np=1:4
  clear all
  exp={'L1','L1sst','L2','L2sst'};
  outpath='./';
  if np < 3
   dtg='2019092721';
  else
   dtg='2019092723';
  end
nest=1;
yy=str2num(dtg(1:4));
mm=str2num(dtg(5:6));
dd=str2num(dtg(7:8));
hh=str2num(dtg(9:10));
[data_grid]=read_datahd_ff(dtg, kpath);
ma=data_grid.nest.nx;
na=data_grid.nest.ny;
kka=data_grid.nz;
mo=num2str(ma,'%04d');
no=num2str(na,'%04d');
Z=flipud(data_grid.sigm);
z2=num2str(data_grid.sigm(1),'%06d');
ztop=data_grid.ztop;
   zup=Z(2:kka);
   zdn=Z(1:kka-1);
   dz=(zup-zdn);
   for i=1:na
     for j=1:ma
       dz3d(i,j,:)=dz;
       z3d(i,j,:)=Z;
       zm3d(i,j,:)=0.5*(zup+zdn);
     end
   end

kka=data_grid.nz;
k1=1;
k2=199;
i1=50;
j1=50;
j2=220;
i2=220;

%% start process
t1=4; % begin model forecast hour
t2=20+55/60;  % end model forecast hour
tint=5/60; % hour
nt=int8((t2-t1)/tint+1);
fld(1:kka,1:nt)=0.0;
datev(1:nt)=0;
cld_base(1:na,1:ma)=0.;
cld_top(1:na,1:ma)=0.;
ic=0;
  cld_en_base(1:nt,1:na,1:ma)=nan;
  cld_de_base(1:nt,1:na,1:ma)=nan;
  cld_en_top(1:nt,1:na,1:ma)=nan;
  cld_de_top(1:nt,1:na,1:ma)=nan;
  entr_base(1:nt,1:na,1:ma)=nan;
  entr_top(1:nt,1:na,1:ma)=nan;
  detr_base(1:nt,1:na,1:ma)=nan;
  detr_top(1:nt,1:na,1:ma)=nan;
  mass(1:na,1:ma,1:kka)=nan;

for t=t1:tint:t2
  ic=ic+1;
  [fy,fm,fd,fh,fmin,fs]=datevec(datenum(yy,mm,dd,0+t,0,0));
  if t >=24
    fh=fh+24;
  end
  str2=[num2str(fh,'%04d') num2str(fmin,'%02d') num2str(fs,'%02d') ]
  datev(ic)=datenum(fy,fm,fd,fh,fmin,fs);
  fname=['cldmix_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];

  if exist([kpath fname]) == 0
   fname
   continue
  end

  cld=readff(kpath,fname,kka);

  fname=['eqvpot_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];
%  epot=readff(kpath,fname,kka);

  fname=['wwwind_sig_' num2str(ztop,'%06d') '_000000_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];

  if exist([kpath fname]) == 0
   continue
  end

  wwind=readff(kpath,fname,kka+1);
  wwind_up=wwind(:,:,2:kka+1);
  wwind_dn=wwind(:,:,1:kka);
  wwindm=0.5*(wwind_up+wwind_dn);

  fname=['thiliq_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];
%  thilq=readff(kpath,fname,kka);

  fname=['airtmp_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];

  if exist([kpath fname]) == 0
   fname
   continue
  end

  airt=readff(kpath,fname,kka);

  fname=['ttlprs_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];

  if exist([kpath fname]) == 0
   fname
   continue
  end

  pres=readff(kpath,fname,kka); % in unit of mb
  pres=pres*100.; % conver to pascale

  fname=['wvapor_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];

  if exist([kpath fname]) == 0
   continue
  end

  qv=readff(kpath,fname,kka);

  fname=['relhum_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];
  if exist([kpath fname]) == 0
   fname
   continue
  end

  rh=readff(kpath,fname,kka);
%% moist static energy temperature=moist static energy
%  Th=T+gz/cp+Lv*qv/cp;
%  T is the air temperature, g: gravity, z: height above surface
%  q: water vapor specific humidity
L=2.501*10^6 ; % specific latent heat of vaporization J/kg
cp=1005.7; % specific heat capacity Jkg-1K-1
grav= 9.80616;
rgas=287.04; % unit of Joule (kg*m-2*s-2)/(kg*K)

density=pres./(rgas*airt); % kg/m3=pascal(kg*m-1*s-2)/ (m-2*s-2)
 
%% moist static energy temperature K
%  thilq=airt+(grav*z3d/cp)+L*qts/cp;
  thilq=airt+(grav*z3d/cp)+L*qv/cp;
   thilq_up=thilq(:,:,2:kka);
   thilq_dn=thilq(:,:,1:kka-1);
   dthl=thilq_up-thilq_dn;
%%
cld_cut=0.0001; % cloud mixing ratio threshold to sample the cloudy region 
cld_cut_e=1e-9; % cloud mixing ratio threshold to sample the environment region 

   index=find(cld>cld_cut);
   index_e=find(cld<=cld_cut_e);
   cld_mask(1:na,1:ma,1:kka)=0;
   smask(1:na,1:ma)=0;
   scmask(ic,1:na,1:ma)=0;
   cld_mask(index)=1;


%% cloud base
   for i=1:na
     for j=1:ma
      index=find(cld(i,j,k1:k2)>=cld_cut);
      ilen=length(index);
      if ilen == 0
        cld_base(i,j)=nan;
        kbase(ic,i,j)=nan;
        ktop(ic,i,j)=nan;
        cld_top(i,j)=nan;
      else
        cld_base(i,j)=Z(index(1));
        cld_top(i,j)=Z(index(ilen));
        kbase(ic,i,j)=index(1);
        ktop(ic,i,j)=index(ilen);
      end
     end
   end 
  
% compute cloud core area fraction and mass flux
   index_u=find(cld_mask==1 & wwindm>0.0); % cloud updraft core
   index2=find(cld_mask==1 & wwindm<=0); % cloud downdraft area
   thilq_u(1:na,1:ma,1:kka)=nan;
% mass unit, kg/m3*m/s=kg m-2 s-1
%   mass(1:na,1:ma,1:kka)=density.*wwindm;
%   mass(1:na,1:ma,1:kka)=nan;
   thilq_u(index_u)=thilq(index_u);
     index3=find(cld_mask==1 & wwindm>0.0);
     index4=find(cld_mask==1 & wwindm<=0.1);
     cld_ucore=length(index3);
     env_cldcore=length(index4);
     clducore_frac=cld_ucore/(env_cldcore+cld_ucore); 
% mass unit kg*m-2*s-1
     mass(index3)=density(index3).*wwindm(index3)*clducore_frac; 

   mass_up=mass(:,:,2:kka);
   mass_dn=mass(:,:,1:kka-1);

   
   dmassdz=(mass_up-mass_dn)./dz3d;
   dmassdz(:,:,kka)=dmassdz(:,:,kka-1);

   thilq_u_up=thilq_u(:,:,2:kka);
   thilq_u_dn=thilq_u(:,:,1:kka-1);
   cld_depth(ic,:,:)=cld_top-cld_base;
% shallow cumulus cloud top, base, and depth sampling
   sindex=find(cld_top<2000 & cld_base >= 500 & (cld_top-cld_base)<= 200);
   smask(sindex)=1;
   slen(ic)=length(sindex);

% compute environment
   thilq_e(1:na,1:ma,1:kka)=nan;
   thilq_e(index_e)=thilq(index_e);

   entrainment(1:na,1:ma,1:kka)=nan;
   detrainment(1:na,1:ma,1:kka)=nan;
 
%% compute entrainment
  for k=1:kka
% k level mean environment moist static energy temperature
      thilq_e_bar(k)=nanmean(nanmean(thilq_e(:,:,k))); 
      thilq_u_bar(k)=nanmean(nanmean(thilq_u(:,:,k))); 
  end

    for i=1:na
     for j=1:ma
        kb=kbase(ic,i,j);
        kt=ktop(ic,i,j);
        if isnan(kb)==0 & kb > 0 & kt>kb
         for k=kb+1:kt
% k level mean environment moist static energy temperature
         thilq_e_bar(k)=nanmean(nanmean(thilq_e(:,:,k)));
         hcb(ic,i,j)=Z(kb);
          dthilq_u_dz(i,j,k)=(thilq_u(i,j,k)-thilq_u(i,j,kb))/(Z(k)-Z(k-1));
          dude2=thilq_u(i,j,k)-thilq_e_bar(k);
          entrainment(i,j,k)=-(dthilq_u_dz(i,j,k)/dude2);
          detrainment(i,j,k)=entrainment(i,j,k)-dmassdz(i,j,k)/mass(i,j,k);
         end
        end
      end
     end

% don't account for negative values

  entrainment2=entrainment;
  detrainment2=detrainment;
  entrainment(find(entrainment<0))=nan;
  detrainment(find(detrainment<0))=nan;

%% find cloud top and base entrainment and detrainment

for i=1:na
  for j=1:ma
  if  cld_top(i,j)<2000 && (cld_top(i,j)-cld_base(i,j))<= 200 && cld_base(i,j)>=500
    cld_en_base(ic,i,j)=entrainment(i,j,kbase(ic,i,j)+1);  
    cld_de_base(ic,i,j)=detrainment(i,j,kbase(ic,i,j)+1);  
    thilq_base(ic,i,j)=thilq(i,j,kbase(ic,i,j)+1);  
    wwindm_base(ic,i,j)=wwindm(i,j,kbase(ic,i,j)+1);  
    cld_en_top(ic,i,j)=entrainment(i,j,ktop(ic,i,j));  
    cld_de_top(ic,i,j)=detrainment(i,j,ktop(ic,i,j));  
    wwindm_top(ic,i,j)=wwindm(i,j,ktop(ic,i,j));  
    entr_base(ic,i,j)=en_rh(i,j,kbase(ic,i,j)+1);
    entr_top(ic,i,j)=en_rh(i,j,ktop(ic,i,j));
    detr_base(ic,i,j)=de_rh(i,j,kbase(ic,i,j)+1);
    detr_top(ic,i,j)=de_rh(i,j,ktop(ic,i,j));
  end

 end
end

for k=1:kka
  en2d=entrainment(:,:,k);
  de2d=detrainment(:,:,k);
  smentr(k,ic)=nanmean(nanmean(en2d(sindex)));
  sdentr(k,ic)=nanmean(nanmean(de2d(sindex)));
  mentr(k,ic)=nanmean(nanmean(entrainment(:,:,k)));
  mdetr(k,ic)=nanmean(nanmean(detrainment(:,:,k)));
  mentr2(k,ic)=nanmean(nanmean(entrainment2(:,:,k)));
  mdetr2(k,ic)=nanmean(nanmean(detrainment2(:,:,k)));
  dmassdz2(k,ic)=nanmean(nanmean(dmassdz(:,:,k)));
  mass2(k,ic)=nanmean(nanmean(mass(:,:,k)));
  massdz2(k,ic)=-nanmean(nanmean(dmassdz(:,:,k)./mass(:,:,k)));
  index_u=find(cld_mask(:,:,k)==1 & wwindm(:,:,k)>0);
  mass2d=mass(:,:,k);
  mmass(k,ic)=nanmean(mass2d(index_u));
  ww2d=wwindm(:,:,k);
  mww(k,ic)=nanmean(ww2d(index_u));
  mthilq_e(k,ic)=nanmean(nanmean(thilq_e(:,:,k)));
  mthilq_u(k,ic)=nanmean(nanmean(thilq_u(:,:,k)));
  mde_rh(k,ic)=nanmean(nanmean(de_rh(:,:,k)));
  men_rh(k,ic)=nanmean(nanmean(en_rh(:,:,k)));
end

mean_cld_base(ic)=mean(mean(cld_base));
mean_cld_top(ic)=mean(mean(cld_top));


end % end t

[nk,nt]=size(mentr) % time-height arrary of mean entrainment

% do hourly average
ic=0;
for i=1:nt/12
    ic=ic+1;
    ib=(i-1)*12+1;
    ie=ib+12-1;
    h_mentr(:,ic)=nanmean(mentr(:,ib:ie),2);
    h_mdetr(:,ic)=nanmean(mdetr(:,ib:ie),2);
    h_mentr2(:,ic)=nanmean(mentr2(:,ib:ie),2);
    h_mdetr2(:,ic)=nanmean(mdetr2(:,ib:ie),2);
    h_smentr(:,ic)=nanmean(smentr(:,ib:ie),2);
    h_smdetr(:,ic)=nanmean(sdentr(:,ib:ie),2);
    h_mdepth(ic,:,:)=nanmean(cld_depth(ib:ie,:,:),1);
    h_mcld_top(ic)=nanmean(cld_top(ib:ie));
end

% do vertical average of every 100 m from 0m to 3500m
 newz=[0:200:7000];
for k=2:length(newz)
 kindex=find(Z<=newz(k) & Z>newz(k-1));
 h_mentr_new(k,1:ic)=nanmean(h_mentr(kindex,1:ic),1);
 h_mdetr_new(k,1:ic)=nanmean(h_mdetr(kindex,1:ic),1);
 h_smentr_new(k,1:ic)=nanmean(h_smentr(kindex,1:ic),1);
 h_smdetr_new(k,1:ic)=nanmean(h_smdetr(kindex,1:ic),1);
end

GMT2LT=8; % conversion from GMT to local time
sunrise_lt=6; % model forecast hour + initial hour + GMT to LT 
sunset_lt=18; % model forecast hour + initial hour + GMT to LT
local_time=mod([t1:1:t2]+GMT2LT+hh,24);
local_tim(find(local_time == 0))= 24;
index_dayt=find(local_time>=sunrise_lt & local_time< sunset_lt);
index_night=find(local_time<sunrise_lt | local_time>= sunset_lt)

for ii=1:length(local_time)
  local_str{ii}=[num2str(local_time(ii),'%02d') 'LT'];
end

for ii=1:length(index_dayt)
  dayt_str{ii}=[num2str(local_time(index_dayt(ii)),'%02d') 'LT'];
end

for ii=1:length(index_night)
  night_str{ii}=[num2str(local_time(index_night(ii)),'%02d') 'LT'];
end


entr_bmean(np)=nanmean(nanmean(nanmean(entr_base)))
entr_tmean(np)=nanmean(nanmean(nanmean(entr_top)))
de_mean_base(np)=nanmean(nanmean(nanmean(detr_base)))
de_mean_top(np)=nanmean(nanmean(nanmean(detr_top)))

sname=[outpath exp{np} '_cldtb_ent_det.mat'];
save(sname,'cld_en_base','cld_de_base','cld_en_top','cld_de_top', ...
'local_time', 'h_mentr', 'h_mdetr','h_mentr2', 'h_mdetr2', 'h_mdepth', ...
'h_mcld_top', ...
'h_mentr_new', 'h_mdetr_new','Z','newz','h_smentr_new', ...
'h_smdetr_new','h_smentr', 'h_smdetr','scmask','slen','entr_base', ...
'entr_top','detr_base','detr_top','mde_rh','men_rh','entr_bmean', ...
'entr_tmean','de_mean_base','de_mean_top','hcb','dmassdz2','mass2','massdz2' );

end % np

