function data=readff(path,ffile,nlevel)
%       PURPOSE
%
%       CALL
%               data=readff(fname)
%       INPUT
%               path,ffile = path to the flatfile (ffile)
%               nlevel = 1 for 2D, nlevel for vars on sigma levels
%       OUTPUT
%               data - flat file contents
%       USES

%       HISTORY
%               Version 1       S. Gabersek 08/22/07
%-----------------------------


fname=[path '/' ffile];
% extract dimensions
%index=strfind(ffile,'x')-4;
% COAMPS flat file name format is fixed; the multiply sign has
% index 28: ('dimx'x'dimy')
index=28;
dimx=str2num(ffile(index:index+3));
dimy=str2num(ffile(index+5:index+8));
clear A B
fid=fopen(fname);
for k=nlevel:-1:1
  A=fread(fid,dimx*dimy,'float32','ieee-be');
  B=reshape(A,dimx,dimy);
  data(:,:,k)=B';
  clear A B
end

fclose(fid);



