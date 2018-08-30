function DATA = readECGsim(dirname)

% modeldir = [dirname  '/model/'];
% ecgdir = [dirname  '/ecgs/'];
% abeatsdir = [dirname  '/atrial_beats/beat1/'];
% vbeatsdir = [dirname  '/ventricular_beats/beat1/'];

modeldir = [dirname  '\model\'];
ecgdir = [dirname  '\ecgs\'];
abeatsdir = [dirname  '\atrial_beats\beat1\'];
vbeatsdir = [dirname  '\ventricular_beats\beat1\'];

[DATA.ATRIA.geom.VER,DATA.ATRIA.geom.ITRI] = loadtri([modeldir  'atria.tri']);
DATA.ATRIA.ADJsurf = loadmat([modeldir  'atria.adj2d']);
DATA.ATRIA.ADJ3D = loadmat([modeldir  'atria.adj3d']);
DATA.ATRIA.DISTsurf = loadmat([modeldir  'atria.dst2d']);
DATA.ATRIA.DIST3D = loadmat([modeldir  'ATRIA.dst3d']);
DATA.ATRIA.ADJANIS = loadmat([modeldir  'ATRIA.adjanis']);
DATA.ATRIA.DISTANIS = loadmat([modeldir  'ATRIA.dstanis']);

DATA.ATRIA.THORAX = loadmat([modeldir  'atria2Thorax.mat']);
DATA.ATRIA.ATRIA = loadmat([modeldir  'atria2atria.mat']);
DATA.ATRIA.VENTRICLES = loadmat([modeldir  'atria2Ventricles.mat']);
DATA.ATRIA.RLUNG = loadmat([modeldir  'atria2RLung.mat']);
DATA.ATRIA.LLUNG = loadmat([modeldir  'atria2LLung.mat']);
DATA.ATRIA.LCAV = loadmat([modeldir  'atria2LCavity.mat']);
DATA.ATRIA.RCAV = loadmat([modeldir  'atria2RCavity.mat']);
DATA.ATRIA.lead12 = loadmat([modeldir  'atria2standard12lead.mat']);
if length(DATA.ATRIA.lead12 ) == 1
    DATA.ATRIA.lead12 = loadmat([modeldir  'atria2standard_12.mat']);
end

[DATA.VENTR.geom.VER,DATA.VENTR.geom.ITRI] = loadtri([modeldir  'ventricle.tri']);
DATA.VENTR.ADJsurf = loadmat([modeldir  'ventricle.adj2d']);
DATA.VENTR.ADJ3D = loadmat([modeldir  'ventricle.adj3d']);
DATA.VENTR.DISTsurf = loadmat([modeldir  'ventricle.dst2d']);
DATA.VENTR.DIST3D = loadmat([modeldir  'ventricle.dst3d']);
DATA.VENTR.ADJANIS = loadmat([modeldir  'ventricle.adjanis']);
DATA.VENTR.DISTANIS = loadmat([modeldir  'ventricle.dstanis']);
DATA.VENTR.THORAX = loadmat([modeldir  'ventricles2Thorax.mat']);
DATA.VENTR.ATRIA = loadmat([modeldir  'ventricles2atria.mat']);
DATA.VENTR.VENTRICLES = loadmat([modeldir  'ventricles2Ventricles.mat']);
DATA.VENTR.RLUNG = loadmat([modeldir  'ventricles2RLung.mat']);
DATA.VENTR.LLUNG = loadmat([modeldir  'ventricles2LLung.mat']);
DATA.VENTR.LCAV = loadmat([modeldir  'ventricles2LCavity.mat']);
DATA.VENTR.RCAV = loadmat([modeldir  'ventricles2RCavity.mat']);
DATA.VENTR.VENTR212lead = loadmat([modeldir  'ventricles2standard12lead.mat']);
if length(DATA.VENTR.VENTR212lead) == 1
    DATA.VENTR.VENTR212lead = loadmat([modeldir  'ventricles2standard_12.mat']);
end



[DATA.GEOM.atria.VER DATA.GEOM.atria.ITRI]     = loadtri([modeldir  'atria.tri']);
[DATA.GEOM.ventr.VER,DATA.GEOM.ventr.ITRI]     = loadtri([modeldir  'ventricle.tri']);
[DATA.GEOM.lcav.VER, DATA.GEOM.lcav.ITRI]      = loadtri([modeldir  'lcav.tri']);
[DATA.GEOM.rcav.VER, DATA.GEOM.rcav.ITRI]      = loadtri([modeldir  'rcav.tri']);
[DATA.GEOM.llung.VER,DATA.GEOM.llung.ITRI]     = loadtri([modeldir  'llung.tri']);
[DATA.GEOM.rlung.VER,DATA.GEOM.rlung.ITRI]     = loadtri([modeldir  'rlung.tri']);
[DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI]   = loadtri([modeldir  'thorax.tri']);


DATA.GEOM.ventr.endoVER=zeros(1,length(DATA.GEOM.ventr.VER));
for i=1:length(DATA.GEOM.ventr.VER)
    if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
           DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
           DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
        DATA.GEOM.ventr.endoVER(i) = 2;
    end
end
for i=1:length(DATA.GEOM.ventr.VER)
    if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
           DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
           DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
        DATA.GEOM.ventr.endoVER(i) = 1;
    end
end

DATA.GEOM.atria.endoVER=zeros(1,length(DATA.GEOM.atria.VER));
for i=1:length(DATA.GEOM.atria.VER)
    if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
           DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
           DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
        DATA.GEOM.atria.endoVER(i) = 2;
    end
end
for i=1:length(DATA.GEOM.atria.VER)
    if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
           DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
           DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
        DATA.GEOM.atria.endoVER(i) = 1;
    end
end



ecgs=dir([ecgdir '*.refECG']);
for i=1:length(ecgs)
    a=['DATA.ECG.' ecgs(i).name(1:end-7) ];%' = loadmat(' [ecgdir ecgs(i).name] ');'];
    a=strrep(a,'(','');a=strrep(a,')','');
    a=strrep(a,' ','_');
    eval([a ' = loadmat(''' [ecgdir ecgs(i).name] ''');']);
    eval([a 'elec = loadmat(''' [ecgdir ecgs(i).name(1:end-6)] 'elec' ''');']);
end

beats=dir([abeatsdir 'user.*']);
for i=1:length(beats)
    a=['DATA.ABEAT.' beats(i).name(6:end) ];%' = loadmat(' [ecgdir ecgs(i).name] ');'];
    a=strrep(a,'(','');a=strrep(a,')','');
    eval([a ' = loadmat(''' [abeatsdir beats(i).name] ''');']);
end

beats=dir([vbeatsdir 'user.*']);
for i=1:length(beats)
    a=['DATA.VBEAT.' beats(i).name(6:end) ];%' = loadmat(' [ecgdir ecgs(i).name] ');'];
    a=strrep(a,'(','');a=strrep(a,')','');
    eval([a ' = loadmat(''' [vbeatsdir beats(i).name] ''');']);
end

if isfield(DATA.ECG,'standard12leadelec')
    ECGelec = DATA.ECG.standard12leadelec(:,2:4);
else
    ECGelec = DATA.ECG.standard_12elec(:,2:4);
end
    
index = zeros(length(ECGelec),1);
for i=1:length(ECGelec)
    len  = DATA.GEOM.thorax.VER - ones(length(DATA.GEOM.thorax.VER),1)*ECGelec(i,:);
    len = sqrt(sum(len.^2,2));
    a= find ( len == min(len));
%     if i==2
%        len = DATA.GEOM.thorax.VER(:,2);
%        a= find(len==max(len));
%     elseif i==1
%        len = DATA.GEOM.thorax.VER(:,2);
%        a= find(len==min(len));
%     else
        if ~isempty(a) && min(len) < 20
        index(i) = a;
    end
end

wct = index(1:3);

DATA.GEOM.wct = wct;
DATA.GEOM.stand12leadsIndex=index;

if size(DATA.ATRIA.THORAX,1)>1
    Awct = calcAwct(DATA.ATRIA.THORAX,wct);
    DATA.ATRIA.THORAX = doWCT(DATA.ATRIA.THORAX, Awct);
    DATA.ATRIA.ATRIA = doWCT(DATA.ATRIA.ATRIA,   Awct);
    DATA.ATRIA.VENTRICLES = doWCT(DATA.ATRIA.VENTRICLES,   Awct);
    DATA.ATRIA.RLUNG = doWCT(DATA.ATRIA.RLUNG,   Awct);
    DATA.ATRIA.LLUNG = doWCT(DATA.ATRIA.LLUNG,   Awct);
    DATA.ATRIA.RCAV = doWCT(DATA.ATRIA.RCAV,   Awct);
    DATA.ATRIA.LCAV = doWCT(DATA.ATRIA.LCAV,   Awct);
end
if size(DATA.VENTR.THORAX,1)>1
    Awct = calcAwct(DATA.VENTR.THORAX,wct);
    DATA.VENTR.THORAX = doWCT(DATA.VENTR.THORAX, Awct);
    DATA.VENTR.ATRIA = doWCT(DATA.VENTR.ATRIA,   Awct);
    DATA.VENTR.VENTRICLES = doWCT(DATA.VENTR.VENTRICLES,   Awct);
    DATA.VENTR.RLUNG = doWCT(DATA.VENTR.RLUNG,   Awct);
    DATA.VENTR.LLUNG = doWCT(DATA.VENTR.LLUNG,   Awct);
    DATA.VENTR.RCAV = doWCT(DATA.VENTR.RCAV,   Awct);
    DATA.VENTR.LCAV = doWCT(DATA.VENTR.LCAV,   Awct);
end


% DATA.VBEAT.SPECS.repslope = DATA.VBEAT.repslope(1);
% DATA.VBEAT.SPECS.repCorrection = 0;
% DATA.VBEAT.SPECS.restpot = DATA.VBEAT.rest;
% DATA.VBEAT.SPECS.ampl = DATA.VBEAT.ampl;
% DATA.VBEAT.SPECS.plateauslope = DATA.VBEAT.platslope(1);
% DATA.VBEAT.SPECS.initialSlope = 0;
% DATA.VBEAT.SPECS.useCumsum = 0;
% DATA.VBEAT.SPECS.useCumsum = 0;
% DATA.VBEAT.SPECS.depSlope = 2;
% 
% DATA.ABEAT.SPECS.repslope = DATA.ABEAT.repslope(1);
% DATA.ABEAT.SPECS.repCorrection = 0;
% DATA.ABEAT.SPECS.restpot = DATA.ABEAT.rest;
% DATA.ABEAT.SPECS.ampl = DATA.ABEAT.ampl;
% DATA.ABEAT.SPECS.plateauslope = DATA.ABEAT.platslope(1);
% DATA.ABEAT.SPECS.initialSlope = 0;
% DATA.ABEAT.SPECS.useCumsum = 0;
% DATA.ABEAT.SPECS.depSlope = 2;


%%=======================================================================

function A = doWCT(Ain,Awct)

if ~isempty(Ain)
    A = Ain - ones(size(Ain,1),1) * Awct;
else
    A= Ain;
end

function Awct = calcAwct(AthorsoIn,wct)

Awct = mean(AthorsoIn(wct,:));

%%=======================================================================
function [M, extraresult]=loadmat(name)

% LOADMAT	Load a MFBF matrix file.
%
%		Usage: m         = loadmat('file');
%		   or  [m,extra] = loadmat('file');
%
%		LOADMAT('file') returns the matrix stored in 'file' and
%		the extra information stored at the bottom of that file.
%		LOADMAT works for binary as well as asci matrix files.
%
%		See also SAVEMAT.
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands
% 20060816; echo (S) switched off

f=fopen(name);
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  M=0;
  extraresult='';
  return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
  fclose(f);
  f=fopen(name);
  [magic ,nr]=fread(f,8,'char');
  if (char(magic')==';;mbfmat')
    fread(f,1,'char');
    hs=fread(f,1,'long');
    fread(f,1,'char');
    fread(f,1,'char');
    fread(f,1,'char');
    N=fread(f,2,'long');
    M=fread(f,[N(2),N(1)],'double');
  else
    fclose(f);
    f=fopen(name);
    N=fread(f,2,'long');
    M=fread(f,[N(2),N(1)],'float');    
  end
else
  M=fscanf(f,'%f',[N(2) N(1)]);
end
[extra,nextra]=fread(f,1000,'char');
fclose(f);
extra = char(extra);

if ~all(isspace(extra))
%   S=sprintf('%s contains the following extra information:\n', name);
%   disp(S);
%   disp(extra');
else
    extra=[];
end
M=M';
extraresult=extra;

%%========================================================================
function [pnt, dhk] = loadtri(fn)

fid = fopen(fn, 'rt');
if fid~=-1

  % read the vertex points
  Npnt = fscanf(fid, '%d', 1);
  pnt  = fscanf(fid, '%f', [4, Npnt]);
  pnt  = pnt(2:4,:)';

  % if present, read the triangles
  if (~(feof(fid)))
    [Ndhk, count] = fscanf(fid, '%d', 1);
    if (count ~= 0)
      dhk = fscanf(fid, '%d', [4, Ndhk]);
      dhk = dhk(2:4,:)';
    end
  else
    dhk = [];
  end
  fclose(fid);

else
  error('unable to open file');
end