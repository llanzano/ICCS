%%   User-fiendly code based on the following article:
%$£   
%$£  M. Oneto, L. Scipioni, M.J. Sarmento, I. Cainero, S. Pelicci, L. Furia, P.G. Pelicci, 
%$£  G.I. Dellino, P. Bianchini, M. Faretta, E. Gratton, A. Diaspro and L. Lanzanò
%$£  Nanoscale distribution of nuclear sites analyzed by superresolution STED-ICCS 
%$£  bioRxiv, 753228
%$£  
%$£  
%$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                              Luca Lanzanò @ Diaspro-Lab             %$£
%$£         Istituto Italiano di Tecnologia - Nanoscopy Department      %$£
%$£                      User-Friendly Version (Oct 2019)               %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£




function [ ]=SR_ICCS_simple(filename );

%open files ... 
[filename,pathname, filterindex] = uigetfile({'*.tif'},'Please select an image ');
filenamefull = [pathname, filename];   
A=simple_ICCS_readfiletif(filenamefull);     
T=size(A,3);
X=size(A,1);
Y=size(A,2);
xch1=A(:,:,1);

% figure
% for i=1:T
% subplot(2,3,i)
% imagesc(A(:,:,i))
% axis image
% end

if T < 2
[filename2,pathname2, filterindex] = uigetfile({'*.tif'},'Select 2nd channel? ');
    if filename2 ~= 0
    filenamefull2 = [pathname2, filename2];   
    A2=simple_ICCS_readfiletif(filenamefull2);     
    T=size(A,3);
    X=size(A,1);
    Y=size(A,2);
    xch2=A2(:,:,1);
    else
        xch2=xch1;
    end
else
    xch2=A(:,:,2);
end
x1=cat(3,xch1,xch2);

%load or define mask
maskch=1;
[filenamemask,pathnamemask, filterindex] = uigetfile({'*.tif'},'Select Mask');
if filenamemask ~= 0
filenamefull3 = [pathnamemask, filenamemask];   
Mask=simple_ICCS_readfiletif(filenamefull3);
% if size(Mask,3)>1
% prompt2 = {'Channel for Mask'}; 
% dlg_title2 = 'Select'; 
% num_lines = 1;
% def2 = {num2str(size(Mask,3))};
% answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
% figcheck=~isempty(answer2); 
% if figcheck==1
% maskch=str2num(answer2{1});
% end
% Mask=squeeze(Mask(:,:,maskch));
% end
 
    Thrmask=0;
%     figcheck=1;
%     while figcheck==1 
    B=simpleICCS_Threshold(Mask,Thrmask);
%     subplot(2,3,3+maskch)
%     imagesc(B)
%     axis image
%     prompt2 = {'Threshold value:'}; 
%     dlg_title2 = 'Input parameters'; 
%     num_lines = 1;
%     def2 = {num2str(Thrmask)};
%     answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
%     figcheck=~isempty(answer2); 
%     if figcheck==1
%     Thrmask=str2num(answer2{1});
%     end
%     end
Mask=B;   
% Mask(Mask>0)=1;
end

%set analysis area
Thr1=0;
Thr2=0;
sm=0.0;
Extra=0;
figcheck=1;
maxlag=ceil(X/4);
figure
while figcheck==1 
% B=simple_ICCS_Threshold(x1,Thr1,Thr2);
if filenamemask ~= 0
[x1pad,Np,Mask,B,Iav]=simple_PadForICS_fromMask(x1,Extra,Mask);
else
[x1pad,Np,Mask,B,Iav]=simple_PadForICS_sm(x1,Extra, [Thr1,Thr2], sm); 
end
b1=x1pad(:,:,1);
[m,n]=size(b1);
b2=x1pad(:,:,2);
subplot(2,2,1)
imagesc(b1)
subplot(2,2,2)
imagesc(b2)
subplot(2,2,3)
imagesc(Mask)

%ICS global analysis with padding G(0) correction
ACF1=simple_ICCS_CCFmean(b1,b1)*m*m*n*n./(m*n-Np(1))./(m*n-Np(1));
ACF2=simple_ICCS_CCFmean(b2,b2)*m*m*n*n./(m*n-Np(2))./(m*n-Np(2));
CCF=simple_ICCS_CCFmean(b1,b2)*m*m*n*n./(m*n-Np(1))./(m*n-Np(2));  %Np(1) and 2 are equal
maxlag=min(maxlag,length(ACF1));
subplot(2,2,4)
plot(ACF1(1:maxlag),'DisplayName','ACF1');hold all;plot(ACF2(1:maxlag),'DisplayName','ACF2');plot(CCF(1:maxlag),'DisplayName','CCF');hold off;

prompt2 = {'Threshold ch1:','Threshold ch2:','smooth','Extra pixels','Points to plot for ACF'}; 
dlg_title2 = 'Set area for analysis'; 
num_lines = 1;
def2 = {num2str(Thr1),num2str(Thr2),num2str(sm),num2str(Extra),num2str(maxlag)};
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck=~isempty(answer2); 
if figcheck==1
Thr1=str2num(answer2{1});
Thr2=str2num(answer2{2});
sm=str2num(answer2{3});
Extra=str2num(answer2{4});
maxlag=str2num(answer2{5});

end
end

%set parameters for fit
lag0=1; %first spatial lag to fit
lagfit=maxlag; %points to fit
FWHM1nm=132;
FWHM2nm=96;
px=0.04; %px size in um
tol=50;
% const=0.19; 
% const=simple_ICCS_Get_const(FWHM1nm/(px*1000),FWHM2nm/(px*1000)); 
figcheck2=1;
figure
subplot(2,2,1);
str = {' f=' , ' d='};
dim = [.30 .15 .7 .76 ];
an=annotation('textbox',dim,'String',str,'FitBoxToText','on');

while figcheck2==1 
subplot(2,2,1)
[param3, fit12, chisq12]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfit-1,CCF(lag0+1:lag0+lagfit),4,1,'Cross');
subplot(2,2,3)
[param1, fit1, chisq1]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfit-1,ACF1(lag0+1:lag0+lagfit),4,1,'auto 1');
subplot(2,2,4)
[param2, fit2, chisq2]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfit-1,ACF2(lag0+1:lag0+lagfit),4,1,'auto 2');
w1=param1(3);
w2=param2(3);
wcc=param3(3);
M1ccg = param3(2)/param2(2) ;  % these are the parameters for coloc global
M2ccg = param3(2)/param1(2) ;  % these are the parameters for coloc global
Size1=w1*px*1.18;
Size2=w2*px*1.18;
Size12=wcc*px*1.18;
wrefg=sqrt(w1*w2);
if chisq12/max(chisq1,chisq2)<tol &&  param3(3) > 0.5*wrefg && param3(3) < 2*wrefg && param3(1)>-0.2*abs(param3(2)) % conditions on global fit
%   M12ccg=sqrt(M1ccg*M2ccg);  
M12ccg=0.5*(M1ccg+M2ccg);    % arithmetic mean 
wrefdelta=sqrt(0.5*( w1^2+w2^2 ));
deltaw=wcc-wrefdelta;
dist=0;
if deltaw>0
const=simple_ICCS_Get_const(FWHM1nm/(px*1000),FWHM2nm/(px*1000));    
dist=sqrt(deltaw/const)*px*1000;
end
str = {[' f=',num2str(M12ccg,2)], [' d=', num2str(dist,3),' nm']};
set(an,'String',str);
else
  M12ccg=0; 
str = {[' f=',num2str(M12ccg,2)]};  
set(an,'String',str);
end

prompt2 = {'Min lag:','Max lag:','tolerance','Pixel size (um)','FWHM1(nm)','FWHM2(nm)'}; 
dlg_title2 = 'Set paramaters for fit'; 
num_lines = 1;
def2 = {num2str(lag0),num2str(lagfit),num2str(tol),num2str(px),num2str(FWHM1nm),num2str(FWHM2nm)}; 
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck2=~isempty(answer2); 
if figcheck2==1
lag0=str2num(answer2{1});
lagfit=str2num(answer2{2});
tol=str2num(answer2{3});
px=str2num(answer2{4});
FWHM1nm=str2num(answer2{5});
FWHM2nm=str2num(answer2{6});
end
end

Results(:,1)=lag0:lag0+lagfit-1;
Results(:,2)=CCF(lag0+1:lag0+lagfit);
Results(:,3)=fit12;
Results(:,4)=ACF1(lag0+1:lag0+lagfit);
Results(:,5)=fit1;
Results(:,6)=ACF2(lag0+1:lag0+lagfit);
Results(:,7)=fit2;


% ResultsPar(i,:)=[N1 w1 N2 w2 M1cc M2cc wcc Area(i) dens1 dens2];
% ResultsPar1000=round(ResultsPar*1000);

filenameout=filenamefull(1:end-4) ;

answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
%save data in Matlab
save([filenameout,'_ICCS.mat'] );
% writing files txt
dlmwrite([filenameout,'.txt'],Results,'delimiter',';','precision',4);
end

MaskP1=Mask;
MaskP2=Mask;
%set parameters for local
BW=1;
tolneg=0.2; %max negative offset %
Gmin=0.001;
Gmax=10*max(param1(2),param2(2));
M=10;  % points in x for local analysis (1 for global only)
Sizemax1=Size1*1.5;
Sizemax2=Size2*1.5;
figcheck3=1;
figure
while figcheck3==1 
Mask2=Mask;
wmax1=Sizemax1/1.18/px;
wmax2=Sizemax2/1.18/px;
%start local if N>1
Dloc=min(m,n);
Mlocal=zeros(m,n);
if M>1
vsampx=ones(1,M,'uint16');   
fos=0.3;  %oversampling ebtween local ROIs
%set size of local moving mask
Dloc=floor( m/(M-(M-1)*fos) );
if iseven(Dloc)
    Dloc=Dloc-1;
end
Dh=uint16((Dloc-1)/2);
vsampx(1)=uint16(1+(Dloc-1)/2);
vsampx(M)=uint16(m-(Dloc-1)/2);
step=uint16(floor((vsampx(M)-vsampx(1))/(M-1)));
for ii=2:M-1
    vsampx(ii)=vsampx(1)+step * (ii-1);
end
My=round(M*n/m);
vsampy=ones(1,My,'uint16');
vsampy(1)=uint16(1+(Dloc-1)/2);
vsampy(My)=uint16(n-(Dloc-1)/2);
stepy=uint16(floor((vsampy(My)-vsampy(1))/(My-1)));
for jj=2:My-1
    vsampy(jj)=vsampy(1)+stepy * (jj-1);
end

%local analysis at points  M x My
Mloc=zeros(M,My);
WW=zeros(M,My);
w1loc=zeros(M,My);
w2loc=zeros(M,My);
wccloc=zeros(M,My);
for iloc=1:M
    for jloc=1:My
        %define the local subimage from padded full image
        Aloc1=b1(vsampx(iloc)-Dh:vsampx(iloc)+Dh, vsampy(jloc)-Dh:vsampy(jloc)+Dh);
        Aloc2=b2(vsampx(iloc)-Dh:vsampx(iloc)+Dh, vsampy(jloc)-Dh:vsampy(jloc)+Dh);
        
        MaskLoc=Mask(vsampx(iloc)-Dh:vsampx(iloc)+Dh, vsampy(jloc)-Dh:vsampy(jloc)+Dh);
        MaskLoc1=MaskP1(vsampx(iloc)-Dh:vsampx(iloc)+Dh, vsampy(jloc)-Dh:vsampy(jloc)+Dh);
        MaskLoc2=MaskP2(vsampx(iloc)-Dh:vsampx(iloc)+Dh, vsampy(jloc)-Dh:vsampy(jloc)+Dh);
        MinPx=(BW-1);
        if ( numel(find(MaskLoc1>0))>MinPx && numel(find(MaskLoc2>0))>MinPx )
        WW(iloc,jloc)=1; 
        end
        
%         WW(iloc,jloc)=numel(find(MaskLoc>0))/numel(MaskLoc);  % weight based on number of significant pixels in Mask
        
        %perform analysis on local subimage and get parameters without plot
        %each time or with a collective plot
        % padding on local image
        [Aloc1,~,~,~,~]=simple_PadForICS_fromMask(Aloc1,0,MaskLoc);
        [Aloc2,~,~,~,~]=simple_PadForICS_fromMask(Aloc2,0,MaskLoc);
        ACFloc1=simple_ICCS_CCFmean(Aloc1,Aloc1);
        ACFloc2=simple_ICCS_CCFmean(Aloc2,Aloc2);
        CCFloc=simple_ICCS_CCFmean(Aloc1,Aloc2);
        disp=0;
        lagfitloc=min(lagfit, length(CCFloc)-lag0);
        [param1loc, ~,chiloc1]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfitloc-1,ACFloc1(lag0+1:lag0+lagfitloc),w1,disp,'');
        [param2loc, ~,chiloc2]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfitloc-1,ACFloc2(lag0+1:lag0+lagfitloc),w2,disp,'');
        [param3loc, ~,chiloc12]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfitloc-1,CCFloc(lag0+1:lag0+lagfitloc),wcc,disp,'');
        if abs(param1loc(2))>Gmin && abs(param1loc(2))<Gmax && abs(param2loc(2))>Gmin && abs(param2loc(2))<Gmax && abs(param3loc(2))<2*max([param1loc(2),param2loc(2)]) && param1loc(3)<wmax1 && param2loc(3)<wmax2  % 1st set of conditions on local fit
        w1loc(iloc,jloc)=param1loc(3);
        w2loc(iloc,jloc)=param2loc(3);
        wref=sqrt(param1loc(3)*param2loc(3));
%         N1=(1/param1(2))*Area(i)/(pi*param1(3)*param1(3));
%         N2=(1/param2(2))*Area(i)/(pi*param2(3)*param2(3));]
        M1cc = param3loc(2)/param2loc(2) ;  % these are the parameters for coloc
        M2cc = param3loc(2)/param1loc(2) ;  % these are the parameters for coloc
                %set some limits for reasonable wcc (to discard values from
                %bad fits) and check relative tolerance
            if  param3loc(3) > 0.5*wref && param3loc(3) < 2*wref && chiloc12/max(chiloc1,chiloc2)<tol  && param3loc(1)>-tolneg*abs(param3loc(2)) % 2nd set of conditions on local fit
%                 Mloc(iloc,jloc)=sqrt(abs(M1cc*M2cc));
                Mloc(iloc,jloc)=mean([M1cc,M2cc]);
                wccloc(iloc,jloc)=param3loc(3);
            end
        end

    end
end
%interpolate from MxN points to the entire image using interp2
[Xs,Ys] = meshgrid(double(vsampx),double(vsampy));
[Xq,Yq] = meshgrid(1:m,1:n);
% if BW==0
%  Vq = interp2(Xs,Ys,Mloc,Xq,Yq);   
% else
%  Vq = interp2(Xs,Ys,((exp(-(BW*(1-WW)/0.5).^2)-  exp(-(BW/0.5).^2))/(1-exp(-(BW/0.5).^2))).*Mloc,Xq,Yq); % formula with weight; weight is 1 for 100% pixels
% end
% Vq = interp2(Xs,Ys,Mloc,Xq,Yq);
Vq = interp2(Xs,Ys,WW.*Mloc,Xq,Yq);
Vq(isnan(Vq)) = 0;  % remove NaN
Mask(isnan(Vq)) = 0; % remove NaN

SizeMap1 = interp2(Xs,Ys,w1loc,Xq,Yq);
SizeMap1=1.18*px*SizeMap1.*Mask;
SizeMap1(isnan(SizeMap1)) = 0;  % remove NaN
SizeMap2 = interp2(Xs,Ys,w2loc,Xq,Yq);
SizeMap2=1.18*px*SizeMap2.*Mask;
SizeMap2(isnan(SizeMap2)) = 0;  % remove NaN
% MaskCh1=simpleICCS_Threshold(SizeMap1,[0,Sizemax]);
% MaskCh2=simpleICCS_Threshold(SizeMap2,[0,Sizemax]);
% Mask2=MaskCh1.*MaskCh2;
Mlocal=Vq.*Mask2;
Mav=mean(Mlocal(Mask2>0));
SizeAv1=mean(SizeMap1(Mask>0));
SizeAv2=mean(SizeMap2(Mask>0));


subplot(2,2,1);
imagesc(Mlocal);
colorbar 
axis image
title(['<Local cross>=',num2str(Mav,2)]);
subplot(2,2,3);
imagesc(SizeMap1);
colorbar 
axis image
title(['<Local size 1>=',num2str(SizeAv1,3),' um']);
subplot(2,2,4);
imagesc(SizeMap2);
colorbar 
axis image
title(['<Local size 2>=',num2str(SizeAv2,3),' um']);


% local test fit
pxtest=[floor(m/2) floor(n/2)];
itest=pxtest(1);  % pixel coord and mask for showing local fit
jtest=pxtest(2);
Aloc1test=b1(itest-Dh:itest+Dh,jtest-Dh:jtest+Dh);
Aloc2test=b2(itest-Dh:itest+Dh,jtest-Dh:jtest+Dh);
CCFloctest=simple_ICCS_CCFmean(Aloc1test,Aloc2test);
subplot(2,2,2);
lagfitloc=min(lagfit, length(CCFloctest)-lag0);
[param3loctest, ~, ~]=simple_ICCS_Fit_ICS_1Dsingle(lag0:lag0+lagfitloc-1,CCFloctest(lag0+1:lag0+lagfitloc),wcc,1,'local test');

end
%end local analysis

prompt2 = {'Points for local ICCS map', 'Max Size 1 (um)','Max Size 2 (um)', 'min G','max G','max neg offset','Min # Pixels in Local'}; 
dlg_title2 = 'Paramaters for Local ICCS'; 
num_lines = 1;
def2 = {num2str(M),num2str(Sizemax1),num2str(Sizemax2),num2str(Gmin),num2str(Gmax),num2str(tolneg),num2str(BW)};
answer3 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck3=~isempty(answer3); 
if figcheck3==1
M=str2num(answer3{1});
Sizemax1=str2num(answer3{2});
Sizemax2=str2num(answer3{3});
Gmin=str2num(answer3{4});
Gmax=str2num(answer3{5});
tolneg=str2num(answer3{6});
BW=str2num(answer3{7});
end
end


answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
%save data in Matlab
save([filenameout,'_locICCS.mat'] );
%export images (using Matlab)
A1=double(Mlocal);
MaxValImg=max(A1,[],'all');
MinValImg=min(A1,[],'all');
Aout=(A1-MinValImg)/(MaxValImg-MinValImg);
outputFileName = [filenameout, '_Local_min',num2str(MinValImg,2),'max',num2str(MaxValImg),'.tiff'];
delete outputFileName ;
imwrite(Aout, outputFileName);

end

end





%required funtions

function A=simple_ICCS_readfiletif(fname)
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end

end

function y=simpleICCS_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end

function [y,Np, varargout]=simple_PadForICS_fromMask(x1,Extra,Mask)
[m,n,p]=size(x1);
Mask=double(Mask);
%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
end
%% adding average on zeros
for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 || isnan(Mask(i,j))
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function [y,Np, varargout]=simple_PadForICS_sm(x1,Extra,Thr, sm)
[m,n,p]=size(x1);
for kk=1:p
x1s(:,:,kk)=simpleICCS_smooth_simple(x1(:,:,kk),sm,2);
end

%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
ys=y;
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
    ys(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1s(:,:,k);
end

%% adding average on zeros
Mask=ys(:,:,1);
Mask2=ys(:,:,2);
Mask(Mask<=Thr(1)& Mask2<=Thr(2))=0;
Mask(Mask2>Thr(2)| Mask>Thr(1) )=1;

for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function  Output=simple_ICCS_CCFmean(x1,x2)

NumberOfAngles=180;
[X,Y]=size(x1);
%ACF=conv2(x1,x2,'same');
F1=fft2(x1);
F2=fft2(x2);
ACF= F1.*conj(F2);
G=((sum(sum(x1)))*(sum(sum(x2)))/X/Y);
ACF= ifft2(ACF);
ACF= fftshift(ACF)./G-1;

[R, C]=size(ACF);
if iseven(R)
r0=R/2+1;
else
r0=(R+1)/2;
end
if iseven(C)
c0=C/2+1;
% Radius=min(R/2,C/2);
else
c0=(C+1)/2;
% Radius=min((R-1)/2,(C-1)/2);
end
Radius=min(r0-1,c0-1);

if NumberOfAngles==1
    Output=ACF(r0,c0:end);
else
ACF1=flipud(ACF(1:r0-1,c0:end));
ACF2=ACF(r0:end,c0:end);
ProfMat=zeros(NumberOfAngles*2,Radius);

for j=1:2
    if j==1
        y=ACF1';
    else
        y=ACF2;
    end
    
% CALCULATION OF ROTATIONAL MEAN
% Definition of angles
t=(pi/NumberOfAngles/2:pi/NumberOfAngles/2:pi/2);
   
% Matrix
y=y(1:Radius,1:Radius);
% Cycle between the 2nd and 2nd to last angles
[~, y1y]=size(y);

for i=1:NumberOfAngles
   rt=ceil(cos(t(i))*(1:Radius));
   ct=ceil(sin(t(i))*(1:Radius));
   profile=y((rt-1).*y1y+ct);

   if j==1
   ProfMat(NumberOfAngles+i,:)=profile;
   else
   ProfMat(i,:)=profile;
   end   
end

end


Output=[double(ACF(r0,c0)) sum(ProfMat)./(2*NumberOfAngles)];
% 
% OrientedProfiles=min_fw_Profile;
% OrientedProfiles(2,:)=max_fw_Profile;
% Angles=[min_th,max_th];

end


end

function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end

function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle(x,y,w0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*exp(-((x-0).*(x-0)./(Param(3)*Param(3)) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My w0]);
param(3)=abs(param(3));
fval=(param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) )));
chisq=sum( (( (param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) ))) -y ).^2)./((param(2)^2)  )) ;

if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));
end
end

function B=simpleICCS_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

end

function const=simple_ICCS_Get_const(FWHM1,FWHM2)

FWHMval=2:0.05:5;

D = [0.2175    0.2203    0.2010    0.1863    0.1755    0.1714    0.1549;
    0.2203    0.1972    0.1850    0.1815    0.1737    0.1545    0.1400;
    0.2010    0.1850    0.1788    0.1772    0.1675    0.1525    0.1435;
    0.1863    0.1815    0.1772    0.1671    0.1539    0.1452    0.1347;
    0.1755    0.1737    0.1675    0.1539    0.1512    0.1389    0.1316;
    0.1714    0.1545    0.1525    0.1452    0.1389    0.1372    0.1309;
    0.15489	  0.13998	0.14346	  0.13471	0.13159   0.13092	0.12654] ;

[Xs,Ys] = meshgrid(double(1:7),double(1:7));
[Xq,Yq] = meshgrid(1:0.1:7,1:0.1:7);
Dint = interp2(Xs,Ys,D,Xq,Yq);
if FWHM1<2
    pos1=1;
elseif FWHM1>5
    pos1=61;
else
    [~, pos1]=min(abs((FWHMval-FWHM1)));
end

if FWHM2<2
    pos2=1;
elseif FWHM2>5
    pos2=61;
else
    [~, pos2]=min(abs((FWHMval-FWHM2)));
end
    
const=Dint(pos1,pos2);

end
