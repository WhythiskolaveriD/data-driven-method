function [ FilteredTS, Recovered, sf, lhat] = Data_drivenApproach(F,col,v,Zbn,Id,frad,dstrp)
% This function implements the data driven approach by vishwakarma et al. 2016
% to repair the damage to the signal due to filtering GRACE products. 
% The details can be found in the relevant publication (doi:10.1002/2016WR018960). 
% Here one can choose either a Gaussian filter or a combination of Gaussian and
% destriping filter. The inputs and outputs are as follows:
%% Input
% F     - field (cell structure, with either stokes coeff or fields).
% col   - which column of the cell data contains field/coeff.
% v     - version of the field, v = 1(GMT in middle) or v = 2(GMT at left extreme).
% Zbn   - Basin map with ID as spatial location. (also available with the bundle)
% Id    - Vector of basin Id under investigation.
% frad  - Half width radius for Gaussian filter.
% dstrp - use destriping filter (1- yes, 0- No).

%% Output
% FilteredTS  - time series from filtered fields.
% Recovered   - corrected time series.
% sf          - scale factors.
% lhat        - estimated true leakage.
%--------------------------------------------------------------------------------
% EXAMPLE: ID = [295, 281, 35, 12]; % 4 catchments: Amazon, Congo, Ob, Yenisie
% [filGRACE,RecGRACE,sfGRACE,lhatGRACE] = Data_drivenApproach(GRACE_GFZcoeffSH,3,2,Zbn,ID,350,1);
% %% GRACE_GFZcoeffSH is a cell structure {157 x 3} with each row denoting
% %% an epoch, first column is for year, second is for month and third
% %% contains the data (SH coeff).
%
%--------------------------------------------------------------------------------
% uses: (Function - fromBundle)
%   gshs - SHBundle
%   gsha - SHBundle
%   resmplbasin - DataDrivenApproach
%   dstrpngmtrx - FilterBundle
%   cs2sc - SHBundle
%   gaussian - FilterBundle
%   globalsptconvgauss - DataDrivenApproach
%   Sptconvgaussbasin - DataDrivenApproach
%   Phase_calc - DataDrivenApproach
%   naninterp - DataDrivenApproach

%--------------------------------------------------------------------------------
% written by Bramha Dutt Vishwakarma, Institute of Geodesy, University of
% Stuttgart.      ----- 14 July 2015 ----
% Citation: Dutt Vishwakarma, B., B. Devaraju, and N. Sneeuw (2016), Minimizing the effects of filtering on catchment scale GRACE solutions, Water Resour. Res., 52, doi:10.1002/2016WR018960.
%--------------------------------------------------------------------------------
%% Checking the input and augmenting initial values
%-----------------------------------------------------
if nargin > 7, error('Too many input arguments!'),          end
if nargin == 6 || isempty(dstrp),        dstrp = 0;          end
if nargin == 5, error('Insufficient input arguments'),       end
%------------------------------------------------------------------------------------------------------------

%% get the size and dimensions of the input elements

[rf, ~] = size(F);% size of field cell array
[~,cid] = size(Id); % to get the number of basins to be analysed
cf = col;
field = cell(rf,1);
% get the resolution of global fields
[lat , long] = size(F{1,cf});

%% check if the given field is spherical harmonic coefficients or grided data, If coefficients, convert it to global grids
if long == lat
    deg = 1; %%%change HERE for different resolution
    for i=1:rf
        [field{i,1}, ~, ~]= gshs(F{i,cf},'water','cell',180/deg,0,0);
    end
    long = 360/deg;
else
    deg = 360/long;
    for i=1:rf
        field{i,1} = F{i,cf};
    end
end
%% Compute the Gaussian filter kernel
Weights = gaussian(90,frad);
filter = repmat(Weights,1,181);
%% Compute the grided Area matrix
[lambda,theta]=meshgrid(0:deg:(359+(1-deg)),0:deg:(179+(1-deg)));
[lambda1,theta1]=meshgrid(1*deg:deg:360,1*deg:deg:180);
Area=(6378.137^2)*(((pi/180)*lambda1)-((pi/180)*lambda)).*(sind(90-theta)-sind(90-theta1));
%% declaration of sizes of the vectors
bfn = cell(cid,1);
ibfn = cell(cid,1);
sf(1:cid,3)=0;
tsleaktotal(1:rf,1:cid)=0;
ftsleaktotal(1:rf,1:cid)=0;
tsleaktotalf(1:rf,1:cid)=0;
fftsleaktotal(1:rf,1:cid)=0;
lhat(1:rf,1:cid)=0;
%UFnoisyTS(1:rf,1:cid)=0;
Recovered(1:rf,1:cid)=0;
FilteredTS(1:rf,1:cid)=0;
leakage(1:rf,1:cid)=0;
leakager(1:rf,1:cid)=0;
%%------------------------------------------------------------------------------------------------------------
%% get the basin functions ready.
for mbasin = 1:cid
    [rbasin2deg,rbasin1deg,rbasinhalfdeg]=resmplbasin( Id(1,mbasin), Zbn );
    if deg==1
        Rb=rbasin1deg;
    elseif deg==2
        Rb=rbasin2deg;
    elseif deg==0.5
        Rb=rbasinhalfdeg;
    end
    
    bfn{mbasin,1}=Rb; % basin function (1 inside, 0 outside)
    ibfn{mbasin,1}=(1-Rb); % inverted basin function (1 outside, 0 inside)
end

%% compute filtered field, leakage for every month, and for every basin
for m=1:rf
    %% get a field, spatially consistent with the basin function
    if v==1
        Fld = field{m,1};
    elseif v==2
        Fld = [field{m,1}(:,(long/2)+1:end) field{m,1}(:,1:(long/2))];
    end
    
    %% Filter the field
    
    if dstrp == 1
        csF = gsha(Fld,'mean','block',long/2);
        [csFd,~] = dstrpngmtrx(csF(1:91,1:91),90);
        csF = cs2sc(csFd(1:91,1:91),0);
        frm = gshs(csF.*filter,'none','cell',long/2,0,0);
    else
        [ frm ] = globalsptconvgauss( Fld,frad,deg ); %  Filtering in spatial domain
        
        %cs = gsha(Fld,'mean','block',long/2);       % uncomment these lines for filtering in spectral domain
        %csF = cs2sc(cs(1:91,1:91),0);
        %frm = gshs(csF.*filter,'none','cell',long/2,0,0);
    end
    %% Filter the field again
    if dstrp == 1
        csF = gsha(frm,'mean','block',long/2);
        [csFd,~] = dstrpngmtrx(csF(1:91,1:91),90);
        csF = cs2sc(csFd(1:91,1:91),0);
        ffrm = gshs(csF.*filter,'none','cell',long/2,0,0);
    else
        [ ffrm ] = globalsptconvgauss( frm,frad,deg ); %  Filtering in spatial domain
        
        %cs = gsha(frm,'mean','block',long/2);       % uncomment these lines for filtering in spectral domain
        %csF = cs2sc(cs(1:91,1:91),0);
        %ffrm = gshs(csF.*filter,'none','cell',long/2,0,0);
    end
    %% compute the relevant time series for each catchment
    for i=1:cid
        
        id=Id(1,i);
        %% leakage from once filtered fields
        Fldi= ibfn{i,1}.*frm;
        if dstrp == 1
            csM = gsha(Fldi,'mean','block',long/2);
            [csMd,~] = dstrpngmtrx(csM(1:91,1:91),90);
            csF = cs2sc(csMd(1:91,1:91),0);
            Fldi = gshs(csF,'none','cell',long/2,0,0);
        end
        [ filF,~ ] = Sptconvgaussbasin(id,Zbn,Fldi,frad,deg); % spatial domain filtering
        tsleaktotal(m,i)=(sum(sum(filF.*bfn{i,1}.*Area))/sum(sum(bfn{i,1}.*Area)));
        %% leakage from twice filtered fields
        Fldi= ibfn{i,1}.*ffrm;
        if dstrp == 1
            csM = gsha(Fldi,'mean','block',long/2);
            [csMd,~] = dstrpngmtrx(csM(1:91,1:91),90);
            csF = cs2sc(csMd(1:91,1:91),0);
            Fldi = gshs(csF,'none','cell',long/2,0,0);
        end
        [ filfF,Area ] = Sptconvgaussbasin(id,Zbn,Fldi,frad,deg); % spatial domain filtering
        tsleaktotalf(m,i)=(sum(sum(filfF.*bfn{i,1}.*Area))/sum(sum(bfn{i,1}.*Area)));
        %% time series from filtered fields and that from unfiltered fields
        FilteredTS(m,i)=sum(sum(frm.*bfn{i,1}.*Area))/sum(sum(bfn{i,1}.*Area));
        %UFnoisyTS(m,i)=(sum(sum(Fld.*bfn{i,1}.*Area))/sum(sum(bfn{i,1}.*Area)));
    end
end

%% calculate the scale factor for every catchment

for i=1:cid
    %% compute the filtered basin function
    if dstrp == 1
        csrb = gsha(bfn{i,1},'mean','block',long/2);
        [dsfld,~] = dstrpngmtrx(csrb(1:91,1:91),90);
        cs = cs2sc(dsfld(1:91,1:91),0);
        fb = gshs(cs.*filter,'none','cell',long/2,0,0);
    else
        fb = globalsptconvgauss( bfn{i,1},frad,deg );
    end
    %% the scale factor:
    sf(i,1)=sum(sum(Area.*bfn{i,1}.*fb)); sf(i,2)=sum(sum(Area.*bfn{i,1}));
    sf(i,3)=sf(i,2)/sf(i,1);
    
    
end
%% Get the phase difference between the time series from once filtered fields and twice filtered fields
[ ps ] = Phase_calc(tsleaktotal,tsleaktotalf); %%%%%


%% Compute the near true leakage
for i=1:cid
    ftsleaktotal(:,i) = naninterp(tsleaktotal(:,i)); % replace the gaps (NaN values) with an interpolated value in leakage time series from once filtered fields
    fftsleaktotal(:,i) = naninterp(tsleaktotalf(:,i)); % replace the gaps (NaN values) with an interpolated value in leakage time series from twice filtered fields
    X = fft(ftsleaktotal(:,i)'); % take fast Fourier transform
    p = -(ps(i))/rf; % compute the fraction of the time period by which the time series is to be shiftes
    Y = exp(1i*pi*p*((0:rf-1)-rf/2)/rf); % compute the Conjugate-Symmetric shift
    Z = X.*Y; % Apply the shift
    a = ifft(Z); % take inverse fast Fourier transform and get back to the time domain
    con = conj(a);
    s = a+con;
    z = s/2;
    leakage(:,i) = z; % shifted time series
    %% shift the time series from once filtered fields in the direction of time series from twice filtered fields, to later compute the amplitude ratio
    p = (ps(i))/rf; % Fraction of a time period to shift data
    Y = exp(1i*pi*p*((0:rf-1)-rf/2)/rf); % Conjugate-Symmetric shift
    Z = X.*Y; % Apply the shift
    a = ifft(Z);
    con = conj(a);
    s = a+con;
    z = s/2;
    leakager(:,i) = z;
end
%% compute the ratio between the amplitude of the shifted leakage from once filtered fields and leakage from twice filtered fields
rfn=(leakage./fftsleaktotal);
rfn((rfn)>=2)=1;
rfn((rfn)<=-2)=-1;
rfn = sum(abs(rfn));
rfn=rfn/rf; %% amplitude ratio

%% get improved leakage estimates using leakage computed from filtered field.
for i=1:rf
    
    lhat(i,:)=leakager(i,:).*(rfn); %% apply the amplitude ratio to the shifted leakage time series from once filtered fields to get the near true leakage
    
end
lhat(isnan(FilteredTS))=nan; %% reintroduce the NaN for data gaps
leakremts=FilteredTS-(lhat); %% Remove near true leakage from time series from filtered fields

%% recover the near true time series
for i=1:rf
    Recovered(i,:)=leakremts(i,:).*sf(:,3)';
end

end

