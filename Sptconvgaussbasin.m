function [ filF,Area ] = Sptconvgaussbasin(Bid,Zbn,Fld,r,deg)

%UNTITLED Summary of this function goes here-- This function convolves the
%gaussian kernel with field in spatial domain for a basin.

%Input.---------------------------------------------------------------------------
%Bid- basin id
%Zbn - matrix denoting location of basins (from file)
%Fld- Field
%r- radius of gaussian kernel (in Km)
%deg -the size of grid of field (1, 2 or 0.5 -- for 1 degree, 2 degree and half degree)
%Area- Area matrix.
%-------------------------------------------------------------------------------------
%output --
%filF- Filtered field

%__ uses --resamplbasin(from leakage scripts)
%-------------------------------------------------------------------------------------
%author- Bramha dutt vishwakarma , date: 03/nov/2014
%-------------------------------------------------------------------------------------
%   Detailed explanation goes here
%---------------------------------------------------------------------------------------

%   Detailed explanation goes here--
%...........------------------------------------------------------------------------
%check if the field is a missing GRACE month, so that we do not have to
%compute convolution in that case since the data gaps can be simply
%allotted NaN

tic
    %Calculate Area:
    [lambda,theta]=meshgrid(0:deg:(359+(1-deg)),0:deg:(179+(1-deg)));
    [lambda1,theta1]=meshgrid(1*deg:deg:360,1*deg:deg:180);
    Area=(6378.137^2)*(((pi/180)*lambda1)-((pi/180)*lambda)).*(sind(90-theta)-sind(90-theta1));
    
    
    %---------------------------------------------------------------------------------------
    
if isnan(Fld(40:70,90:150))==1
    filF(1:180/deg,1:360/deg)=NaN;
else
    
    %   otherwise we generate the basin function for right grid size to be used further using
    %   resmplbasin function
    [rbasin2deg,rbasin1deg,rbasinhalfdeg]=resmplbasin( Bid, Zbn );
    if deg==1
        Rb=rbasin1deg;
    elseif deg==2
        Rb=rbasin2deg;
    elseif deg==0.5
        Rb=rbasinhalfdeg;
    end
    
    
    %------------------------------------------------------------------
    [i, j]=find(Rb==1); % find the pixels on which convolution is to be done.
    
    %------------------------------------------------------------------
    
    
    %------------------------------------------------------------------
    % get the size of field and create meshgrid for lat and long
    [row col]=size(Fld);
    [long,lat]=meshgrid(1*deg:deg:360,1*deg:deg:180);
    
    %------------------------------------------------------------------
    
    %now calculate inverse of the second moment of the gaussian distribution
    a=(log(2)/(1-cos(r/6378.137)));
    
    %------------------------------------------------------------------
    %declare the size of output field
    filF(1:row,1:col)=0;
    
    %------------------------------------------------------------------
    
    % generate loop to traverse each pixel of the basin
    for loop=1:size(i)
        th=i(loop,1)*deg;
        lm=j(loop,1)*deg;
        %-get cosine of the central angle between two coordinates
        cang=(cos(th*pi/180)*cos(lat*pi/180)+ sin(th*pi/180)*sin(lat*pi/180).*cos((long*pi/180)-(lm*pi/180)));
        %-weight of gaussian kernel at that point
        W=(a*(exp(-a*(1-cang)))/(2*pi*(1-exp(-2*a))));
        W=(W)/(max(max(W)));%/sum(sum(Area)))
        
        SumA=sum(sum(W.*Area));
        %summation of weight times the field along with area average
        filF(th/deg,(lm/deg))=sum(sum(Fld.*W.*Area))/SumA;
        %filF(th/deg,(lm/deg))=filF(th/deg,(lm/deg))/(sum(sum(W.*Area))/sum(sum(Area)));
        
    end
end
toc
end


