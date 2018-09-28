function [ filF ] = globalsptconvgauss( Fld,r,deg )
%UNTITLED Summary of this function goes here-
%This function convolves the gaussian kernel with field in spatial domain for whole globe.

%Input.---------------------------------------------------------------------------

%Fld- Field
%r- radius of gaussian kernel (in Km)
%deg -the size of grid of field (1, 2 or 0.5 -- for 1 degree, 2 degree and half degree)

%-------------------------------------------------------------------------------------
%output --
%filF- Filtered field
%-------------------------------------------------------------------------------------
%author- Bramha dutt vishwakarma , date: 03/nov/2014
%-------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------
%
tic
%Calculate Area:
[lambda,theta]=meshgrid(0:deg:(359+(1-deg)),0:deg:(179+(1-deg)));
[lambda1,theta1]=meshgrid(1*deg:deg:360,1*deg:deg:180);
Area=(6378.137^2)*(((pi/180)*lambda1)-((pi/180)*lambda)).*(sind(90-theta)-sind(90-theta1));


[row col]=size(Fld); % get size of the field


%create meshgrid for colatitide and longitude,
[long,lat]=meshgrid(1*deg:deg:360,1*deg:deg:180);


%now calculate inverse of the second moment of the gaussian distribution
a=log(2)/(1-cos(r/6378.137));
%declare the size of output field
filF(1:row,1:col)=0;

FldA=Fld.*Area;

%start loop to traverse each pixel on the globe.
for th=1*deg:deg:90
    %lm=1;
    
    
    for lm=1*deg:deg:180
        
        %-get cosine of the central angle between kernel and exery point on
        %globe
        cang=(cosd(th)*cosd(lat)+ sind(th)*sind(lat).*cosd(long-lm));
        
        %------------------------------------------------------------------
        %- calculate weight of gaussian kernel
        W=a*(exp(-a*(1-cang))/(2*pi*(1-exp(-2*a))));
        W=(W)/(max(max(W)));
        D=sum(sum(W.*Area));
        
        %------------------------------------------------------------------
        % transfer the weight to a point in opposite hemisphere
        W4=flipud(W);
        
        
        %------------------------------------------------------------------
        % transfer the weight on a latitude in opposite hemisphere but 180 degrees away in
        % terms of longitude
        W2=[W4(:,(col/2)+1:end) W4(:,1:col/2)];
        
        %------------------------------------------------------------------
        % transfer the weight on same latitude but 180 degrees away in
        % terms of longitude
        W3=[W(:,(col/2)+1:end) W(:,1:col/2)];
        DD=sum(sum(W2.*Area));
        
        %------------------------------------------------------------------
        %summation of weight times the field along with area average to
        %yield the value wich should be alloted to kernel.
        filF(th/deg,(lm/deg))=sum(sum(FldA.*W))/D;
        
        filF((180-th+deg)/deg,(lm/deg)+col/2)=sum(sum(FldA.*W2))/DD;
        
        filF((th)/deg,(lm/deg)+col/2)=sum(sum(FldA.*W3))/D;
        
        filF((180-th+deg)/deg,lm/deg)=sum(sum(FldA.*W4))/DD;
           
        
    end
end

toc

end

