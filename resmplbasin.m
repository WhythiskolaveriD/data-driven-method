function [ b2deg,b1deg,bhalf ] = resmplbasin( id, Zbn )
%UNTITLED2 Summary of this function goes here
% This function creates mask of 1 and 2 degree grid resolution from 0.5 degree grid
% resolution mask file Zbn. It samples down the catchment mask to the resolution required
%
%--------------------------------------------------------------------------------
% written by Bramha Dutt Vishwakarma, Institute of Geodesy, University of
% Stuttgart. 17 Oct 2014

%--------------------------------------------------------------------------------

[lat,long]=size(Zbn);
d(1:lat,1:long)=0;
d(Zbn==id)=1;

if lat == 360
    bhalf = d;
else
    bhalf = zeros(360,720);
end
b2deg = zeros(90,180);
b1deg = zeros(180,360);
%
if lat==360
for r=1:90
    for c=1:180
        s=(sum(sum(d(4*(r-1)+1:4*(r-1)+4,4*(c-1)+1:4*(c-1)+4))));
        if s>=8
            b2deg(r,c)=1;
            
        end
    end
end
end

if lat==360
for r=1:180
    for c=1:360
        s=(sum(sum(d(2*(r-1)+1:2*(r-1)+2,2*(c-1)+1:2*(c-1)+2))));
        if s>=2
            b1deg(r,c)=1;
            
        end
    end
end
end

if lat==180
    b1deg = d;
elseif lat == 90
    b2deg = d;
end

end

