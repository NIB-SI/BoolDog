function [ u ] = projection( u)
%Projectes input u into [0,1] such that output u fulfils 0<=u<=1

u=min(max(0,u),1);
end

