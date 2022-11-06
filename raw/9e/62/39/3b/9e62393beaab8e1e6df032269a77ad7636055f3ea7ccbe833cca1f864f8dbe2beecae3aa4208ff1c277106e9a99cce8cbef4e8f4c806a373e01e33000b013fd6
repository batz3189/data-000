function [p,t]=plate(X1, X2, Y1, Y2, Z1, Z2, Size);
%STRIP_ Creates strip (or plate) mesh in an arbitrary plane (xy or xz or yz)
%
%   Essentially equivalent to the script plate.m of Chapter 2 or to the 
%   function plate.m from Chapter 6 but is more flexible 
%   
%   Copyright 2002 AEMM. Revision 2002/04/23 
%   AppendixA 

if nargin < 7
  error('function strip_ requires seven input arguments.')
end
warning off

%   W   Plate width (along the x-axis)
%   L   Plate length (along the y-axis)
%   Nx  Discretization parameter (width)
%   Ny  Discretization parameter (length)

SizeX=X2-X1;
SizeY=Y2-Y1;
SizeZ=Z2-Z1;

if (SizeZ==0) %plate mesh is in the xy-plane
    W=SizeX;
    L=SizeY;
    Nx=round(SizeX/Size);
    Ny=round(SizeY/Size);
end
if (SizeY==0) %plate mesh is in the xz-plane
    W=SizeX;
    L=SizeZ;
    Nx=round(SizeX/Size);
    Ny=round(SizeZ/Size);    
end
if (SizeX==0) %plate mesh is in the yz-plane
    W=SizeY;
    L=SizeZ;
    Nx=round(SizeY/Size);
    Ny=round(SizeZ/Size);
end

%Set vertexes
epsilon=1e-9;
M=1;
for i=1:Nx+1
    for j=1:Ny+1
        X(M)=-W/2+(i-1)/Nx*W;
        Y(M)=-L/2+(j-1)/Ny*L-epsilon*X(M);
        M=M+1;
    end
end

%Use Delaunay triangulation
TRI = delaunay(X,Y,0);
t=TRI';
t(4,:)=0;

%Scale plate
if (SizeZ==0) %plate mesh is in the xy-plane
    X=X1+SizeX*(X-min(X))/(max(X)-min(X));
    Y=Y1+SizeY*(Y-min(Y))/(max(Y)-min(Y));
    p=[X; Y; ones(1,length(X))*Z1];end
if (SizeY==0) %plate mesh is in the xz-plane
    X=X1+SizeX*(X-min(X))/(max(X)-min(X));
    Z=Z1+SizeZ*(Y-min(Y))/(max(Y)-min(Y));
    p=[X; ones(1,length(X))*Y1; Z];end
end
if (SizeX==0) %plate mesh is in the yz-plane
    Z=Z1+SizeZ*(Y-min(Y))/(max(Y)-min(Y));
    Y=Y1+SizeY*(X-min(X))/(max(X)-min(X));
    p=[ones(1,length(X))*X1; Y; Z];end
end

