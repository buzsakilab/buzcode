function SpatBits=SpatialInfo(pf,occ)

%  USAGE
%  	SpatBits=CellsSpatialInfo(S,XS,YS,epoch)
%  	
%  compute for each cells the number of spatial bits per spikes
%  
%  INPUT:
%  	pf: placeField (In Hz per pixel)
%  	occ: Occupancy probability
%  OUTPUT:
%  	SpatBits: the spatial information (in bits/spikes)
%  
%  Adrien Peyrache 2008

pf = pf(:);
occ = occ(:);
occ = occ/sum(occ);
f = sum(occ.*pf);
pf = pf/f;
ix = pf~=0;
SB = (occ(ix).*pf(ix)).*log2(pf(ix));
SpatBits = sum(SB);

