function basis_set=RotateBasis(xgrid,ygrid, degrees,rfSize,xx,yy)

[thgrid, rgrid] = cart2pol(xgrid,ygrid); % CHECK THIS FOR up/down
thgrid = thgrid + deg2rad(degrees); % NEED TO KNOW WHY ANGLE OFFSET? -TONY
[xgrid_adj, ygrid_adj] = pol2cart(thgrid,rgrid);
rfGridX = reshape(xgrid_adj,numel(xgrid_adj),1);rfGridY = reshape(ygrid_adj,numel(ygrid_adj),1);
basis_set = nan(size(xx,1),size(rfGridX,1)); % initialize
for bb = 1:size(basis_set,2)
    basis_set(:,bb) = make2dcos(xx,yy,rfGridX(bb),rfGridY(bb),2*rfSize*2.5166,7);
end
end