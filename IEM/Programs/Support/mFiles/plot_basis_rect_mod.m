% plot_basis_rect.m

% adapted from plot_basis by TCS 10/25/13
%
% plot_basis_rect(b),n_rfX,n_rfY;

function plot_basis_rect_mod(b,n_rfX,n_rfY,resX,resY)

figure; clf; 

nr = n_rfY;
nc = n_rfX;

%res = sqrt(size(b,2));

 cidx = [2.5 3.5 4.5 5.5 9 10 11 12 13 15.5 16.6 17.5 18.5 19.5 20.5 22 23 24 25 26 27 28 29.5 30.5 31.5 32.5 33.5 34.5 37 38 39 40 41 44.5 45.5 46.4 47.5];
c=1;
for bb = 1:size(b,1);
    
    % goes down each column first
    
    subplot(nr,nc,cidx(c));
    %subplot(nr,nc,bb);
    imagesc(reshape(b(bb,:),resY,resX));
    set(gca,'XTick',[],'YTick',[]);
    axis equal;
    axis tight;
%     if ridx == nr
%         ridx = 1;
%         cidx = cidx+1;
%     else
%         ridx = ridx+1;
%     end
    c=c+1;
end

figure;
%imagesc(reshape(sum(b,1),resY,resX));
surf(reshape(sum(b,1),resY,resX));
%set(gca,'XTick',[],'YTick',[]); axis square;


return