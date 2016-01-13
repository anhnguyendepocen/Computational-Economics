clc
close all; clear all;

global H F S nh nw nb n rep Tempvec thres ccrit lcrit zcrit

nh = 225; nw = 110; nb = 110; rep = 45; thres = input('The discrimination rate:');

n = sqrt(nh);

H = [zeros(1,nh-nw-nb), ones(1,nb), -1*ones(1,nw)]; % '1' = black, '-1' = white
H = reshape(H(randperm(numel(H))), n, n); % Initial random distribution of households

F = [1 1 1; 1 0 1; 1 1 1]; % Filter

ccrit = max(-3, 3-2*ceil(8*thres));
lcrit = max(-5, 5-2*ceil(8*thres));
zcrit = 8-2*ceil(8*thres);

m = 1;
figure
colormap(flipud(gray))
set(gcf,'color','w');
subplot(2,2,m);
imagesc(H)
title('Initial distribution')

for l=1:rep
    
    Tempvec = []; % Temporary storage of moving households
    Crit = zeros(n); % Matrix storing moving decisions
    S = conv2(H, F, 'same'); % The sum of the surrounding elements
    
    for i=1:n
        for j=1:n
            if ((i==1)||(i==n))&&((j==1)||(j==n)) % corners
                if ((H(i,j)==1)&&(S(i,j)<=ccrit))||((H(i,j)==-1)&&(S(i,j)>=-ccrit)) Crit(i,j) = 1; Tempvec = [Tempvec H(i,j)]; end;
            elseif (i==1)||(i==n)||(j==1)||(j==n) % lines
                if ((H(i,j)==1)&&(S(i,j)<=lcrit))||((H(i,j)==-1)&&(S(i,j)>=-lcrit)) Crit(i,j) = 1; Tempvec = [Tempvec H(i,j)]; end;
            else % centre
                if ((H(i,j)==1)&&(S(i,j)<=zcrit))||((H(i,j)==-1)&&(S(i,j)>=-zcrit)) Crit(i,j) = 1; Tempvec = [Tempvec H(i,j)]; end;
            end;
        end;
    end;
    
    Tempvec = Tempvec(randperm(numel(Tempvec)));
    k = 1;

    for i=1:n
        for j=1:n
            if (Crit(i,j)==1)
                H(i,j) = Tempvec(k);
                k = k+1;
            end;
        end;
    end;
    
    if mod(l,15)==0
        m = m+1;
        subplot(2,2,m);
        imagesc(H)
        title(['After ', num2str(l), ' periods'])
    end

end;
