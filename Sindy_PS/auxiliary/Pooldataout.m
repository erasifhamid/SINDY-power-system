function yout = Pooldataout(yin,polyorder,delta,gamma,K)
nVars=size(yin,2);
n = size(yin,1);
ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

% poly order 2
if(polyorder>=2)
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

% poly order 3
if(polyorder>=3)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end
 for i = 1:nVars
     for j=1:nVars 
         if i==j
             continue
         else
        yout(:,ind) = K(i,j) * sin(delta(:,i) - delta(:,j) - gamma(i,j));
        ind = ind+1;
         end
     end
 end
 
%   yout = [yout1 yout2 yout3]; 
% poly order 1
% for i=1:nVars
%     yout(:,ind) = yin(:,i);
%     ind = ind+1;
% end



% if(polyorder>=4)
%     % poly order 4
%     for i=1:nVars
%         for j=i:nVars
%             for k=j:nVars
%                 for l=k:nVars
%                     yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
%                     ind = ind+1;
%                 end
%             end
%         end
%     end
% end

% if(polyorder>=5)
%     % poly order 5
%     for i=1:nVars
%         for j=i:nVars
%             for k=j:nVars
%                 for l=k:nVars
%                     for m=l:nVars
%                         yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
%                         ind = ind+1;
%                     end
%                 end
%             end
%         end
%     end
% end

% if(usesine)
%     for k=1:10;

%         yout = [yout yout2 sin(yin)];
%     end
% end