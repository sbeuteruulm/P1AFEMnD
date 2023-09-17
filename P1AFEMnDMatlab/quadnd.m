function [lambda, weights] = quadnd(nD,order)
max_order = 5;
if order > max_order
  disp(['Warning: Gauss type rule of order ',num2str(order),...
        ' not supported. '])
  disp(['Using Gauss quadrature of order ',num2str(max_order),'.'])  
end
%*** Table of quadrature formulae on triangle of different order
wdata = cell(3,1);
if order <= 1       % Reference: Grundmann&Moeller, Stroud (1971)
  wdata{1} = 1;
elseif  order <= 2  % Reference: Stroud (1971) 
  xdata{2} =  1.0/(nD+1) * (1.0-1.0/sqrt(nD+2));    
  wdata{2} =  1.0/(nD+1);
elseif  order <= 3  % Reference: Grundmann&Moeller, Stroud (1971)
  wdata{1} = -0.25 * ( nD + 1 ) * ( nD + 1 ) / ( nD + 2 );
  xdata{2} =  1.0 / ( nD + 3 );    
  wdata{2} =  0.25 * ( nD + 3 ) * ( nD + 3 ) / ( ( nD + 1 ) * ( nD + 2 ) );
elseif order <= 5  % Reference: Stroud (1969) 
  wdata{1} = (39*nD^2-307*nD+626)*(nD+1)^5/(72*prod(nD+(1:5)));
  xdata{2} = [((nD+4)-sqrt(15.0))/((nD+8)*nD+1); ...
              ((nD+4)+sqrt(15.0))/((nD+8)*nD+1)];
  wdata{2} = [(585-151*sqrt(15))*(93-7*nD+16*sqrt(15))*(nD+4+sqrt(15))^5; ...
              (585+151*sqrt(15))*(93-7*nD-16*sqrt(15))*(nD+4-sqrt(15))^5] ...
              ./(7560*prod(nD+(1:5)));
  xdata{3} =  [(nD+7+2*sqrt(15))/((nD+14)*nD-11); ...
               (nD+7-2*sqrt(15))/((nD+14)*nD-11)];
  wdata{3} =  [(585+151*sqrt(15))*(nD+7-2*sqrt(15))^5; ...
               (585-151*sqrt(15))*(nD+7+2*sqrt(15))^5]./(1080*prod(nD+(1:5)));
end
%*** Generate quadrature points and weights from data above
lambda = zeros(0,nD+1); weights=[];
if ~isempty(wdata{1})
  lambda = ones(1,nD+1)/(nD+1); 
  weights = wdata{1}; 
end
if ~isempty(wdata{2})
   perm = ones(nD+1)+eye(nD+1);
   points = [xdata{2},1-nD*xdata{2}];
   for i=1:nD+1
     lambda = [lambda;points(:,perm(i,:))]; 
     weights = [weights;wdata{2}]; 
   end
end
if ~isempty(wdata{3})
  idx = reshape([(1:nD*(nD+1)/2)'*[1,1],nchoosek(1:nD+1,2)],[],2);
  perm = 1 + accumarray(idx,1);
  points = [xdata{3},(1-(nD-1)*xdata{3})/2];
  for j=1:size(perm,1)
    lambda = [lambda;points(:,perm(j,:))]; 
    weights = [weights;wdata{3}]; 
  end
end
