function As=SmoothScene(A,n)
% As=SmoothScene(A,n)
%
% n is horizontal extent in pixel coords.

%[~,M]=size(A);

%A0=A;

P=11; % make sure this is odd.
B=ones(P,P);

i0 = (P+1)/2;
j0 = (P+1)/2;

for i=1:P
    for j=1:P

        dist=sqrt((i-i0).^2 + (j-j0).^2);
        
        B(i,j)=exp(-(dist/n).^2);

    end;
end
B=B/sum(sum(B));

As=conv2(A,B);

P2=(P-1)/2;
As = As((P2+1):(end-P2),(P2+1):(end-P2) );

% back fill original values where it was nan'ed owing to smoothing
ib=find(isnan(As) & ~isnan(A));
As(ib)=A(ib);