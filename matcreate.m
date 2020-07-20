function K=matcreate(N)
w= zeros(N);
[row, col] = size(w);
for i=1:row
    for j=1:col
        if rem(i+j,2)==0
            w(i,j)=0;
        else
            w(i,j)=1;
        end
    end
end
K=w;
end