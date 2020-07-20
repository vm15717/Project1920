function [S_x,S_z] =hamiltonmat(N)
    s =[(N-1)/2:-1:-(N-1)/2];
    s_max= max(s);
    S_z = diag(s);
    S_z(S_z>0)=1;
    S_z(S_z<0)=-1;
    A= ones(N);
    B=(triu(A)-triu(A,2));
    C=(tril(A)-tril(A,-2));
    S_x = abs(B-C);
end
