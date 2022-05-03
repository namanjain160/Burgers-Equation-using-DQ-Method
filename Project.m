%Set domain
a=0;
b=1;
c=0;
d=1;

%Set number of grid points, parameters, number of iterations
N=80;
M=50;
alp=1;
mu=0.001;
iterations=10;

t=a+(b-a)*(1-cos(pi*([1:N]-1)/(N-1)))*0.5;
x=c+(d-c)*(1-cos(pi*([1:M]-1)/(M-1)))*0.5;

A=W(t, 1);
B=W(x, 2);

A_k = kron(eye(M-2), squeeze(A(1, 2:N, 2:N)));
B_k1 = kron(squeeze(B(1, 2:M-1, 2:M-1)), eye(N-1));
B_k2 = kron(squeeze(B(2, 2:M-1, 2:M-1)), eye(N-1));

U=zeros(N,M);

exact=zeros(N,M);

%save the exact solution
for n=1:N
    for m=1:M
        exact(n,m)=f(x(m), t(n), mu, alp);
    end
end

%assign the initial values
for m=1:M
    U(1,m)=f(x(m), t(1), mu, alp);
end

%assign boundary conditions
for n=1:N
    U(n,M)=f(x(M),t(n),mu, alp);
    U(n,1)=f(x(1), t(n), mu, alp);
end

for k = 1:iterations
    F=zeros(N,M);
    
    for n=2:N
        for m=2:M-1
            F(n,m)=-A(1,n,1)*U(1,m)+mu*(B(2,m,1)*U(n,1) + B(2,m,M)*U(n,M))-alp*U(n,m)*(B(1,m,1)*U(n,1) + B(1,m,M)*U(n,M));
        end
    end
    
    F=F(2:N,2:M-1);
    F=squeeze(transpose(reshape(F, 1, [])));
    
    Ck=A_k+alp*diag(vec(U(2:N,2:M-1)))*B_k1-mu*B_k2;
    
    V=Ck\F;
    U(2:N, 2:M-1)=reshape(V, N-1, M-2);
    
    disp(norm(abs(exact(N,:)-U(N,:))', 'inf') + "    " + norm(abs(exact(N,:)-U(N,:))', 2));
end