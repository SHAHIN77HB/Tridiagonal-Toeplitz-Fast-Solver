function x = FastSolver (alpha,beta,gamma,n,b)
%% Fast solvers for tridiagonal Toeplitz linear systems

% This Function uses the idea of article " Fast solvers for tridiagonal
% Toeplitz linear systems " by Zhongyun Liu , Shan Li , Yi Yin and Yulin Zhang
% DOI : https://doi.org/10.1007/s40314-020-01369-3
% For solving linear system of equations with tridiagonal Toeplitz coefficient matrix.
%% Inputs : 
%           alpha , beta and gamma are the digonal elements of tridigonal
%           Toeplitz matrix A = [ alpha gamma
%                                 beta  alpha gammma
%                                 0     beta  alpha gammma ...
%                                .       .     .     .    .  .   .
%                                                                   ]
%           n : matrix dimention
%           b : right handside vector 

%This MATLAB function is written by Shahin Hasanbeigi :
%Github.com/shahin77hb
%%
A=toeplitz([alpha beta zeros(1,n-2)],[alpha gamma zeros(1,n-2)])
c=zeros(n);
c(n,1)=1;
for i=1:n
    for j=1:n
        if j==i+1
            c(i,j)=1;
        end
    end
end
A_new=c*A;
%%
A11=A_new(1:n-1,1:n-1);
w=A_new(n,1:n-1);
p=A_new(1:n-1,n);
b_2=b(1:n-1);
b_1=b(end);
%%
u=BS(alpha,beta,gamma,n,p)';
v=BS(alpha,beta,gamma,n,b_2)';
xn=(w*v-b_1)/(w*u);
x_1=v-xn*u;
x=[x_1;xn];