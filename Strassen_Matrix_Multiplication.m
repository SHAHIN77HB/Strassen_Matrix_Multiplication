function C = Strassen_Matrix_Multiplication(A,B)
% This function computes multiplication of two N x N matrix A and B with
% help of Strassen Algorithm which needs less matrix multiplications than
% normal matrix multiplication A*B

n=length(A);
A11=A(1:n/2,1:n/2);
A12=A(1:n/2,n/2+1:end);
A21=A(n/2+1:end,1:n/2);
A22=A(n/2+1:end,n/2+1:end);
B11=B(1:n/2,1:n/2);
B12=B(1:n/2,n/2+1:end);
B21=B(n/2+1:end,1:n/2);
B22=B(n/2+1:end,n/2+1:end);

M1=(A11+A22)*(B11+B22);
M2=(A21+A22)*B11;
M3=A11*(B12-B22);
M4=A22*(B21-B11);
M5=(A11+A12)*B22;
M6=(A21-A11)*(B11+B12);
M7=(A12-A22)*(B21+B22);

C11=M1+M4-M5+M7;
C12=M3+M5;
C21=M2+M4;
C22=M1-M2+M3+M6;

C=[C11 C12;C21 C22];