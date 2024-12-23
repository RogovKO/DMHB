function [C] = kronsum(A,B)

C = kron(A, eye(length(B))) +  kron(eye(length(A)),B);
end

