function B=CirShift(A)
l=length(A);
B = [A(end) A(1:end-1)];