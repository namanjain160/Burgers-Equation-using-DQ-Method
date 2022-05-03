function A = W(x, k)
M = length(x);
A = zeros(k, M, M);
for r = 1:k
    if r == 1
        for i = 1:M
            for j = 1:M
                A(r, i, j) = 1;
                if i ~= j
                    for l = 1:M
                        if (l ~= i) && (l ~= j) 
                            A(r, i, j) = A(r, i, j)*(x(i) - x(l))/(x(j) - x(l));
                        end
                    end
                    A(r, i, j) = A(r, i, j)/(x(j) - x(i));                        
                else
                    A(r, i, j) = 0;
                    for l = 1:M
                        if l ~= i
                            A(r, i, j) = A(r, i, j) + 1/(x(i) - x(l));
                        end
                    end
                end
            end
        end
    else
        for i=1:M
            for j=1:M
                if i~= j
                    A(r, i, j) = 2*(A(1,i,j)*A(r-1,i,i) - A(r-1,i,j)/(x(i) - x(j)));
                end
            end
            
            for j=1:M
                if i~=j
                    A(r, i, i) = A(r, i, i) - A(r, i, j);
                end
            end
        end
    end
end
end