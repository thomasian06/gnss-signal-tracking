function str = latex_matrix(A)
    
    [r, c] = size(A);
    str = "\begin{bmatrix}";
    for i = 1:r
        str = str + " & ";
        str = str + sprintf("%0.2f", A(i,1));
        if c > 1
            for j = 2:c
                str = str + " & ";
                str = str + sprintf("%0.2f", A(i,j));
            end
        end
        str = str + " \\ \hline";
    end
    str = str + " \end{bmatrix}";

end