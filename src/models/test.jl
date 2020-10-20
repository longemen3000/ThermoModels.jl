using CodeTracking, Revise, Latexify


f(x,y)=x^2+y^2
latexify(@code_expr(f(1,2)))
