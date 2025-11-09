#
# Função objetivo
#
function Esfera1(x) 
     (x[1]-2)^2 + (x[2]-3)^2
end

#
# Derivada da função objetivo
#
function dEsfera1(x) 
    [ 2*(x[1]-2) ;
      2*(x[2]-3) ]
end


#
# Matriz Hessiana da função objetivo
#
function HEsfera1(x)
    [2.0 0.0 ; 
     0.0 2.0]
end

#
# Função objetivo
#
function Rosenbrook(x) 
     100*(x[2]-x[1]^2)^2 + (1-x[1])^2
end

#
# Derivada da função objetivo
#
function dRosenbrook(x) 
    [ 200*(x[2]-x[1]^2)*(-2*x[1]) + 2*(1-x[1])*(-1) ;
      200*(x[2]-x[1]^2) ]
end

#
# Matriz Hessiana da função objetivo
#
function HRosenbrook(x)
    H11 =   -(400*(x[2]-x[1]^2))+800*x[1]^2+2
    H12 = -400*x[1]
    H21 = H12 
    H22 = 200.0

    [H11 H12 ; 
     H21 H22]
end

#
# Função apostila
#
function apostila(x) 
     exp(x) -4*x +2
end

function dapostila(x) 
     exp(x) -4
end