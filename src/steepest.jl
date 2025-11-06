

#
# Steepest descent
#
function Steepest(x0::Vector,xmin::Vector,xmax::Vector,f::Function,df::Function,tol=1E-6,nmaxiter=1000)


    # Verifica se o x0 está dentro dos limites
    if !all(xmin .<= x0 .<= xmax) 
        error("Steepest:: valores iniciais das variáveis de projeto não respeitam as restrições laterais")
    end

    # Copia o ponto inicial para outro vetor 
    x = copy(x0)

    # Loop até que a condição de parada seja obtida
    for iter=1:nmaxiter

        # Gradiente no ponto
        ∇f = df(x)

        # Norma da derivada no ponto
        norma = norm(∇f) 

        # Direção de busca
        d = -∇f/norma

        # Etapa de bloqueio da direção de descida
        for i in LinearIndices(d)

            # Se d[i]<0 e já estamos em xmin 
            if (d[i]<0.0 && x[i] ≈ xmin[i]) ||  (d[i]>0.0 && x[i] ≈ xmax[i] )
               d[i] = 0.0
            end

        end

        # Critério de parada agora é direto no d
        if norm(d) <= tol 

            # Anuncia a convergência
            println("Tolerância atingida $norma e $iter iterações")

            # Testes de consistência 
            flag_consistencia_livre = true
            flag_consistencia_min   = true
            flag_consistencia_max   = true
            

            # Loop pelas componentes da solução 
            for i in LinearIndices(d)

                # Teste de consistência para variável livre
                if (xmin[i] <  x[i] < xmax[i]) && abs(d[i])> tol
                    @show x[i], d[i]
                   flag_consistencia_livre = false
                end
                
                # Teste de consistência para variável bloqueada para baixo
                if (x[i]==xmin[i] && -∇f[i]>0.0)
                   flag_consistencia_min = false
                end
                
                # Teste de consistência para variável bloqueada para cima
                if (x[i]==xmax[i] && -∇f[i]<0.0)
                   flag_consistencia_max = false
                end

            end
        
            println("Testes de consistência indicam") 
            println("Variáveis livres               :: $(flag_consistencia_livre)")
            println("Variáveis bloqueadas para baixo:: $(flag_consistencia_min)")
            println("Variáveis bloqueadas para cima :: $(flag_consistencia_max)")

            break
        end

        #
        # Calcula o α_lim para o line-search
        #
        α_lim = maxintfloat(1.0)

        # Loop por cada componente 
        for i in LinearIndices(d)
            
            # Se d[i]>0 temos que verificar o quanto podemos andar
            # até chegar no limite máximo
            # x + αd == xmax  -> α = (xmax - x)/d
            if d[i]>0.0 
               α_lim = min(α_lim, (xmax[i]-x[i])/d[i])

            # Se d[i]<0 temos que verificar o quanto podemos andar
            # até chegar no limite máximo
            # x + αd == xmin  -> α = (xmin - x)/d
            elseif d[i]<0.0
               α_lim = min(α_lim, (xmin[i]-x[i])/d[i])

            end   

        end

        # Line-search para encontrar o passo ótimo nessa direção
        α =  LineSearch_NR_constr(x, d, f, α_lim)

        # Atualiza o ponto 
        x .= x .+ α*d

    end

    #
    # Retorna as variáveis ótimas do problema
    #
    return x

end

