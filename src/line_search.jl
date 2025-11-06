#
# Função auxiliar para calcular a primeira e a segundas derivadas em relação a α
# utilizando diferenças finitas com uma perturbação δ
#
function d1d2α(x::Vector,d::Vector,f::Function,α::Float64, δ, δ_min=1E-12)

     # Vamos TENTAR calcular uma boa perturbação 
     δ =  max(δ_min, δ*max(abs(α),1.0))

     # Função no ponto atual 
     f0 = f(x .+ α*d)

     # Perturba para frente
     ff = f(x .+ (α+δ)*d)

     # Perturba para tras
     ft = f(x .+ (α-δ)*d)

     # Estimativa da primeira derivada
     d1 = (ff-ft)/(2*δ)

     # Estimativa da segunda derivada
     d2 = (ff -2*f0 + ft)/(δ^2)

     # Retorna as derivadas
     return d1, d2

end


#
#
# Line-Search por Newton-Raphson 
#
# α = arg min f(x + α*d)
#
# com as derivadas calculadas por DF de passo δ
#
function LineSearch_NR(x::Vector,d::Vector,f::Function;δ=1E-6,tol=1E-10,nmaxiter=100,verbose=true)

    # Estimativa inicial do α
    α = 0.1

    #
    # Loop para encontrar a raiz 
    #
    for iter=1:nmaxiter

        # Calcula a primeira e a segunda derivada utilizando DF
        d1, d2 = d1d2α(x,d,f,α, δ)

        # Se tivermos uma segunda derivada nula, não podemos calcular o próximo ponto e devemos 
        # retornar com um erro (α negativo)
        if abs(d2)<tol
            println("LineSearch_NR:: segunda derivada é nula no ponto α=$α")
            return -1.0
        end

        # Nova estimativa do α
        αn = α - d1/d2

        # Diferença na estimativa
        difalpha = abs(αn-α)

        # Testa se atingimos a tolerância 
        if  verbose && difalpha < tol 
            println("LineSearch_NR:: atingimos a tolerância $(difalpha) em $(iter) iterações")
            α = αn
            break
        end

        # Atualiza a estimativa do alpha 
        α = αn  

    end

    # Retorna a estimativa do α 
    return α

end



#############################################################################################################
#
# Soluciona o problema de Line-Search com restrição lateral no passo α
#
#  α = arg min f(α)
#      S.t
#          0 < α < α_lim
#
# Para isso, montamos a função Lagrangiana 
#
# L(α,λ,s) = f(α) + 0.5 λ (α - α_lim + s^2)
#
# que tem condições de mínimo em 
#
#
#  F1 = dL/dα = df/dα + λ = 0 
#  
#  F2 = dL/dλ = 0.5*(α - α_lim + s^2) = 0  --> α - α_lim + s^2 = 0
#
#  F3 = dL/ds = λs - μ = 0
# 
# onde μ > 0 é só um artifício numérico para evitar problemas nas iterações de NR (é como se 
# ao invés de igualarmos a zero, igualamos a um número que tende a zero).
#
# Assim, temos um problema definido por 3 equações não lineares em relação a α, μ e s. 
#
#  Linearizando cada equação em relação a α, s e μ (nessa ordem)
#
#  F1(α,s,μ)  = F1(α_k,s_k,μ_k) + dF1/dα Δα + dF1/ds Δs + dF1/dλ Δλ ≈ 0
#  F2(α,s,μ)  = F2(α_k,s_k,μ_k) + dF2/dα Δα + dF2/ds Δs + dF2/dλ Δλ ≈ 0
#  F3(α,s,μ)  = F3(α_k,s_k,μ_k) + dF3/dα Δα + dF3/ds Δs + dF3/dλ Δλ ≈ 0
#
# tal que obtemos 
#
#  d2f/dα2 Δα +   0  Δs   + 1   Δλ ≈ -F1(α_k,s_k,μ_k)
#       1  Δα + 2s_k Δs   + 0   Δλ ≈ -F2(α_k,s_k,μ_k)
#       0  Δα + μ_k  Δs   + s_k Δλ ≈ -F2(α_k,s_k,μ_k)
#
# ou seja, a cada iteração temos um sistema 3 × 3 definido pela 
# matriz Jacobiana (J), por um vetor de incrementos (Δ) e por um 
# vetor do lado direito (F)
#
#
#   [ d2f/dα2  0    1     [Δα       [F1
#        1    2s_k  0      Δs   = -  F2 (α_k,s_k,μ_k)
#        0    μ_k  s_k]    Δλ]       F3]
#
#
#  Esse sistema pode ser solucionado tranquilamente, mas pode ocorrer 
#  alguma instabilidade numérica se s+Δs e/ou λ+Δλ não forem >= 0. Isso 
#  pode ser "melhorado" fazendo um line-search interno para achar 
#  δ ∈ (0,1] (1 seria o passo da própria iteração do NR) tal que 
#
#  α = α + δΔα
#  s = s + δΔs
#  λ = λ + δΔλ       
#
# Por fim, temos que observar que utilizei (para evitar problemas mesmo) uma 
# variável auxiliar μ na definição do F3. Esse cara pode ser ajustado a cada 
# iteração do NR por algo do tipo μ = τμ com τ ∈ (0,1) para ir reduzindo esse valor.
# Como não queremos que esse cara vá a zero (senão nem precisava dele), podemos
# colocar uma condição 
#
# μ = max(μ_min, τμ) 
#
# com μ_min ≈ 1E-8 
#
# Se tudo der certo, ou temos a restrição ativa (α vai bater na parede) e μ>0 com s=0
# ou o problema é irrestrito e μ=0 e s>0.
#   
#
# Um problema que pode ocorrer com a matriz é que s→0 vai fazer dois termos da diagonal 
# serem nulos. Calculando a inversa da matriz e fazendo o limite pela direita obtemos 
#
#  [ 0  1   0 
#    0  0  1/λ
#    1 -d2  0]
#
# tal que a inversa existe 
#
function LineSearch_NR_constr(x::Vector, d::Vector, f::Function, α_lim=0.0; tol=1E-9, maxiter=100, δ_DF=1E-6,
                              μ0=1E-3, μ_tol=1E-12, τ=0.1, α_min=1E-8)::Float64
    
        
    # Se α_lim for nulo, então o problema é irrestrito e podemos chamar o NR tradicional 
    if α_lim == 0.0
       return LineSearch_NR(x,d,f,δ_DF,tol)::Float64
    end
    
    ########################################################################################
    #                       Variáveis do problema restrito
    #
    #                        Com as estimativas iniciais
    ########################################################################################

    # Estimativa inicial do α 
    α = max(α_min, α_lim / 2)

    # Variável de folga. Aqui podemos calcular diretamente quanto de folga temos
    s = α_lim - α

    # Multiplicador de KT. Aqui podemos fazer uma bruxaria usando a definição de F3
    # e isolando o λ. Como realmente não sabemos o valor desta penalização, vamos 
    # deixar umas salvaguardas para não ficar nem muito pequeno e nem muito grande 
    λ = max(1e-6, μ0 / max(1e-6, s))

    # Perturbação (que deve tender a zero na condição de complementariedade)
    μ = μ0

    ########################################################################################
    #                          Iterações do NR com restrições 
    ########################################################################################
    for iter in 1:maxiter

        # Calcula as derivadas em relação a α, usando DF
        d1, d2 = d1d2α(x,d,f,α,δ_DF)

        # Resíduos
        F1 = d1 + λ
        F2 = α - α_lim + s^2
        F3 = s*λ - μ

        # Verifica a convergência (primal + estacionário)
        if abs(F1) < tol && abs(F2) < tol && abs(F3) < max(tol, μ_tol)
            return α
        end

        # Matriz Jacobiana
        J = [ d2   0     1;
              1    2*s   0;
              0    λ     s ]

        # R.H.S do sistema 
        F = [F1; F2; F3]

        # Calcula as variações nesta iteração
        Δ = -J \ F   

        # Para facilitar
        Δα, Δs, Δλ = Δ[1], Δ[2], Δ[3]

        # Aqui teríamos a busca em linha interna,  para garantir que 
        # s e λ sejam positivos. O objetivo aqui é minimizar o resíduo 
        # na direção de Δ
        # Por hora, sem isso
        δ = 1.0
 
        # Agora podemos atualizar os valores
        α += δ*Δα
        s += δ*Δs
        λ += δ*Δλ

        # Só garantindo α ∈ (0,α_lim] e que s e λ sejam positivos
        α = clamp(α, 1e-12, α_lim)
        s = max(s, 1e-12)
        λ = max(λ, 1e-12)

        # Reduz μ 
        μ = max(μ_tol, τ*μ)

    end

    #
    # Se chegamos aqui, então não atingimos a tolerância...snif
    #
    println("LineSearch_NR_constr:: não foi possível atingir a tolerância.")
    return -1.0
end
