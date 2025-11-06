module RotinasOtim

    using LinearAlgebra

    # Carrega as rotinas definidas em 
    include("line_search.jl")

    # Funções de teste
    include("funcoes_teste.jl")

    # Steepest com bloqueio
    include("steepest.jl")

    # Exporta as rotinas 
    export  Steepest


end
