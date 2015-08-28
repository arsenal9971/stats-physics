#Script de modelo de Ising 

module Metropolis
#T=float(ARGS[1]) #Temperatura del sistema
#N=int(ARGS[2]) #Tamaño en x
#M=int(ARGS[2]) #Tamaño en y

export β,energia,magnet,metropolis,config,energiaprom,magnetprom
#Temperatura y beta
kB=1
β(T)=1/(kB*T);


#Creamos el tipo configuracion
type Config
    N::Int
    M::Int
    conf::Array{Float64,2}
end

#Creamos un constructor de ese tipo 
function config(a,b)
    Config(a,b,2*(int(rand(a,b)))-1)
end

#Función módulo arreglada
function mod(i,L)
    if i>L
        k=1
    elseif i==0
        k=L
    else
        k=i
    end
    k
end

#Calculo de Energiawse
function energia(σ::Config)
    L1=σ.N
    L2=σ.M
    #Calculamos la energía total
    E=0
        for i in 1:L1
            for j in 1:L2
                E+=-σ.conf[i,j]*(σ.conf[mod(i-1,L1),j]
                    +σ.conf[mod(i+1,L1),j]+σ.conf[i,mod(j-1,L2)]
                    +σ.conf[i,mod(j+1,L2)])
            end
        end
        E
end

function magnet(σ::Config)
    sum=0
    for i in 1:σ.N
        for j in 1:σ.M
            sum+=σ.conf[i,j]
        end
    end
    sum
end

#Proceso estocastico
function metropolis(σ::Config,T)
    L1=σ.N
    L2=σ.M
    i=int((L1-1)*rand())+1
    j=int((L2-1)*rand())+1 #Elegimos entrada aleatoria
    ΔE=float(2*σ.conf[i,j]*(σ.conf[mod(i-1,L1),j]+σ.conf[i,mod(j-1,L2)]
        +σ.conf[mod(i+1,L1),j]+σ.conf[i,mod(j+1,L2)]))
    α=exp(-β(T)*ΔE)
    r=rand()
    if r<α
        σ.conf[i,j]=-σ.conf[i,j]
    end
    σ
end

function energiaprom(N,M,T,iter)
    σ=config(N,M)
    Energias=Float64[]
    for i in 1:int(iter)
        push!(Energias, energia(σ))
        σ=metropolis(σ,T)
    end
    sum(Energias[int(iter)/2:(iter)])/length(Energias[int(iter)/2:(iter)])
end

function magnetprom(N,M,T,iter)
    σ=config(N,M)
    Magnets=Float64[]
    for i in 1:int(iter)
        push!(Magnets, magnet(σ))
        σ=metropolis(σ,T)
    end
    sum(Magnets[int(iter)/2:(iter)])/length(Magnets[int(iter)/2:(iter)])
end

end




