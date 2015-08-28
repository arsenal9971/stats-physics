#Algoritmo de Wang y Landau
module wanglandau
#T=float(ARGS[1]) #Temperatura del sistema
#N=int(ARGS[2]) #Tamaño en x
#M=int(ARGS[2]) #Tamaño en y


export β,config,mod,energia, magnet,flip_one,energia_ij

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

#Calculo de Energias y magnetización 
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
        E/2
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

#Funciones de flip y ΔE
function flip_one(A::Array{Float64,2},i::Int64,j::Int64)
    A[i,j]*=-1
    A
end

function energia_ij(configuracion::Array{Float64,2},n::Int64,m::Int64,i::Int64,j::Int64)
    2*configuracion[i,j]*(configuracion[mod1(i-1,n),j]+configuracion[mod1(i+1,n),j]+
        configuracion[i,mod1(j-1,m)]+configuracion[i,mod1(j+1,m)])/2
end

end




