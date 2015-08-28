
using PyPlot

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

### histograma y esplano

function histograma(A::Array, Cajas::Int=10, normalizado::Int=1) # 1 -> true, 0 -> false
    num_datos = length(A)
    if normalizado == 1
        fact_norm = num_datos
    else 
        fact_norm = 1
    end
    out = zeros(Int, Cajas)
    limitescajas = linspace(minimum(A),maximum(A),Cajas+1)
    for i in 1:Cajas
        for a in 1:num_datos
            if A[a] >= limitescajas[i] && A[a] < limitescajas[i+1]
                out[i] += 1
            end
        end
    end
    return out/fact_norm
end

function esplano(histo::Array) # 100 -> max, 0 -> min, el arreglo que se le dá debe ser un histograma
    prom = mean(histo)
    desviacion = 0
    for i in histo
        desviacion += abs(prom - i)
    end
    return (1 - desviacion/(length(histo)*sum(histo)))*100 # desviación total porcentual
end

function flip_one(A::Array{Float64,2},i::Int64,j::Int64)
    A[i,j]*=-1
    A
end

function energia_ij(configuracion::Array{Float64,2},n::Int64,m::Int64,i::Int64,j::Int64)
    -configuracion[i,j]*(configuracion[mod1(i-1,n),j]+configuracion[mod1(i+1,n),j]+
        configuracion[i,mod1(j-1,m)]+configuracion[i,mod1(j+1,m)])/2
end

#Algoritmo de wang-landau

function wanglandau(N,M)
    σ=config(N,M)
    Emin, Emax = -2*N*M, 2*N*M
    Energias = [Emin:Emax]      #[i for i in Emin:Emax] #Vector de energías
    H = zeros(length(Energias)) #Histograma
    S = zeros(length(Energias)) #Vector de entropías
    E = energia(σ)
    H[E-Emin+1] += 1 #Aumentamos en 1 el histograma
    f = 1 #factor de modificación inicia
    while f >= 1e-6
        #LLenamos el histograma al principio
        for k in 1:(2*length(H))
            i,j=rand(1:N),rand(1:M)
            ΔE=-2*energia_ij(σ.conf,N,M,i,j)
            η = min(1, exp(S[E-Emin+1])-exp(S[E+ΔE-Emin+1]))
            if rand() < η
                E = E + ΔE
                σ.conf = flip_one(σ.conf,i,j)
            end
            H[E-Emin+1] += 1
            S[E-Emin+1] += f
        end
        #Ahora empezamos el algoritmo
        contador = 1
        while minimum(H[H.!=0])<0.8*mean(H[H.!=0])
            i,j=rand(1:N),rand(1:M)
            ΔE=-2*energia_ij(σ.conf,N,M,i,j)
            η = min(1, exp(S[E-Emin+1])-exp(S[E+ΔE-Emin+1]))
            if rand() < η
                E = E + ΔE
                σ.conf = flip_one(σ.conf,i,j)
            end
            H[E-Emin+1] += 1
            S[E-Emin+1] += f
            contador += 1
        end
        H = zeros(4*N*M+1) #Histograma
        f = f/2
    end
    S
end

#FUnción de particion
function Z(T)
    ZZ=0
    for E in Energias
        ZZ+=g(SS)[E-Emin+1]*exp(-β(T)*E)
    end
    ZZ
end


#Entropia final
N=5
M=5

SS=wanglandau(N,N)
g(SS)=[exp(i) for i in SS] 
Emin, Emax = -2*N*M , 2*N*M
Energias = [Emin:Emax] 






