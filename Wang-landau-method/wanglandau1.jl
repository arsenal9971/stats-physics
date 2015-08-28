
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

#Algoritmo de wangalandau

type Wanglandau
    S::Array{Float64,2}
    σf::Config
end



Nitt=10000000 #Número de iteraciones del algoritmo
planitud=0.9 

function wanglandau(Nitt,N,M, planitud)
    σ=config(N,M)
    L1, L2 = σ.N, σ.M
    Emin, Emax = -2*L1*L2 , 2*L1*L2
    Energias = [Emin:Emax]      #[i for i in Emin:Emax] #Vector de energías
    H = zeros(4*L1*L2+1) #Histograma
    S = zeros(4*L1*L2-1) #Vector de entropías
    E = energia(σ)
    H[E-Emin+1] += 1 #Aumentamos en 1 el histograma
    f = 1 #factor de modificación inicial
    for k in 1:Nitt
        i,j=rand(1:L1),rand(1:L2)
        ΔE=-2*energia_ij(σ.conf,L1,L2,i,j)
        η = exp(S[E-Emin+1])-exp(S[E+ΔE-Emin+1])
        if rand() < η
            E = E + ΔE
            σ.conf = flip_one(σ.conf,i,j)
        end
        H[E-Emin+1] += 1
        S[E-Emin+1] += f
        if k%100==0
            aH=mean(H)
            mH=minimum(H)
            if mH>aH*planitud
                H = zeros(4*L1*L2+1) #Se reinicia el histograma
                f=f/2 #Reducimos el factor de modificació
            end
        end
    end
    S,H
end

N,M=5,5
wanglandau(Nitt,N,M, planitud)










