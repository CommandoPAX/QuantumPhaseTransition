using ITensors
using Plots
using Random
using ITensorMPS
using LaTeXStrings

print("begin\n")
#https://www.overleaf.com/2663516136vbjstqbfdvgk#c9d482

function true_rng(N,max_per_site)
    state = [rand(0:max_per_site) for j in 1:N]
    tot = sum(state)
    while tot > N
        j = rand(1:N)
        if state[j] != 0
            sub = rand(1:state[j])
            if tot-sub < N
                sub += N-tot
            end
            state[j] -= sub
            tot -= sub
        end
    end
    tot = sum(state)
    while tot < N
        r = rand(1:N)
        if state[r] == 0
            state[r] = N-tot
        end
    end
    return state
end

function Googoogaga()
    # Initializes N bosons sites
    N = 10
    sites = siteinds("Qudit", N, dim=N+1;conserve_qns=true)
    U = 5
    J = 0.001

    # Trying to build the Hamiltonian
    os = OpSum()
    for j=1:N-1
        os += -J,"A",j,"Adag",j+1
        os += - J,"Adag",j,"A",j+1
        os += U/2,"n",j,"n",j 
        os += -U/2,"n",j
    end
    H = MPO(os,sites)
    TestState = ["1","1","1","1","1","1","1","1","1","1"]
    psi_test = MPS(sites, TestState)
    psi0 = random_mps(sites, TestState;linkdims=10)


    nsweeps = 50# number of sweeps : 5
    maxdim = [20i+50 for i in 1:nsweeps] # bonds dimension
    cutoff = [1E-10] # truncation error

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    sing_density=zeros(N,N)
    for j in 1:N
        for k in 1:N
            temp=deepcopy(psi)
            temp=apply(op("A",sites,k),temp)
            temp=apply(op("Adag",sites,j),temp)
            sing_density[j,k]= abs(inner(psi,temp))^2
           
        end
    end
    diag = [sing_density[j,j] for j in 1:N]
    plt = Plots.plot(1:N,diag)
    display(plt)
    @show(sum(diag))
    
    return sing_density
end

function plot(N,proba)
    col_grad = cgrad([:blue, :orange], [0.1, 0.5, 1.5])
    Plots.surface(1:N,1:N,proba,xlabel="i",ylabel="j",zlabel="Proba",color=col_grad)
    
end

function HITMAN(N,proba)
    xs = 1:N
    ys = 1:N
    plt = Plots.heatmap(xs,ys,proba)
    display(plt)
end
    
let 
    Googoogaga()
end
