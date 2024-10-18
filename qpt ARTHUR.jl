using ITensors
using Plots
using Random
using ITensorMPS
using LaTeXStrings

print("begin\n")
#https://www.overleaf.com/2663516136vbjstqbfdvgk#c9d482

function true_rng(size,nb_bosons,max_per_site)
    state = [rand(0:max_per_site) for j in 1:size] #creation of random boson patches
    tot = sum(state)
    while tot > nb_bosons #Suppression of bosons if we have too many
        j = rand(1:size)
        if state[j] != 0
            tot = sum(state)
            sub = rand(1:state[j]) #subtract a random value
            if tot-sub < nb_bosons #check if said value doesn't bring the total under the N total of boson required
                sub = tot-nb_bosons
            end
            state[j] -= sub
        end
    end
    while tot < nb_bosons #Addition of bosons if we have not enough
        j = rand(1:size)
        if state[j] < max_per_site
            add = rand(state[j]:max_per_site)
            tot = sum(state)
            if tot+add > nb_bosons #check if said value doesn't bring the total above the N total of boson required
                add = nb_bosons-tot
            end
            state[j] += add
        end
    end
    result = map(string,state)
    return result
end

function Googoogaga(N,U,J)
    # Initializes N bosons sites
    sites = siteinds("Qudit", N, dim=N+1;conserve_number=true, conserve_qns = true)

    # Trying to build the Hamiltonian
    os = OpSum()
    for j=1:N-1
        os += -J,"A",j,"Adag",j+1
        os += - J,"Adag",j,"A",j+1
    end
    for j=1:N 
        os += U/2,"n",j,"n",j 
        os += -U/2,"n",j
    end
    H = MPO(os,sites)
    #TestState = true_rng(N, 2)
    #@show(TestState)
    TestState = ["2", "0", "1", "1", "1", "1", "1", "1", "1", "1"]
    #sum(map(Int, TestState))
    psi_test = MPS(sites, TestState)
    #@show(psi_test)
    psi0 = random_mps(sites, TestState;linkdims=10)
    #psi = psi_test
    nsweeps = 50 # number of sweeps
    maxdim = [10,20,100,100,200] # bonds dimension
    cutoff = [1E-10] # truncation error

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

    sing_density=zeros(N,N)
    for j in 1:N
        for k in 1:N
            temp=deepcopy(psi)
            temp=apply(op("A",sites,k),temp)
            temp=apply(op("Adag",sites,j),temp)
            sing_density[j,k]= inner(psi,temp)
           
        end
    end
    diag = [sing_density[j, j] for j in 1:N]
    print(sum(diag))
    return sing_density
end

function plot(N,proba)
    col_grad = cgrad([:orange, :blue], [0.1, 0.3, 0.8])
    Plots.surface(1:N,1:N,proba,xlabel="i",ylabel="j",zlabel="Proba",color=col_grad)
    
end

function HITMAN(N,proba)
    xs = 1:N
    ys = 1:N
    plt = Plots.heatmap(xs,ys,proba)
    display(plt)
end
    
let 
    N=10
    U=10
    J=1
    plot(N, Googoogaga(N,U,J))
end

function MarieAntoinette(single_density,middle)
    list = map(log10,single_density[:,middle])
    list_ = list[1:Int(length(list)/2)]
    lslope = []
    for j in 1:(length(list_)-1)
        push!(lslope,list_[j]-list_[j+1])
    end
    slope = sum(lslope)/length(lslope)
    print("\n")
    print(slope)
    print("\n")
    print(lslope)
    print("\n")
    s = 0
    for j in 1:middle
        if slope-0.1 < lslope[j] < slope+0.1
            s+=1
        end
    end
    if s > N/2
        print("linear")
    else 
        print("not linear")
    end
end
