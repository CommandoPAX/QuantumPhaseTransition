using ITensors
using Plots
using Random
using ITensorMPS
using LaTeXStrings
using ForwardDiff
using EasyFit

#https://www.overleaf.com/2663516136vbjstqbfdvgk#c9d482

function true_rng(size,nb_bosons,max_per_site)
    max_per_site = floor(Int, max_per_site)
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

function SPDM(sites, psi, N)
    sing_density=zeros(N,N)
    for j in 1:N
        for k in 1:N
            temp=deepcopy(psi) # deepcopy of psi to avoid accidentaly messing with the files
            temp=apply(op("A",sites,k),temp) # applies the annihilation operator on site k
            temp=apply(op("Adag",sites,j),temp) # applies the creation operator on site j 
            sing_density[j,k]= inner(psi,temp) # computes the inner product of psi and the single particle density operator
        end
    end
    Partic_Numb = sum([sing_density[j, j] for j in 1:N])
    return sing_density, Partic_Numb
end

function isLinear(x,y,thresh)
    fit = fitlinear(x,y)
    if(fit.R < thresh)
        return false
    else 
        return true
    end
end

function Plot_3D(N,DATA)
    col_grad = cgrad([:orange, :blue], [0.1, 0.3, 0.8])
    Plots.surface(1:N,1:N,DATA,xlabel="i",ylabel="j",zlabel="Proba",color=col_grad)
end

#Plot the density for the site j in log-log scale
function Plot_one_site_density(single_density,j)
    one_site_density = log.(single_density[:,j])
    N=length(one_site_density)
    one_site_density=one_site_density[Int(N/2)+1:N]
    plt=Plots.plot(log.(Int(N/2)+1:N),one_site_density,xlabel="site #",ylabel="one-site density",title="One site density for site " * string(j),
    legend=false, linewidth=2,linecolor=[:black])
    display(plt)
end 

function Hitmap(N,DATA)
    xs = 1:N
    ys = 1:N
    plt = Plots.heatmap(xs,ys,DATA)
    display(plt)
end

function Delta(Enp1, En, N)
    return Enp1 - (N+1)*(En/N) 
end

function Plot_Energy(delta, U, J)
    Plots.plot(U/J, delta, xlabel="U/J", ylabel="Delta", title="Delta as a function of U/J")
end

function size_increasing_behavior(filenames,loglogscale=true)

    plt = Plots.plot()
   
    
        for f in filenames
            
            @show(f)
            single_densities = import_density(f)
            N= length(single_densities[:,1])
            log_sp = [log(abs(single_densities[i])) for i in Int(N/2)+1:N ]
            if(loglogscale)
                log_x = [log(abs(i-N/2)) for i in Int(N/2)+1:N]
            else
                log_x = [abs(i-N/2) for i in Int(N/2)+1:N]
            
            end
            
            plot!(log_x,log_sp,label="N= "*string(N))
        end
        plot!(xlabel="log(distance)",ylabel="log(density)")
    display(plt)
    
end 

#Export data in a .txt file. Data could be reimported using 'import_density' function.
function export_density(single_matrix_density,U,J)
    N= length(single_matrix_density[:,1])

    file = open("simulation_data_size_"*string(N)*"_"*string(Dates.format(now(), "yyyy-mm-dd-HH_MM_SS"))*".txt","w") do f

        write(f,"U="*string(U)*" J="*string(J)*" N="*string(N))

        for i in 1:N
            write(f,"\n")

            for j in 1:N
                write(f, string(single_matrix_density[i,j])*" ")
            end

        end
        print("Data saved in file")
        close(f)
    end
end

#Import data from a .txt file, cf 'export_density' function.
function import_density(filepath)
    single_density_matrix = []
    open(filepath,"r") do f
        data = readlines(f)
        U,J,N = split(data[1]," ")
        U= parse(Float64,U[3:end])
        J= parse(Float64,J[3:end])
        N= parse(Int,N[3:end])
       
        single_density_matrix = zeros(N,N)
        
        for i in 2:N+1
            line = split(data[i]," ")
            for j in 1:N
            single_density_matrix[i-1,j]=parse.(Float64,line[j])
            end
            
        end
  
    print("Data retrieved : "* "U="*string(U)*" J="*string(J)*" N="*string(N))
    close(f)
    
    end
    return single_density_matrix
    
end

function Run_Simulation(N, U, J)
    #=
    This function handles creating the sites with a random boson distribution and running the simulation using the DMRG algorithm.
    N -- Number of bosons and bosons sites
    U -- Value of the pairwise interaction
    J -- Value of the tunneling rate
    =#

    # Initializes N bosons sites
    print("Initializing...\n")
    sites = siteinds("Qudit", N, dim=N+2;conserve_number=true, conserve_qns = true)
    
    # Variables needed for the dmrg algorithm
    nsweeps = 40 # number of sweeps
    maxdim = [10,20,100,100,200] # bonds dimension
    cutoff = [1E-10] # truncation error

    # Creates the Bose-Hubbard model Hamiltonian
    print("Computing Hamiltonian...\n")
    os = OpSum()
    for j=1:N-1
        os += -J,"A",j,"Adag",j+1
        os += - J,"Adag",j,"A",j+1
    end
    for j=1:N 
        os += U/2,"n",j,"n",j 
        os += -U/2,"n",j
    end
    H = MPO(os,sites) # Transforms the Hamiltonian into a matrix product operator

    # Intialises the random state from a random distribution of states
    print("Computing initial state...\n")
    Init_State = true_rng(N, N, N/4)
    Init_State2 = true_rng(N, N+1, (N+1)/4)
    psi0 = ITensors.ITensorMPS.MPS(sites, Init_State)
    psi02 = ITensors.ITensorMPS.MPS(sites, Init_State2)
    #psi02 = MPS(sites2, Init_State2; linkdims=10)

    # Executes the DMRG algorithm
    print("Applying DMRG...\n")
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff, outputlevel=1)
    energy2,psi2 = dmrg(H,psi02;nsweeps,maxdim,cutoff, outputlevel=1)

    return energy, energy2 # This is an int

    # Gets the single particle density matrix from the DMRG results
    #print("Getting single particle densities...\n")
    #Single_Particle_Density, Particle_Number = SPDM(sites, psi, N)

    # Creates the plots and display them (with the particle number at the end to verify)
    #print("Final particle number : ", Particle_Number)
    #Plot_3D(N, Single_Particle_Density)
    #Hitmap(N, Single_Particle_Density)
    #Plot_one_site_density(Single_Particle_Density, Int(N/2))

    # Check for the phase we're in 
    #return isLinear(log.(Single_Particle_Density[Int(N/2):N,Int(N/2)]) , [1.0*i for i in Int(N/2):N],0.995)

end
    
let 
    #U = [round(2.5+0.1*j, digits = 3) for j in 0:15]
        U = [2, 3, 4, 5, 6]
    N = 50
    Index = []
    Values = []
    Old = 0
    for k in U 
        print("########################## Iteration : ", k, " ##########################\n")
        En, Enp1 = Run_Simulation(N, k, 1)
        push!(Index, k)
        New=Delta(Enp1, En, N)
        push!(Values, New)
        if abs(New-Old) > 0.02 && Old !=0
            print("########################## Phase transition found at U/J : ", k, "##########################\n")
        end 
        Old = New
    end
    Plots.plot(Index, Values, xlabel="U/J", ylabel="Delta", title="Energy variation as a function of the ratio U/J", legend=false, linewidth=2,linecolor=[:black])
end
# Valeur bizarre, refaire la simu entre 3 et 4.5 pour vérifier position transition de phase