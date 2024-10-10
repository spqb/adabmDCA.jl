# module eaDCAsrc

# using Base.Threads
# using StatsBase
# using Distances
# using Logging
# using Random

# include("utils.jl")
# # using .utils 


# COMPUTE SWITCHING TIME ########################################################################################################################################################################

    function compute_sampling_switch_time_ea(J, vbias, v_model, nsweeps, method)
        println("computing sampling switch time...")
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        Nq, Nv, Ns = size(v_model) 
        useless1, useless2 = copy(v_model), copy(v_model)
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        @threads :dynamic for t in 1:thread_count
                start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                useless1[:, :, start_idx:end_idx] = sampling_function(J, vbias, useless1[:, :, start_idx:end_idx], nsweeps)
        end
        switch_time = @elapsed begin
            @threads :dynamic for t in 1:thread_count
                start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                useless2[:, :, start_idx:end_idx] = sampling_function(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps)
            end
        end
        return switch_time
    end


# SAMPLING FUNCTIONS ########################################################################################################################################################################

    function sampling_ea(J, vbias, contact_list, site_degree, v_model::BitArray{3}, nsweeps, switch_time, switch_flag, method)
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        sampling_function_couplingwise = method == "metropolis" ? metropolis_sampling_couplingwise : gibbs_sampling_couplingwise
        Nq, Nv, Ns = size(v_model) 
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        if switch_flag == false 
            sampling_time = @elapsed begin
                @threads :dynamic for t in 1:thread_count
                    start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                    v_model[:, :, start_idx:end_idx] = sampling_function_couplingwise(J, vbias, v_model[:, :, start_idx:end_idx], contact_list, site_degree, nsweeps)
                end
            end
            switch_flag = (switch_flag == false && sampling_time >= switch_time) ? true : false
        elseif switch_flag == true
            @threads :dynamic for t in 1:thread_count
                    start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                    v_model[:, :, start_idx:end_idx] = sampling_function(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps)
            end
        end
        return v_model, switch_flag
    end


# COUPLINGS ACTIVATION  ########################################################################################################################################################################

    function compute_dkl_matrix(pij_model, fij_target, Nq, Nv)
        ijs = [(0, 0) for _ in 1:Nv*(Nv-1)/2]; index = 0
        for i in 1:Nv, j in i+1:Nv
            index += 1; ijs[index] = (i, j)
        end
        dkl_matrix = zeros(Float32, 3, Int64(length(ijs) * Nq*Nq))
        # multi-threads computation over the number of edges
        thread_count = nthreads(); chunk_size = div(length(ijs), thread_count)
        @threads :dynamic for t in 1:thread_count        
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, length(ijs))
            for (ij, idx) in zip(ijs[start_idx:end_idx], start_idx:end_idx)
                i, j, f_tar, f_mod, f_mod1, f_mod2 = ij[1], ij[2], zeros(Float32, 2), zeros(Float32, 2), zeros(Float32, 2), zeros(Float32, 2)
                count = 1
                for iq in id(i, 1:Nq, Nq), jq in id(j, 1:Nq, Nq)
                    f_tar[1], f_tar[2] = fij_target[iq, jq], 1 - fij_target[iq, jq]
                    f_mod[1], f_mod[2] = pij_model[iq, jq], 1 - pij_model[iq, jq]
                    dkl_matrix[:,id(idx, count, Nq*Nq)] = [iq, jq, Distances.kl_divergence(f_tar, f_mod)]
                    count += 1
                end
            end
        end
        return dkl_matrix
    end

    function select_n_couplings(pij_model, fij_target, n_couplings, Nq, Nv) 
        dkl_matrix = compute_dkl_matrix(pij_model, fij_target, Nq, Nv)
        dkl_matrix = dkl_matrix[:, sortperm(dkl_matrix[3, :], rev=true)] # sort the dkl_matrix
        chosen_couplings = dkl_matrix[1:2, 1:n_couplings] # chose the first n_couplings
        return chosen_couplings
    end

    function update_graph(chosen_couplings::Matrix{Float32}, filter, contact_list, site_degree)
        new, updated = 0, 0
        for idx in 1:size(chosen_couplings, 2)
            iq, jq  = Int64(chosen_couplings[1, idx]), Int64(chosen_couplings[2, idx])
            if filter[iq, jq] == 0 && filter[jq, iq] == 0 
                site_degree[iq] += 1; site_degree[jq] += 1
                contact_list[iq, site_degree[iq]] = jq; contact_list[jq, site_degree[jq]] = iq
                filter[iq, jq], filter[jq, iq] = true, true
                new += 1
            else
                updated += 1
            end
        end       
        println("new coupligs: ", new ,",  updated couplings: ", updated)
        return filter, contact_list, site_degree
    end



# GRADIENT UPDATE ########################################################################################################################################################################

    function gradient_update_ea(J::Matrix{Float32}, filter, v_model::BitArray{3}, fij_natural, pij_model, lr)
        J .+= lr * ((fij_natural .- pij_model) .* filter)
        return J
    end
    

# DO ONE EPOCH ########################################################################################################################################################################

    function do_convergence_ea(J, vbias, filter, contact_list, site_degree, v_model, nsweeps, fij_natural, cij_natural, target_cij, lr, pseudo_count, max_conervgence_step, method, switch_time=0, switch_flag=true)
        slope_err = 100
        cij_model = oneHotCijFast(v_model, ones(size(v_model, 3)), 0) 
        cij_model_act, cij_natural_act = filter .* cij_model, filter .* cij_natural
        p, slope = cor(vec(cij_model_act), vec(cij_natural_act)), compute_slope(cij_model_act, cij_natural_act)
        (p <= target_cij) ? println("initial cij: ", p, " --> starting convergence...") : nothing; flush(stdout) # (th_eps < exp_eps)
        (abs(slope - 1) > slope_err) ? println("initial slope: ", slope, " --> starting convergence...") : nothing; flush(stdout)
        Nq, Nv, Ns = size(v_model)
        pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), 0)
        count = 0
        
        while (count <= max_conervgence_step) && (!(p <= target_cij) && (abs(slope - 1) > slope_err) || (p <= target_cij))
            count += 1
            v_model, _ = sampling_ea(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method)
            pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), 0)
            cij_model = oneHotCijFast(v_model, ones(size(v_model, 3)), 0) 
            cij_model_act, cij_natural_act = filter .* cij_model, filter .* cij_natural
            p, slope = cor(vec(cij_model_act), vec(cij_natural_act)), compute_slope(cij_model_act, cij_natural_act)
            println("step: ", count, "; slope: ", slope, ", pearson Cij: ", p); flush(stdout)
            J = gradient_update_ea(J, filter, v_model, fij_natural, pij_model, lr) # update parameters
        end
        return J, v_model, pij_model
    end

    function do_epoch(J::Matrix{Float32}, vbias::Matrix{Float32}, n_couplings::Float64, filter::BitMatrix, contact_list::Matrix{Int64}, site_degree::Vector{Int64}, v_model::BitArray{3}, fij_natural::Matrix{Float32}, lr::Float64, nsweeps::Int64, n_gradsteps::Int64, pseudo_count, switch_time::Float64, switch_flag::Bool, method)
        Nq, Nv, Ns = size(v_model)
        v_model, switch_flag = sampling_ea(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling
        pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), pseudo_count) # compute 2-point frequencies


        Ntot = div(Nv * (Nv-1), 2) * Nq*Nq 
        n_couplings = round(Int64, (n_couplings) * (Ntot - sum(filter)))
        # println("n_couplings: ", n_couplings)

        chosen_couplings = select_n_couplings(pij_model, fij_natural, n_couplings, Nq, Nv) # chose couplings
        filter, contact_list, site_degree = update_graph(chosen_couplings, filter, contact_list, site_degree) # update filter & contact list & site degree
        J = gradient_update_ea(J, filter, v_model, fij_natural, pij_model, lr) # update parameters
        for step in 2:n_gradsteps
            v_model, switch_flag = sampling_ea(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
            pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), 0) # compute 2-point frequencies
            J = gradient_update_ea(J, filter, v_model, fij_natural, pij_model, lr) # update parameters
        end
        GC.gc()
        return J, v_model, filter, contact_list, site_degree, switch_flag
    end
    

    function do_epoch_with_convergence(J::Matrix{Float32}, vbias::Matrix{Float32}, n_couplings::Float64, filter::BitMatrix, contact_list::Matrix{Int64}, site_degree::Vector{Int64}, v_model::BitArray{3}, fij_natural::Matrix{Float32}, cij_natural::Matrix{Float32}, target_cij, lr::Float64, nsweeps::Int64, n_gradsteps::Int64, max_conervgence_step::Int64, pseudo_count, switch_time::Float64, switch_flag::Bool, method)
        Nq, Nv, Ns = size(v_model)
        v_model, switch_flag = sampling_ea(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling
        pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), pseudo_count) # compute 2-point frequencies

        Ntot = div(Nv * (Nv-1), 2) * Nq*Nq 
        n_couplings = round(Int64, (n_couplings) * (Ntot - sum(filter)))

        chosen_couplings = select_n_couplings(pij_model, fij_natural, n_couplings, Nq, Nv) # chose couplings
        filter, contact_list, site_degree = update_graph(chosen_couplings, filter, contact_list, site_degree) # update filter & contact list & site degree
        J = gradient_update_ea(J, filter, v_model, fij_natural, pij_model, lr) # update parameters
        J, v_model, pij_model = do_convergence_ea(J, vbias, filter, contact_list, site_degree, v_model, nsweeps, fij_natural, cij_natural, target_cij, lr, pseudo_count, max_conervgence_step, method, switch_time, switch_flag)
       
        return J, v_model, filter, contact_list, site_degree, switch_flag
    end
 

# SAVE AND RESTORE MODEL ########################################################################################################################################################################

    # function save_model(J, vbias, v_model, alphabet, save_list, pearsonCij, density, outputpath)
    #     Nq, Nv, Ns = size(v_model)
    #     v_cat = oneHot2Categorical(v_model, Nq)
    #     # save chains
    #     code_pearson = string(round(save_list[save_list .<= pearsonCij][end], digits=3))
    #     code_density = string(round(density, digits=3))
    #     file_chains = open(outputpath * "/trainingchains_density" * code_density * "_pearson" * code_pearson * ".fasta", "w")
    #     for m in 1:Ns-1
    #         head = ">chain $m\n"
    #         line = "$(alphabet[v_cat[:, m]])\n"
    #         write(file_chains, head); write(file_chains, line)
    #     end
    #     head = ">chain $Ns\n"; line = "$(alphabet[v_cat[:, Ns]])"
    #     write(file_chains, head); write(file_chains, line)
    #     close(file_chains)
    #     # save model
    #     file_model = open(outputpath * "/model_density" * code_density * "_pearson" * code_pearson * ".dat", "w")
    #     for i in 1:Nv, j in i+1:Nv
    #         for iq in 1:Nq, jq in 1:Nq
    #         line = "J $(i-1) $(j-1) $(iq-1) $(jq-1) $(J[id(i, iq, Nq), id(j, jq, Nq)])\n"
    #             write(file_model, line)
    #         end
    #     end
    #     for i in 1:Nv, iq in 1:Nq
    #         line = "h $(i-1) $(iq-1) $(vbias[iq, i])\n"
    #         write(file_model, line)
    #     end
    #     close(file_model)
    #     save_list = save_list[save_list .> pearsonCij]
    #     return save_list
    # end

    function save_new_ea(J, vbias, filter, v_model, alphabet, save_list, nsave, pearsonCij, outputpath, label, n_saved)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # code_pearson = string(round(save_list[save_list .<= pearsonCij][end], digits=3))
        # code_density = string(round(density, digits=3))      
        if nsave > 1
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params"*string(n_saved)*".dat" : outputpath*"/params"*string(n_saved)*".dat"
        else
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params.dat" : outputpath*"/params.dat"
        end
        chains_path = (label != nothing) ? outputpath*"/"*label*"_"*"chains.fasta" : outputpath*"/chains.fasta"
        # chains_path, model_path = outputpath * "/trainingchains_density" * code_density * "_pearson" * code_pearson * ".fasta", outputpath * "/model_density" * code_density * "_pearson" * code_pearson * ".dat"
        save_chains_new(J, vbias, v_model, alphabet, chains_path); save_model_new(J, vbias, filter, alphabet, model_path)
        save_list = save_list[save_list .> pearsonCij]
        return save_list
    end

    
    function save_model_chains_ea(J, vbias, filter, v_model, alphabet, outputpath, label)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params.dat" : outputpath*"/params.dat"
        chains_path = (label != nothing) ? outputpath*"/"*label*"_"*"chains.fasta" : outputpath*"/chains.fasta"
        save_chains_new(J, vbias, v_model, alphabet, chains_path); save_model_new(J, vbias, filter, alphabet, model_path)
        return 
    end
    



# FIT MODEL ########################################################################################################################################################################

    function fit_eaDCA(datapath, alphabet, weights, nchains, n_couplings, pseudo_count, lr, nepochs, nsweeps, n_gradsteps, outputpath, target, path_params, path_chains, label, method, seed) # graph
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)
        restore = (path_params != nothing && path_chains != nothing) ? true : false
        (restore == true) && (alphabet !== nothing) ? println("Warning: please make sure that the alphabet and parameters given are coherent.") : nothing
        (restore == false) && (alphabet == nothing) ? error("Provide an alphabet") : nothing

        Random.seed!(seed)
        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        # save_list = linear_division(target, nsave)


        # initialize model and chains        
        if !restore 
            println("reading natural data..."); flush(stdout)
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
            _, fi = oneHotFreqs(v_natural, natural_weights, pseudo_count)
            J, vbias = zeros(Float32, Nq*Nv, Nq*Nv),  Float32.(log.(fi) .- 1/Nq * sum(log.(fi), dims=1))
            v_model = sample_from_profile(vbias, nchains)
        elseif restore 
            J, vbias, v_model, alphabet = restore_model_new(path_params, path_chains) #restore_model(path_params, path_chains, alphabet, Nv, Nq)
            nchains = size(v_model, 3)
            println("reading natural data..."); flush(stdout)
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        end
        model_weights = ones(Float32, nchains)

        # initialize graph
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv) # (graph != nothing) ? initialize_graph_couplingwise(read_graph(graph, Nv, Nq)[1], Nq, Nv) :
        
        # epoch-0 statistics
        tot_params = (Nv*Nq * (Nv*Nq - 1))/2
        fij_natural = oneHotFijSymmFast(v_natural, natural_weights, pseudo_count) 
        fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count); fi_model, _ =  oneHotFreqs(v_model, model_weights, 0)
        cij_natural, cij_model = oneHotCijFast(v_natural, natural_weights, pseudo_count), oneHotCijFast(v_model, model_weights, 0) 
        pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
        println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", perasonFi, "\n"); flush(stdout)
        # save_list = (restore == true) ? save_list[save_list .> pearsonCij] : save_list
        
        # conpute switch time
        switch_time, switch_flag = compute_sampling_switch_time_ea(J, vbias, v_model, nsweeps, method), false

        # training
        n_saved = 1
        for epoch in 1:nepochs
            epoch_time = @elapsed begin         
                # if convergence == true
                #     max_conervgence_step = 1_000
                #     target_couplingwise = 0.80
                #     J, v_model, filter, contact_list, site_degree, switch_flag = do_epoch_with_convergence(J, vbias, n_couplings, filter, contact_list, site_degree, v_model, fij_natural, cij_natural, target_couplingwise, lr, nsweeps, n_gradsteps, max_conervgence_step, pseudo_count, switch_time, switch_flag, method)
                # else
                    J, v_model, filter, contact_list, site_degree, switch_flag = do_epoch(J, vbias, n_couplings, filter, contact_list, site_degree, v_model, fij_natural, lr, nsweeps, n_gradsteps, pseudo_count, switch_time, switch_flag, method)
                # end
            end
            cij_model = oneHotCijFast(v_model, model_weights, 0) 
            density = sum(filter) / (2*tot_params)
            pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
            # if pearsonCij >= save_list[1]
            #     # save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, pearsonCij, outputpath, label, n_saved)
            #     save_list = save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, pearsonCij, outputpath, label, n_saved) 
            #     n_saved += 1
            # end
            # save_list = (pearsonCij >= save_list[1]) ? save_new(J, vbias, filter, v_model, alphabet, save_list, pearsonCij, density, outputpath, label) : save_list  # save_model(J, vbias, v_model, alphabet, save_list, pearsonCij, density, outputpath) : save_list
            println("epoch: ", epoch, " time: ", epoch_time); flush(stdout)
            println("pearson Cij: ", pearsonCij, ", pearson Fi: ", perasonFi); flush(stdout)
            println("model density: ", density, "\n"); flush(stdout)

            (epoch % 50 == 0) ? save_model_chains_ea(J, vbias, filter, v_model, alphabet, outputpath, label) : nothing
            (pearsonCij >= target) ? break : nothing
        end
        save_model_chains_ea(J, vbias, filter, v_model, alphabet, outputpath, label)
        println("training time: ", training_time); flush(stdout)
        close(logfile)
    end

#     export fit_eaDCA
# end