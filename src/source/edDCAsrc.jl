# module edDCAsrc

using Base.Threads
using LinearAlgebra
using StatsBase
using Distances
using Random

include("utils.jl")
# using .utils 




# COMPUTE SWITCHING TIME ########################################################################################################################################################################

    function compute_sampling_switch_time(J, vbias, v_model, nsweeps, method)
        println("computing sampling switch time..."); flush(stdout)
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

    function sampling(J, vbias, contact_list, site_degree, v_model::BitArray{3}, nsweeps, switch_time, switch_flag, method)
        Nq, Nv, Ns = size(v_model) 
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        sampling_function_couplingwise = method == "metropolis" ? metropolis_sampling_couplingwise : gibbs_sampling_couplingwise
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
                    v_model[:, :, start_idx:end_idx] = sampling_function(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps) #gibbs_sampling
            end
        end
        return v_model, switch_flag
    end


# COUPLINGS ACTIVATION  ########################################################################################################################################################################

    function update_graph_dec(J, chosen_couplings::Matrix{Float32}, filter, contact_list, site_degree)
        for idx in 1:size(chosen_couplings, 2)
            iq, jq  = Int64(chosen_couplings[1, idx]), Int64(chosen_couplings[2, idx])
            if filter[iq, jq] == true && filter[jq, iq] == true
                idx_iq, idx_jq = findfirst(isequal(jq), contact_list[iq, :]), findfirst(isequal(iq), contact_list[jq, :]) # find position in contact list
                contact_list[iq, idx_iq], contact_list[jq, idx_jq] = 0, 0 # set it to zero
                contact_list[iq, 1:site_degree[iq]], contact_list[jq, 1:site_degree[jq]] = sort(contact_list[iq, 1:site_degree[iq]], rev=true), sort(contact_list[jq, 1:site_degree[jq]], rev=true) # re-order the contact list
                site_degree[iq] -= 1; site_degree[jq] -= 1
                filter[iq, jq], filter[jq, iq] = false, false
                J[iq, jq], J[jq, iq] = 0, 0
            else
                error("Trying to remove a coupling that is not present.")
            end
        end       
        return J, filter, contact_list, site_degree
    end

    function compute_dkl_matrix_dec(pij_model, fij_target, filter)
        nonzero_iqjqs = findall(x -> x == 1, triu(filter))
        dkl_matrix = zeros(Float32, 3, length(nonzero_iqjqs))
        thread_count = nthreads(); chunk_size = div(length(nonzero_iqjqs), thread_count)
        @threads :dynamic for t in 1:thread_count       
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, length(nonzero_iqjqs))
            for (idx, iqjq) in zip(start_idx:end_idx, nonzero_iqjqs[start_idx:end_idx])
                iq, jq, f_tar, f_mod = iqjq[1], iqjq[2], zeros(Float32, 2), zeros(Float32, 2)
                f_tar[1], f_tar[2] = fij_target[iq, jq], 1 - fij_target[iq, jq]
                f_mod[1], f_mod[2] = pij_model[iq, jq], 1 - pij_model[iq, jq]
                dkl_matrix[:,idx] = [iq, jq, Distances.kl_divergence(f_tar, f_mod)]
            end
        end
        return dkl_matrix
    end

    function compute_dkl_matrix_dec2(pij_model, J, filter)
        nonzero_iqjqs = findall(x -> x == 1, triu(filter))
        dkl_matrix = zeros(Float32, 3, length(nonzero_iqjqs))
        thread_count = nthreads(); chunk_size = div(length(nonzero_iqjqs), thread_count)
        @threads :dynamic for t in 1:thread_count       
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, length(nonzero_iqjqs))
            for (idx, iqjq) in zip(start_idx:end_idx, nonzero_iqjqs[start_idx:end_idx])
                iq, jq = iqjq[1], iqjq[2]
                d = J[iq, jq]* pij_model[iq, jq] + log(pij_model[iq, jq] * exp(- J[iq, jq]) + 1 - pij_model[iq, jq])
                dkl_matrix[:,idx] = [iq, jq, d]
            end
        end
        return dkl_matrix
    end


# GRADIENT UPDATE ########################################################################################################################################################################

    function gradient_update(J::Matrix{Float32}, filter, v_model::BitArray{3}, fij_natural, pij_model, lr)
        J .+= lr * ((fij_natural .- pij_model) .* filter)
        return J
    end
    

# DO ONE EPOCH ########################################################################################################################################################################
    
    function do_convergence(J, vbias, filter, contact_list, site_degree, v_model, nsweeps, fij_natural, cij_natural, target_cij, lr, pseudo_count, max_conervgence_step, method, switch_time=0, switch_flag=true)
        slope_err = 0.5
        cij_model = oneHotCijFast(v_model, ones(size(v_model, 3)), 0) 
        p, slope = cor(vec(cij_model), vec(cij_natural)), compute_slope(cij_model, cij_natural)
        (p <= target_cij) ? println("initial cij: ", p, " --> starting convergence...") : nothing; flush(stdout)
        (abs(slope - 1) > slope_err) ? println("initial slope: ", slope, " --> starting convergence...") : nothing; flush(stdout)
        Nq, Nv, Ns = size(v_model)
        pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), 0)
        count = 0
        
        while (count <= max_conervgence_step) && (!(p <= target_cij) && (abs(slope - 1) > slope_err) || (p <= target_cij))
            count += 1
            v_model, _ = sampling(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method)
            pij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), 0)
            cij_model = oneHotCijFast(v_model, ones(size(v_model, 3)), 0) 
            exp_eps = maximum(abs.(cij_natural .- cij_model))
            p, slope = cor(vec(cij_model), vec(cij_natural)), compute_slope(cij_model, cij_natural)
            p_active = cor(vec(cij_model[filter]), vec(cij_natural[filter]))
            println("step: ", count, "; slope: ", slope, ", pearson Cij: ", p, ", pearson Cij active: ", p_active); flush(stdout)
            J = gradient_update(J, filter, v_model, fij_natural, pij_model, lr) # update parameters
        end
        return J, v_model, pij_model
    end

    function do_decimation(J, pij_model, fij_target, filter, contact_list, site_degree, th_d)
        dkl_matrix = compute_dkl_matrix_dec2(pij_model, J, filter)
        dkl_matrix = dkl_matrix[:, sortperm(dkl_matrix[3, :])] # sort the dkl_matrix
        n_cps = Int64(floor(th_d * size(dkl_matrix, 2)))
        chosen_couplings = dkl_matrix[1:2, 1:n_cps] # chose the first th_d % of activated coupling
        println("number of decimated couplings: ", n_cps), flush(stdout)
        J, filter, contact_list, site_degree = update_graph_dec(J, chosen_couplings, filter, contact_list, site_degree)
        return J, filter, contact_list, site_degree
    end 

    function do_epoch(J::Matrix{Float32}, vbias::Matrix{Float32}, filter::BitMatrix, contact_list::Matrix{Int64}, site_degree::Vector{Int64}, v_model::BitArray{3}, fij_natural::Matrix{Float32}, cij_natural::Matrix{Float32}, target_cij, lr::Float64, nsweeps::Int64, pseudo_count, th_d, max_conervgence_step, method)
        # convergence
        J, v_model, pij_model = do_convergence(J, vbias, filter, contact_list, site_degree, v_model, nsweeps, fij_natural, cij_natural, target_cij, lr, pseudo_count, max_conervgence_step, method)
        # decimation
        J, filter, contact_list, site_degree = do_decimation(J, pij_model, fij_natural, filter, contact_list, site_degree, th_d)
        # sampling with decimated model
        v_model, _ = sampling(J, vbias, contact_list, site_degree, v_model, nsweeps, 0, true, method)
        GC.gc()
        return J, v_model, filter, contact_list, site_degree
    end
 

# SAVE AND RESTORE MODEL ########################################################################################################################################################################

    function save_model(J, vbias, v_model, alphabet, save_list, sparsity, pearsonCij, outputpath)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # save chains
        code_density = string(round(sparsity, digits=3))
        code_pearson = string(round(pearsonCij, digits=3))
        file_chains = open(outputpath * "/trainingchains_density" * code_density * "_pearson" * code_pearson * ".fasta", "w")
        for m in 1:Ns-1
            head = ">chain $m\n"
            line = "$(alphabet[v_cat[:, m]])\n"
            write(file_chains, head); write(file_chains, line)
        end
        head = ">chain $Ns\n"; line = "$(alphabet[v_cat[:, Ns]])"
        write(file_chains, head); write(file_chains, line)
        close(file_chains)
        # save model
        file_model = open(outputpath * "/model_density" * code_density * "_pearson" * code_pearson * ".dat", "w")
        for i in 1:Nv, j in i+1:Nv
            for iq in 1:Nq, jq in 1:Nq
            line = "J $(i-1) $(j-1) $(iq-1) $(jq-1) $(J[id(i, iq, Nq), id(j, jq, Nq)])\n"
                write(file_model, line)
            end
        end
        for i in 1:Nv, iq in 1:Nq
            line = "h $(i-1) $(iq-1) $(vbias[iq, i])\n"
            write(file_model, line)
        end
        close(file_model)
        save_list = save_list[save_list .> 1-sparsity]
        return save_list
    end

    function save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, density, outputpath, label, n_saved)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # code_pearson = string(round(pearsonCij, digits=3))
        # code_density = string(round(density, digits=3))      
        if nsave > 1
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params"*string(n_saved)*".dat" : outputpath*"/params"*string(n_saved)*".dat"
        else
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params.dat" : outputpath*"/params.dat"
        end
        chains_path = (label != nothing) ? outputpath*"/"*label*"_"*"chains.fasta" : outputpath*"/chains.fasta"


        # chains_path, model_path = outputpath * "/trainingchains_density" * code_density * "_pearson" * code_pearson * ".fasta", outputpath * "/model_density" * code_density * "_pearson" * code_pearson * ".dat"
        save_chains_new(J, vbias, v_model, alphabet, chains_path); save_model_new(J, vbias, filter, alphabet, model_path)
        # println("save_list: ", save_list, ", density:", density)
        save_list = save_list[save_list .> 1-density]
        return save_list
    end

    function save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params.dat" : outputpath*"/params.dat"
        chains_path = (label != nothing) ? outputpath*"/"*label*"_"*"chains.fasta" : outputpath*"/chains.fasta"
        save_chains_new(J, vbias, v_model, alphabet, chains_path); save_model_new(J, vbias, filter, alphabet, model_path)
        return 
    end
 
    

# FIT MODEL ########################################################################################################################################################################

    function fit_edDCA(datapath, path_params, path_chains, target_density, target_Cij, drate, lr, nsweeps, nchains, alphabet, weights, pseudo_count, nepochs, outputpath, max_conervgence_step, label, method, seed)
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)
        restore = (path_params != nothing && path_chains != nothing) ? true : false
        (restore == false && alphabet == nothing) ? error("Provide an alphabet") : nothing
        (restore == false && nchains == nothing) ? error("Provide nchains") : nothing

        Random.seed!(seed)
        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        # save_list = linear_division(1-target_density, nsave)
        
        if !restore 
            println("reading natural data..."); flush(stdout)
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
            _, fi = oneHotFreqs(v_natural, natural_weights, pseudo_count)
            J, vbias = zeros(Float32, Nq*Nv, Nq*Nv),  Float32.(log.(fi) .- 1/Nq * sum(log.(fi), dims=1))
            v_model = sample_from_profile(vbias, nchains)
            # initialize graph
            filter, contact_list, site_degree = initialize_graph_couplingwise(ones(Float32, Nq*Nv, Nq*Nv), Nq, Nv)
        elseif restore 
            J, vbias, v_model, alphabet = restore_model_new(path_params, path_chains) #restore_model(path_params, path_chains, alphabet, Nv, Nq)
            (nchains == nothing) ? nchains = size(v_model, 3) : v_model = v_model[:, :, rand(1:size(v_model, 3), nchains)]
            println("reading natural data..."); flush(stdout)
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
            # initialize graph
            filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
        end
        model_weights = ones(Float32, nchains)
        tot_params = (Nv*Nq * (Nv*Nq - 1))/2
        
        # epoch-0 statistics
        fij_natural = oneHotFijSymmFast(v_natural, natural_weights, pseudo_count) 
        fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count); fi_model, _ =  oneHotFreqs(v_model, model_weights, 0)
        cij_natural, cij_model = oneHotCijFast(v_natural, natural_weights, pseudo_count), oneHotCijFast(v_model, model_weights, 0) 
        pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
        density = sum(filter) / (2*tot_params)
        println("\nt = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", perasonFi, ", Model Density: ", round(density, digits=2), "\n"); flush(stdout)

        # training
        n_saved = 1
        training_time = @elapsed begin
            for epoch in 1:nepochs
                epoch_time = @elapsed begin                
                    J, v_model, filter, contact_list, site_degree = do_epoch(J, vbias, filter, contact_list, site_degree, v_model, fij_natural, cij_natural, target_Cij, lr, nsweeps, pseudo_count, drate, max_conervgence_step, method) 
                end
                density = sum(filter) / (2*tot_params)
                cij_model = oneHotCijFast(v_model, model_weights, 0) 
                pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                # if 1-density >= save_list[1]
                #     save_list = save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, density, outputpath, label, n_saved) 
                #     n_saved += 1
                # end
                # save_list = (1-sparsity >= save_list[1]) ? save_new(J, vbias, filter, v_model, alphabet, save_list, pearsonCij, sparsity, outputpath, label) : save_list #save_model(J, vbias, v_model, alphabet, save_list, sparsity, pearsonCij, outputpath) : save_list
                println("fraction of parameters: ", density); flush(stdout)
                println("pearson Cij: ", pearsonCij, ", pearson Fi: ", perasonFi); flush(stdout)
                println("epoch: ", epoch, " time: ", epoch_time, "\n"); flush(stdout)
                (epoch % 50 == 0) ? save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label) : nothing
                (density <= target_density) ? break : nothing
            end
        end 
        println("training time: ", training_time); flush(stdout)
        save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label)
        close(logfile)
    end

#     export fit_edDCA

# end