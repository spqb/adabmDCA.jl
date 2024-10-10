# module bmDCAsrc

# using Base.Threads
# using StatsBase
# using Logging
# using Random



# include("utils.jl")
# using .utils 




# COMPUTE SWITCHING TIME ########################################################################################################################################################################

    function compute_sampling_switch_time_bm(J, vbias, v_model, nsweeps, method)
        println("computing sampling switch time..."); flush(stdout)
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        Nq, Nv, Ns = size(v_model) 
        useless1, useless2= copy(v_model), copy(v_model)
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

    function sampling_bm(J, vbias, contact_list, site_degree, v_model::BitArray{3}, nsweeps, switch_time, switch_flag, method)
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        sampling_function_edgewise = gibbs_sampling_edgewise

        Nq, Nv, Ns = size(v_model) 
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        if switch_flag == false 
            sampling_time = @elapsed begin
                @threads :dynamic for t in 1:thread_count
                    start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                    v_model[:, :, start_idx:end_idx] = sampling_function_edgewise(J, vbias, v_model[:, :, start_idx:end_idx], contact_list, site_degree, nsweeps)
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


# GRADIENT UPDATE ########################################################################################################################################################################

    function gradient_update(J::Matrix{Float32}, filter, v_model::BitArray{3}, fij_natural, lr, pseudo_count=0)
        Nq, Nv, Ns = size(v_model)
        fij_model = oneHotFijSymmFast(v_model, ones(Float32, Ns), pseudo_count)
        J .+= lr * ((fij_natural .- fij_model) .* filter)
        return J
    end
    

# DO ONE EPOCH ########################################################################################################################################################################

    function do_epoch(J::Matrix{Float32}, vbias::Matrix{Float32}, filter::BitMatrix, contact_list::Matrix{Int64}, site_degree::Vector{Int64}, v_model::BitArray{3}, fij_natural::Matrix{Float32}, lr::Float64, nsweeps::Int64, switch_time::Float64, switch_flag::Bool, method)
        v_model, switch_flag = sampling_bm(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method)
        J = gradient_update(J, filter, v_model, fij_natural, lr)
        GC.gc()
        return J, v_model, switch_flag
    end
    
 

# SAVE AND RESTORE MODEL ########################################################################################################################################################################

    function save_model(J, vbias, v_model, alphabet, save_list, pearsonCij, outputpath)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # save chains
        code = string(round(save_list[save_list .<= pearsonCij][end], digits=3))
        file_chains = open(outputpath * "/trainingchains_pearson" * code * ".fasta", "w")
        for m in 1:Ns-1
            head = ">chain $m\n"
            line = "$(alphabet[v_cat[:, m]])\n"
            write(file_chains, head); write(file_chains, line)
        end
        head = ">chain $Ns\n"; line = "$(alphabet[v_cat[:, Ns]])"
        write(file_chains, head); write(file_chains, line)
        close(file_chains)


        # save model
        file_model = open(outputpath * "/model_pearson" * code * ".dat", "w")
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
        save_list = save_list[save_list .> pearsonCij]
        return save_list
    end


    function save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, pearsonCij, outputpath, label, n_saved)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # code = string(round(save_list[save_list .<= pearsonCij][end], digits=3))
        if nsave > 1
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params"*string(n_saved)*".dat" : outputpath*"/params"*string(n_saved)*".dat"
        else
            model_path = (label != nothing) ? outputpath*"/"*label*"_"*"params.dat" : outputpath*"/params.dat"
        end
        chains_path = (label != nothing) ? outputpath*"/"*label*"_"*"chains.fasta" : outputpath*"/chains.fasta"

        
        # chains_path, model_path = outputpath * "/trainingchains_pearson" * code * ".fasta", outputpath * "/model_pearson" * code * ".dat"
        save_chains_new(J, vbias, v_model, alphabet, chains_path); save_model_new(J, vbias, filter, alphabet, model_path)
        save_list = save_list[save_list .> pearsonCij]
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

     function fit_bmDCA(datapath, alphabet, weights, nchains, pseudo_count, lr, nepochs, nsweeps, outputpath, target, graph, path_params, path_chains, label, method, seed) # nsave
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
        filter, contact_list, site_degree = (graph != nothing) ? initialize_graph(read_graph(graph, Nv, Nq)[1], Nq, Nv) : initialize_graph(ones(Float32, Nq*Nv, Nq*Nv), Nq, Nv)
        
        # epoch-0 statistics
        fij_natural = oneHotFijSymmFast(v_natural, natural_weights, pseudo_count) 
        fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count); fi_model, _ =  oneHotFreqs(v_model, model_weights, 0)
        cij_natural, cij_model = oneHotCijFast(v_natural, natural_weights, pseudo_count), oneHotCijFast(v_model, model_weights, 0) 
        pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
        println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", perasonFi); flush(stdout)
        # save_list = (restore == true) ? save_list[save_list .> pearsonCij] : save_list
       
        # compute switch time
        switch_time, switch_flag = (graph != nothing) ? (compute_sampling_switch_time_bm(J, vbias, v_model, nsweeps, method), false) : 0.0, true
       
        # training
        n_saved = 1
        training_time = @elapsed begin
            for epoch in 1:nepochs
                epoch_time = @elapsed begin         
                    J, v_model, switch_flag = do_epoch(J, vbias, filter, contact_list, site_degree, v_model, fij_natural, lr, nsweeps, switch_time, switch_flag, method)
                end
                cij_model = oneHotCijFast(v_model, model_weights, 0) 
                pearsonCij, perasonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                
                # if pearsonCij >= save_list[1]
                #     save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label)
                #     # save_list = save_new(J, vbias, filter, v_model, alphabet, save_list, nsave, pearsonCij, outputpath, label, n_saved) 
                #     # n_saved += 1
                # end

                println("epoch: ", epoch, " time: ", epoch_time); flush(stdout)
                println("pearson Cij: ", pearsonCij, ", pearson Fi: ", perasonFi); flush(stdout)
                (epoch % 50 == 0) ? save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label) : nothing
                (pearsonCij >= target) ? break : nothing
            end
        end
        save_model_chains(J, vbias, filter, v_model, alphabet, outputpath, label)
        println("training time: ", training_time); flush(stdout)
        close(logfile)
    end

    # export fit_bmDCA
# end