# COMPUTE SWITCHING TIME ########################################################################################################################################################################

    function compute_sampling_switch_time_re(J, vbias, v_model, nsweeps)
        println("computing sampling switch time..."); flush(stdout)
        Nq, Nv, Ns = size(v_model) 
        useless1, useless2 = copy(v_model), copy(v_model)
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        @threads :dynamic for t in 1:thread_count
                start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                useless1[:, :, start_idx:end_idx] .= gibbs_sampling(J, vbias, useless1[:, :, start_idx:end_idx], nsweeps)
        end
        switch_time = @elapsed begin
            @threads :dynamic for t in 1:thread_count
                start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                useless2[:, :, start_idx:end_idx] .= gibbs_sampling(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps)
            end
        end
        return switch_time
    end


# SAMPLING FUNCTIONS ########################################################################################################################################################################

    function sampling_sa(J, vbias, contact_list, site_degree, v_model::BitArray{3}, nsweeps, switch_time, switch_flag, method)
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        sampling_function_couplingwise = method == "metropolis" ? metropolis_sampling_couplingwise : gibbs_sampling_couplingwise
        Nq, Nv, Ns = size(v_model) 
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        if switch_flag == false 
            sampling_time = @elapsed begin
                @threads :dynamic for t in 1:thread_count
                    start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                    v_model[:, :, start_idx:end_idx] .= sampling_function_couplingwise(J, vbias, v_model[:, :, start_idx:end_idx], contact_list, site_degree, nsweeps)
                end
            end
            if switch_flag == false && sampling_time >= switch_time
                switch_flag = true
                println("\nSwitched sampling to full!\n"); flush(stdout)
            else
                switch_flag = false
            end
        elseif switch_flag == true
            @threads :dynamic for t in 1:thread_count
                    start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
                    v_model[:, :, start_idx:end_idx] .= sampling_function(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps) 
            end
        end
        return v_model, switch_flag
    end


    function sampling_TD(J, vbias, contact_list, site_degree, v_model::BitArray{3}, nsweeps, method)
        sampling_function = method == "metropolis" ? metropolis_sampling : gibbs_sampling
        sampling_function_couplingwise = method == "metropolis" ? metropolis_sampling_couplingwise : gibbs_sampling_couplingwise
        Nq, Nv, Ns = size(v_model) 
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        @threads :dynamic for t in 1:thread_count
            start_idx, end_idx, upper_s = index_interval(t, thread_count, chunk_size, Ns)
            v_model[:, :, start_idx:end_idx] .= sampling_function(J, vbias, v_model[:, :, start_idx:end_idx], nsweeps) 
        end
        return v_model
    end




# SAVE AND RESTORE MODEL ########################################################################################################################################################################

    function save_fasta(v_model, alphabet, label, outputpath)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        # save chains
        path = (label != nothing) ? outputpath*"/"*label*"_resampled_chains.fasta" : outputpath * "/resampled_chains.fasta"
        file_chains = open(path, "w")
        for m in 1:Ns-1
            head = ">chain $m\n"
            line = "$(alphabet[v_cat[:, m]])\n"
            write(file_chains, head); write(file_chains, line)
        end
        head = ">chain $Ns\n"; line = "$(alphabet[v_cat[:, Ns]])"
        write(file_chains, head); write(file_chains, line)
        close(file_chains)
    end

    function restore_params(params_path, Nv, Nq)
        J, vbias = read_graph(params_path, Nv, Nq)
        return Float32.(J), Float32.(vbias) 
    end

    function restore_params_new(params_path)
        J, vbias, alphabet = read_graph_new(params_path)
        return Float32.(J), Float32.(vbias), alphabet
    end

    function restore_chains(path_chains, alphabet)
        data = read_fasta2(path_chains, alphabet)
        v_model = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
        return v_model
    end


# AVERAGE DISTANCE ########################################################################################################################################################################

    # function oneHotHammingDistance(seq1, seq2)
    #     length(seq1) == length(seq2) || throw(ArgumentError("Sequences must have same length"))
    #     return count(seq1 .!= seq2) / 2
    # end

    function averageDistance(v1, v2)
        size(v1) == size(v2) || throw(ArgumentError("Samples must have same size"))
        Nq, Nv, Ns = size(v1)
        v1, v2 = reshape(v1, (Nq*Nv, Ns)), reshape(v2, (Nq*Nv, Ns))
        ave_distance = 0
        for i_s in 1:Ns
            ave_distance += oneHotHammingDistance(v1[:, i_s], v2[:, i_s])
        end
        return ave_distance / Ns
    end

    function sampleDecorrelation_andSTD(v1, v2)
        size(v1) == size(v2) || throw(ArgumentError("Samples must have same size"))
        Nq, Nv, Ns = size(v1)
        v1, v2 = reshape(v1, (Nq*Nv, Ns)), reshape(v2, (Nq*Nv, Ns))
        distance = zeros(Int64, Ns)
        for i_s in 1:Ns
            distance[i_s] = oneHotHammingDistance(v1[:, i_s], v2[:, i_s])
        end
        decor = 1 .- (distance / Nv)
        ave_decor, var_decor = mean_and_var(decor)
        return ave_decor, var_decor
    end

    function sampleDecorrelation(v1, v2)
        Nq, Nv, Ns = size(v1)
        return 1 - averageDistance(v1, v2) / Nv
    end

    

# PLOT DECORRELATION ########################################################################################################################################################################

    function plot_decorrelation(decorrelation_compare, decorrelation_back, outputpath, label, epoch, nsweeps)
        # Define the formatter to use two decimal places
        decimal_formatter = x -> @sprintf("%.2f", x)
        t = length(decorrelation_compare)
        # Create the plot with the specified formatters
        plot(decorrelation_compare, label="SeqID(t)", xlabel="t", ylabel="1 - average_distance", xscale=:log10, yscale=:log10, xformatter=decimal_formatter, yformatter=decimal_formatter, yticks=0:0.1:1, linewidth=:3) #, marker = :o)
        plot!(decorrelation_back, label="SeqID(t,t/2)", xscale=:log10, yformatter=decimal_formatter, yticks=0:0.1:1, linewidth=:3) #, marker = :o)
        title!("Mixing Time")
        xlabel!("epoch x " *string(nsweeps)* " sweeps")
        dec_path = (label != nothing) ? outputpath*"/"*label*"_correlation.png" : outputpath * "/correlation.png"
        savefig(dec_path)
        return 
    end

    function plot_hamming(hamm_dist, outputpath, label)

        histogram(hamm_dist, label="")
        title!("Generate Sequences")
        xlabel!("hamming distance from target sequence")
        
        hamm_path = (label != nothing) ? outputpath*"/"*label*"hamming.png" : outputpath * "/hamming.png"
        savefig(hamm_path)
        return 
    end


# RESAMPLE MODEL ########################################################################################################################################################################

    function sample_DCA(datapath, alphabet, weights, nchains, pseudo_count, nepochs, nsweeps, outputpath, path_params, nmeasure, nmix, mixing_time, label, showplot, seed, method)
        println("used threads: ", Threads.nthreads())
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        println("\nInput arguments")
        println("data: ", datapath)
        println("weights: ", weights)
        println("nchains: ", nchains)
        println("nchains for mixing time: ", nmeasure)
        println("max number of epochs: ", nepochs)
        println("pseudo_count: ", pseudo_count)
        println("nmix: ", nmix)
        println("nsweeps: ", nsweeps)
        (path_params != nothing) ? println("parameters path: ", path_params) : nothing
        println("sampling method: ", method)
        println("random seed: ", seed); flush(stdout)


        path_dec = (label != nothing) ? model_dir*"/"*label*"_autocorrelation.txt" : model_dir * "/autocorrelation.txt"
        decorr_file = open(path_dec, "w") 
        path_Cij = (label != nothing) ? model_dir*"/"*label*"_pearsonCij.txt" : model_dir * "/pearsonCij.txt"
        Cij_file =  open(path_Cij, "w")  
        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)
        (Nq, Nv) = size(vbias)
        data = read_fasta2(datapath, alphabet)
        v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
        (Nq, Nv, Ns) = size(v_natural)
        natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
        natural_weights = natural_weights / sum(natural_weights)
        Random.seed!(0);  
        extracted = sample(1:length(natural_weights), Weights(natural_weights), nmeasure; replace=true)  # rand(Categorical(natural_weights), nmeasure); 
        v_model = copy(v_natural[:, :, extracted]) # v_model = sample_from_profile(vbias, div(nchains, 2))
        model_weights = ones(Float32, size(v_model, 3))
        
        fi_model, _ =  oneHotFreqs(v_model, model_weights, 0); fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count)
        cij_model, cij_natural = oneHotCijFast(v_model, model_weights, 0), oneHotCijFast(v_natural, natural_weights, pseudo_count)
        pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))    
        println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", pearsonFi, "\n"); flush(stdout)
        write(Cij_file, "0 $pearsonCij\n"); flush(Cij_file)

        switch_time, switch_flag = 0, true

        if mixing_time == true
            v_compare, v_back = copy(v_model[:, :, shuffle(1:size(v_model, 3))]), copy(v_model)
            ave1, sigma1 = sampleDecorrelation_andSTD(v_model, v_compare)
            ave2, sigma2 = sampleDecorrelation_andSTD(v_model, v_back)
            decorrelation_compare, decorrelation_back = [ave1], [ave2]
            write(decorr_file, "0 $(decorrelation_compare[1]) $(decorrelation_back[1])\n"); flush(decorr_file)
            
            # compute switch time
            
            t_mix = 0
            for epoch in 1:nepochs
                Random.seed!(epoch); v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
                v_compare = copy(v_model[:, :, shuffle(1:size(v_model, 3))])
                if epoch % 2 == 0
                    println("epoch n: ", div(epoch, 2), ", sweep n: ", div(epoch, 2)*nsweeps); flush(stdout)
                    Random.seed!(div(epoch, 2)); v_back, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_back, nsweeps, switch_time, switch_flag, method) # sampling 
                    ave1, sigma1 = sampleDecorrelation_andSTD(v_model, v_compare)
                    ave2, sigma2 = sampleDecorrelation_andSTD(v_model, v_back)
                    push!(decorrelation_compare, ave1); push!(decorrelation_back, ave2)
                    cij_model = oneHotCijFast(v_model, model_weights, 0) 
                    pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                    println("SeqID independent chains ", round(decorrelation_compare[end], digits=4) , "  chains after t ", round(decorrelation_back[end], digits=4)); flush(stdout)
                    println("pearson Cij ", pearsonCij, " pearson Fi ", pearsonFi); flush(stdout)
                    write(decorr_file, "$epoch (sweeps: $(epoch*nsweeps)) $(decorrelation_compare[end]) $(decorrelation_back[end])\n"); flush(decorr_file)
                    write(Cij_file, "$epoch (sweeps: $(epoch*nsweeps)) $pearsonCij\n"); flush(Cij_file)
                    # (showplot == true) ? plot_decorrelation(decorrelation_compare, decorrelation_back, outputpath) : nothing
                    plot_decorrelation(decorrelation_compare, decorrelation_back, outputpath, label, div(epoch, 2), nsweeps)
                    if abs(ave1 - ave2)  / sqrt(sigma1 + sigma2) < 0.01
                        t_mix = div(epoch, 2) 
                        println("Chains are at equilibrium! mixing time is: ", t_mix * nsweeps, " sweeps \n"); flush(stdout)
                        break
                    end
                end
                GC.gc()
            end
        end

        Random.seed!(seed)
        println("\nInitializing in profile model and sampling...", "\n"); flush(stdout)
        ntot = (mixing_time == true) ? nmix*t_mix : nepochs
        v = sample_from_profile(vbias, nchains, inv_temp)
        model_weights = ones(Float32, size(v, 3))

        for i in 1:ntot
            println("\nepoch n: ", i); flush(stdout)
            v, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v, nsweeps, switch_time, switch_flag, method)
            cij_model = oneHotCijFast(v, ones(Float32, size(v, 3)), 0) 
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
            println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi); flush(stdout)
            write(Cij_file, "$i $pearsonCij\n"); flush(Cij_file)
            if i % 10 == 0
                save_fasta(v, alphabet, label, outputpath)
            end
            GC.gc()
        end
        save_fasta(v, alphabet, label, outputpath)
        close(logfile)
        close(Cij_file)
        close(decorr_file)
    end

# IMPORTANCE SAMPLING

    function importance_sample_DCA(datapath, alphabet, weights, nchains, pseudo_count, nepochs, nsweeps, outputpath, path_params, nmeasure, nmix, mixing_time, label, showplot, seed, method, target_seq_path, theta)
        println("used threads: ", Threads.nthreads())
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        println("\nInput arguments")
        println("data: ", datapath)
        println("weights: ", weights)
        println("target sequence: ", target_seq_path)
        println("theta integration: ", theta)
        println("nchains: ", nchains)
        println("nchains for mixing time: ", nmeasure)
        println("max number of epochs: ", nepochs)
        println("nmix: ", nmix)
        println("nsweeps: ", nsweeps)
        println("pseudo_count: ", pseudo_count)
        (path_params != nothing) ? println("parameters path: ", path_params) : nothing
        println("sampling method: ", method)
        println("random seed: ", seed); flush(stdout)


        path_dec = (label != nothing) ? model_dir*"/"*label*"_autocorrelation.txt" : model_dir * "/autocorrelation.txt"
        decorr_file = open(path_dec, "w") 
        path_Cij = (label != nothing) ? model_dir*"/"*label*"_pearsonCij.txt" : model_dir * "/pearsonCij.txt"
        Cij_file =  open(path_Cij, "w")  
        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)
        (Nq, Nv) = size(vbias)
        data = read_fasta2(datapath, alphabet)


        target_seq = oneHotEncoding(permutedims(read_fasta2(target_seq_path, alphabet), [2, 1]), length(alphabet))
        println(size(target_seq), " ", size(vbias))
        vbias .+= theta .* target_seq 
        target_seq = reshape(target_seq, (Nq*Nv))

        v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
        (Nq, Nv, Ns) = size(v_natural)
        natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
        natural_weights = natural_weights / sum(natural_weights)
        Random.seed!(0);  
        #extracted = sample(1:length(natural_weights), Weights(natural_weights), nmeasure; replace=true)  # rand(Categorical(natural_weights), nmeasure); 
        v_model = sample_from_profile(vbias, nmeasure, inv_temp) #copy(v_natural[:, :, extracted]) 
        model_weights = ones(Float32, size(v_model, 3))
        
        fi_model, _ =  oneHotFreqs(v_model, model_weights, 0); fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count)
        cij_model, cij_natural = oneHotCijFast(v_model, model_weights, 0), oneHotCijFast(v_natural, natural_weights, pseudo_count)
        pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))    
        println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", pearsonFi, "\n"); flush(stdout)
        write(Cij_file, "0 $pearsonCij\n"); flush(Cij_file)

        switch_time, switch_flag = 0, true

        hamm_dist = zeros(size(v_model, 3))

        if mixing_time == true
            v_compare, v_back = copy(v_model[:, :, shuffle(1:size(v_model, 3))]), copy(v_model)
            ave1, sigma1 = sampleDecorrelation_andSTD(v_model, v_compare)
            ave2, sigma2 = sampleDecorrelation_andSTD(v_model, v_back)
            decorrelation_compare, decorrelation_back = [ave1], [ave2]
            write(decorr_file, "0 $(decorrelation_compare[1]) $(decorrelation_back[1])\n"); flush(decorr_file)
            
            # compute switch time
            
            t_mix = 0
            for epoch in 1:nepochs
                Random.seed!(epoch); v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
                v_compare = copy(v_model[:, :, shuffle(1:size(v_model, 3))])
                if epoch % 2 == 0
                    println("epoch n: ", div(epoch, 2), ", sweep n: ", div(epoch, 2)*nsweeps); flush(stdout)
                    Random.seed!(div(epoch, 2)); v_back, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_back, nsweeps, switch_time, switch_flag, method) # sampling 
                    ave1, sigma1 = sampleDecorrelation_andSTD(v_model, v_compare)
                    ave2, sigma2 = sampleDecorrelation_andSTD(v_model, v_back)
                    push!(decorrelation_compare, ave1); push!(decorrelation_back, ave2)
                    cij_model = oneHotCijFast(v_model, model_weights, 0) 
                    pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                    println("SeqID independent chains ", round(decorrelation_compare[end], digits=4) , "  chains after t ", round(decorrelation_back[end], digits=4)); flush(stdout)
                    println("pearson Cij ", pearsonCij, " pearson Fi ", pearsonFi); flush(stdout)
                    write(decorr_file, "$epoch (sweeps: $(epoch*nsweeps)) $(decorrelation_compare[end]) $(decorrelation_back[end])\n"); flush(decorr_file)
                    write(Cij_file, "$epoch (sweeps: $(epoch*nsweeps)) $pearsonCij\n"); flush(Cij_file)
                    # (showplot == true) ? plot_decorrelation(decorrelation_compare, decorrelation_back, outputpath) : nothing
                    plot_decorrelation(decorrelation_compare, decorrelation_back, outputpath, label, div(epoch, 2), nsweeps)

                    v = reshape(v_back, (Nq*Nv, size(v_back, 3)))
                    for i in 1:size(v_model, 3)
                        hamm_dist[i] = oneHotHammingDistance(v[:, i], target_seq)
                    end
                    plot_hamming(hamm_dist, outputpath, label)
                    

                    # if abs(ave1 - ave2)  / sqrt(sigma1 + sigma2) < 0.01
                    #     # if length(decorrelation_compare) < 51
                    #     #     m1, s1 = mean_and_var(decorrelation_compare[1:end])
                    #     #     m2, s2 = mean_and_var(decorrelation_back[1:end])
                    #     # else
                    #     #     m1, s1 = mean_and_var(decorrelation_compare[end-50:end])
                    #     #     m2, s2 = mean_and_var(decorrelation_back[end-50:end])
                    #     # end
                    #     # println("s1: ", s1, ", m1: ", m1, ", s2: ", s2, ", m2: ", m2 )
                       
                    #     m1, m11 = mean(decorrelation_compare), mean(decorrelation_compare[1:end-1])
                    #     m2, m21 = mean(decorrelation_back), mean(decorrelation_back[1:end-1])
                    #     println("delta1: ", abs(m1-m11), ", delta2: ", abs(m2-m21) )
                    #     if abs(m1-m11) < 0.00001 &&  abs(m2-m21) < 0.00001 #(s1 <= m1/1000) && (s2 <= m2/1000) #abs(m1 - m2) < 0.005 #(
                    #         t_mix = div(epoch, 2) 
                    #         println("Chains are at equilibrium! mixing time is: ", t_mix * nsweeps, " sweeps \n"); flush(stdout)
                    #         break
                    #     end
                    # end
                end
                GC.gc()
            end
        end

        Random.seed!(seed)
        println("\nInitializing in profile model and sampling...", "\n"); flush(stdout)
        ntot = (mixing_time == true) ? nmix*t_mix : nepochs
        v = sample_from_profile(vbias, nchains, 1)
        model_weights = ones(Float32, size(v, 3))

        hamm_dist = zeros(size(v, 3))
        for i in 1:ntot
            println("\nepoch n: ", i); flush(stdout)
            v, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v, nsweeps, switch_time, switch_flag, method)
            cij_model = oneHotCijFast(v, ones(Float32, size(v, 3)), 0) 
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
            println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi); flush(stdout)
            write(Cij_file, "$i $pearsonCij\n"); flush(Cij_file)
            v_hamm = reshape(v, (Nq*Nv, size(v, 3)))
            for i in 1:size(v, 3)
                hamm_dist[i] = oneHotHammingDistance(v_hamm[:, i], target_seq)
            end
            plot_hamming(hamm_dist, outputpath, label)
            GC.gc()
            if i % 10 == 0
                save_fasta(v, alphabet, label, outputpath)
            end
        end
       
        save_fasta(v, alphabet, label, outputpath)
        close(logfile)
        close(Cij_file)
        close(decorr_file)
    end







    function TD_integration(datapath, alphabet, weights, nchains, nsweeps, outputpath, path_params, path_chains, label, method, target_seq_path, intstep, theta_max)
        println("used threads: ", Threads.nthreads())

        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        println("\nInput arguments")
        println("data: ", datapath)
        println("weights: ", weights)
        println("target sequence: ", target_seq_path)
        println("max theta integration: ", theta_max)
        println("nchains: ", nchains)
        println("nsweeps: ", nsweeps)
        (path_params != nothing) ? println("parameters path: ", path_params) : nothing
        println("sampling method: ", method)

        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)
        (Nq, Nv) = size(vbias)
        data = read_fasta2(datapath, alphabet)

        target_seq = oneHotEncoding(permutedims(read_fasta2(target_seq_path, alphabet), [2, 1]), length(alphabet))[:,:,1]
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
       
        initial_sweeps = 1_000
        println("sampling to thermalize at theta = 0..."); flush(stdout)
        println("(nsweeps = ", initial_sweeps, ")"); flush(stdout)
        println("initializing chains for theta = 0 ..."); flush(stdout)

        chains_0 = (path_chains != nothing) ? oneHotEncoding(permutedims(read_fasta2(path_chains, alphabet), [2, 1]), length(alphabet))[:, :, 1:nchains] : sample_from_profile(vbias, nchains, inv_temp)
        println("size of chains_0 :", size(chains_0)); flush(stdout)
        v = sampling_TD(J, vbias, contact_list, site_degree, chains_0, initial_sweeps, method)
        ave_ene = mean(compute_energy(J, vbias, v))
        println("average energy sample theta 0: ", ave_ene)

        println("sampling to thermalize at theta_max = ", theta_max); flush(stdout)
        vbias_theta_max = vbias + theta_max .* target_seq
        v_max = sampling_TD(J, vbias_theta_max, contact_list, site_degree, sample_from_profile(vbias_theta_max, nchains, inv_temp), 100, method)
        v2_max = reshape(v_max, (Nq*Nv, nchains))
        hamm_dist_max = zeros(nchains)
        for i in 1:nchains
            hamm_dist_max[i] = oneHotHammingDistance(v2_max[:, i], reshape(target_seq, Nq*Nv))
        end
        n_wt = count(hamm_dist_max .== 0)
        p_wt = n_wt / nchains
        println(p_wt * 100, "% sequences collapse to wt"); flush(stdout)
        while p_wt <= 0.1
            theta_max = theta_max + 0.1 * theta_max
            println("Number of sequences collapsed to WT is less then 10%. Increasing theta max to: ", theta_max); flush(stdout)
            vbias_theta_max = Float32.(vbias + theta_max .* target_seq)
            v_max = sampling_TD(J, vbias_theta_max, contact_list, site_degree, sample_from_profile(vbias_theta_max, nchains, inv_temp), 100, method)
            v2_max = reshape(v_max, (Nq*Nv, nchains))
            hamm_dist_max = zeros(nchains)
            for i in 1:nchains
                hamm_dist_max[i] = oneHotHammingDistance(v2_max[:, i], reshape(target_seq, Nq*Nv))
            end
            n_wt = count(hamm_dist_max .== 0)
            p_wt = n_wt / nchains
            println(p_wt * 100, "% sequences collapse to wt"); flush(stdout)
        end

        F_max = log(p_wt) + mean(compute_energy(J, vbias_theta_max, v_max[:, :, hamm_dist_max .== 0]))
        seqID, F, S, integral = zeros(size(v, 3)), F_max, 0, 0
        thetas = collect(range(0, theta_max, length=intstep))  
        println("F_theta_max = ", F_max)
        factor = theta_max / (2*intstep)

        t = @elapsed begin
            i = 1
            theta = thetas[i]
            println("\n", i," theta : ", theta)
            vbias_theta = Float32.(vbias +  theta .* target_seq)
            v = sampling_TD(J, vbias_theta, contact_list, site_degree, v, nsweeps, method)
            v2 = reshape(v, (Nq*Nv, nchains))
            for i in 1:nchains
                seqID[i] = Nv .- oneHotHammingDistance(v2[:, i], reshape(target_seq, Nq*Nv))
            end
                
            integral += factor * mean(seqID) 
            F += factor * mean(seqID) 
            S = ave_ene - F
            println("mean(seqID): ", mean(seqID))  
            println("integral: ", integral)
            println("F: ", F)
            println("in itinere entropy: ", S)
            GC.gc(); flush(stdout)
        end
        println("elapsed time: ", t)
    
        for theta in thetas[2:end-1]
            t = @elapsed begin
                i += 1
                println("\n", i," theta : ", theta) 
                vbias_theta = Float32.(vbias +  theta .* target_seq)
                v = sampling_TD(J, vbias_theta, contact_list, site_degree, v, nsweeps, method)
                v2 = reshape(v, (Nq*Nv, nchains))
                for i in 1:nchains
                    seqID[i] = Nv .- oneHotHammingDistance(v2[:, i], reshape(target_seq, Nq*Nv))
                end
                F +=  2 * factor * mean(seqID) 
                integral +=  2 * factor * mean(seqID) 
                S = ave_ene - F
                println("mean(seqID): ", mean(seqID))
                println("integral: ", integral)
                println("F: ", F)
                println("in itinere entropy: ", S)
                GC.gc(); flush(stdout)
            end
            println("elapsed time: ", t)
        end

        t = @elapsed begin
            i = intstep
            println("\n", i," theta : ", thetas[i]) 
            theta = thetas[i]
            vbias_theta = Float32.(vbias +  theta .* target_seq)
            v = sampling_TD(J, vbias_theta, contact_list, site_degree, v, nsweeps, method)
            v2 = reshape(v, (Nq*Nv, nchains))
            for i in 1:nchains
                seqID[i] = Nv .- oneHotHammingDistance(v2[:, i], reshape(target_seq, Nq*Nv))
            end
            F += factor * mean(seqID) 
            integral += factor * mean(seqID) 
            S = ave_ene - F
            println("mean(seqID): ", mean(seqID))
            println("integral: ", integral) 
            println("F: ", F) 
            println("in itinere entropy: ", S) 
            GC.gc(); flush(stdout)
        end
        println("elapsed time: ", t)
        println("\nFinal Entropy: ", S) 
    end








































function resampling_old(datapath, alphabet, weights, nchains, pseudo_count, nepochs, nsweeps, outputpath, target, path_params, path_chains, mixing_time, label, method)
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        path_dec = (label != nothing) ? model_dir*"/"*label*"_autocorrelation.txt" : model_dir * "/autocorrelation.txt"
        decorr_file = open(path_dec, "w") 

        path_Cij = (label != nothing) ? model_dir*"/"*label*"_pearsonCij.txt" : model_dir * "/pearsonCij.txt"
        Cij_file = (datapath != nothing) ? open(path_Cij, "w") : nothing

        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)

        (Nq, Nv) = size(vbias)

        if datapath != nothing
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        end
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
        Random.seed!(0); v_model = (path_chains != nothing) ? restore_chains(path_chains, alphabet) : sample_from_profile(vbias, div(nchains, 2))
        model_weights = ones(Float32, size(v_model, 3))
        

        v_compare, v_back = copy(v_model[:, :, shuffle(1:size(v_model, 3))]), copy(v_model)
        decorrelation_compare, decorrelation_back = [sampleDecorrelation(v_model, v_compare)], [sampleDecorrelation(v_model, v_back)]
        write(decorr_file, "0 $(decorrelation_compare[1]) $(decorrelation_back[1])\n"); flush(decorr_file)
    
        if datapath != nothing
            fi_model, _ =  oneHotFreqs(v_model, model_weights, 0); fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count)
            cij_model, cij_natural = oneHotCijFast(v_model, model_weights, 0), oneHotCijFast(v_natural, natural_weights, pseudo_count)
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))    
            println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", pearsonFi, "\n"); flush(stdout)
            write(Cij_file, "0 $pearsonCij\n"); flush(Cij_file)
        end
    
        # compute switch time
        switch_time, switch_flag = 0, true
        decorr_time = 0
        for epoch in 1:nepochs
            Random.seed!(epoch); v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
            
            
            v_compare = copy(v_model[:, :, shuffle(1:size(v_model, 3))])
            if epoch % 2 == 0
                println("sweep n: ", div(epoch, 2)*nsweeps); flush(stdout)
                Random.seed!(div(epoch, 2)); v_back, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_back, nsweeps, switch_time, switch_flag, method) # sampling 
                push!(decorrelation_compare, sampleDecorrelation(v_model, v_compare)); push!(decorrelation_back, sampleDecorrelation(v_model, v_back))
                println("1-average_distance - independent chains: ", round(decorrelation_compare[end], digits=4) , ", chains after t: ", round(decorrelation_back[end], digits=4)); flush(stdout)
                plot_decorrelation(decorrelation_compare, decorrelation_back, epoch, outputpath)
                write(decorr_file, "$epoch $(decorrelation_compare[end]) $(decorrelation_back[end])\n"); flush(decorr_file)
                if (decorrelation_compare[end] - decorrelation_back[end]) * (decorrelation_compare[end-1] - decorrelation_back[end-1]) < 0 
                    decorr_time = div(epoch, 2) 
                    println("Chains are at equilibrium!\n"); flush(stdout)
                    break
                end
                if datapath != nothing
                    cij_model = oneHotCijFast(v_model, model_weights, 0) 
                    pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                    println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi, "\n"); flush(stdout)
                    write(Cij_file, "$epoch $pearsonCij\n"); flush(Cij_file)
                end
            end
            
        end

        v = cat(v_model, v_back, dims=3)
        for i in 1:decorr_time
            v, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v, nsweeps, switch_time, switch_flag, method)
            if datapath != nothing
                cij_model = oneHotCijFast(v, ones(Float32, size(v, 3)), 0) 
                pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi, "\n"); flush(stdout)
                write(Cij_file, "$i $pearsonCij\n"); flush(Cij_file)
            end
        end
        save_fasta(v, alphabet, outputpath)
        close(logfile)
        close(Cij_file)
        close(decorr_file)
    end

    function resampling(datapath, alphabet, weights, nchains, pseudo_count, nepochs, nsweeps, outputpath, target, path_params, path_chains, mixing_time, label, method)
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        path_dec = (label != nothing) ? model_dir*"/"*label*"_autocorrelation.txt" : model_dir * "/autocorrelation.txt"
        decorr_file = (mixing_time == true) ? open(path_dec, "w") : nothing

        path_Cij = (label != nothing) ? model_dir*"/"*label*"_pearsonCij.txt" : model_dir * "/pearsonCij.txt"
        Cij_file = (datapath != nothing) ? open(path_Cij, "w") : nothing

        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)

        (Nq, Nv) = size(vbias)

        if datapath != nothing
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        end
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)
        v_model = (path_chains != nothing) ? restore_chains(path_chains, alphabet) : sample_from_profile(vbias, div(nchains, 2))
        sample_size = size(v_model, 3)
        model_weights = ones(Float32, size(v_model, 3))
        if mixing_time == true
            v_shuffled = copy(v_model[:, :, shuffle(1:sample_size)])
            v_compare, v_1, v_2 = copy(v_model), copy(v_shuffled[:, :, 1:div(sample_size, 2)]), copy(v_shuffled[:, :, div(sample_size, 2)+1:end])
            decorrelation_compare, decorrelation_back = [sampleDecorrelation(v_model, v_compare)], [sampleDecorrelation(v_1, v_2)]
            write(decorr_file, "0 $(decorrelation_compare[1]) $(decorrelation_back[1])\n"); flush(decorr_file)
        end
        if datapath != nothing
            fi_model, _ =  oneHotFreqs(v_model, model_weights, 0); fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count)
            cij_model, cij_natural = oneHotCijFast(v_model, model_weights, 0), oneHotCijFast(v_natural, natural_weights, pseudo_count)
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))    
            println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", pearsonFi, "\n"); flush(stdout)
            write(Cij_file, "0 $pearsonCij\n"); flush(Cij_file)
        end
    
        # compute switch time
        switch_time, switch_flag = compute_sampling_switch_time_re(J, vbias, v_model, nsweeps), false

        cross_epoch, cross_flag = 0, false
        epoch, correlated = 1, true
        while correlated
            println("epoch: ", epoch); flush(stdout)
            if mixing_time == true
                v = cat(v_model, v_compare, dims=3)
                time = @elapsed begin
                v, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v, nsweeps, switch_time, switch_flag, method)
                end
                println("sampled time: ", time); flush(stdout)
                v_model, v_compare = v[:, :, 1:sample_size], v[:, :, sample_size+1:end]
                v_shuffled = copy(v_model[:, :, shuffle(1:sample_size)])
                v_1, v_2 = copy(v_shuffled[:, :, 1:div(sample_size, 2)]), copy(v_shuffled[:, :, div(sample_size, 2)+1:end])
                push!(decorrelation_compare, sampleDecorrelation(v_model, v_compare)); push!(decorrelation_back, sampleDecorrelation(v_1, v_2))
                println("1-average_distance - chains after t: ", round(decorrelation_compare[end], digits=4) , ", independent chains: ", round(decorrelation_back[end], digits=4)); flush(stdout)
                write(decorr_file, "$epoch $(decorrelation_compare[end]) $(decorrelation_back[end])\n"); flush(decorr_file)
                if (decorrelation_compare[end] - decorrelation_back[end]) * (decorrelation_compare[end-1] - decorrelation_back[end-1]) < 0 && cross_flag == false
                    println("Chains are at equilibrium!"); flush(stdout)
                    cross_epoch, cross_flag = epoch, true
                end
                plot_decorrelation(decorrelation_compare, decorrelation_back, epoch, outputpath)
                (epoch == 2*cross_epoch) ? correlated = false : nothing
                epoch += 1
            else
                v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
            end
            if datapath != nothing
                cij_model = oneHotCijFast(v_model, model_weights, pseudo_count) 
                pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi, "\n"); flush(stdout)
                write(Cij_file, "$epoch $pearsonCij\n"); flush(Cij_file)
            end
        end
        v_save = (mixing_time == true) ? cat(v_model, v_compare, dims=3) : v_model 
        save_fasta(v_save, alphabet, outputpath)
        close(logfile)
        close(Cij_file)
        close(decorr_file)
    end

    function resampling_nat(datapath, alphabet, weights, nchains, pseudo_count, nepochs, nsweeps, outputpath, target, path_params, path_chains ,mixing_time, label, method)
        model_dir = outputpath; (!isdir(model_dir)) ? mkdir(model_dir) : nothing
        path_log = (label != nothing) ? model_dir*"/"*label*"_adabmDCA.log" : model_dir * "/adabmDCA.log"
        logfile = open(path_log, "w"); redirect_stdout(logfile)

        path_dec = (label != nothing) ? model_dir*"/"*label*"_autocorrelation.txt" : model_dir * "/autocorrelation.txt"
        decorr_file = (mixing_time == true) ? open(path_dec, "w") : nothing

        path_Cij = (label != nothing) ? model_dir*"/"*label*"_pearsonCij.txt" : model_dir * "/pearsonCij.txt"
        Cij_file = open(path_Cij, "w") 

        inv_temp = 1
        alphabet = set_alphabet(alphabet)
        J, vbias, alphabet = restore_params_new(path_params)

        (Nq, Nv) = size(vbias)

        if datapath != nothing
            data = read_fasta2(datapath, alphabet)
            v_natural = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
            (Nq, Nv, Ns) = size(v_natural)
            natural_weights, Meff, pseudo_count = assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        end
        filter, contact_list, site_degree = initialize_graph_couplingwise(J, Nq, Nv)

        # mixing time
        v_nat_mt, v_model_mt = copy(v_natural), sample_from_profile(vbias, 10_000)
        sample_size = size(v_model_mt, 3)
        model_weights = ones(Float32, size(v_model_mt, 3))
        if mixing_time == true
            v_model_shuffled, v_nat_shuffled = copy(v_model_mt[:, :, shuffle(1:sample_size)]), copy(v_nat_mt[:, :, shuffle(1:Ns)])
            decorrelation_model, decorrelation_nat = [sampleDecorrelation(v_model_mt, v_model_shuffled)], [sampleDecorrelation(v_nat_mt, v_nat_shuffled)]
            write(decorr_file, "0 $(decorrelation_model[1]) $(decorrelation_nat[1])\n"); flush(decorr_file)
        end
        if datapath != nothing
            fi_model, _ =  oneHotFreqs(v_model_mt, model_weights, 0); fi_natural, _ = oneHotFreqs(v_natural, natural_weights, pseudo_count)
            cij_model, cij_natural = oneHotCijFast(v_model_mt, model_weights, 0), oneHotCijFast(v_natural, natural_weights, pseudo_count)
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))    
            println("t = 0 - Pearson Cij: ", pearsonCij, ", Pearson Fi: ", pearsonFi, "\n"); flush(stdout)
            write(Cij_file, "0 $pearsonCij\n"); flush(Cij_file)
        end
    
        # compute switch time
        # switch_time, switch_flag = compute_sampling_switch_time_re(J, vbias, v_model, nsweeps), false
        switch_flag, switch_time = true, 0

        cross_epoch, cross_flag = 0, false
        epoch, correlated = 1, true
        while correlated
            println("epoch: ", epoch); flush(stdout)
            if mixing_time == true
                v_mt = cat(v_model_mt, v_nat_mt, dims=3)
                time = @elapsed begin
                    v_mt, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_mt, nsweeps, switch_time, switch_flag, method)
                end
                println("sampling time: ", time); flush(stdout)

                v_model_mt, v_nat_mt = v_mt[:, :, 1:sample_size], v_mt[:, :, sample_size+1:end]
                v_model_shuffled, v_nat_shuffled = copy(v_model_mt[:, :, shuffle(1:sample_size)]), copy(v_nat_mt[:, :, shuffle(1:Ns)])

                push!(decorrelation_model, sampleDecorrelation(v_model_mt, v_model_shuffled)); push!(decorrelation_nat, sampleDecorrelation(v_nat_mt, v_nat_shuffled))
                println("1-average_distance - chains after t: ", round(decorrelation_model[end], digits=4) , ", naturals after t: ", round(decorrelation_nat[end], digits=4)); flush(stdout)
                write(decorr_file, "$epoch $(decorrelation_model[end]) $(decorrelation_nat[end])\n"); flush(decorr_file)
                if (decorrelation_model[end] - decorrelation_nat[end]) * (decorrelation_model[end-1] - decorrelation_nat[end-1]) < 0 && cross_flag == false
                    println("Chains are at equilibrium!"); flush(stdout)
                    cross_epoch, correlated = epoch, false
                end
                correlated = true
                plot_decorrelation(decorrelation_model, decorrelation_nat, epoch, outputpath)
                epoch += 1
            else
                v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method) # sampling 
            end
            if datapath != nothing
                cij_model = oneHotCijFast(v_model_mt, model_weights, 0) 
                pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
                println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi, "\n"); flush(stdout)
                write(Cij_file, "$epoch $pearsonCij\n"); flush(Cij_file)
            end
        end

        println("\nMIXING TIME COMPUTED: ", cross_epoch ,"\n")
        v_model = sample_from_profile(vbias, nchains)
        for epoch in 1:(2*cross_epoch)
            v_model, switch_flag = sampling_sa(J, vbias, contact_list, site_degree, v_model, nsweeps, switch_time, switch_flag, method)
            cij_model = oneHotCijFast(v_model, ones(Float32, nchains), 0) 
            pearsonCij, pearsonFi = cor(vec(cij_model), vec(cij_natural)), cor(vec(fi_natural), vec(fi_model))
            println("pearson Cij: ", pearsonCij, ", pearson Fi: ", pearsonFi, "\n"); flush(stdout)
        end

        save_fasta(v_model, alphabet, outputpath)
        close(logfile)
        close(Cij_file)
        close(decorr_file)
    end

#     export sampleDecorrelation, plot_decorrelation, save_fasta, sample_DCA
# end