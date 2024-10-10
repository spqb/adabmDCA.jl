# module utils 

    # using Base.Threads
    # using Random
    # using StatsBase




# INDEXING FUNCTION ########################################################################################################################################################################

    id(i, a, q) = (i .- 1).*q .+ a  

    function index_interval(t, n_threads, chunk_size, total_size)
        start_idx = (t - 1) * chunk_size + 1
        end_idx = t == n_threads ? total_size : t * chunk_size
        upper_idx = end_idx - start_idx + 1
        return start_idx, end_idx, upper_idx
    end

# UTILITY FUNCTIONS ########################################################################################################################################################################

    function linear_division(x, n)
        step = x / n
        intervals = [i * step for i in 1:n]
        return intervals
    end

    function exponential_division_euler(x, n)
        intervals = Float64[]
        for i in 0:n
            val = x * (exp(i/n) - 1) / (exp(1) - 1)
            push!(intervals, val)
        end
        return intervals
    end

    function compute_slope(x, y)
        x, y = vec(x), vec(y)
        N = length(x)
        return (N * sum(x.*y) - sum(x)*sum(y)) / (N * sum(x.*x) - sum(x)*sum(x))
    end

# DATA FORMAT ########################################################################################################################################################################

    function oneHotEncoding(V, Nq)
        Nv, Ns = size(V)
        oneHotV = BitArray(undef, Nq, Nv, Ns)
        all_v = collect(1:Nq)
        for i_s in 1:Ns, i_v in 1:Nv
            oneHotV[:, i_v, i_s] = (all_v .== V[i_v, i_s])
        end
        return oneHotV
    end

    function oneHot2Categorical(oneHotV, Nq)
        Nq, Nv, Ns = size(oneHotV)
        V = zeros(Int16, Nv, Ns)
        all_v = collect(1:Nq)
        for i_s in 1:Ns, i_v in 1:Nv
            V[i_v, i_s] = findfirst(x -> x==1, oneHotV[:, i_v, i_s])
        end
        return V
    end

    function oneHot2Categorical2D(oneHotdata)
        Nf, Ns = size(oneHotdata)
        data = zeros(Int16, Ns)
        all_features = collect(1:Nf)
        for i_s in 1:Ns
            data[i_s] = findfirst(x -> x==1, oneHotdata[:, i_s])
        end
        return data
    end


# DATA STATISTICS ########################################################################################################################################################################

    function remove_autocorrelation(Matij, Nq, Nv)
        for i in 1:(Nv*Nq)
            Matij[i, 1:i] .= 0
        end
        thread_count = nthreads()
        chunk_size = div(Nv, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv)
            for i in start_idx:end_idx, q1 in 1:Nq, q2 in 1:Nq       
                Matij[id(i, q1, Nq), id(i, q2, Nq)] = 0
            end
        end
        return Matij
    end

    # Fi

    function oneHotFreqs(sample, weights, pseudo_count)
        Nq, Nv, Ns = size(sample)
        sample, fi = reshape(sample, (Nq*Nv, Ns)), zeros(Float32, Nv*Nq)
        fi .= sum(weights' .* sample, dims=2) ./ sum(weights)
        fi = (1 - pseudo_count) .* fi .+ pseudo_count / Nq
        return fi, reshape(fi, (Nq, Nv))
    end

    # Fij

    function oneHotFij(sample, weights, pseudo_count)
        Nq, Nv, Ns = size(sample)
        fij = zeros(Float32, Nv*Nq, Nv*Nq)

        sample = reshape(sample, (Nv*Nq, Ns))
        weighted_sample = (weights .* sample')
        
        thread_count = nthreads(); chunk_size = div(Nv*Nq, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv*Nq)
            fij[start_idx:end_idx, :] .= sample[start_idx:end_idx, :] * weighted_sample 
        end
        fij ./= sum(weights)
        fij = remove_autocorrelation(fij, Nq, Nv)
        return (1 - pseudo_count) .* fij .+ pseudo_count / (Nq*Nq)
    end 

    function oneHotFijFast(sample, weights, pseudo_count) 
        Nq, Nv, Ns = size(sample)
        sample_cat, fij = oneHot2Categorical(sample, Nq), zeros(Float32, Nv*Nq, Nv*Nq)     
        for s in 1:Ns, i in 1:Nv, j in i+1:Nv  
            fij[id(i, sample_cat[i, s], Nq), id(j, sample_cat[j, s], Nq)] += weights[s]
        end
        fij ./= sum(weights)
        thread_count = nthreads(); chunk_size = div(Nv, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv)
            for i in start_idx:end_idx, q1 in 1:Nq, q2 in 1:Nq       
                fij[id(i, q1, Nq), id(i, q2, Nq)] = 0
            end
        end
        return (1 - pseudo_count) .* fij .+ pseudo_count / (Nq*Nq)
    end

    function oneHotFijSymm(sample, weights, pseudo_count)
        Nq, Nv, Ns = size(sample)
        fij, sample = zeros(Float32, Nv*Nq, Nv*Nq), reshape(sample, (Nv*Nq, Ns))
        weighted_sample = (weights .* sample')

        thread_count = nthreads(); chunk_size = div(Nv*Nq, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv*Nq)
            fij[start_idx:end_idx, :] .= sample[start_idx:end_idx, :] * weighted_sample 
        end
        fij ./= sum(weights)
        thread_count = nthreads(); chunk_size = div(Nv, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv)
            for i in start_idx:end_idx, q1 in 1:Nq, q2 in 1:Nq       
                fij[id(i, q1, Nq), id(i, q2, Nq)] = 0
            end
        end
        return (1 - pseudo_count) .* fij .+ pseudo_count / (Nq*Nq)
    end

    function oneHotFijSymmFast(sample, weights, pseudo_count=0) 
        Nq, Nv, Ns = size(sample)
        sample_cat, fij = oneHot2Categorical(sample, Nq), zeros(Float32, Nv*Nq, Nv*Nq)     
        for s in 1:Ns, i in 1:Nv, j in i+1:Nv  
            fij[id(i, sample_cat[i, s], Nq), id(j, sample_cat[j, s], Nq)] += weights[s]
            fij[id(j, sample_cat[j, s], Nq), id(i, sample_cat[i, s], Nq)] += weights[s]
        end
        fij ./= sum(weights)

        thread_count = nthreads(); chunk_size = div(Nv, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Nv)
            for i in start_idx:end_idx, q1 in 1:Nq, q2 in 1:Nq       
                fij[id(i, q1, Nq), id(i, q2, Nq)] = 0
            end
        end
        return (1 - pseudo_count) .* fij .+ pseudo_count / (Nq*Nq)
    end

    # Cij

    function oneHotCij(sample, weights, pseudo_count)
        Nq, Nv, Ns = size(sample)
        fi, _ = oneHotFreqs(sample, weights, pseudo_count)
        fij = oneHotFij(sample, weights, pseudo_count)
        cij = zeros(Float32, Nv*Nq, Nv*Nq)
        cij .= fij .- (fi .* transpose(fi))
        cij = remove_autocorrelation(cij, Nq, Nv)
        return cij
    end

    function oneHotCijFast(sample, weights, pseudo_count)
        Nq, Nv, Ns = size(sample)
        fi, _ = oneHotFreqs(sample, weights, pseudo_count)
        fij = oneHotFijFast(sample, weights, pseudo_count)
        cij = zeros(Float32, Nv*Nq, Nv*Nq)
        cij .= fij .- (fi .* transpose(fi))
        cij = remove_autocorrelation(cij, Nq, Nv)
        return cij
    end

# DATA AND MODEL INITIALIZATION ########################################################################################################################################################################

    function read_fasta(data_path, vocab)
        translation = Dict(zip(vocab, 1:length(vocab)))
        fasta_content = read(data_path, String)
        lines = split(fasta_content, '\n')
        lines = [s for s in lines if !startswith(s, '>')]
        data = zeros(Int, length(lines), length(lines[1]))
        for (m, line) in enumerate(lines)
            for (i, letter) in enumerate(line)
                data[m, i] = get(translation, line[i], 0) 
            end
        end
        data_final = copy(data)
        m = 1
        while m <= size(data_final, 1)
            if any(x -> x == 0, data_final[m, :])
                println("Warning: An aminoacid not in the vocabulary is present in the sequence number ", m, " and was removed"); flush(stdout)
                data_final = data_final[1:end .!= m, :]
            else
                m += 1
            end
        end
        return data_final
    end


    function read_fasta2(data_path, vocab)
        translation = Dict(zip(vocab, 1:length(vocab)))
        fasta_content = read(data_path, String)
        lines = split(fasta_content, '\n')
        
        # Unisci le linee che appartengono alla stessa sequenza
        sequences = []
        current_sequence = ""
        for line in lines
            if startswith(line, '>')
                if current_sequence != ""
                    push!(sequences, current_sequence)
                    current_sequence = ""
                end
            else
                current_sequence *= line
            end
        end
        if current_sequence != ""
            push!(sequences, current_sequence)
        end
    
        data = zeros(Int, length(sequences), length(sequences[1]))
        for (m, sequence) in enumerate(sequences)
            for (i, letter) in enumerate(sequence)
                data[m, i] = get(translation, letter, 0)
            end
        end
    
        data_final = copy(data)
        m = 1
        while m <= size(data_final, 1)
            if any(x -> x == 0, data_final[m, :])
                println("Warning: An aminoacid not in the vocabulary is present in the sequence number ", m, " and was removed")
                flush(stdout)
                data_final = data_final[1:end .!= m, :]
            else
                m += 1
            end
        end
        
        return data_final
    end
    

    function modify_fasta_headers(filename, alphabet, v_model, energies)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        fasta = open(filename, "r") do file
            readlines(file)
        end
        modified_fasta, idx, header = [], 1, ""
        for line in fasta
            line_labeled = "$(alphabet[v_cat[:, idx]])"
            if startswith(line, ">")
                header = line 
            elseif line_labeled == line
                push!(modified_fasta, header * " DCAenergy $(energies[idx])\n"); push!(modified_fasta, line*"\n")
                idx += 1
            end
        end
        return modified_fasta
    end

    function compute_weights(sample::BitArray{3}, outputpath, label, threshold=0.8)
        """
        threshold : sequence identity threshold
        """
        Nq, Nv, Ns = size(sample)
        sample_flat = reshape(sample, (Nq*Nv, Ns))
        seqID_matrix = zeros(Float32, Ns, Ns)
        weights = ones(Float32, Ns)
        # compute sequence identity matrix
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Ns)
            seqID_matrix[start_idx:end_idx, :] .= (sample_flat[:, start_idx:end_idx]' * sample_flat) ./ Nv
        end
        # compute weights
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Ns)
            for i in start_idx:end_idx
                weights[i] = sum(1 != seqID_matrix[i, :] .>= threshold)
            end
        end
        weights = weights.^(-1)
        Meff = sum(weights)

        weights_file = (label != nothing) ? outputpath*"/"*label*"_"*"weights.dat" : outputpath*"/weights.dat"
        file = open(weights_file, "w")
        for w in weights
            write(file, "$(w)\n")
        end
        close(file)
        return weights, Meff
    end

    function set_alphabet(alphabet)
        if alphabet == "protein"
            alphabet = "-ACDEFGHIKLMNPQRSTVWY"
        elseif alphabet == "rna"
            alphabet = "-ACGU"
        elseif alphabet == "dna"
            alphabet = "-ACGT"
        else
            alphabet = alphabet
            # error("alphabet can be chosen between: 'dna', 'rna' and 'protein'.")
        end
        return alphabet
    end

    function read_graph(filepath, Nv, Nq)
        file = open(filepath, "r"); content = read(file, String)[1:end-1]; close(file)
        content = split(content, "\n")
        J, vbias = zeros(Float32, Nv*Nq, Nv*Nq), zeros(Float32, Nq, Nv)
        for line ∈ content
            if string(line[1]) == "J"
                _, i, j, a, b, value = split(line, " ")
                i, j, a, b = parse(Int64, i) + 1, parse(Int64, j) + 1, parse(Int64, a) + 1, parse(Int64, b) + 1
                # a = (a != 0) ? a : Nq 
                # b = (b != 0) ? b : Nq 
                id1, id2 = id(i, a, Nq), id(j, b, Nq)
                J[id1, id2] = parse(Float32, value)
                J[id2, id1] = parse(Float32, value)
            elseif string(line[1]) == "h"
                _, i, a, value = split(line, " ")
                i, a = parse(Int64, i) + 1, parse(Int64, a) + 1
                # a = (a != 0) ? a : Nq 
                vbias[a, i] = parse(Float32, value)
            end
        end
        return J, vbias
    end


    function read_graph_new(filepath)
        file = open(filepath, "r"); content = read(file, String)[1:end-1]; close(file)
            content = split(content, "\n")
            alphabet, idx = "", 0
            for line ∈ content
                if string(line[1]) == "h"
                    _, i, a, value = split(line, " ")
                    if alphabet == "" || !occursin(a, string(alphabet)) 
                        alphabet = alphabet * a
                    end
                    idx += 1
                end
            end
            Nq = length(alphabet)
            Nv = div(idx, Nq)
            J, vbias = zeros(Float32, Nv*Nq, Nv*Nq), zeros(Float32, Nq, Nv)
            for line ∈ content
                if string(line[1]) == "J"
                    _, i, j, a, b, value = split(line, " ")
                    i, j = parse(Int64, i) + 1, parse(Int64, j) + 1
                    id1, id2 = id(i, findfirst(a, string(alphabet)), Nq)[1], id(j, findfirst(b, string(alphabet)), Nq)[1]
                    J[id1, id2], J[id2, id1] = parse(Float32, value), parse(Float32, value)
                elseif string(line[1]) == "h"
                    _, i, a, value = split(line, " ")
                    i = parse(Int64, i) + 1
                    vbias[findfirst(a, string(alphabet))[1], i] = parse(Float32, value)
                end
            end
        return J, vbias, alphabet
    end


    function initialize_graph(J, Nq, Nv)
        filter, contact_list, site_degree = (J .!= 0), zeros(Int64, Nv, Nv), zeros(Int64, Nv)
        for i in 1:Nv, j in (i+1):Nv
            if sum(filter[id(i, 1:Nq, Nq), id(j, 1:Nq, Nq)]) != 0 
                site_degree[i] += 1; site_degree[j] += 1
                contact_list[i, site_degree[i]] = j; contact_list[j, site_degree[j]] = i
            end
        end
        return filter, contact_list, site_degree
    end

    function initialize_graph_couplingwise(J, Nq, Nv)
        filter, contact_list, site_degree = (J .!= 0), zeros(Int64, Nv*Nq, Nv*Nq), zeros(Int64, Nv*Nq)
        for iq in 1:Nv*Nq, jq in iq+1:Nv*Nq
            if filter[iq, jq] != 0 
                site_degree[iq] += 1; site_degree[jq] += 1
                contact_list[iq, site_degree[iq]] = jq; contact_list[jq, site_degree[jq]] = iq
            end
        end
        return filter, contact_list, site_degree
    end


# SAMPLING FUNCTIONS ########################################################################################################################################################################

    function sample_from_profile(vbias, Ns, inv_temp=1)
        Nq, Nv = size(vbias)
        v_model_cat = zeros(Nv, Ns)
        for i_v in 1:Nv, i_s in 1:Ns
            v_model_cat[i_v, i_s] = sample(Weights(exp.(inv_temp * vbias[:, i_v])))
        end
        return oneHotEncoding(v_model_cat, Nq)
    end

    function gibbs_sampling_edgewise(J::Matrix{Float32}, vbias::Matrix{Float32}, V::BitArray{3}, contact_list, site_degree, nsweeps, inv_temp=1)
        Nq, Nv, Ns = size(V)
        V, energy0, all_v = reshape(V, Nq*Nv, Ns), zeros(Float32, Nq, Ns), collect(1:Nq)
        for sweep in 1:nsweeps, i_v in randperm(Nv)
            indices = id(i_v, 1:Nq, Nq) # Nq indices of flat encoding corresponding to position i_v
            energy = energy0 .+ vbias[:, i_v]
            for i_s in 1:Ns # loop over the sample size  
                if site_degree[i_v] != 0
                    for contact in contact_list[i_v, 1:site_degree[i_v]]
                        contacts_jq = id(contact, 1:Nq, Nq)
                        energy[:, i_s] += sum(J[indices, contacts_jq] * V[contacts_jq, i_s], dims=2) 
                    end
                end
                V[indices, i_s] .= (all_v .== sample(Weights(exp.(inv_temp * energy[:, i_s]))))
            end      
        end
        return reshape(V, (Nq, Nv, Ns))
    end

    function gibbs_sampling_couplingwise(J::Matrix{Float32}, vbias::Matrix{Float32}, V::BitArray{3}, contact_list, site_degree, nsweeps, inv_temp=1)
        Nq, Nv, Ns = size(V)
        V, energy, all_v = reshape(V, Nq*Nv, Ns), zeros(Float32, Nq, Ns), collect(1:Nq)
        indices = zeros(Int32, Nq)
        for sweep in 1:nsweeps, i_v in randperm(Nv)
            indices .= id(i_v, 1:Nq, Nq) # Nq indices of flat encoding corresponding to position i_v
            energy .=  vbias[:,i_v]
            for i_s in 1:Ns # loop over the sample size                
                for (iq, idx) in zip(indices, 1:Nq) # loop over the iq corresponding to i_v
                    energy[idx, i_s] += (site_degree[iq] != 0) ? sum(J[iq, contact_list[iq, 1:site_degree[iq]]] .* V[contact_list[iq, 1:site_degree[iq]], i_s]) : 0
                end
                V[indices, i_s] .= (all_v .== sample(Weights(exp.(inv_temp * energy[:, i_s]))))
            end      
        end
        return reshape(V, (Nq, Nv, Ns))
    end

    function  gibbs_sampling(J::Matrix{Float32}, vbias::Matrix{Float32}, V::BitArray{3}, nsweeps, inv_temp=1)
        Nq, Nv, Ns = size(V)
        V, energy, all_v = reshape(V, Nq*Nv, Ns), zeros(Float32, Nq, Ns), collect(1:Nq)
        w = similar(energy)
        indices = zeros(Int32, Nq)
        for sweep in 1:nsweeps
            for i_v in randperm(Nv)
                indices .= id(i_v, 1:Nq, Nq) # Nq indices of flat encoding corresponding to position i_v
                energy .= vbias[:, i_v] .+ J[indices, :] * V  
                w = exp.(inv_temp * energy)
                
                for i_s in 1:Ns
                    V[indices, i_s] .= (all_v .== sample(Weights(w[:, i_s]))) # exp.(inv_temp * energy[:, i_s])
                end
            end
        end
        return reshape(V, (Nq, Nv, Ns))
    end





    function  metropolis_sampling(J::Matrix{Float32}, vbias::Matrix{Float32}, V::BitArray{3}, nsweeps, inv_temp=1)
        Nq, Nv, Ns = size(V)
        V, deltaE, all_v, old_res = reshape(V, Nq*Nv, Ns), zeros(Float32, Ns), collect(1:Nq), zeros(Int, Ns)
        for sweep in 1:nsweeps, i_v in randperm(Nv)
            indices = id(i_v, 1:Nq, Nq)
            for i_s in 1:Ns
                old_res[i_s] = findfirst(==(1), V[id(i_v, 1:Nq, Nq), i_s])
            end
            propositions = rand(1:Nq, Ns)
            deltaE = vbias[propositions, i_v] .- vbias[old_res, i_v]  .+ sum((J[id(i_v, propositions, Nq), :] .- J[id(i_v, old_res, Nq), :] ) .* V', dims=2)  #* V[:, i_s]
            for i_s in 1:Ns
                (rand() <= exp(inv_temp * deltaE[i_s])) ?  V[indices, i_s] .= (all_v .== propositions[i_s]) : nothing
            end
        end
        return reshape(V, (Nq, Nv, Ns))
    end
    
    function metropolis_sampling_couplingwise(J::Matrix{Float32}, vbias::Matrix{Float32}, V::BitArray{3}, contact_list, site_degree, nsweeps, inv_temp=1)
        Nq, Nv, Ns = size(V)
        V, deltaE, all_v, old_res = reshape(V, Nq*Nv, Ns), 0, collect(1:Nq), zeros(Int, Ns)
        for sweep in 1:nsweeps, i_v in randperm(Nv)
            indices = id(i_v, 1:Nq, Nq) # Nq indices of flat encoding corresponding to position i_v
            for i_s in 1:Ns
                old_res[i_s] = findfirst(==(1), V[id(i_v, 1:Nq, Nq), i_s])
            end
            propositions = rand(1:Nq, Ns)
            for i_s in 1:Ns # loop over the sample size    
                deltaE = vbias[propositions[i_s], i_v] .- vbias[old_res[i_s], i_v]             
                #for (iq, idx) in zip(indices, 1:Nq) # loop over the iq corresponding to i_v
                    iq_old, iq_prop = id(i_v, old_res[i_s], Nq), id(i_v, propositions[i_s], Nq)
                    deltaE += (site_degree[iq_prop]!= 0) ? sum(J[iq_prop, contact_list[iq_prop, 1:site_degree[iq_prop]]] .* V[contact_list[iq_prop, 1:site_degree[iq_prop]], i_s]) : 0
                    deltaE += (site_degree[iq_old]!= 0) ? - sum(J[iq_old, contact_list[iq_old, 1:site_degree[iq_old]]] .* V[contact_list[iq_old, 1:site_degree[iq_old]], i_s]) : 0
                #end
                (rand() <= exp(inv_temp * deltaE)) ?  V[indices, i_s] .= (all_v .== propositions[i_s]) : nothing
            end      
        end
        return reshape(V, (Nq, Nv, Ns))
    end


# COMPUTE ENERGY ########################################################################################################################################################################
   
    function compute_energy(J::Matrix{Float32}, vbias::Matrix{Float32}, sequences::Union{BitArray{3}, BitMatrix})
        Ns = size(sequences, 3); Nv = size(sequences, 2); Nq = size(sequences, 1)
        sequences, vbias, energy = reshape(sequences, Nq*Nv, Ns), reshape(vbias, Nq*Nv), zeros(Float32, Ns)
        J_perm = permutedims(J, [2, 1])
        thread_count = nthreads(); chunk_size = div(Ns, thread_count)
        @threads :static for t in 1:thread_count
            start_idx, end_idx, _ = index_interval(t, thread_count, chunk_size, Ns)
            for idx in start_idx:end_idx    
                J_ene = sum((J * sequences[:, idx]) .* sequences[:, idx])
                energy[idx] = - sum(vbias .* sequences[:, idx]) - 0.5 * J_ene
            end
        end
        return energy
    end

    function compute_energy_1thread(J::Matrix{Float32}, vbias::Matrix{Float32}, sequences::Union{BitArray{3}, BitMatrix})
        Ns = size(sequences, 3); Nv = size(sequences, 2); Nq = size(sequences, 1)
        sequences, vbias, energy = reshape(sequences, Nq*Nv, Ns), reshape(vbias, Nq*Nv), zeros(Float32, Ns)
        for idx in 1:Ns
            energy[idx] = - sum(vbias .* sequences[:, idx])  .+ (1/2) * sum((J * sequences[:, idx]) .* sequences[:, idx])
        end
        
        return energy
    end

    function compute_energy_from_fasta(path_params, path_chains, outputpath) 
        # logfile = open(model_dir * "/adabmDCA.log", "w"); redirect_stdout(logfile)
        inv_temp = 1
        J, vbias, v_model, alphabet = restore_model_new(path_params, path_chains)
        energies = compute_energy(J, vbias, v_model)
        modified_fasta = modify_fasta_headers(path_chains, alphabet, v_model, energies)
        open(outputpath, "w") do file
            for line in modified_fasta
                write(file, line)
            end
        end
        # close(logfile)
    end

# DMS ENERGY ########################################################################################################################################################################

    
    function compute_DMS_energies(path_params, wt_path, outputpath)
        J, vbias, v_model, alphabet = restore_model_new(path_params, wt_path)
        Nq, Nv = size(vbias)
        WT_seq = v_model[:, :, 1]
        DMS, headers = build_DMS(WT_seq, alphabet)
        energies = compute_energy(J, vbias, oneHotEncoding(DMS, Nq))
        DCA_score = energies .- energies[1]
        save_DMS_chains(DMS, headers, DCA_score, alphabet, outputpath)
    end


    function save_DMS_chains(DMS, headers, energies, alphabet, outputpath)
        Nv, Ns = size(DMS)
        file_chains = open(outputpath, "w")
            for m in 1:Ns-1
                header = headers[m] * " DCAscore " * string(energies[m]) * "\n"
                line = "$(alphabet[DMS[:, m]])\n"
                write(file_chains, header); write(file_chains, line)
            end
            head = headers[Ns] * " DCAscore " * string(energies[Ns]) * "\n"; line = "$(alphabet[DMS[:, Ns]])"
            write(file_chains, head); write(file_chains, line)
        close(file_chains)
    end

    function build_DMS(WT_seq, alphabet)
        (Nq, Nv) = size(WT_seq)
        WT_seq_cat, all_v = zeros(Int16, Nv), collect(1:Nq)
        for i_v in 1:Nv
            WT_seq_cat[i_v] = findfirst(x -> x==1, WT_seq[:, i_v])
        end
        DMS, headers = zeros(Int32, Nv, (Nq-1)*Nv + 1), fill("", (Nq-1)*Nv + 1)
        DMS[:, 1], headers[1] = WT_seq_cat, ">wt sequence"
        idx = 2
        for v in 1:Nv, q in 1:Nq
            if q != WT_seq_cat[v] && string(alphabet[q]) != "-"
                headers[idx] = ">" * string(alphabet[WT_seq_cat[v]]) * string(v) * string(alphabet[q])
                mut_seq = copy(WT_seq_cat); mut_seq[v] = q
                DMS[:, idx] = mut_seq
                idx += 1
            end
        end
        return DMS[:, 1:idx-1], headers[1:idx-1]
    end




# SAVE AND RESTORE MODEL ########################################################################################################################################################################

    function restore_model(params_path, path_chains, alphabet)
        data = read_fasta(path_chains, alphabet)
        v_model = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
        Nq, Nv, _ = size(v_model)
        J, vbias = read_graph(params_path, Nv, Nq)
        return Float32.(J), Float32.(vbias), v_model
    end

    function restore_model_new(params_path, path_chains)
        J, vbias, alphabet = read_graph_new(params_path)
        data = read_fasta2(path_chains, alphabet)
        v_model = oneHotEncoding(permutedims(data, [2, 1]), length(alphabet))
        Nq, Nv, _ = size(v_model)
        (Nq != size(vbias, 1) || Nv != size(vbias, 2)) ? error("model and data do not correspond!") : nothing  
        return Float32.(J), Float32.(vbias), v_model, alphabet
    end

    function save_model_new(J, vbias, filter, alphabet, model_path)
        Nq, Nv = size(vbias)
        file_model = open(model_path, "w")
        for i in 1:Nv, j in i+1:Nv
            for iq in 1:Nq, jq in 1:Nq
                if filter[id(i, iq, Nq), id(j, jq, Nq)] != 0
                    line = "J $(i-1) $(j-1) $(alphabet[iq]) $(alphabet[jq]) $(J[id(i, iq, Nq), id(j, jq, Nq)])\n"
                    write(file_model, line)
                end
            end
        end
        for i in 1:Nv, iq in 1:Nq
            line = "h $(i-1) $(alphabet[iq]) $(vbias[iq, i])\n"
            write(file_model, line)
        end
        close(file_model)
    end

    function save_chains_new(J, vbias, v_model, alphabet, chains_path)
        Nq, Nv, Ns = size(v_model)
        v_cat = oneHot2Categorical(v_model, Nq)
        file_chains = open(chains_path, "w")
            for m in 1:Ns-1
                head = ">chain $m\n"
                line = "$(alphabet[v_cat[:, m]])\n"
                write(file_chains, head); write(file_chains, line)
            end
            head = ">chain $Ns\n"; line = "$(alphabet[v_cat[:, Ns]])"
            write(file_chains, head); write(file_chains, line)
        close(file_chains)
    end


    function assign_natural_weights(weights, v_natural, pseudo_count, outputpath, label)
        (Nq, Nv, Ns) = size(v_natural)
        if weights == nothing
            println("computing weights..."); flush(stdout)
            (natural_weights, Meff) = compute_weights(v_natural, outputpath, label)
        else
            natural_weights = parse.(Float32, readlines(weights))
            Meff = sum(natural_weights)
            (length(natural_weights) !== Ns) ? error("number of weights in the file != number of sequences in the data") : error
        end
        println("effective number of sequences in the dataset: ", Meff)
        (pseudo_count == nothing) ? pseudo_count = 1/Meff : nothing
        return natural_weights, Meff, pseudo_count
    end



    # OTHERS 

    function old2new_model(old_pathparams, new_pathparams, vocab, Nv)
        vocab = set_alphabet(vocab)
        Nq = length(vocab)
        J, vbias = read_graph(old_pathparams, Nv, Nq)
        save_model_new(J, vbias, (J.!==0), vocab, new_pathparams)
        return 
    end

    # function distance_from_MSA(seq, V, Nv)
    #     Ns = size(V, 1)
    #     length(seq) == Ns || throw(ArgumentError("Samples must have same size"))
    #     distances = zeros(Int32, Ns)
    #     M = 0
    #     for i in 1:Ns
    #         distances[i] = Nv - count(seq .!= V[:, i]) / 2 #count(seq .== V[:, i]) 
    #         M = (distances[i] >= M) ? distances[i] : M
    #     end
    #     return M, distances
    # end

    function sample_sid_from_MSA(sample, MSA, Nv)
        Ns, L_msa = size(sample, 2), size(MSA, 2)
        sids = zeros(Int32, Ns)
        for i_s in 1:Ns, i_msa in 1:L_msa
            m = Nv - count(sample[:, i_s] .!= MSA[:, i_msa]) / 2 
            sids[i_s] = (m >= sids[i_s]) ? m : sids[i_s]
        end
        return sids
    end

    function oneHotHammingDistance(seq1, seq2)
        length(seq1) == length(seq2) || throw(ArgumentError("Sequences must have same length"))
        return count(seq1 .!= seq2) / 2
    end






#     export id, index_interval
#     export linear_division, exponential_division_euler, compute_slope
#     export oneHotEncoding, oneHot2Categorical, oneHot2Categorical2D
#     export remove_autocorrelation, oneHotFreqs, oneHotFij, oneHotFijFast, oneHotFijSymm, oneHotFijSymmFast, oneHotCij, oneHotCijFast
#     export read_fasta, read_fasta2, modify_fasta_headers, compute_weights, set_alphabet, read_graph, read_graph_new, initialize_graph, initialize_graph_couplingwise
#     export sample_from_profile, gibbs_sampling_edgewise, gibbs_sampling_couplingwise, gibbs_sampling, metropolis_sampling, metropolis_sampling_couplingwise
#     export restore_model, restore_model_new, save_model_new, save_chains_new, old2new_model
#     export compute_energy, compute_energy_1thread, compute_energy_from_fasta
#     export assign_natural_weights
#     export wt2dms, distance_from_MSA, oneHotHammingDistance, sample_sid_from_MSA
# end