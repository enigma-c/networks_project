using LinearAlgebra
using SparseArrays
using ControlSystems
using LightGraphs
using StatsBase
using InvertedIndices

function regularize(A)
    # Not Tested

    # Prune dead nodes
    # i.e those node that have nod edges to and from.them.
    to_delete = []
    for i in size(A,1)
        if sum(A[i,:]) == 0 == sum(A[:,i])
            push!(to_delete, i)
        end
    end

    A = A[ Not(to_delete) , Not(to_delete) ]

    # If Σ_{j != i} w_ij = w_kj, set w_kj to 0
    for i in size(A,1)
        for j in size(A,1)
            total = sum(A[i,:])
            if total == A[j,i]
                A[i,j] = 1
                break
            end
        end
    end

    # A_ij -> A_ij / Σ_i A_ij
    rowSums = sum(A, dims=2) # dims=1 sums across columns
    for (i, x) in enumerate(rowSums)
        if x == 0 continue end
        A[i,:] = A[i,:] ./ x
    end

    return A
end

# Make a canonical basis vector
e = function(i, n)
    out = zeros(n)
    out[i] = 1
    return out
end

function trMi(A, node_index)
    n = size(A,1)
    Di = e(node_index, n) * e(node_index, n)'
    acc = 0
    C = A
    for t=1:n
        if t != 1
            C *= A
        end
        acc += tr( C' * Di * C )
    end
    acc
end

function get_node_metrics(A, node_index)
    B = reshape(e(node_index, size(A,1)), (size(A,1), 1))
    system = ss(A,B, zeros(size(A)), zeros(size(B)), size(A,1)) # Some digit, don't think it applies here
    Ci = ctrb(A,B)
    Wi = Ci * Ci'
    tW = node_to_network = tr(Wi)
    tM = trMi(A, node_index)
    return ( tW = tW, tM = tM )
end

function get_metrics(A)
    Ws = Vector{Float64}(undef, size(A,1))
    Ms = Vector{Float64}(undef, size(A,1))
    for i=1:(size(A,1))
        tW, tM = get_node_metrics(A, i)
        Ws[i] = tW
        Ms[i] = tM
    end
    return Ws, Ms
end

n=8
k=22
A = erdos_renyi(n, k) |> adjacency_matrix
A = convert(Matrix{Float64}, A)
A = regularize(A)
# A = A[ sum(A, dims=1) .= 0 , :]
B = begin
    x = zeros(size(A,1))
    x[2] = 1
    reshape(x, (size(A,1),1) )
    end


