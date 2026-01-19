# This code is a specialized version (with skew symmetric foot) of the Pade approach in computing Dexp proposed in computing the Frechet derivative of the matrix exponential ... by Al-Mohy, Higham in https://eprints.maths.manchester.ac.uk/1218/1/alhi09.pdf


#!/usr/bin/env julia
import Pkg

# list only non-stdlib deps
const DEPS = ["Printf", "DelimitedFiles", "Random"]

function ensure_deps()
    # use a local env next to the script (better than polluting global env)
    script_dir = @__DIR__
    Pkg.activate(script_dir; shared=false)

    # install if needed
    for pkg in DEPS
        if Base.find_package(pkg) === nothing
            @info "Installing $pkg ..."
            Pkg.add(pkg)
        end
    end

    Pkg.instantiate()
    Pkg.precompile()
end

ensure_deps()

using LinearAlgebra, Printf, DelimitedFiles, Random

## Workspace Class ##

import Base: getindex, setindex!

mutable struct WSP
    # wsp: the type that stores references of workspace 
    vec::Vector{Any}
    function WSP(c...)
        n = length(c)
        if n == 1 && typeof(c[1]) <: Int
            vec = Vector{Any}(undef, c)
        else
            vec = Vector{Any}(undef, n)
            for ii = 1:n
                vec[ii] = Ref(c[ii])
            end
        end
        return new(vec)
    end
end

(workspace::WSP)(i::Int) = @inbounds workspace.vec[i];
getindex(workspace::WSP, i) = @inbounds workspace.vec[i][];
setindex!(workspace::WSP, val, i) = workspace.vec[i][] = val;

## Workspace Class ##


retrieve(workspace::WSP, i) = workspace.vec[i][]
safe_access(workspace::WSP, i) = workspace.vec[i][]


const REAL_TYPE = BigFloat

# No rescale is included

@inline get_wsp_expm_Pade(n::Int) = WSP(Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n),
    Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n),
    Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n), Matrix{REAL_TYPE}(undef, n, n))


mutable struct expm_Pade_system
    A::Ref{Matrix{REAL_TYPE}}
    A2::Ref{Matrix{REAL_TYPE}}
    A4::Ref{Matrix{REAL_TYPE}}
    A6::Ref{Matrix{REAL_TYPE}}
    W::Ref{Matrix{REAL_TYPE}}
    W1::Ref{Matrix{REAL_TYPE}}
    W2::Ref{Matrix{REAL_TYPE}}
    Z1::Ref{Matrix{REAL_TYPE}}
    Z2::Ref{Matrix{REAL_TYPE}}
    U::Ref{Matrix{REAL_TYPE}}
    V::Ref{Matrix{REAL_TYPE}}
    b::Ref{Vector{Int}}
    # lu::LU{REAL_TYPE,Matrix{REAL_TYPE},Vector{Int64}}
    scale::Int

    expm_Pade_system(n::Int) = new(Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)),
        Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(REAL_TYPE, n, n)), Ref(zeros(Int64, 14)), 0)
end


function Expm_Pade(S_sys::expm_Pade_system, R::Ref{Matrix{REAL_TYPE}}, S::Ref{Matrix{REAL_TYPE}}, M::Ref{Matrix{REAL_TYPE}})

    MatR = R[]
    MatS = S[]
    MatM = M[]

    MatA = S_sys.A[]
    MatA2 = S_sys.A2[]
    MatA4 = S_sys.A4[]
    MatA6 = S_sys.A6[]
    MatW = S_sys.W[]
    MatW1 = S_sys.W1[]
    MatW2 = S_sys.W2[]
    MatZ1 = S_sys.Z1[]
    MatZ2 = S_sys.Z2[]
    MatU = S_sys.U[]
    MatV = S_sys.V[]
    VecB = S_sys.b[]

    lm_13::REAL_TYPE = 4.74
    VecB .= [64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600,
        670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1]

    S_sys.scale = max(0, ceil(log2(opnorm(MatS, 1) / lm_13)))

    copy!(MatA, MatS)
    ldiv!(2^S_sys.scale, MatA)
    # Scale S to get A = S/2^s

    mul!(MatA2, MatA, MatA)
    mul!(MatA4, MatA2, MatA2)
    mul!(MatA6, MatA2, MatA4)

    fill!(MatW1, 0.0)
    axpy!(VecB[14], MatA6, MatW1)
    axpy!(VecB[12], MatA4, MatW1)
    axpy!(VecB[10], MatA2, MatW1)

    fill!(MatW2, 0.0)
    axpy!(VecB[8], MatA6, MatW2)
    axpy!(VecB[6], MatA4, MatW2)
    axpy!(VecB[4], MatA2, MatW2)
    @inbounds for ind in axes(MatW2, 1)
        MatW2[ind, ind] += VecB[2]
    end

    fill!(MatZ1, 0.0)
    axpy!(VecB[13], MatA6, MatZ1)
    axpy!(VecB[11], MatA4, MatZ1)
    axpy!(VecB[9], MatA2, MatZ1)

    fill!(MatZ2, 0.0)
    axpy!(VecB[7], MatA6, MatZ2)
    axpy!(VecB[5], MatA4, MatZ2)
    axpy!(VecB[3], MatA2, MatZ2)
    @inbounds for ind in axes(MatZ2, 1)
        MatZ2[ind, ind] += VecB[1]
    end

    copy!(MatW, MatW2)
    mul!(MatW, MatA6, MatW1, 1, 1)

    mul!(MatU, MatA, MatW)

    copy!(MatV, MatZ2)
    mul!(MatV, MatA6, MatZ1, 1, 1)

    copy!(MatR, MatV)
    axpy!(-1, MatU, MatR)

    copy!(MatM, MatU)
    axpy!(1, MatV, MatM)
    MatR .= MatR \ MatM

    for ind = 1:S_sys.scale
        copy!(MatM, MatR)
        mul!(MatR, MatM, MatM)
    end
end

function _expm_data_mesh(dim_beg::Int, dim_end::Int)
    mesh = Int[]
    
    # From 3 up to 100, step 1
    append!(mesh, 3:1:min(100, dim_end))
    dim_end <= 100 && return mesh[mesh .>= dim_beg]
    
    # From 100 to 400, step 10
    append!(mesh, 100:10:min(400, dim_end))
    dim_end <= 400 && return mesh[mesh .>= dim_beg]
    
    # From 400 to 600, step 20
    append!(mesh, 400:20:min(600, dim_end))
    dim_end <= 600 && return mesh[mesh .>= dim_beg]
    
    # From 600 to 1000, step 50
    append!(mesh, 600:50:min(1000, dim_end))
    dim_end <= 1000 && return mesh[mesh .>= dim_beg]
    
    # Above 1000, step 100
    append!(mesh, 1000:100:dim_end)

    return mesh[mesh .>= dim_beg]
end

function expm_accurate_data(dim=(3, 100); seed = 9527, ioformat = "w")
    dim_beg, dim_end = dim;
    rand_eng = MersenneTwister(seed)
    @printf "Creating accurate exp(A) for skew A. Random seed: %i.\n\n" seed


    file_skewsymm = open("expm_accurate_skewsymm_s$(seed).out", ioformat)
    file_specorth = open("expm_accurate_specorth_s$(seed).out", ioformat)
    file_xforward = open("expm_accurate_xforward_s$(seed).out", ioformat)
    file_yinverse = open("expm_accurate_yinverse_s$(seed).out", ioformat)

    try
        for n in _expm_data_mesh(dim_beg,dim_end)
            @printf "\rComputing %i x %i case" n n

            varep::REAL_TYPE = 1e-20

            l = div(n * (n - 1), 2)
            Vec = 2 .* rand(rand_eng, l) .- 1

            MatA1 = zeros(Float64, n, n)
            MatR1 = zeros(Float64, n, n)
            MatX1 = zeros(Float64, n, n)
            MatY1 = zeros(Float64, n, n)

            vec_ind = 1
            for ind in 1:(n-1)
                pos = (ind - 1) * n + ind + 1
                len = n - ind
                MatA1[pos:(pos+len-1)] .= Vec[vec_ind:(vec_ind+len-1)]
                vec_ind += len
            end
            MatA1 .-= MatA1'

            Vec .= 2 .* rand(rand_eng, l) .- 1

            vec_ind = 1
            for ind in 1:(n-1)
                pos = (ind - 1) * n + ind + 1
                len = n - ind
                MatX1[pos:(pos+len-1)] .= Vec[vec_ind:(vec_ind+len-1)]
                vec_ind += len
            end
            MatX1 .-= MatX1'

            MatA2 = zeros(REAL_TYPE, n, n)
            MatR2 = zeros(REAL_TYPE, n, n)
            MatM2 = zeros(REAL_TYPE, n, n)
            MatX2 = zeros(REAL_TYPE, n, n)
            MatY2 = zeros(REAL_TYPE, n, n)

            MatA2 .= MatA1
            MatX2 .= MatX1
            sys = expm_Pade_system(n)
            Expm_Pade(sys, Ref(MatR2), Ref(MatA2), Ref(MatM2));
            # display(MatR2 .- exp(MatA1))

            MatY2 .= MatR2
            MatR1 .= MatR2

            MatA2 .= MatA1
            MatA2 .+= varep .* MatX2
            Expm_Pade(sys, Ref(MatR2), Ref(MatA2), Ref(MatM2));

            MatY2 .= (MatR2 .- MatY2) ./ varep
            MatY1 .= MatY2

            # append!(VecA, n)
            # append!(VecA, vec(MatA1))
            # append!(VecQ, n)
            # append!(VecQ, vec(MatR1))

            write(file_skewsymm, Float64(n))
            write(file_skewsymm, vec(MatA1))

            write(file_specorth, Float64(n))
            write(file_specorth, vec(MatR1))

            write(file_xforward, Float64(n))
            write(file_xforward, vec(MatX1))

            write(file_yinverse, Float64(n))
            write(file_yinverse, vec(MatY1))
        end
        @printf "\r\nExit the routine.\n\n"
    finally

        close(file_skewsymm)
        close(file_specorth)
        close(file_xforward)
        close(file_yinverse)
    end
end