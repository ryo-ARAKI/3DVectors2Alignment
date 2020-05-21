#=
Julia program to compute PDF of alignment of vectors in  distributed in 3D space uniform-random-angle
=#

"""
Module for parameters & variables
"""
module ParamVar
    """
    Global parameters
    """
    struct Parameters
        x_lim::Float64  # Boundary of cubic region [-x_lim, x_lim]³
        num_vectors::Int64  # Number of vectors
    end

    """
    Property of vector
    """
    mutable struct Vector
        x::Float64  # x coordinate of origin
        y::Float64
        z::Float64
        vx::Float64  # Component in x direction
        vy::Float64
        vz::Float64
    end

    """
    Alignment in terms of 3D polar coordinates
    """
    mutable struct Alignment_Polar
        cos_θ::Float64
        ϕ::Float64
        cos_ϕ::Float64
    end
end


"""
Module for defining origin & components of vectors
"""
module DefineVectors
    """
    Define 3d random points (origin of vectors)
    """
    function distribute_points(param, vectors)

    end

    """
    Define vectors in uniform random direction
    """
    function distribute_angles(param, vectors)

    end
end


"""
Module for computing alignments
"""
module ComputeAlignments
    """
    Compute PDF of vector alignment in 3D polar coordinates
    """
    function compute_alignment_polar(param, vectors, aligns)

    end
end



# ========================================
# Main function
# ========================================

## Declare modules
using Printf
using Plots
gr(
    markerstrokewidth = 0,
    markersize = 10
)

using .ParamVar
using .DefineVectors:
    distribute_points,
    distribute_angles
using .ComputeAlignments:
    compute_alignment_polar


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
x_lim = 1.0
num_vectors = 100
param = ParamVar.Parameters(
    x_lim, num_vectors
)

# Define mutable struct for vectors
vectors = Array{ParamVar.Vector}(undef, param.num_vectors)


# ----------------------------------------
## Define 3d random points (origin of vectors)
# ----------------------------------------
distribute_points(param, vectors)


# ----------------------------------------
## Define vectors in uniform random direction
# ----------------------------------------
distribute_angles(param, vectors)


# Define mutable structs for alignments
aligns_polar = Array{ParamVar.Alignment_Polar}(undef, param.num_vectors)

# ----------------------------------------
## Compute PDF of vector alignment in 3D polar coordinates
# ----------------------------------------
compute_alignment_polar(param, vectors, aligns_polar)

#=
# ----------------------------------------
## Compute PDF of vector alignment in 3D cylindrical coordinates
# ----------------------------------------
compute_alignment_polar(param, vectors, align_cylindrical)
=#
