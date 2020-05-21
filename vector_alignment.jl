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
        # Constructor
        Vector() = new()
    end

    """
    Alignment in terms of 3D polar coordinates
    """
    mutable struct Alignment_Polar
        cos_θ::Float64
        ϕ::Float64
        cos_ϕ::Float64
        # Constructor
        Alignment_Polar() = new()
    end
end


"""
Module for defining origin & components of vectors
"""
module DefineVectors
    using Distributions

    """
    Define 3d random points (origin of vectors)
    """
    function distribute_points(param, vectors)
        for itr_vec in 1:param.num_vectors
            vectors[itr_vec].x = rand(Uniform(-param.x_lim, param.x_lim))
            vectors[itr_vec].y = rand(Uniform(-param.x_lim, param.x_lim))
            vectors[itr_vec].z = rand(Uniform(-param.x_lim, param.x_lim))
        end

        #=
        x = getfield.(vectors, :x)
        println(x)
        =#

    end

    """
    Define vectors in uniform random direction
    """
    function distribute_angles(param, vectors)
        θ = rand(Uniform(0.0, π), param.num_vectors)  # 0≤θ≤π
        ϕ = rand(Uniform(0.0, 2.0*π), param.num_vectors)  # 0≤ϕ≤2π

        for itr_vec in 1:param.num_vectors
            vectors[itr_vec].vx = sin(θ[itr_vec]) * cos(ϕ[itr_vec])
            vectors[itr_vec].vy = sin(θ[itr_vec]) * sin(ϕ[itr_vec])
            vectors[itr_vec].vz = cos(θ[itr_vec])
        end

        #=
        println(θ)
        vx = getfield.(vectors, :vx)
        println(vx)
        =#

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


"""
Module for plotting
"""
module PlotFigures
    using PyPlot
    # Font setting
    fs=12
    rc("font",family ="Times New Roman",size=fs)
    rc("font",serif ="Times New Roman",size=fs)
    rc("text",usetex ="true")

    """
    Plot 3D vectors
    """
    function plot_3d_vectors(param, vectors)
        # Figure setting
        pygui(true)
        fig = figure()
        ax = gca(
            projection="3d",
            xlabel=L"$x$", ylabel=L"$y$", zlabel=L"$z$",
            xlim=[-1.0, 1.0], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0],
            xticks=[-1.0, -0.5, 0.0, 0.5, 1.0], yticks=[-1.0, -0.5, 0.0, 0.5, 1.0], zticks=[-1.0, -0.5, 0.0, 0.5, 1.0]
        )

        # Obtain fields from array of mutable struct
        x = getfield.(vectors, :x)
        y = getfield.(vectors, :y)
        z = getfield.(vectors, :z)
        vx = getfield.(vectors, :vx)
        vy = getfield.(vectors, :vy)
        vz = getfield.(vectors, :vz)

        # Vector plot
        quiver_coef = 0.2
        ax.quiver(x,y,z, quiver_coef*vx, quiver_coef*vy, quiver_coef*vz)

        savefig("./tmp/vectors.png", bbox_inches="tight", pad_inches=0.1)
    end
end



# ========================================
# Main function
# ========================================

## Declare modules
using Printf
using Plots
# gr(
#     markerstrokewidth = 0,
#     markersize = 10
# )

using .ParamVar
using .DefineVectors:
    distribute_points,
    distribute_angles
using .ComputeAlignments:
    compute_alignment_polar
using .PlotFigures:
    plot_3d_vectors


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
x_lim = 1.0
num_vectors = 10
param = ParamVar.Parameters(
    x_lim, num_vectors
)

# Define mutable struct for vectors
vectors = Array{ParamVar.Vector}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    vectors[itr_vec] = ParamVar.Vector()
end


# ----------------------------------------
## Define 3d random points (origin of vectors)
# ----------------------------------------
distribute_points(param, vectors)


# ----------------------------------------
## Define vectors in uniform random direction
# ----------------------------------------
distribute_angles(param, vectors)

plot_3d_vectors(param, vectors)

# Define mutable structs for alignments
aligns_polar = Array{ParamVar.Alignment_Polar}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    aligns_polar[itr_vec] = ParamVar.Alignment_Polar()
end

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
