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
        cosθ::Float64
        ϕ::Float64
        cosϕ::Float64
        # Constructor
        Alignment_Polar() = new()
    end

    """
    Alignment in terms of 3D cylindrical coordinates
    """
    mutable struct Alignment_Cylindrical
        cosθ::Float64
        ϕ::Float64
        cosϕ::Float64
        # Constructor
        Alignment_Cylindrical() = new()
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
    Compute vector alignment in 3D polar coordinates
    """
    function compute_alignment_polar(param, vectors, aligns)

        for itr_vec in 1:param.num_vectors
            vx = vectors[itr_vec].vx
            vy = vectors[itr_vec].vy
            vz = vectors[itr_vec].vz

            # θ: angle between vector & z axis
            # ϕ: angle between vector projected onto xy plane & x axis
            cosθ = vz / sqrt(vx^2 + vy^2 + vz^2)
            cosϕ = vx / sqrt(vx^2 + vy^2)
            ϕ = acos(cosϕ)

            aligns[itr_vec].cosθ = cosθ
            aligns[itr_vec].cosϕ = cosϕ
            aligns[itr_vec].ϕ = ϕ
        end

        #=
        cosθ = getfield.(aligns, :cosθ)
        println(cosθ)
        =#

    end

    """
    Compute vector alignment in 3D cylindrical coordinates
    """
    function compute_alignment_cylindrical(param, vectors, aligns)

        for itr_vec in 1:param.num_vectors
            x = vectors[itr_vec].x
            y = vectors[itr_vec].y
            vx = vectors[itr_vec].vx
            vy = vectors[itr_vec].vy
            vz = vectors[itr_vec].vz

            # θ: angle between vector & z axis
            cosθ = vz / sqrt(vx^2 + vy^2 + vz^2)

            # ψ: angle between vector origin on xy plane & coordinate origin
            cosψ = x / sqrt(x^2 + y^2)
            sinψ = y / sqrt(x^2 + y^2)

            # ϕ: angle between vector projected onto xy plane & ψ vector
            vϕ = -sinψ * vx + cosψ * vy
            cosϕ = vϕ / sqrt(vx^2 + vy^2)
            ϕ = acos(cosϕ)

            aligns[itr_vec].cosθ = cosθ
            aligns[itr_vec].cosϕ = cosϕ
            aligns[itr_vec].ϕ = ϕ
        end

    end
end


"""
Module for plotting
"""
module PlotFigures
    using PyPlot
    using PyCall
    sns = pyimport("seaborn")
    sns.set(
        context="talk",
        style="white",
        palette="plasma",
        font="sans-serif",
        font_scale=1,
        color_codes=false,
    )
    rc("text", usetex ="true")

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
            xlim=[-param.x_lim, param.x_lim], ylim=[-param.x_lim, param.x_lim], zlim=[-param.x_lim, param.x_lim],
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

    """
    Plot PDF of vector alignment in polar/cylindrical coordinate
    """
    function plot_alignment(param, dat_x, dat_y, file_prefix)
        # Figure setting
        pygui(true)
        joint_plot = sns.JointGrid(
            x=dat_x,
            y=dat_y,
        )

        # Define range of kernel density estimation
        if file_prefix[2] == "phi_costheta"  # x:ϕ y:cosθ
            clip_range = [(0.0, π), (-1.0, 1.0)]
        elseif file_prefix[2] == "cosphi_costheta"  # x:cosϕ y:cosθ
            clip_range = [(-1.0, 1.0), (-1.0, 1.0)]
        end

        # Plot KDE & scatter
        joint_plot = joint_plot.plot_joint(
            sns.kdeplot,
            # cbar=True,
            cmap="plasma",
            clip=clip_range
        )
        joint_plot = joint_plot.plot_joint(
            sns.scatterplot,
            marker="+",
            color="black",
            linewidth=2
        )

        # Plot attribution
        if file_prefix[1] == "polar"
            joint_plot.ax_joint.set_title("Alignment in 3D polar coordinates")
        elseif file_prefix[1] == "cylindrical"
            joint_plot.ax_joint.set_title("Alignment in 3D cylindrical coordinates")
        end
        if file_prefix[2] == "phi_costheta"  # x:ϕ y:cosθ
            joint_plot.ax_joint.set_xlabel(L"$\phi$")
            joint_plot.ax_joint.set_ylabel(L"$\cos \theta$")
            joint_plot.ax_joint.set_xlim([-0.1, π+0.1])
            joint_plot.ax_joint.set_xticks([0.0, π/4.0, π/2.0, π*3.0/4.0, π])
            joint_plot.ax_joint.set_xticklabels([L"$0$", L"$\pi/4$", L"$\pi/2$", L"$3\pi/4$", L"$\pi$"])
        elseif file_prefix[2] == "cosphi_costheta"  # x:cosϕ y:cosθ
            joint_plot.ax_joint.set_xlabel(L"$\cos \phi$")
            joint_plot.ax_joint.set_ylabel(L"$\cos \theta$")
            joint_plot.ax_joint.set_xlim([-1.1, 1.1])
            joint_plot.ax_joint.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        end
        joint_plot.ax_joint.set_ylim([-1.1, 1.1])
        joint_plot.ax_joint.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        joint_plot.ax_marg_x.set_axis_off()
        joint_plot.ax_marg_y.set_axis_off()

        savefig(string("./tmp/", file_prefix[1], "_", file_prefix[2], ".png"), bbox_inches="tight", pad_inches=0.1)
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
    plot_3d_vectors,
    plot_alignment


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
x_lim = 1.0
num_vectors = 500
param = ParamVar.Parameters(
    x_lim, num_vectors
)


# ----------------------------------------
## Define 3d random points (origin of vectors)
# ----------------------------------------

# Define mutable struct for vectors
vectors = Array{ParamVar.Vector}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    vectors[itr_vec] = ParamVar.Vector()
end

distribute_points(param, vectors)


# ----------------------------------------
## Define vectors in uniform random direction
# ----------------------------------------
distribute_angles(param, vectors)

# plot_3d_vectors(param, vectors)


# ----------------------------------------
## Compute vector alignment in 3D polar coordinates
## Output its PDF
# ----------------------------------------

# Define mutable structs for alignments
aligns_polar = Array{ParamVar.Alignment_Polar}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    aligns_polar[itr_vec] = ParamVar.Alignment_Polar()
end

compute_alignment_polar(param, vectors, aligns_polar)

plot_alignment(
    param,
    getfield.(aligns_polar, :ϕ),
    getfield.(aligns_polar, :cosθ),
    ["polar", "phi_costheta"]
)

plot_alignment(
    param,
    getfield.(aligns_polar, :cosϕ),
    getfield.(aligns_polar, :cosθ),
    ["polar", "cosphi_costheta"]
)

# ----------------------------------------
## Compute vector alignment in 3D cylindrical coordinates
# ----------------------------------------

# Define mutable structs for alignments
align_cylindrical = Array{ParamVar.Alignment_Cylindrical}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    align_cylindrical[itr_vec] = ParamVar.Alignment_Cylindrical()
end

compute_alignment_polar(param, vectors, align_cylindrical)

plot_alignment(
    param,
    getfield.(aligns_polar, :ϕ),
    getfield.(aligns_polar, :cosθ),
    ["cylindrical", "phi_costheta"]
)

plot_alignment(
    param,
    getfield.(aligns_polar, :cosϕ),
    getfield.(aligns_polar, :cosθ),
    ["cylindrical", "cosphi_costheta"]
)
