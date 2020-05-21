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
        flow_swirling::Bool  # true: vectors swirling in xy plane. false: uniformly random vectors
        file_postfix::String
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
        cosθ::Float64  # angle between vector & z axis
        ϕ::Float64  # angle between vector projected onto xy plane & x axis
        cosϕ::Float64
        # Constructor
        Alignment_Polar() = new()
    end

    """
    Alignment in terms of 3D cylindrical coordinates
    """
    mutable struct Alignment_Cylindrical
        cosθ::Float64  # θ: angle between vector & z axis
        ϕ::Float64  # ϕ: angle between vector projected onto xy plane & ψ vector
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
    function distribute_random_angles(param, vectors)
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

    """
    Define vectors in swirling (in xy plane) direction
    """
    function distribute_swirling_angles(param, vectors)
        θ = rand(Uniform(π/2-0.05*π/2, π/2+0.05*π/2), param.num_vectors)  # 0≤θ≤π

        for itr_vec in 1:param.num_vectors
            x = vectors[itr_vec].x
            y = vectors[itr_vec].y

            # ψ: angle between vector origin on xy plane & coordinate origin
            cosψ = x / sqrt(x^2 + y^2)
            sinψ = y / sqrt(x^2 + y^2)

            # vr: vector in r direction (perpendicular to ψ)
            # vψ: vector in ψ direction (perpendicular to r)
            vr = 0.0
            vψ = 0.5 * sqrt(x^2 + y^2)

            vectors[itr_vec].vx = cosψ * vr - sinψ * vψ
            vectors[itr_vec].vy = sinψ * vr + cosψ * vψ
            vectors[itr_vec].vz = cos(θ[itr_vec])
        end
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
            if 1.0<cosϕ<1.001  # Handle error with swirling flow
                cosϕ = 1.0
            end
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
        savefig(
            string("./tmp/vectors", param.file_postfix),
            bbox_inches="tight", pad_inches=0.1
        )
    end

    """
    Plot PDF of vector alignment in polar/cylindrical coordinate
    """
    function plot_alignment(param, dat_x, dat_y, file_prefix)
        pygui(true)

        # Define range of kernel density estimation
        if file_prefix[2] == "phi_costheta"  # x:ϕ y:cosθ
            clip_range = [(0.0, π), (-1.0, 1.0)]
        elseif file_prefix[2] == "cosphi_costheta"  # x:cosϕ y:cosθ
            clip_range = [(-1.0, 1.0), (-1.0, 1.0)]
        end

        # Plot kde & scatter
        g = (sns.jointplot(
            x=dat_x,
            y=dat_y,
            kind="kde",
            # cbar=True,
            cmap="plasma",
            clip=clip_range,
            shade=false,
            color="gold"
        ).plot_joint(
            sns.scatterplot,
            marker="+",
            color="black",
            linewidth=1
        ))

        # Plot attribution
        if file_prefix[1] == "polar"
            g.ax_joint.text(
                0.25, 1.05,
                "polar"
            )
        elseif file_prefix[1] == "cylindrical"
            g.ax_joint.text(
                0.25, 1.05,
                "cylindrical"
            )
        end
        if file_prefix[2] == "phi_costheta"  # x:ϕ y:cosθ
            g.ax_joint.set_xlabel(L"$\phi$")
            g.ax_joint.set_ylabel(L"$\cos \theta$")
            g.ax_joint.set_xlim([-0.1, π+0.1])
            g.ax_joint.set_xticks([0.0, π/4.0, π/2.0, π*3.0/4.0, π])
            g.ax_joint.set_xticklabels([L"$0$", L"$\pi/4$", L"$\pi/2$", L"$3\pi/4$", L"$\pi$"])
        elseif file_prefix[2] == "cosphi_costheta"  # x:cosϕ y:cosθ
            g.ax_joint.set_xlabel(L"$\cos \phi$")
            g.ax_joint.set_ylabel(L"$\cos \theta$")
            g.ax_joint.set_xlim([-1.1, 1.1])
            g.ax_joint.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        end
        g.ax_joint.set_ylim([-1.1, 1.1])
        g.ax_joint.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

        savefig(
            string("./tmp/", file_prefix[1], "_", file_prefix[2], param.file_postfix),
            bbox_inches="tight", pad_inches=0.1
        )
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
    distribute_random_angles,
    distribute_swirling_angles
using .ComputeAlignments:
    compute_alignment_polar,
    compute_alignment_cylindrical
using .PlotFigures:
    plot_3d_vectors,
    plot_alignment


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
x_lim = 1.0
num_vectors = 1000
flag_siwrling = true
if flag_siwrling
    file_postfix = "_swirling.png"
else
    file_postfix = "_random.png"
end
param = ParamVar.Parameters(
    x_lim, num_vectors, flag_siwrling, file_postfix
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
if flag_siwrling
    distribute_swirling_angles(param, vectors)
else
    distribute_random_angles(param, vectors)
end

plot_3d_vectors(param, vectors)


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
## Output its PDF
# ----------------------------------------

# Define mutable structs for alignments
align_cylindrical = Array{ParamVar.Alignment_Cylindrical}(undef, param.num_vectors)
for itr_vec = 1:param.num_vectors
    align_cylindrical[itr_vec] = ParamVar.Alignment_Cylindrical()
end

compute_alignment_cylindrical(param, vectors, align_cylindrical)

plot_alignment(
    param,
    getfield.(align_cylindrical, :ϕ),
    getfield.(align_cylindrical, :cosθ),
    ["cylindrical", "phi_costheta"]
)

plot_alignment(
    param,
    getfield.(align_cylindrical, :cosϕ),
    getfield.(align_cylindrical, :cosθ),
    ["cylindrical", "cosphi_costheta"]
)
