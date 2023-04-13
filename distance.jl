import Pkg; Pkg.activate("ComputeDistances",shared=true)
#First time only:
Pkg.add(["Chemfiles", "PDBTools", "Plots","CellListMap"])

# loading packages
using LinearAlgebra: norm
using StaticArrays
using PDBTools
using Plots
import CellListMap: wrap_relative_to
import Chemfiles

# function that computes the distance between atoms for each frame
function compute_distance(;pdbname, trajectory_name, selection1, selection2)
    pdb = readPDB(pdbname)
    dcd = Chemfiles.Trajectory(trajectory_name)
    atom1 = select(pdb, selection1)
    atom2 = select(pdb, selection2)
    distances = Float64[]
    for frame in dcd
        # read unit cell matrix from trajectory
        matrix_read = Chemfiles.matrix(Chemfiles.UnitCell(frame))
        unit_cell = SMatrix{3,3}(transpose(matrix_read))
        # read coordinates, convert them to small static vectors
        coordinates = Chemfiles.positions(frame)
        x = SVector{3}(@view(coordinates[:,atom1[1].index]))
        y = SVector{3}(@view(coordinates[:,atom2[1].index]))
        # wrap coordinates according to PBC. Note that here we use an
        # internal function of `CellListMap`. 
        y_wrapped = wrap_relative_to(y,x,unit_cell)
        # compute distance
        d = norm(y_wrapped-x)
        # add data to distance array
        push!(distances, d)
    end
    return distances
end

#
# Examples
#

# NAMD trajectory file
function run_namd()
    distances = compute_distance(
        pdbname = "./data/structure.pdb", 
        trajectory_name ="./data/trajectory.dcd",
        selection1 = "residue 1 and name CA",
        selection2 = "residue 8000 and name OH2",
    )
    Plots.default(fontfamily="Computer Modern")
    plt = plot(distances,
        xlabel = "frame",
        ylabel = "distance / Å",
        framestyle=:box,
        linewidth=2,
        label=nothing
    )
    # savefig("./plot.svg") # to save the figure
    display(plt)
    return plt, distances
end
        
        
# Gromacs trajectory file
function run_gromacs()
    distances = compute_distance(
        pdbname = "./data/system.pdb", 
        trajectory_name ="./data/trajectory.xtc",
        selection1 = "index 81044",
        selection2 = "index 78821",
    )
    Plots.default(fontfamily="Computer Modern")
    plt = plot(distances,
        xlabel = "frame",
        ylabel = "distance / Å",
        framestyle=:box,
        linewidth=2,
        label=nothing
    )
    # savefig("./plot.svg") # to save the figure
    display(plt)
    return plt, distances
end
        
