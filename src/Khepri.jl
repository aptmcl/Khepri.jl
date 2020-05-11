module Khepri
using Base.Iterators
using Sockets
using Dates
using ColorTypes
using StaticArrays
using LinearAlgebra
using IntervalSets
using Interpolations
using PyCall
import Base.show, Base.zero, Base.iterate
import Base.convert
import Base.getindex, Base.firstindex, Base.lastindex, Base.broadcastable

include("Parameters.jl")
include("Utils.jl")
include("Coords.jl")
include("Paths.jl")
include("Shapes.jl")
include("BIM.jl")
include("Geometry.jl")
include("Encoders.jl")
include("Primitives.jl")
include("Tikz.jl")
include("AutoCAD.jl")
include("Rhinoceros.jl")
include("Revit.jl")
include("Grasshopper.jl")
#include("ArchiCAD.jl")
include("Radiance.jl")
include("POVRay.jl")
#include("Jupyter.jl")
#include("Plot.jl")
include("Unity.jl")
include("Unreal.jl")
include("Camera.jl")
const com = PyNULL()
include("Robot.jl")
include("BIMFamilies.jl")
end
