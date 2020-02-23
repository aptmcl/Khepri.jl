module Khepri
using Base.Iterators
using Sockets
using Dates
using ColorTypes
#using PyCall

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
include("Unity.jl")
include("Unreal.jl")
include("Camera.jl")
#const com = PyNULL()
#include("Robot.jl")
end
