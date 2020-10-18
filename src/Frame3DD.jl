export frame3dd,
       frame3DD_circular_tube_truss_bar_family

Base.@kwdef struct Frame3DDBackend{K,T} <: LazyBackend{K,T}
  realized::Parameter{Bool}=Parameter(false)
  truss_nodes::Vector{<:TrussNode}=TrussNode[]
  truss_bars::Vector{<:TrussBar}=TrussBar[]
  shapes::Vector{<:Shape}=Shape[] # This contains all the rest that is not treated yet
  truss_node_data::Vector{TrussNodeData}=TrussNodeData[]
  truss_bar_data::Vector{TrussBarData}=TrussBarData[]
end

abstract type FR3DDKey end
const FR3DDId = Any
const FR3DDIds = Vector{FR3DDId}
const FR3DDRef = GenericRef{FR3DDKey, FR3DDId}
const FR3DDRefs = Vector{FR3DDRef}
const FR3DDEmptyRef = EmptyRef{FR3DDKey, FR3DDId}
const FR3DDUniversalRef = UniversalRef{FR3DDKey, FR3DDId}
const FR3DDNativeRef = NativeRef{FR3DDKey, FR3DDId}
const FR3DDUnionRef = UnionRef{FR3DDKey, FR3DDId}
const FR3DDSubtractionRef = SubtractionRef{FR3DDKey, FR3DDId}
const FR3DD = Frame3DDBackend{FR3DDKey, FR3DDId}

void_ref(b::FR3DD) = FR3DDNativeRef(-1)

const frame3dd = FR3DD()

# Frame3DD needs to merge nodes and bars
save_shape!(b::FR3DD, s::TrussNode) = maybe_merged_node(b, s)
save_shape!(b::FR3DD, s::TrussBar) = maybe_merged_bar(b, s)

# Frame3DD does not need layers
with_family_in_layer(f::Function, backend::FR3DD, family::Family) = f()

# Frame3DD Families
abstract type Frame3DDFamily <: BackendFamily{Any} end

Base.@kwdef struct Frame3DDTrussBarFamily <: Frame3DDFamily
  Ax  # cross section area
  Asy # shear area y-direction
  Asz # shear area z-direction
  Jxx # torsional moment of inertia - x axis
  Iyy # bending moment of inertia - y axis
  Izz # bending moment of inertia - z axis
  E   # elastic modulus
  G   # shear   modulus
  p   # roll angle
  d   # mass density
end

truss_bar_family_cross_section_area(f::Frame3DDTrussBarFamily) = f.Ax

#=
Circular Tube (outer radius= Ro, inner radius = Ri):

Ax = π ( Ro2 - Ri2 )
Asy = Asz = Ax / ( 0.54414 + 2.97294(Ri/Ro) - 1.51899(Ri/Ro)2 ) ± 0.05%
Jxx = (1/2) π ( Ro4 - Ri4 )
Ixx = Iyy = (1/4) π ( Ro4 - Ri4 )
=#

frame3DD_circular_tube_truss_bar_family(rₒ, rᵢ; E, G, p, d) =
  let Ax = annulus_area(rₒ, rᵢ),
      Asyz = Ax/(0.54414 + 2.97294*(rᵢ/rₒ) - 1.51899*(rᵢ/rₒ)^2),
      Jxx = π/2*(rₒ^4 - rᵢ^4),
      Ixxyy = Jxx/2
    Frame3DDTrussBarFamily(
      Ax  =Ax,
      Asy =Asyz,
      Asz =Asyz,
      Jxx =Jxx,
      Iyy =Ixxyy,
      Izz =Ixxyy,
      E   =E,
      G   =G,
      p   =p,
      d   =d)
 end

#=
 @rad1 = 25.0                   # Section diameter (mm)
 @thk1 = 1.6                    # Section thickness (mm)

 @section = "CHS_25x1.6"        # Tube name
 @material = "Steel C450"       # Duragal DualGrade RHS Grade C450Lo (AS 1163)
 @ftype = "Tube"                # Tube type (Tube, Beam, Channel, Angle)
 @emod = "210.0"                # Young's Modulus (GN per square metre) for steel
 @gmod = "79.0"                 # Modulus of Rigidity (GN per square metre) for steel
 @pang = "0.0"                  # section axis rotation angle (deg)
 @ax = "568.3"                  # cross-sectional area (mm^2)
 @asy = "3.2"                   # 2 x thickness (mm) shear area
 @asz = "3.2"                   # 2 x thickness (mm) shear area
 @jxx = "0.0585E+06"            # J (mm^4)
 @iyy = "0.0702E+06"            # Iyy (mm^4)
 @izz = "0.0237E+06"            # Izz (mm^4)
=#

set_backend_family(
   default_truss_bar_family(),
   frame3dd,
   frame3DD_circular_tube_truss_bar_family(0.0213/2, 0.0213/2-0.0026,
     E=210000000000.0, # (Young's modulus)
     G=81000000000.0, # (Kirchoff's or Shear modulus)
     p=0.0, # Roll angle
     d=77010.0)) # Density

#
delete_all_shapes(b::FR3DD) =
  begin
    empty!(b.truss_nodes)
    empty!(b.truss_bars)
    b.realized(false)
    b
  end

realize(b::FR3DD, s::TrussNode) =
  error("BUM")

realize(b::FR3DD, s::TrussBar) =
  error("BUM")

realize(b::FR3DD, s::Panel) =
  error("BUM")


#new_truss_analysis(v=nothing; self_weight=false, backend=frame3dd) =
displacements_from_frame3dd(b::FR3DD, filename, load) =
  let nodes = process_nodes(b.truss_nodes, load),
      bars = process_bars(b.truss_bars, nodes),
      supports = unique(filter(s -> s != false, map(n -> n.family.support, nodes))),
      loaded_nodes = filter(!truss_node_is_supported, nodes),
      supported_nodes = filter(truss_node_is_supported, nodes),
      shear = 0,		     # 1: include shear deformation effects, 0: don't
      geom  = 0,		     # 1: include geometric stiffness effects, 0: don't
      exagg = 1,		     # exaggeration factor for Gnuplot output
      scale = 1.0,		   # zoom scale factor for Gnuplot plotting
      dx = 1.0		       # x-axis increment for internal force calc's
    empty!(b.truss_node_data)
    append!(b.truss_node_data, nodes)
    empty!(b.truss_bar_data)
    append!(b.truss_bar_data, bars)
    open(filename, "w") do io
      println(io, "Frame analysis for Khepri");
      println(io, "%% node data ...");
      println(io, "$(length(nodes))\t\t%% number of nodes");
      println(io, "%% J\t\tX\t\tY\t\tZ\t\tr\t\tnode coordinates");
      for n in nodes
        i = n.id
        p = n.loc
        println(io, "$(i)\t$(p.x)\t$(p.y)\t$(p.z)\t0")
      end
      println(io, "");
      println(io, "%% reaction data ...");
      #nR = count(has_support, nodes);
      println(io, "$(length(supported_nodes))    %% number of nodes with reaction forces");
      println(io, "%% j\tRx\tRy\tRz\tRxx\tRyy\tRzz");
      for n in supported_nodes
        s = n.family.support
        print(io, "$(n.id),\t$(Int(s.ux)),\t$(Int(s.uy)),\t$(Int(s.uz))\t")
        println(io, "$(Int(s.rx)),\t$(Int(s.ry)),\t$(Int(s.rz))")
      end
      println(io, "")
      println(io, "%% frame element section property data ...")
      println(io, "$(length(bars))\t\t%% number of frame elements")
      println(io, "%% m\tn1\tn2\t\tAx\t\tAsy\t\tAsz\t\tJxx\t\tIyy\t\tIzz\t\tE\t\tG\t\tp\tdensity");
      for bar in bars
        f = family_ref(b, bar.family)
        print(io, "$(bar.id),\t$(bar.node1.id),\t$(bar.node2.id),")
        print(io, "\t$(f.Ax),\t$(f.Asy),\t$(f.Asz),")
        print(io, "\t$(f.Jxx),\t$(f.Iyy),\t$(f.Izz),")
        println(io, "\t$(f.E),\t$(f.G),\t$(f.p),\t$(f.d)")
      end
      println(io, "\n");
      println(io, "$(shear)\t\t%% 1: include shear deformation, 0: do not");
      println(io, "$(geom)\t\t%% 1: include geometric stiffness, 0: do not");
      println(io, "$(exagg)\t%% exagerate deformations in plotting");
      println(io, "$(scale)\t%% zoom scale factor for 3D plotting");
      println(io, "$(dx)\t%% x-axis increment for internal forces calc");
      println(io, "\n");
      println(io, "%% static load data ...");
      println(io, "1\t\t%% number of static load cases");
      println(io, "\t\t%% begin static load case 1 of 1");
      println(io, "%% gravitational acceleration for self-weight loading");
      println(io, "%% gX         gY         gZ");
      println(io, "  0.0        0.0        0.0\n");
      println(io, "$(length(loaded_nodes))\t\t%% number of loaded nodes");
      println(io, "%% j\t\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz");
      for n in loaded_nodes
        println(io, "$(n.id),\t$(load.x),\t$(load.y),\t$(load.z),\t0,\t0,\t0")
      end
      println(io, "\n");
      println(io, "0\t\t%% number of members with uniform distributed loads")
      println(io, "0\t\t%% number of members with trapezoidal loads");
      println(io, "0\t\t%% number of members with internal point loads");
      println(io, "0\t\t%% number of members with temperature loads");
      println(io, "0\t\t%% number of nodes with prescribed displacements");
      println(io, "\t\t%% end   static load case 1 of 1\n");
      println(io, "%% inertial load data ...");
      println(io, "0\t\t%% number of dynamic modes to analyze");
    end
  end

frame3dd_plugin = joinpath(dirname(dirname(abspath(@__FILE__))), "Plugins", "Frame3DD", "frame3dd.exe")

frame3dd_simulation_path() =
  mktempdir(tempdir(), prefix="Frame3DD_")

##########################################
backend_truss_analysis(b::FR3DD, load::Vec) =
  # Ensure extension is FMM to force Matlab mode
  let simulation_folder = frame3dd_simulation_path(),
#      input_path = joinpath(simulation_folder, "IOdata.FMM"),
#      output_path = joinpath(simulation_folder, "IOdata.OUT") ,
#      output_m_path = joinpath(simulation_folder, "IOdata_out.m")
       input_path = joinpath(simulation_folder, "IOdata.IN"),
       output_path = joinpath(simulation_folder, "IOdata.OUT")
    @info input_path
    displacements_from_frame3dd(b, input_path, load)
    try
      #run(`$frame3dd_plugin -w -i $input_path -o $output_path`)
      #run(`$frame3dd_plugin -i $input_path -o $output_path`)
      withenv("FRAME3DD_OUTDIR"=>simulation_folder) do
        run(`$frame3dd_plugin -q -i $input_path -o $output_path`)
      end
    catch e
    end
    lines = readlines(output_path)
    idx1 = findfirst(l->startswith(l, "N O D E   D I S P L A C E M E N T S"), lines) + 2
    idx2 = findnext(l->startswith(l, "F R A M E   E L E M E N T   E N D   F O R C E S"), lines, idx1)
    Dict([let strs = split(line)
            parse(Int, strs[1])=>vxyz(parse.(Float64, strs[2:4])...)
          end
          for line in lines[idx1:idx2-1]])
  end

node_displacement_function(b::FR3DD, results) =
  n -> get(results, n.id, vx(0))
