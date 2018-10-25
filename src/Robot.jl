export robot,
       create_node_support,
       new_robot_analysis,
       node_displacement,
       displacements,
       nodes


# generate code to test whether expr is in the given set of values
function membershiptest(expr, values)
    lo, hi = extrema(values)
    if length(values) == hi - lo + 1
        :($lo <= $expr <= $hi)
    elseif length(values) < 20
        foldl((x1,x2)->:($x1 || ($expr == $x2)), values[2:end]; init=:($expr == $(values[1])))
    else
        :($expr in $(Set(values)))
    end
end


import Core.Intrinsics.bitcast

macro enums(T, syms...)
    if isempty(syms)
        throw(ArgumentError("no arguments given for Enum $T"))
    end
    basetype = Int32
    typename = T
    if isa(T, Expr) && T.head == :(::) && length(T.args) == 2 && isa(T.args[1], Symbol)
        typename = T.args[1]
        basetype = Core.eval(__module__, T.args[2])
        if !isa(basetype, DataType) || !(basetype <: Integer) || !isbitstype(basetype)
            throw(ArgumentError("invalid base type for Enum $typename, $T=::$basetype; base type must be an integer primitive type"))
        end
    elseif !isa(T, Symbol)
        throw(ArgumentError("invalid type expression for enum $T"))
    end
    vals = Vector{Tuple{Symbol,Integer}}()
    lo = hi = 0
    i = zero(basetype)
    hasexpr = false

    if length(syms) == 1 && syms[1] isa Expr && syms[1].head == :block
        syms = syms[1].args
    end
    for s in syms
        s isa LineNumberNode && continue
        if isa(s, Symbol)
            if i == typemin(basetype) && !isempty(vals)
                throw(ArgumentError("overflow in value \"$s\" of Enum $typename"))
            end
        elseif isa(s, Expr) &&
               (s.head == :(=) || s.head == :kw) &&
               length(s.args) == 2 && isa(s.args[1], Symbol)
            i = Core.eval(__module__, s.args[2]) # allow exprs, e.g. uint128"1"
            if !isa(i, Integer)
                throw(ArgumentError("invalid value for Enum $typename, $s; values must be integers"))
            end
            i = convert(basetype, i)
            s = s.args[1]
            hasexpr = true
        else
            throw(ArgumentError(string("invalid argument for Enum ", typename, ": ", s)))
        end
        if !Base.isidentifier(s)
            throw(ArgumentError("invalid name for Enum $typename; \"$s\" is not a valid identifier."))
        end
        push!(vals, (s,i))
        if length(vals) == 1
            lo = hi = i
        else
            lo = min(lo, i)
            hi = max(hi, i)
        end
        i += oneunit(i)
    end
    values = basetype[i[2] for i in vals]
#    if hasexpr && values != unique(values)
#        throw(ArgumentError("values for Enum $typename are not unique"))
#    end
    blk = quote
        # enum definition
        Base.@__doc__(primitive type $(esc(typename)) <: Enum{$(basetype)} $(sizeof(basetype) * 8) end)
        function $(esc(typename))(x::Integer)
            $(membershiptest(:x, values)) || enum_argument_error($(Expr(:quote, typename)), x)
            return bitcast($(esc(typename)), convert($(basetype), x))
        end
        Base.Enums.basetype(::Type{$(esc(typename))}) = $(esc(basetype))
        Base.typemin(x::Type{$(esc(typename))}) = $(esc(typename))($lo)
        Base.typemax(x::Type{$(esc(typename))}) = $(esc(typename))($hi)
        Base.isless(x::$(esc(typename)), y::$(esc(typename))) = isless($basetype(x), $basetype(y))
        let insts = ntuple(i->$(esc(typename))($values[i]), $(length(vals)))
            Base.instances(::Type{$(esc(typename))}) = insts
        end
        function Base.print(io::IO, x::$(esc(typename)))
            for (sym, i) in $vals
                if i == $(basetype)(x)
                    print(io, sym); break
                end
            end
        end
        function Base.show(io::IO, x::$(esc(typename)))
            if get(io, :compact, false)
                print(io, x)
            else
                print(io, x, "::")
                show(IOContext(io, :compact => true), typeof(x))
                print(io, " = ")
                show(io, $basetype(x))
            end
        end
        function Base.show(io::IO, ::MIME"text/plain", t::Type{$(esc(typename))})
            print(io, "Enum ")
            Base.show_datatype(io, t)
            print(io, ":")
            for (sym, i) in $vals
                print(io, "\n", sym, " = ")
                show(io, i)
            end
        end
    end
    if isa(typename, Symbol)
        for (sym,i) in vals
            push!(blk.args, :(const $(esc(sym)) = $(esc(typename))($i)))
        end
    end
    push!(blk.args, :nothing)
    blk.head = :toplevel
    return blk
end


@enums IRobotProjectType begin
  I_PT_FRAME_2D=1
  I_PT_TRUSS_2D
  I_PT_GRILLAGE
  I_PT_FRAME_3D
  I_PT_TRUSS_3D
  I_PT_PLATE
  I_PT_SHELL
  I_PT_AXISYMMETRIC
  I_PT_VOLUMETRIC
  I_PT_CONCRETE_BEAM
  I_PT_CONCRETE_COLUMN
  I_PT_FOUNDATION
  I_PT_PARAMETRIZED
  I_PT_STEEL_CONNECTION
  I_PT_SECTION
  I_PT_PLANE_STRESS
  I_PT_PLANE_DEFORMATION
  I_PT_DEEP_BEAM
  I_PT_CONCRETE_DEEP_BEAM=18
  I_PT_BUILDING
end

@enums IRobotObjectType begin
  I_OT_NODE = 0
  I_OT_BAR
  I_OT_CASE
  I_OT_FAMILY
  I_OT_PANEL
  I_OT_FINITE_ELEMENT
  I_OT_GEOMETRY
  I_OT_VOLUME
  I_OT_OBJECT = -2
  I_OT_UNDEFINED = -1
end

@enums IRobotLabelType begin
  I_LT_NODE_SUPPORT = 0
  I_LT_NODE_RELEASE
  I_LT_NODE_COMPATIBILITY
  I_LT_BAR_SECTION
  I_LT_BAR_RELEASE
  I_LT_BAR_OFFSET
  I_LT_BAR_CABLE
  I_LT_BAREND_BRACKET
  I_LT_BAR_MATERIAL
  I_LT_EDGE_SUPPORT = 10
  I_LT_PANEL_THICKNESS
  I_LT_PANEL_REINFORCEMENT
  I_LT_UNKNOWN = -1
  #I_LT_SUPPORT
  I_LT_MATERIAL = 8
  I_LT_LINEAR_RELEASE = 13
  I_LT_BAR_ELASTIC_GROUND
  I_LT_NODE_RIGID_LINK
  I_LT_MEMBER_TYPE
  I_LT_VEHICLE
  I_LT_SOLID_PROPERTIES
  I_LT_BAR_GEO_IMPERFECTIONS
  I_LT_BAR_NONLINEAR_HINGE
  I_LT_CLADDING
  I_LT_PANEL_CALC_MODEL
  I_LT_MEMBER_REINFORCEMENT_PARAMS = 25
  I_LT_COUNT
end

@enums IRobotMaterialType begin
  I_MT_ALL
  I_MT_STEEL
  I_MT_CONCRETE
  I_MT_ALUMINIUM
  I_MT_TIMBER
  I_MT_OTHER
end

@enums IRobotBarSectionType begin
  I_BST_STANDARD
  I_BST_NS_BOX
  I_BST_NS_I
  I_BST_NS_II
  I_BST_NS_TUBE
  I_BST_NS_RECT
  I_BST_NS_C
  I_BST_NS_L
  I_BST_NS_LP
  I_BST_NS_Z
  I_BST_NS_ZP
  I_BST_NS_T
  I_BST_NS_H
  I_BST_NS_XT
  I_BST_NS_XI
  I_BST_NS_DRECT
  I_BST_NS_CROSS
  I_BST_NS_HOLE = 99
  I_BST_COMPLEX = 1000
  I_BST_NS_BOX_2 = 17
  I_BST_NS_POLYGONAL
  I_BST_SPECIAL = 1001
  I_BST_JOIST
  I_BST_NS_BOX_3 = 20
end

@enums IRobotBarSectionShapeType begin
  I_BSST_UNKNOWN
  I_BSST_CAE
  I_BSST_CAEP
  I_BSST_CAI
  I_BSST_CAIP
  I_BSST_DCEC
  I_BSST_DCED
  I_BSST_DCEP
  I_BSST_DCIG
  I_BSST_DCIP
  I_BSST_HEA
  I_BSST_HEAA
  I_BSST_HEB
  I_BSST_HEC
  I_BSST_HEM
  I_BSST_HER
  I_BSST_HHEA
  I_BSST_HHEB
  I_BSST_HHEM
  I_BSST_IIPE
  I_BSST_IPE
  I_BSST_IPEA
  I_BSST_IPEO
  I_BSST_IPER
  I_BSST_IPEV
  I_BSST_IPN
  I_BSST_MHEA
  I_BSST_MHEB
  I_BSST_MHEM
  I_BSST_MIPE
  I_BSST_PRS
  I_BSST_TCAR
  I_BSST_TEAE
  I_BSST_TEAI
  I_BSST_THEX
  I_BSST_TREC
  I_BSST_TRON
  I_BSST_UAP
  I_BSST_UPN
  I_BSST_UUAP
  I_BSST_UUPN
  I_BSST_FRTG
  I_BSST_UPAF
  I_BSST_BOX = 91
  I_BSST_RECT
  I_BSST_TUBE
  I_BSST_ISYM
  I_BSST_INSYM
  I_BSST_TUSER
  I_BSST_CUSER
  I_BSST_TBETC
  I_BSST_WELD_CROSS = 201
  I_BSST_DRECT = 99
  I_BSST_COLD_SIGMA1 = 1001
  I_BSST_COLD_SIGMA2
  I_BSST_COLD_ZED1
  I_BSST_COLD_U
  I_BSST_COLD_CE1
  I_BSST_COLD_ANGL
  I_BSST_COLD_OMEGA
  I_BSST_COLD_SO1
  I_BSST_COLD_RIVE1
  I_BSST_USER_BOX = 91
  I_BSST_USER_RECT
  I_BSST_USER_TUBE
  I_BSST_USER_I_BISYM
  I_BSST_USER_I_MONOSYM
  I_BSST_USER_T_SHAPE
  I_BSST_USER_C_SHAPE
  I_BSST_USER_CROSS = 201
  I_BSST_WOOD_RECT = 41
  I_BSST_WOOD_DRECT = 99
  I_BSST_RECT_FILLED = 43
  I_BSST_CIRC_FILLED
  I_BSST_CONCR_COL_R = -108
  I_BSST_CONCR_COL_T
  I_BSST_CONCR_COL_L
  I_BSST_CONCR_COL_Z
  I_BSST_CONCR_COL_P
  I_BSST_CONCR_COL_C
  I_BSST_CONCR_COL_CH
  I_BSST_CONCR_COL_CQ
  I_BSST_CONCR_BEAM_RECT = -3
  I_BSST_CONCR_BEAM_T
  I_BSST_CONCR_BEAM
  I_BSST_COMP_2C_FACE = 1101
  I_BSST_COMP_2C_BACK
  I_BSST_COMP_2I
  I_BSST_COMP_CI
  I_BSST_COMP_2LI
  I_BSST_COMP_4L_FACE
  I_BSST_COMP_4L_BACK
  I_BSST_COMP_2L_SHORT
  I_BSST_COMP_2L_LONG
  I_BSST_COMP_2L_CROSS
  I_BSST_USER_BOX_2 = 102
  I_BSST_CCL = 45
  I_BSST_URND
  I_BSST_TRND
  I_BSST_CUAP
  I_BSST_WOOD_CIRC = 100
  I_BSST_USER_CIRC_FILLED
  I_BSST_USER_POLYGONAL = 103
  I_BSST_COLD_C_PLUS = 1010
  I_BSST_COLD_SIGMA_SL
  I_BSST_COLD_SIGMA
  I_BSST_COLD_Z
  I_BSST_COLD_L_LIPS
  I_BSST_COLD_Z_ROT
  I_BSST_COMP_2L_FACE_SHORT = 1111
  I_BSST_COMP_2L_FACE_LONG
  I_BSST_COMP_CI_BACK
  I_BSST_COMP_2C_FACE_WELD = 1201
  I_BSST_COMP_2C_BACK_WELD
  I_BSST_COMP_2I_WELD
  I_BSST_COMP_CI_WELD
  I_BSST_COMP_2LI_WELD
  I_BSST_COMP_4L_FACE_WELD
  I_BSST_COMP_4L_BACK_WELD
  I_BSST_COMP_2L_SHORT_WELD
  I_BSST_COMP_2L_LONG_WELD
  I_BSST_COMP_2L_CROSS_WELD
  I_BSST_COMP_2L_FACE_SHORT_WELD
  I_BSST_COMP_2L_FACE_LONG_WELD
  I_BSST_COMP_CI_BACK_WELD
  I_BSST_SPEC_CORRUGATED_WEB = 104
  I_BSST_SPEC_CASTELLATED_WEB_HEXAGONAL_OPENINGS = 50
  I_BSST_SPEC_CASTELLATED_WEB_HEXAGONAL_OPENINGS_SHIFTED = 52
  I_BSST_SPEC_CASTELLATED_WEB_ROUND_OPENINGS = 51
  I_BSST_CONCR_BEAM_I = -4
  I_BSST_SPEC_SFB = 53
  I_BSST_SPEC_IFBA
  I_BSST_SPEC_IFBB
  I_BSST_JOIST_K = 500
  I_BSST_JOIST_LH
  I_BSST_JOIST_KCS
  I_BSST_JOIST_DLH
  I_BSST_JOIST_SLH
  I_BSST_USER_BOX_3 = 105
  I_BSST_JOIST_G = 510
  I_BSST_JOIST_VG
  I_BSST_JOIST_BG
end

@enums IRobotBarEndReleaseValue begin
  I_BERV_NONE
  I_BERV_STD
  I_BERV_PLUS
  I_BERV_MINUS
  I_BERV_FIXED = 1
  I_BERV_ELASTIC = 4
  I_BERV_ELASTIC_PLUS
  I_BERV_ELASTIC_MINUS
  I_BERV_NONLINEAR
  I_BERV_ELASTIC_REDUCED
  I_BERV_ELASTIC_REDUCED_PLUS
  I_BERV_ELASTIC_REDUCED_MINUS
end

@enums IRobotBarSectionDataValue begin
  I_BSDV_SURFACE = 10
  I_BSDV_WEIGHT = 11
  I_BSDV_D = 12
  I_BSDV_BF = 13
  I_BSDV_TW = 14
  I_BSDV_TF = 15
  I_BSDV_RA = 16
  I_BSDV_RI = 17
  I_BSDV_S = 18
  I_BSDV_ZY = 19
  I_BSDV_ZZ = 20
  I_BSDV_WX = 21
  I_BSDV_WY = 22
  I_BSDV_WZ = 23
  I_BSDV_GAMMA = 24
  I_BSDV_IOMEGA = 25
  I_BSDV_P1_LENGTH = 26
  I_BSDV_P1_THICKNESS = 27
  I_BSDV_P2_LENGTH = 28
  I_BSDV_P2_THICKNESS = 29
  I_BSDV_P3_LENGTH = 30
  I_BSDV_P3_THICKNESS = 31
  I_BSDV_P4_LENGTH = 32
  I_BSDV_P4_THICKNESS = 33
  I_BSDV_BF2 = 34
  I_BSDV_TF2 = 35
  I_BSDV_DIM1 = 36
  I_BSDV_DIM2 = 37
  I_BSDV_DIM3 = 38
  I_BSDV_ANGLE1 = 39
  I_BSDV_ANGLE2 = 40
end

@enums IRobotBarSectionNonstdDataValue begin
  I_BSNDV_BOX_H
  I_BSNDV_BOX_B
  I_BSNDV_BOX_TF
  I_BSNDV_BOX_TW
  I_BSNDV_I_H = 0
  I_BSNDV_I_B
  I_BSNDV_I_TF
  I_BSNDV_I_TW
  I_BSNDV_II_H = 0
  I_BSNDV_II_TW
  I_BSNDV_II_TF1
  I_BSNDV_II_TF2
  I_BSNDV_II_B1
  I_BSNDV_II_B2
  I_BSNDV_TUBE_D = 0
  I_BSNDV_TUBE_T
  I_BSNDV_RECT_B = 0
  I_BSNDV_RECT_H
  I_BSNDV_RECT_T
  I_BSNDV_C_H = 0
  I_BSNDV_C_B
  I_BSNDV_C_TF
  I_BSNDV_C_TW
  I_BSNDV_L_H = 0
  I_BSNDV_L_B
  I_BSNDV_L_TF
  I_BSNDV_L_TW
  I_BSNDV_Z_H = 0
  I_BSNDV_Z_B
  I_BSNDV_Z_TF
  I_BSNDV_Z_TW
  I_BSNDV_T_H = 0
  I_BSNDV_T_B
  I_BSNDV_T_TF
  I_BSNDV_T_TW
  I_BSNDV_H_H = 0
  I_BSNDV_H_B
  I_BSNDV_H_B1
  I_BSNDV_H_TF
  I_BSNDV_H_TW
  I_BSNDV_XT_H = 0
  I_BSNDV_XT_H1
  I_BSNDV_XT_B
  I_BSNDV_XT_TF
  I_BSNDV_XT_TW
  I_BSNDV_XI_H = 0
  I_BSNDV_XI_H1
  I_BSNDV_XI_B
  I_BSNDV_XI_TF
  I_BSNDV_XI_TW
  I_BSNDV_DRECT_B = 0
  I_BSNDV_DRECT_H
  I_BSNDV_DRECT_D
  I_BSNDV_CROSS_P1_L = 0
  I_BSNDV_CROSS_P2_L
  I_BSNDV_CROSS_P3_L
  I_BSNDV_CROSS_P4_L
  I_BSNDV_CROSS_P1_T
  I_BSNDV_CROSS_P2_T
  I_BSNDV_CROSS_P3_T
  I_BSNDV_CROSS_P4_T
  I_BSNDV_HOLE_B = 0
  I_BSNDV_HOLE_H
  I_BSNDV_HOLE_X
  I_BSNDV_HOLE_Z
  I_BSNDV_BOX_2_H = 0
  I_BSNDV_BOX_2_B
  I_BSNDV_BOX_2_B1
  I_BSNDV_BOX_2_TF
  I_BSNDV_BOX_2_TW
  I_BSNDV_POLYGONAL_N = 0
  I_BSNDV_POLYGONAL_D
  I_BSNDV_POLYGONAL_T
  I_BSNDV_POLYGONAL_D_IS_INT
  I_BSNDV_BOX_3_H = 0
  I_BSNDV_BOX_3_B
  I_BSNDV_BOX_3_TW
  I_BSNDV_BOX_3_TF
  I_BSNDV_BOX_3_B1
  I_BSNDV_BOX_3_B2
  I_BSNDV_BOX_3_TF2
end

@enums IRobotCaseNature begin
  I_CN_PERMANENT
  I_CN_EXPLOATATION
  I_CN_WIND
  I_CN_SNOW
  I_CN_TEMPERATURE
  I_CN_ACCIDENTAL
  I_CN_SEISMIC
end

@enums IRobotCaseAnalizeType begin
  I_CAT_COMB_NONLINEAR_INCREMENTAL = -5
  I_CAT_COMB_BUCKLING
  I_CAT_COMB_NONLINEAR_BUCKLING
  I_CAT_COMB_INCREMENTAL
  I_CAT_COMB_NONLINEAR
  I_CAT_COMB
  I_CAT_STATIC_LINEAR
  I_CAT_STATIC_NONLINEAR
  I_CAT_STATIC_INCREMENTAL
  I_CAT_STATIC_FLAMBEMENT
  I_CAT_STATIC_LINEAR_AUXILIARY
  I_CAT_STATIC_NONLINEAR_INCREMENTAL
  I_CAT_STATIC_NONLINEAR_FLAMBEMENT
  I_CAT_STATIC_NONLINEAR_MODAL
  I_CAT_DYNAMIC_MODAL = 11
  I_CAT_DYNAMIC_SPECTRAL
  I_CAT_DYNAMIC_SEISMIC
  I_CAT_DYNAMIC_HARMONIC
  I_CAT_TEMPORAR = 20
  I_CAT_MOBILE_MAIN = 30
  I_CAT_MOBILE_DERIVED
  I_CAT_NULL = 99
  I_CAT_COMB_CODE = -12
  I_CAT_STATIC_BUCKLING = 4
  I_CAT_STATIC_NONLINEAR_BUCKLING = 7
  I_CAT_DYNAMIC_NONLINEAR_MODAL_WITH_STATIC_FORCE
  I_CAT_TIME_HISTORY = 20
  I_CAT_PUSH_OVER = 25
  I_CAT_DYNAMIC_FRF = 15
  I_CAT_DYNAMIC_FOOTFALL
  I_CAT_STATIC_ELF_SEISMIC = 9
end

@enums IRobotLoadRecordType begin
  I_LRT_NODE_FORCE
  I_LRT_NODE_DISPLACEMENT
  I_LRT_BAR_DILATATION
  I_LRT_BAR_FORCE_CONCENTRATED
  I_LRT_BAR_MOMENT_DISTRIBUTED
  I_LRT_BAR_UNIFORM
  I_LRT_BAR_TRAPEZOIDALE
  I_LRT_BAR_DEAD
  I_LRT_BAR_THERMAL
  I_LRT_LINEAR = 21
  I_LRT_LINEAR_3D = 19
  I_LRT_NODE_AUXILIARY
  I_LRT_POINT_AUXILIARY = 23
  I_LRT_IN_3_POINTS = 22
  I_LRT_PRESSURE = 24
  I_LRT_THERMAL_IN_3_POINTS
  I_LRT_IN_CONTOUR = 28
  I_LRT_NODE_FORCE_MASS = 30
  I_LRT_BAR_FORCE_CONCENTRATED_MASS = 33
  I_LRT_BAR_UNIFORM_MASS = 35
  I_LRT_BAR_TRAPEZOIDALE_MASS
  I_LRT_MASS_ACTIVATION = 39
  I_LRT_SPECTRUM_VALUE
  I_LRT_UNIFORM = 26
  I_LRT_THERMAL = 8
  I_LRT_NODE_FORCE_IN_POINT = 23
  I_LRT_LINEAR_ON_EDGES = 69
  I_LRT_DEAD = 7
  I_LRT_SURFACE_ON_OBJECT = 70
  I_LRT_MOBILE_POINT_FORCE = 53
  I_LRT_MOBILE_DISTRIBUTED = 55
  I_LRT_NODE_VELOCITY = 10
  I_LRT_NODE_ACCELERATION
end

@enums IRobotNodeForceRecordValues begin
  I_NFRV_FX
  I_NFRV_FY
  I_NFRV_FZ
  I_NFRV_CX
  I_NFRV_CY
  I_NFRV_CZ
  I_NFRV_ALPHA = 8
  I_NFRV_BETA
  I_NFRV_GAMMA
end

#=
using PyCall
@pyimport array
@pyimport win32com.client as com
rh = com.Dispatch("Rhino.Interface.6")
rh[:Visible] = true
rh[:IsInitialized]
script = rh[:GetScriptObject]
pycall(rh, PyObject, "GetScriptObject")
for k in keys(script) println(k) end
script.AddCircle()
get(script, PyAny, "AddCircle")
script[:AddCircle]
([1,2,3],4)


app = com.GetActiveObject("AutoCAD.Application")
doc = app[:ActiveDocument]
db = doc[:ModelSpace]
util = doc[:Utility]
db["AddCircle"](PyVector([0.0, 0.0, 0.0]), float(4))

=#

#=
r = com.Dispatch("Robot.Application")
r[:Visible] = 1
p = r[:Project]
p[:New](5)
s = p[:Structure]
ns = s[:Nodes]
bs = s[:Bars]

ns[:Create](1, 2.0, 3.0, 4.0)
ns[:Create](2, 4.0, 5.0, 6.0)
bs[:Create](1, 1, 2)
=#

## Properties
name_method(name) =
  isa(name, Symbol) ?
    (name, name |> string |> titlecase |> s -> replace(s, r"_" => "") |> Symbol) :
    (name.args[1], name.args[2])

macro def_rw_property(name, InType, OutType)
    name, method = name_method(name)
    quote
        function $(esc(name))(in::$(esc(InType)))::$(esc(OutType)) $(esc(OutType))(in[$(QuoteNode(method))]) end
        function $(esc(name))(in::$(esc(InType)), val::$(esc(OutType)))::Nothing
            in[$(QuoteNode(method))] = $(Core.eval(__module__, OutType) <: Enum ? :(Int(val)) : :(val))
            nothing
        end
    end
end

macro def_ro_property(name, InType, OutType)
    name, method = name_method(name)
    quote
        function $(esc(name))(in::$(esc(InType)))::$(esc(OutType)) $(esc(OutType))(in[$(QuoteNode(method))]) end
    end
end

#=
@macroexpand(@def_rw_property(visible, Any, Bool))
=#

macro def_com_type(name)
    quote
        $(esc(name)) = $(esc(PyObject))
    end
end

@def_com_type IRobotApplication
@def_com_type IRobotProject
#@def_com_type IRobotProjectType
@def_com_type IRobotStructure
@def_com_type IRobotLabel
@def_com_type IRobotNodeServer
@def_com_type IRobotBarServer
@def_com_type IRobotLabelServer
@def_com_type IRobotNode
@def_com_type IRobotBar
@def_com_type IRobotSimpleCase
@def_com_type IRobotLoadRecordMngr
@def_com_type IRobotLoadRecord
@def_com_type IRobotNodeSuportData
@def_com_type IRobotBarSectionData
@def_com_type RobotBarSectionConcreteData
@def_com_type IRobotBarSectionNonstdData
@def_com_type IRobotBarReleaseData
@def_com_type StartNode
@def_com_type EndNode
@def_com_type IRobotMaterialData
@def_com_type IRobotMaterialTimberType
@def_com_type IRobotCaseServer
@def_com_type IRobotResultsServer
@def_com_type IRobotNodeResultServer
@def_com_type IRobotBarResultServer
@def_com_type IRobotCalcEngine
@def_com_type IRobotBarSectionConcreteData
@def_com_type IRobotNodeDisplacementServer
@def_com_type IRobotBarDisplacementServer
@def_com_type IRobotBarStressServer
@def_com_type IRobotBarStressData
@def_com_type IRobotSelectionFactory
@def_com_type IRobotSelection

Double = Float64
Long = Int64
Void = Nothing
Boolean = Bool

@def_ro_property project IRobotApplication IRobotProject
@def_ro_property nodes IRobotStructure IRobotNodeServer
@def_ro_property bars IRobotStructure IRobotBarServer

@def_ro_property structure IRobotProject IRobotStructure
@def_ro_property calc_engine IRobotProject IRobotCalcEngine
@def_ro_property results IRobotStructure IRobotResultsServer
@def_ro_property labels IRobotStructure IRobotLabelServer
@def_ro_property cases IRobotStructure IRobotCaseServer
#ERROR @def_ro_property selections IRobotStructure IRobotSelection
@def_ro_property selections IRobotStructure IRobotSelectionFactory
@def_ro_property records IRobotSimpleCase IRobotLoadRecordMngr
@def_ro_property objects IRobotLoadRecord IRobotSelection
@def_ro_property data IRobotLabel IRobotNodeSuportData
#@def_ro_property data IRobotLabel IRobotBarSectionData
@def_rw_property UX IRobotNodeSuportData Integer
@def_rw_property UY IRobotNodeSuportData Integer
@def_rw_property UZ IRobotNodeSuportData Integer
@def_rw_property RX IRobotNodeSuportData Integer
@def_rw_property RY IRobotNodeSuportData Integer
@def_rw_property RZ IRobotNodeSuportData Integer
@def_rw_property Gamma IRobotBar Double #Angle of rotation

@def_rw_property shape_type IRobotBarSectionData IRobotBarSectionShapeType
@def_rw_property concrete IRobotBarSectionData RobotBarSectionConcreteData
@def_rw_property MaterialName IRobotBarSectionData String
@def_rw_property start_node IRobotBarReleaseData StartNode
@def_rw_property end_node IRobotBarReleaseData EndNode

@def_rw_property CB71_Category IRobotMaterialData Integer
@def_rw_property CB71_Humidity IRobotMaterialData Double
@def_rw_property CB71_Nature IRobotMaterialData Integer
@def_rw_property CB71_Retreat IRobotMaterialData Double
@def_rw_property CS IRobotMaterialData Double
@def_rw_property Default IRobotMaterialData Bool
@def_rw_property DumpCoef IRobotMaterialData Double
@def_rw_property E IRobotMaterialData Double
@def_rw_property E_5 IRobotMaterialData Double
@def_rw_property E_Trans IRobotMaterialData Double
@def_rw_property EC_Deformation IRobotMaterialData Double
@def_rw_property GMean IRobotMaterialData Double
@def_rw_property Kirchoff IRobotMaterialData Double
@def_rw_property LX IRobotMaterialData Double
@def_rw_property Name IRobotMaterialData String
@def_rw_property NU IRobotMaterialData Double
@def_rw_property Nuance IRobotMaterialData String
@def_rw_property PN_Deformation IRobotMaterialData Double
@def_rw_property PN_E_Additional IRobotMaterialData Double
@def_rw_property PN_E_Trans IRobotMaterialData Double
@def_rw_property RE IRobotMaterialData Double
@def_rw_property RE_AxCompr IRobotMaterialData Double
@def_rw_property RE_AxTens IRobotMaterialData Double
@def_rw_property RE_Bending IRobotMaterialData Double
@def_rw_property RE_Shear IRobotMaterialData Double
@def_rw_property RE_TrCompr IRobotMaterialData Double
@def_rw_property RE_TrTens IRobotMaterialData Double
@def_rw_property RO IRobotMaterialData Double
@def_rw_property RT IRobotMaterialData Double
@def_rw_property SecondName IRobotMaterialData Double
@def_rw_property Steel_Thermal IRobotMaterialData Bool
@def_rw_property Timber_Type IRobotMaterialData IRobotMaterialTimberType
@def_rw_property (MaterialType, Type) IRobotMaterialData IRobotMaterialType
@def_rw_property (SectionType, Type) IRobotBarSectionData IRobotBarSectionType
@def_ro_property Nodes IRobotResultsServer IRobotNodeResultServer
@def_ro_property Bars IRobotResultsServer IRobotBarResultServer
@def_ro_property Displacements IRobotNodeResultServer IRobotNodeDisplacementServer
@def_ro_property Stresses IRobotBarResultServer IRobotBarStressServer

macro def_com(name, Type, params...)
  out = last(params)
  params = params[1:end-1]
  name, method = name_method(name)
  param_names = map(param -> param.args[1], params)
  param_types = map(param -> param.args[2], params)
  args = map((name, typ) -> Core.eval(__module__, typ) <: Enum ? :(Int($name)) : name, param_names, param_types)
  quote
    function $(esc(name))(receiver :: $(Type), $(params...)) :: $(esc(out))
        receiver[$(QuoteNode(method))]($(args...))
    end
  end
end

@def_com calculate IRobotCalcEngine Int
@def_com new IRobotProject typ::IRobotProjectType Void
@def_com new IRobotLoadRecordMngr typ::IRobotLoadRecordType Long
@def_com (create_node, Create) IRobotNodeServer node_number::Long x::Double y::Double z::Double Void
@def_com (create_bar, Create) IRobotBarServer bar_number::Long start_node::Long end_node::Long Void
@def_com (create_label, Create) IRobotLabelServer typ::IRobotLabelType name::String IRobotLabel
@def_com is_available IRobotLabelServer typ::IRobotLabelType name::String Boolean
@def_com delete IRobotLabelServer typ::IRobotLabelType name::String Void
@def_com store IRobotLabelServer label::IRobotLabel Void
@def_com (get_node, Get) IRobotNodeServer idx::Int IRobotNode
@def_com (get_bar, Get) IRobotBarServer idx::Int IRobotBar
@def_com (get_record, Get) IRobotLoadRecordMngr idx::Int IRobotLoadRecord
@def_com set_label IRobotNode typ::IRobotLabelType name::String Void
#@def_com set_label IRobotBar typ::IRobotLabelType name::String Void
@def_com (get_selection, Get) IRobotSelectionFactory typ::IRobotObjectType IRobotSelection
@def_com (set_selection_label, SetLabel) IRobotBarServer selection::IRobotSelection typ::IRobotLabelType name::String Void
#@def_com set_value IRobotBarSectionConcreteData attr::IRobotBarSectionConcreteDataValue value::Double Void
@def_com set_value IRobotBarSectionData attr::IRobotBarSectionDataValue value::Double Void
@def_com set_value IRobotBarSectionNonstdData attr::IRobotBarSectionNonstdDataValue value::Double Void
@def_com set_value IRobotLoadRecord value_id::IRobotNodeForceRecordValues value::Double Void
#@def_com set_value IRobotBarSectionSpecialData attr::IRobotBarSectionSpecialDataValue value::Double Void
@def_com CreateNonstd IRobotBarSectionData rel_pos::Double IRobotBarSectionNonstdData
@def_com CalcNonstdGeometry IRobotBarSectionData Void
@def_com create_simple IRobotCaseServer number::Int name::String nature::IRobotCaseNature analize_type::IRobotCaseAnalizeType IRobotSimpleCase
@def_com SaveToDBase IRobotMaterialData Void
@def_com LoadFromDBase IRobotMaterialData name::String Boolean
@def_com add_one IRobotSelection id::Int Void
@def_com (node_displacement, Value) IRobotNodeDisplacementServer node::Int case::Int IRobotDisplacementData
@def_com (bar_displacement, Value) IRobotBarDisplacementServer bar::Int pos::Double case::Int IRobotDisplacementData
@def_com (bar_stress, Value) IRobotBarStressServer bar::Int case::Int pos::Double IRobotBarStressData
@def_com from_text IRobotSelection ids::String Void
#;Stress
@def_rw_property Smin IRobotBarStressData Double
@def_rw_property Smax IRobotBarStressData Double
@def_rw_property Torsion IRobotBarStressData Double

robot_app = nothing
function application()
    global robot_app
    if robot_app == nothing
        robot_app = let r = com[:Dispatch]("Robot.Application")
                        r[:Visible] = 1
                        r
                    end
    else
        robot_app
    end
end

# Labels
new_label(typ, name, fn) =
  let labels = labels(structure(project(application())))
    if is_available(labels, typ, name)
      delete(labels, typ, name)
    end
    label = create_label(labels, typ, name)
    fn(data(label))
    store(labels, label)
    name
  end

# Node support
create_node_support_label(name, ux, uy, uz, rx, ry, rz) =
  new_label(I_LT_NODE_SUPPORT, name,
    support_data -> begin
      UX(support_data, ux ? 1 : 0)
      UY(support_data, uy ? 1 : 0)
      UZ(support_data, uz ? 1 : 0)
      RX(support_data, rx ? 1 : 0)
      RY(support_data, ry ? 1 : 0)
      RZ(support_data, rz ? 1 : 0)
  end)

# Nodes
mutable struct node_support
    name::String
    ux::Bool
    uy::Bool
    uz::Bool
    rx::Bool
    ry::Bool
    rz::Bool
    created::Bool
end

create_node_support(
    name::String;
    ux::Bool=false,
    uy::Bool=false,
    uz::Bool=false,
    rx::Bool=false,
    ry::Bool=false,
    rz::Bool=false) = node_support(name, ux, uy, uz, rx, ry, rz, false)

struct truss_node_data
    id::Int
    loc::Loc
    family::Any
    load::Any
end

node_counter = Parameter(0)
added_nodes = Parameter(Dict())
case_counter = Parameter(0)
bar_counter = Parameter(0)
added_bars = Parameter(Dict())

add_node!(p, family, load=false, reuse=false) =
    begin
        #isreuse && for/or(k(v)(in_hash(added_nodes()))(),
        #                  distance(k, p) < isreuse && k) || #We should check that the families are the same. What about the loads?
        #                   ||
        node_counter(node_counter()+1)
        added_nodes()[p] = truss_node_data(node_counter(), p, family, load) # Should we check for collisions here? (nodes at the same location);
    end

current_nodes_ids() = 0:node_counter()

# Bars

struct truss_bar_data
    id::Int
    node0::truss_node_data
    node1::truss_node_data
    rotation::Double
    family::Any
end

add_bar!(p0, p1, rotation, family) =
    begin
        bar_counter(bar_counter()+1)
        added_bars()[bar_counter()] =
            truss_bar_data(bar_counter(),
                           added_nodes()[p0],
                           added_nodes()[p1],
                           rotation,
                           family)
    end

current_bars_ids() = 0:bar_counter()

new_robot_analysis(process_results, create_truss, v=nothing) =
    with(node_counter, 0,
         added_nodes, Dict(),
         added_bars, Dict(),
         case_counter, 0) do
        create_truss()
        let struc = structure(project(application()))
            nds = nodes(struc)
            brs = bars(struc)
            node_loads = Dict(v==nothing ? [] : [v => map(n -> n.id, values(added_nodes()))])
          for node_data in values(added_nodes())
            let (node_id, p, node_family, node_load) = (node_data.id, node_data.loc, node_data.family, node_data.load)
                create_node(nds, node_id, p.x, p.y, p.z)
                support = node_family.support
                if support != false
                    if ! support.created
                        create_node_support_label(support.name,
                                                  support.ux, support.uy, support.uz,
                                                  support.rx, support.ry, support.rz)
                        support.created = true
                    end
                    set_label(get_node(nds, node_id), I_LT_NODE_SUPPORT, support.name)
                end
                if node_load != 0
                    node_loads[node_load] = [node_id, get_node(node_loads, node_load, [])...]
                end
            end
          end
            family_bars = Dict()
            for bar_data in values(added_bars())
                let (bar_id, node_id0, node_id1, rotation, bar_family) = (bar_data.id, bar_data.node0.id, bar_data.node1.id, bar_data.rotation, bar_data.family)
                    create_bar(brs, bar_id, node_id0, node_id1)
                    if abs(rotation) > 1e-16 #fix this
                      Gamma(get_bar(brs, bar_id), rotation)
                    end
                    family_bars[bar_family] = [bar_id, get(family_bars, bar_family, [])...]
                end
            end
            for (bar_family, bars_ids) in family_bars
                if ! bar_family.created()
                    create_bar_material_label(bar_family.material...)
                    create_bar_tube_section_label(bar_family.section...)
                    bar_family.created(true)
                end
                let selection = get_selection(selections(struc), I_OT_BAR)
                    ids = IOBuffer()
                    for bar_id in bars_ids
                        print(ids, bar_id, " ")
                    end
                    str = String(take!(ids))
                    from_text(selection, str)
                    let (name, material_name, wood, specs) = bar_family.section
                        set_selection_label(brs, selection, I_LT_BAR_SECTION, name)
                    end
                end
            end
            case_counter(case_counter()+1)
            new_case(case_counter(),
                     "Test-$(case_counter())",
                     I_CN_PERMANENT, # I_CN_EXPLOATATION I_CN_WIND I_CN_SNOW I_CN_TEMPERATURE I_CN_ACCIDENTAL I_CN_SEISMIC,
                     I_CAT_STATIC_LINEAR, #I_CAT_STATIC_NONLINEAR I_CAT_STATIC_FLAMBEMENT,
                     records -> new_node_loads(records, node_loads),
                     process_results)
      end
end
#

# Bar release

create_bar_release_label(name, sux, suy, suz, srx, sry, srz, eux, euy, euz, erx, ery, erz) =
  new_label(I_LT_BAR_RELEASE,
            name,
            bar_data -> let (start_node, end_node) = (bar_data.start_node, bar_data.end_node);
                        UX(start_node, sux)
                        UY(start_node, suy)
                        UZ(start_node, suz)
                        RX(start_node, srx)
                        RY(start_node, sry)
                        RZ(start_node, srz)
                        UX(end_node, eux)
                        UY(end_node, euy)
                        UZ(end_node, euz)
                        RX(end_node, erx)
                        RY(end_node, ery)
                        RZ(end_node, erz)
                     end)

set_bar_release!(bar, label) =
  set_label(get_bar(bars(structure(project(application()))), bar),
            I_LT_BAR_RELEASE,
            label)


# Bar material

create_bar_timber_material_label(name, _Type, _Timber_Type, _Name, _Nuance, _E, _NU, _GMean, _RO, _LX, _DumpCoef, _RE_Bending, _RE_AxTens, _RE_TrTens, _RE_AxCompr, _RE_TrCompr, _RE_Shear, _E_5, _E_Trans) =
  new_label(I_LT_BAR_MATERIAL,
            name,
            bar_data -> begin
                        #bar_data["Type"] = Int(_Type)
                        MaterialType(bar_data, _Type)
                        Timber_Type(bar_data, _Timber_Type)
                        Name(bar_data, _Name)
                        Nuance(bar_data, _Nuance)
                        E(bar_data, _E)
                        NU(bar_data, _NU)
                        GMean(bar_data, _GMean)
                        RO(bar_data, _RO)
                        LX(bar_data, _LX)
                        DumpCoef(bar_data, _DumpCoef)
                        RE_Bending(bar_data, _RE_Bending)
                        RE_AxTens(bar_data, _RE_AxTens)
                        RE_TrTens(bar_data, _RE_TrTens)
                        RE_AxCompr(bar_data, _RE_AxCompr)
                        RE_TrCompr(bar_data, _RE_TrCompr)
                        RE_Shear(bar_data, _RE_Shear)
                        E_5(bar_data, _E_5)
                        E_Trans(bar_data, _E_Trans)
                        SaveToDBase(bar_data)
                    end)

create_bar_material_label(name, _Type, _Name, _Nuance, _E, _NU, _Kirchoff, _RO, _LX, _DumpCoef, _RE, _RT) =
  new_label(I_LT_BAR_MATERIAL,
            name,
            bar_data -> begin
                        MaterialType(bar_data, _Type)
                        #bar_data["Type"] = Int(_Type)
                        Name(bar_data, _Name)
                        Nuance(bar_data, _Nuance)
                        E(bar_data, _E)
                        NU(bar_data, _NU)
                        Kirchoff(bar_data, _Kirchoff)
                        RO(bar_data, _RO)
                        LX(bar_data, _LX)
                        DumpCoef(bar_data, _DumpCoef)
                        RE(bar_data, _RE)
                        RT(bar_data, _RT)
                        SaveToDBase(bar_data)
                    end)

# Bar section

create_bar_tube_section_label(name, material_name, iswood, specs) =
  new_label(I_LT_BAR_SECTION,
            name,
            bar_data -> begin
                        #bar_data["Type"] = Int(I_BST_NS_TUBE)
                        SectionType(bar_data, I_BST_NS_TUBE)
                        #bar_data["ShapeType"] = 93
                        shape_type(bar_data, iswood ? I_BSST_WOOD_CIRC : I_BSST_TUBE)
                        MaterialName(bar_data, material_name)
                        for (spec, relative) in zip(specs, division(0.0, 1.0, length(specs)))
                            let (issolid, diameter, thickness) = spec
                                robotBarSectionNonstdData = CreateNonstd(bar_data, relative)
                                set_value(robotBarSectionNonstdData, I_BSNDV_BOX_H, diameter)
                                if ! issolid
                                    set_value(robotBarSectionNonstdData, I_BSNDV_BOX_B, thickness)
                                end
                            end
                        end
                        CalcNonstdGeometry(bar_data)
                    end)

create_bar_rectangle_section_label(name, material_name, iswood, specs) =
  new_label(I_LT_BAR_SECTION,
            name,
            bar_data -> begin
                        SectionType(bar_data, I_BST_NS_RECT)
                        shape_type(bar_data, iswood ? I_BSST_FRTG : I_BSST_RECT)
                        MaterialName(bar_data, material_name)
                        for (spec, relative) in zip(specs, division(0.0, 1.0, length(specs)))
                            let (issolid, width, height, thickness) = spec
                                robotBarSectionNonstdData = CreateNonstd(bar_data, relative)
                                set_value(robotBarSectionNonstdData, I_BSNDV_BOX_H, width)
                                set_value(robotBarSectionNonstdData, I_BSNDV_BOX_B, height)
                                if ! issolid
                                    set_value(robotBarSectionNonstdData, I_BSNDV_BOX_TF, thickness)
                                end
                            end
                        end
                        CalcNonstdGeometry(bar_data)
                    end)


set_bar_section!(bar, label) =
  set_label(get_bar(bars(structure(project(application()))), bar),
            I_LT_BAR_SECTION,
            label)

new_case(number, name, nature, analize_type, setup, process_results) =
  let case = create_simple(cases(structure(project(application()))),
                          number,
                          name,
                          nature,
                          analize_type)
    @time(setup(records(case)))
    @time(calculate(calc_engine(project(application()))))
    process_results(results(structure(project(application()))))
  end

new_node_loads(records, loads) =
  for (vec, ids) in loads
    let idx = new(records, I_LRT_NODE_FORCE)
        record = get_record(records, idx)
        objs = objects(record)
      for node_id in ids
        add_one(objs, node_id)
      end
      set_value(record, I_NFRV_FX, vec.x)
      set_value(record, I_NFRV_FY, vec.y)
      set_value(record, I_NFRV_FZ, vec.z)
    end
  end

node_displacement_vector(results, id, case_id) =
  let d = node_displacement(Displacements(nodes(results)),
                            id,
                            case_id)
    vxyz(UX(d), UY(d), UZ(d))
  end

bar_max_stress(results, id, case_id) =
  Smax(bar_stress(Stresses(bars(results)), id, case_id, 0.0)) ##The position should be changeable



#
struct COM_Backend{K,T} <: Backend{K,T}
  com::Any
end
connection(backend::COM_Backend) = backend.com

abstract type ROBOTKey end
const ROBOTId = Any
const ROBOTIds = Vector{ROBOTId}
const ROBOTRef = GenericRef{ROBOTKey, ROBOTId}
const ROBOTRefs = Vector{ROBOTRef}
const ROBOTEmptyRef = EmptyRef{ROBOTKey, ROBOTId}
const ROBOTUniversalRef = UniversalRef{ROBOTKey, ROBOTId}
const ROBOTNativeRef = NativeRef{ROBOTKey, ROBOTId}
const ROBOTUnionRef = UnionRef{ROBOTKey, ROBOTId}
const ROBOTSubtractionRef = SubtractionRef{ROBOTKey, ROBOTId}
const ROBOT = COM_Backend{ROBOTKey, ROBOTId}

void_ref(b::ROBOT) = ROBOTNativeRef(-1)

create_ROBOT_connection() =
    begin
        application()
    end

const robot = ROBOT(LazyParameter(Any, create_ROBOT_connection))

realize(b::ROBOT, s::TrussNode) =
    add_node!(s.p, s.family)

realize(b::ROBOT, s::TrussBar) =
    add_bar!(s.p0, s.p1, s.angle, s.family)
