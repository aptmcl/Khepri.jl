export robot,
       project_kind,
       new_project!,
       new_3d_project!,
       new_robot_analysis,
       node_displacement,
       node_displacement_vector,
       displacements,
       nodes,
       added_nodes,
       added_bars,
       robot_truss_bar_family


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

@enums IRobotDeadRecordValues begin
  I_DRV_X = 0
  I_DRV_Y = 1
  I_DRV_Z = 2
  I_DRV_ENTIRE_STRUCTURE = 15
  I_DRV_COEFF = 3
end

@enums IRobotComponentType begin
  I_CT_VALUES_ARRAY = 1
  I_CT_NUMBERS_ARRAY = 2
  I_CT_NAMES_ARRAY = 3
  I_CT_OBJECTS_ARRAY = 4
  I_CT_GEO_POINT_2D = 5
  I_CT_GEO_POINT_3D = 6
  I_CT_GEO_CURVE_DIV = 7
  I_CT_GEO_SEGMENT = 8
  I_CT_GEO_SEGMENT_LINE = 9
  I_CT_GEO_SEGMENT_ARC = 10
  I_CT_GEO_OBJECT = 11
  I_CT_GEO_POLYLINE = 12
  I_CT_GEO_CONTOUR = 13
  I_CT_GEO_ARC = 14
  I_CT_GEO_CIRCLE = 15
  I_CT_GEO_SEGMENT_COLLECTION = 16
  I_CT_GEO_POINT_3D_COLLECTION = 17
  I_CT_GEO_LAYER = 18
  I_CT_EMITTER = 19
  I_CT_JOINT_LOAD = 20
  I_CT_JOINT_KNEE_LOAD = 20
  I_CT_JOINT_ANGLE_LOAD = 20
  I_CT_JOINT_TUBE_LOAD = 20
  I_CT_JOINT_PINNED_LOAD = 20
  I_CT_JOINT_FIXED_LOAD = 20
  I_CT_JOINT_GUSSET_SIMPLE_LOAD = 20
  I_CT_JOINT_GUSSET_CROSS_LOAD = 20
  I_CT_JOINT_GUSSET_FLANGE_LOAD = 20
  I_CT_MODIF_EXTRUSION = 29
  I_CT_MODIF_LATHE = 30
  I_CT_MODIF_PYRAMID = 31
  I_CT_OPER_TRANSLATION = 32
  I_CT_OPER_SCALING = 33
  I_CT_OPER_ROTATION = 34
  I_CT_OPER_MESHING = 35
  I_CT_SPECTRAL_ANALYSIS_POINTS_COLLECTION = 36
  I_CT_SPECTRAL_ANALYSIS_SPECTRUM = 37
  I_CT_CASE_ANALYSIS_MODES_FILTER = 38
  I_CT_FE_RESULT_PARAMS = 39
  I_CT_RTF_VIEW = 40
  I_CT_POINTS_ARRAY = 41
  I_CT_EXTREME_PARAMS = 42
  I_CT_STRUCTURE_GEO_ANALYSER = 43
  I_CT_FE_EXTREME_PARAMS = 44
  I_CT_FE_MULTI_RESULT_TYPE = 45
  I_CT_SW_STRUCT3D_GEN_PARAMS = 46
  I_CT_SW_STRUCT3D_PURLIN_GEN_PARAMS = 47
  I_CT_SW_STRUCT3D = 48
  I_CT_SW_STRUCT3D_FRAME = 49
  I_CT_SW_STRUCT3D_ELEMENT = 50
  I_CT_HTML_VIEW = 52
  I_CT_STRUCTURE_MERGE_DATA = 53
  I_CT_NUMBERS_DICTIONARY = 54
  I_CT_TABLE_SCREEN_CAPTURE_PARAMS = 55
  I_CT_PARAM_DEF = 56
  I_CT_RESULT_QUERY_PARAMS = 57
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
    esc(quote
          $(name)(in :: $(InType)) :: $(OutType) = $(OutType)(getproperty(in, $(QuoteNode(method))))
          $(name)(in :: $(InType), val :: $(OutType)) :: Nothing = begin
            setproperty!(in, $(QuoteNode(method)), $(Core.eval(__module__, OutType) <: Enum ? :(Int(val)) : :(val)))
            nothing
          end
      end)
end

macro def_ro_property(name, InType, OutType)
    name, method = name_method(name)
    esc(:($(name)(in :: $(InType)) :: $(OutType) = $(OutType)(getproperty(in, $(QuoteNode(method))))))
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
@def_com_type IRobotComponentFactory
@def_com_type IRobotSelection
@def_com_type IRobotDisplacementData
@def_com_type IRobotGeoSegmentLine
@def_com_type RobotGeoSegmentLine
@def_com_type IRobotPointsArray
@def_com_type IRobotObjObjectServer
@def_com_type IRobotObjObject
@def_com_type IRobotObjPartMain
@def_com_type IRobotObjAttributes
@def_com_type IRobotGeoObject
@def_com_type IRobotGeoContour
#Generic types that does not actually exists in the Robot API
@def_com_type IRobotComponent


Double = Float64
Long = Int64
Void = Nothing
Boolean = Bool

@def_ro_property project IRobotApplication IRobotProject
@def_ro_property CmpntFactory IRobotApplication IRobotComponentFactory
@def_ro_property nodes IRobotStructure IRobotNodeServer
@def_ro_property bars IRobotStructure IRobotBarServer

@def_ro_property structure IRobotProject IRobotStructure
@def_ro_property calc_engine IRobotProject IRobotCalcEngine
@def_ro_property results IRobotStructure IRobotResultsServer
@def_ro_property labels IRobotStructure IRobotLabelServer
@def_ro_property cases IRobotStructure IRobotCaseServer
#ERROR @def_ro_property selections IRobotStructure IRobotSelection
@def_ro_property selections IRobotStructure IRobotSelectionFactory
@def_rw_property results_freeze IRobotStructure Boolean
@def_ro_property records IRobotSimpleCase IRobotLoadRecordMngr
@def_ro_property objects IRobotLoadRecord IRobotSelection
@def_ro_property data IRobotLabel IRobotNodeSuportData
#@def_ro_property data IRobotLabel IRobotBarSectionData
@def_rw_property UX IRobotNodeSuportData Real
@def_rw_property UY IRobotNodeSuportData Real
@def_rw_property UZ IRobotNodeSuportData Real
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
@def_ro_property displacements IRobotNodeResultServer IRobotNodeDisplacementServer
@def_ro_property Stresses IRobotBarResultServer IRobotBarStressServer
@def_ro_property (MainPart, Main) IRobotObjObject IRobotObjPartMain
@def_rw_property Geometry IRobotObjPartMain IRobotGeoObject
@def_ro_property Attribs IRobotObjPartMain IRobotObjAttributes
@def_rw_property Meshed IRobotObjAttributes Boolean
#@def_rw_property P1 IRobotGeoSegmentLine


macro def_com(name, Type, params...)
  out = last(params)
  params = params[1:end-1]
  name, method = name_method(name)
  param_names = map(param -> param.args[1], params)
  param_types = map(param -> param.args[2], params)
  args = map(param_names, param_types) do name, typ
           Core.eval(__module__, typ) <: Enum ?
             :(Int($name)) :
             typ == :Double ? :(float($name)) : name
         end
  types = map(param_types) do t
            t == :Double ? :Real : t
          end
  params = map(param_names, types) do name, type
             :($(name) :: $(type))
           end
  esc(
    quote
      function $name(receiver :: $(Type), $(params...)) :: $(out)
          getproperty(receiver, $(QuoteNode(method)))($(args...))
      end
  end)
end

@def_com calculate IRobotCalcEngine Int
@def_com new IRobotProject typ::IRobotProjectType Void
@def_com new IRobotLoadRecordMngr typ::IRobotLoadRecordType Long
@def_com (close_project, Close) IRobotProject Void
@def_com (create_node, Create) IRobotNodeServer node_number::Long x::Double y::Double z::Double Void
@def_com (create_bar, Create) IRobotBarServer bar_number::Long start_node::Long end_node::Long Void
@def_com (create_label, Create) IRobotLabelServer typ::IRobotLabelType name::String IRobotLabel
@def_com (create_component, Create) IRobotComponentFactory typ::IRobotComponentType IRobotComponent
@def_com is_available IRobotLabelServer typ::IRobotLabelType name::String Boolean
@def_com (delete_label, Delete) IRobotLabelServer typ::IRobotLabelType name::String Void
@def_com store IRobotLabelServer label::IRobotLabel Void
@def_com (get_node, Get) IRobotNodeServer idx::Int IRobotNode
@def_com (get_bar, Get) IRobotBarServer idx::Int IRobotBar
@def_com (get_record, Get) IRobotLoadRecordMngr idx::Int IRobotLoadRecord
@def_com (get_contour, Get) IRobotObjObjectServer idx::Int IRobotGeoContour
@def_com set_label IRobotNode typ::IRobotLabelType name::String Void
#@def_com set_label IRobotBar typ::IRobotLabelType name::String Void
@def_com (get_selection, Get) IRobotSelectionFactory typ::IRobotObjectType IRobotSelection
@def_com (set_selection_label, SetLabel) IRobotBarServer selection::IRobotSelection typ::IRobotLabelType name::String Void
#@def_com set_value IRobotBarSectionConcreteData attr::IRobotBarSectionConcreteDataValue value::Double Void
@def_com set_value IRobotBarSectionData attr::IRobotBarSectionDataValue value::Double Void
@def_com set_value IRobotBarSectionNonstdData attr::IRobotBarSectionNonstdDataValue value::Double Void
@def_com set_value IRobotLoadRecord value_id::IRobotNodeForceRecordValues value::Double Void
@def_com set_value IRobotLoadRecord value_id::IRobotDeadRecordValues value::Double Void
#@def_com set_value IRobotBarSectionSpecialData attr::IRobotBarSectionSpecialDataValue value::Double Void
@def_com set_size IRobotPointsArray size::Long Void
#@def_com (get_position, Get) IRobotPointsArray idx::Long x::Double y::Double z::Double Void
@def_com (set_position, Set) IRobotPointsArray idx::Long x::Double y::Double z::Double Void

@def_com CreateNonstd IRobotBarSectionData rel_pos::Double IRobotBarSectionNonstdData
@def_com CalcNonstdGeometry IRobotBarSectionData Void
@def_com create_simple IRobotCaseServer number::Int name::String nature::IRobotCaseNature analize_type::IRobotCaseAnalizeType IRobotSimpleCase
@def_com (delete_case, Delete) IRobotCaseServer number::Int Void
@def_com (get_case, Get) IRobotCaseServer idx::Int IRobotSimpleCase
@def_com SaveToDBase IRobotMaterialData Void
@def_com LoadFromDBase IRobotMaterialData name::String Boolean
@def_com add_one IRobotSelection id::Int Void
@def_com (node_displacement, Value) IRobotNodeDisplacementServer node::Int case::IRobotLoadRecordType IRobotDisplacementData
@def_com (bar_displacement, Value) IRobotBarDisplacementServer bar::Int pos::Double case::Int IRobotDisplacementData
@def_com (bar_stress, Value) IRobotBarStressServer bar::Int case::Int pos::Double IRobotBarStressData
@def_com from_text IRobotSelection ids::String Void

#Again, these operate on generic types
@def_com initialize IRobotComponent Void
@def_com update IRobotComponent Void

#@def_com create_arc IRobotObjObjectServer number::Long pts::IRobotPointsArray creation_type::IRobotGeoArcDefinitionMethod=I_GADM_CENTER_BEGIN_END Void
#@def_com create_circle IRobotObjObjectServer number::Long pts::IRobotPointsArray Void
#@def_com CreateCone
@def_com create_contour IRobotObjObjectServer number::Long pts::IRobotPointsArray Void
@def_com (create_object, Create) IRobotObjObjectServer number::Long IRobotObjObject
#@def_com CreateCube
#@def_com CreateCylinder
#@def_comCreateOnFiniteElems
#@def_com CreatePolyline
#@def_com CreateSolid


#;Stress
@def_rw_property Smin IRobotBarStressData Double
@def_rw_property Smax IRobotBarStressData Double
@def_rw_property Torsion IRobotBarStressData Double

robot_app = nothing
function application()
    global robot_app
    if robot_app == nothing
        copy!(com, pyimport("win32com.client"))
        robot_app = let r = com.Dispatch("Robot.Application")
                        r.Visible = 1
                        r.Interactive = 1
                        r
                    end
    else
        robot_app
    end
end

new_project!(typ::IRobotProjectType) =
  let proj = project(application())
    close_project(proj)
    new(proj, typ)
  end

# Labels
new_label(typ, name, fn) =
  let labels = labels(structure(project(application())))
    if is_available(labels, typ, name)
      delete_label(labels, typ, name)
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

truss_node_support_label(s) =
  let k(v) = v ? "x" : "f"
    "Support$(k(s.ux))$(k(s.uy))$(k(s.uz))$(k(s.rx))$(k(s.ry))$(k(s.rz))"
  end

create_node_support(s) =
  create_node_support_label(truss_node_support_label(s),
                            s.ux, s.uy, s.uz,
                            s.rx, s.ry, s.rz)

# Cladding

struct cladding_data
    id::Int
    pts::Locs
    family::Any
end

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

create_cladding(id, pts) =
  let pts = in_world.(pts),
      rpts = create_component(CmpntFactory(application()), I_CT_POINTS_ARRAY),
      objserver = objects(structure(project(application())))
    println("Creating cladding for $(id) and points $(pts)")
    set_size(rpts, length(pts))
    for i in 1:length(pts)
        set_position(rpts, i, pts[i].x, pts[i].y, pts[i].z)
    end
    create_contour(objserver, id, rpts)
    let contour = get_contour(objserver, id) #,
#        m = Main(contour)
#      Geometry(m, contour)
#      Meshed(Attribs(m), false)
      initialize(contour)
      MainPart(contour)
      Attribs(MainPart(contour))
      Meshed(Attribs(MainPart(contour)))
      set_label(contour, I_LT_CLADDING, "Two-way")
      update(contour)
      contour
    end
  end

#=
mypts = [xyz(-4.0,2.0,0.0), xyz(1.0,-3.0,0.0), xyz(1.0,2.0,4.0)]
rpts = create_component(CmpntFactory(application()), I_CT_POINTS_ARRAY)
objserver = objects(structure(project(application())))
set_size(rpts, length(mypts))
for i in 1:length(mypts)
    set_position(rpts, i, mypts[i].x, mypts[i].y, mypts[i].z)
end
create_contour(objserver, 1234, rpts)
contour = get_contour(objserver, 1234)

Main(contour)
Attribs(Main(contour))
Meshed(Attribs(Main(contour)))
initialize(contour)
set_label(contour, I_LT_CLADDING, "Two-way")
update(contour)
=#

analyze_case(number, name, nature, analize_type, setup) =
  let project = project(application()),
      structure = structure(project),
      cases = cases(structure)
    results_freeze(structure, false)
    delete_case(cases, 1)
    let case = create_simple(cases,
                             1, "KhepriTest", #number,
                             #name,
                             nature,
                             analize_type)
      @time(setup(records(case)))
      @time(calculate(calc_engine(project)))
      results(structure)
    end
  end

new_node_loads(records, loads, self_weight) =
  begin
    for (vec, ids) in loads
      let idx = new(records, I_LRT_NODE_FORCE),
          record = get_record(records, idx),
          objs = objects(record)
        for node_id in ids
          add_one(objs, node_id)
        end
        set_value(record, I_NFRV_FX, vec.x)
        set_value(record, I_NFRV_FY, vec.y)
        set_value(record, I_NFRV_FZ, vec.z)
      end
    end
    if self_weight != false
      let idx = new(records, I_LRT_DEAD),
          record = get_record(records, idx)
        set_value(record, I_DRV_COEFF, 1.0)
        set_value(record, I_DRV_ENTIRE_STRUCTURE, 1.0)
        set_value(record, I_DRV_Z, -1.0)
      end
    end
  end

node_displacement_vector(displacements, id, case_id) =
  let d = node_displacement(displacements, id, case_id)
    vxyz(UX(d), UY(d), UZ(d))
  end

bar_max_stress(results, id, case_id) =
  Smax(bar_stress(Stresses(bars(results)), id, case_id, 0.0)) ##The position should be changeable

########################################################################
const project_kind = Parameter(I_PT_SHELL)
create_ROBOT_connection() = new_project!(project_kind())

#
Base.@kwdef struct RobotBackend{K,T} <: LazyBackend{K,T}
  realized::Parameter{Bool}=Parameter(false)
  com::Any=LazyParameter(Any, create_ROBOT_connection)
  truss_nodes::Vector{<:TrussNode}=TrussNode[]
  truss_bars::Vector{<:TrussBar}=TrussBar[]
  shapes::Vector{<:Shape}=Shape[] # This contains all the rest that is not treated yet
  truss_node_data::Vector{TrussNodeData}=TrussNodeData[]
  truss_bar_data::Vector{TrussBarData}=TrussBarData[]
end

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
const ROBOT = RobotBackend{ROBOTKey, ROBOTId}

connection(b::ROBOT) = b.com
void_ref(b::ROBOT) = ROBOTNativeRef(-1)

const robot = ROBOT()

# Robot needs to merge nodes and bars
save_shape!(b::ROBOT, s::TrussNode) = maybe_merged_node(b, s)
save_shape!(b::ROBOT, s::TrussBar) = maybe_merged_bar(b, s)

# Robot does not need layers
with_family_in_layer(f::Function, backend::ROBOT, family::Family) = f()

# Robot Families
abstract type RobotFamily <: Family end

struct RobotTrussBarFamily <: RobotFamily
    section::Any
    material::Any
    ref::Parameter{Any}
end

robot_truss_bar_family(;section, material) = RobotTrussBarFamily(section, material, Parameter{Any}(nothing))
backend_get_family_ref(b::ROBOT, f::Family, rf::RobotTrussBarFamily) =
  begin
    create_bar_material_label(rf.material...)
    create_bar_tube_section_label(rf.section...)
    rf
  end

truss_bar_family_cross_section_area(f::RobotTrussBarFamily) =
  let section_data = f.section[4][1],
      rₒ = section_data[2]/2,
      rᵢ = rₒ - section_data[3]
    annulus_area(rₒ, rᵢ)
  end

set_backend_family(
  default_truss_bar_family(),
  robot,
  robot_truss_bar_family(
  material=[
    "ElasticIsotropic",   # name
    Khepri.I_MT_STEEL,    # Type
    "Steel",              # Name
    "I'm really steel",   # Nuance
    210000000000.0,       # E (Young's modulus)
    0.3,                  # NU (Poisson's ratio)
    81000000000.0,        # G (Kirchoff's or Shear modulus)
    77010.0,              # RO (Density)
    1.2e-05,              # LX (Thermal expansion)
    0.04,                 # DUMPCOEF (Damping coefficient)
    235000000.0,          # RE (Design resistence)
    360000000.0],         # RT (Limit tension resistance)
  section=[
    "Tube",               #name
    "ElasticIsotropic",   #material_name
    false,                #iswood
    [(true,               #solid?
      0.0213,             #diameter
      0.0026)]]))         #thickness

delete_all_shapes(b::ROBOT) =
  begin
    empty!(b.truss_nodes)
    empty!(b.truss_bars)
    b.realized(false)
    new_project!(project_kind())
    b
  end

realize(b::ROBOT, s::TrussNode) =
  error("BUM")

realize(b::ROBOT, s::TrussBar) =
  error("BUM")

realize(b::ROBOT, s::Panel) =
  error("BUM")

##################################################################
# Robot analysis
struct RobotAnalysis <: StructuralAnalysis
  v::Vec
  vs::Vecs
end

robot_analysis(; node_load, nodes_load) =
  RobotAnalysis(node_load, nodes_load)

ensure_realized_structure(b::ROBOT) =
  if ! b.realized()
    realize_structure(b)
    b.realized(true)
  end

realize_structure(b::ROBOT) =
  let ns = process_nodes(b.truss_nodes),
      bs = process_bars(b.truss_bars, ns),
      struc = structure(project(application())),
      nds = nodes(struc),
      brs = bars(struc),
      supports = unique(map(n -> n.family.support, filter(truss_node_is_supported, ns))),
      family_bars = Dict()
    empty!(b.truss_node_data)
    append!(b.truss_node_data, ns)
    empty!(b.truss_bar_data)
    append!(b.truss_bar_data, bs)
    for support in supports
      create_node_support(support)
    end
    for node in ns
      let p = node.loc,
          support = node.family.support
        create_node(nds, node.id, p.x, p.y, p.z)
        if truss_node_is_supported(node)
          set_label(get_node(nds, node.id), I_LT_NODE_SUPPORT, truss_node_support_label(support))
        end
      end
    end
    for bar in bs
      create_bar(brs, bar.id, bar.node1.id, bar.node2.id)
      if abs(bar.rotation) > 1e-16 #fix this
        Gamma(get_bar(brs, bar.id), bar.rotation)
      end
      family_bars[bar.family] = [bar.id, get(family_bars, bar.family, [])...]
    end
    for (bar_family, bars_ids) in family_bars
      let robot_bar_family = realize(b, bar_family),
          selection = get_selection(selections(struc), I_OT_BAR),
          ids = IOBuffer()
        for bar_id in bars_ids
          print(ids, bar_id, " ")
        end
        str = String(take!(ids))
        from_text(selection, str)
        let (name, material_name, wood, specs) = robot_bar_family.section
          set_selection_label(brs, selection, I_LT_BAR_SECTION, name)
        end
      end
    end
  end

case_counter = Parameter(0)

new_robot_analysis(v=nothing; self_weight=false, backend=robot) =
  let node_loads = Dict(v==nothing ?
                    [] :
                    [v => findall(! truss_node_is_supported, backend.truss_nodes)])
    ensure_realized_structure(backend)
    case_counter(case_counter()+1)
    analyze_case(case_counter(),
                 "KhepriTest-$(case_counter())",
                 I_CN_PERMANENT, # I_CN_EXPLOATATION I_CN_WIND I_CN_SNOW I_CN_TEMPERATURE I_CN_ACCIDENTAL I_CN_SEISMIC,
                 I_CAT_STATIC_LINEAR, #I_CAT_STATIC_NONLINEAR I_CAT_STATIC_FLAMBEMENT,
                 records -> new_node_loads(records, node_loads, self_weight))
  end

#
backend_truss_analysis(b::ROBOT, load::Vec) =
  new_robot_analysis(load, backend=b)

node_displacement_function(b::ROBOT, results) =
  let disps = displacements(nodes(results))
    n -> node_displacement_vector(disps, n.id, I_LRT_NODE_DISPLACEMENT)
  end
#=
Change view, save image
https://forums.autodesk.com/t5/robot-structural-analysis-forum/api-command-vba-projection-capture-of-a-view/m-p/3188566/highlight/true#M937

Public Sub CreeVueRM()
' définition d'une vue

robapp.Interactive = True
robapp.Visible = True
robapp.Window.Activate


robapp.Project.CalcEngine.Calculate ' calculate model



Dim mavueRobot As IRobotView3 ' this is important to set IRobotView3 if you want to make screen capture of this view

Set mavueRobot = robapp.Project.ViewMngr.GetView(1) ' it seems CreateView makes this strange affect, use GetView instead

mavueRobot.Selection.Get(I_OT_CASE).FromText ("3") ' selecting case for results display
mavueRobot.Redraw (True)

mavueRobot.Projection = I_VP_3DXYZ

' displaying map
mavueRobot.ParamsFeMap.CurrentResult = I_VFMRT_GLOBAL_DISPLACEMENT_Z

mavueRobot.Visible = True
mavueRobot.Redraw (True)

robapp.Project.ViewMngr.Refresh

' making screen capture
Dim ScPar As RobotViewScreenCaptureParams
Set ScPar = robapp.CmpntFactory.Create(I_CT_VIEW_SCREEN_CAPTURE_PARAMS)

ScPar.Name = "My screen capture"
ScPar.UpdateType = I_SCUT_UPDATED_UPON_PRINTING
mavueRobot.MakeScreenCapture ScPar



robapp.Project.PrintEngine.SaveReportToOrganizer


End Sub




Dim ScPar As RobotViewScreenCaptureParams
Set ScPar = robapp.CmpntFactory.Create(I_CT_VIEW_SCREEN_CAPTURE_PARAMS)



ScPar.Name = "My screen capture"
ScPar.UpdateType = I_SCUT_CURRENT_VIEW
ScPar.Resolution = I_VSCR_4096
mavueRobot.MakeScreenCapture ScPar



robapp.Project.PrintEngine.AddScToReport "My screen capture"

robapp.Project.PrintEngine.SaveReportToOrganizer



' saving printout \ report to file

robapp.Project.PrintEngine.SaveReportToFile "c:\Autodesk\output\rr.rtf", I_OFF_RTF_JPEG



'or directly opening printout in Word

robapp.Project.PrintEngine.ExternalPreviewReport EPF_MS_OFFICE

=#
