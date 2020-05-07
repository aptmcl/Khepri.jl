
families_registry = Dict{DataType, Vector{Family}}()

available_families(element::BIMShape) = available_families(element.family)
available_families(family::Family) = available_families(typeof(family))
available_families(family_type::DataType) = families_registry[family_type]

register_family(family::Family) = push!(get!(families_registry, typeof(family), Family[]), family)
export available_families, register_family

#=
We define standard family types
=#

for thickness in [11, 15, 30]
  register_family(
    wall_family_element(
      default_wall_family(),
      name="Wall_$(thickness)",
      thickness=thickness))
end

for thickness in [20, 22, 25, 30]
  register_family(
    slab_family_element(
      default_slab_family(),
      name="Slab_$(thickness)",
      thickness=thickness))
end
