#=
We want to support these two syntaxes

a < Number
a < Number("Number")
a < Number("Number", "N")
a < Number("Number", "N", "Parameter N")
a < Number("Number", "N", "Parameter N", 0.0)

a < Number(0.0)
a < Number("Number", 0.0)
a < Number("Number", "N", 0.0)
a < Number("Number", "N", "Parameter N", 0.0)

a < String()
a < String("String")
a < String("String")
a < Number(0.0, "a")

=#

#in_gh(sym) = Expr(:., :GH, QuoteNode(sym))
in_gh(sym) = Symbol("GH$(sym)")
kgh_io_function_names = Symbol[]
is_kgh_io_function_name(sym) =
  sym in kgh_io_function_names

macro ghdef(name, init)
  let str = string(name),
      ghname = esc(in_gh(name))
    quote
      push!(kgh_io_function_names, $(QuoteNode(name)))
      $(ghname)(description::Base.String, short_description=description[1:1], message=description*" parameter", value=$init) =
        [$(str), description, short_description, message, value]
      $(ghname)(value::Base.Any, description=$(str), short_description=description[1:1], message=description*" parameter") =
        [$(str), description, short_description, message, value]
    end
  end
end
@ghdef(String, "")
@ghdef(Path, "")
@ghdef(Boolean, false)
@ghdef(Number, 0.0)
@ghdef(Integer, 0)
@ghdef(Point, u0())
@ghdef(Vector, vx(1))
@ghdef(Any, nothing)
@ghdef(Eval, nothing)
@ghdef(JL, nothing)
@ghdef(Strings, [])
@ghdef(Paths, [])
@ghdef(Booleans, [])
@ghdef(Numbers, [])
@ghdef(Integers, [])
@ghdef(Points, [])
@ghdef(Vectors, [])
@ghdef(Many, [])
@ghdef(Evals, [])
@ghdef(JLs, [])
@ghdef(Stringss, [])
@ghdef(Pathss, [])
@ghdef(Booleanss, [])
@ghdef(Numberss, [])
@ghdef(Integerss, [])
@ghdef(Pointss, [])
@ghdef(Vectorss, [])
@ghdef(Manies, [])
@ghdef(Evalss, [])
@ghdef(JLss, [])

export define_kgh_function

kgh_forms(text, idx=1) =
  let (expr, idx) = Meta.parse(text, idx, greedy=true, depwarn=false)
    isnothing(expr) ?
      [] :
      [expr, kgh_forms(text, idx)...]
  end

is_kgh_io_function_call(e) =
  e isa Expr && e.head === :call && is_kgh_io_function_name(e.args[1])

match_expr(e1, e2) =
  e2 === :_ ||
  e1 == e2 ||
  e1 isa Expr && e2 isa Expr &&
  match_expr(e1.head, e2.head) &&
  length(e1.args) == length(e2.args) &&
  all([match_expr(e1, e2) for (e1, e2) in zip(e1.args, e2.args)])

is_kgh_input(form) =
  match_expr(form, :(_ < _)) &&
  is_kgh_io_function_call(form.args[3])
is_kgh_output(form) =
  (match_expr(form, :(_ > _)) &&
   is_kgh_io_function_call(form.args[3])) ||
  is_kgh_io_function_call(form)

kgh_io_param(form) =
  form.args[2] == :_ ? :__result : form.args[2]
kgh_io_call(form) =
  let param = form.args[2],
      form = form.args[3],
      func = Expr(:., :Khepri, QuoteNode(in_gh(form.args[1])))
    match_expr(form, :(_())) ||
    (match_expr(form, :(_(_))) && !(form.args[2] isa String)) ?
        :($(func)($(form.args[2:end]...), $(string(param)))) :
        :($(func)($(form.args[2:end]...)))
  end

kgh_inputs(forms) = filter(is_kgh_input, forms)
kgh_outputs(forms) = filter(is_kgh_output, forms)

create_kgh_function(name::String, body::String) =
  let forms = kgh_forms(body),
      inputs = filter(is_kgh_input, forms),
      outputs = filter(is_kgh_output, forms),
      forms = filter(f -> !(f in inputs || f in outputs), forms),
      inp_params = map(kgh_io_param, inputs),
      out_params = map(kgh_io_param, outputs),
      inp_forms = map(kgh_io_call, inputs),
      out_forms = map(kgh_io_call, outputs),
      inps = Symbol("__inps_$name"),
      inp_inits = [:($(inp) = $(inps)[$(i)]) for (i, inp) in enumerate(inp_params)],
      outs = Symbol("__outs_$name"),
      out_inits = [:($(outs)[$(i)] = $(out)) for (i, out) in enumerate(out_params)],
      inp_docs = Symbol("__doc_inps_$name"),
      out_docs = Symbol("__doc_outs_$name"),
      shapes = Symbol("__shapes_$name")
    if isempty(inp_forms) && isempty(out_forms)
      quote
        $(inp_docs) = Any[$(inp_forms...)]
        $(out_docs) = Any[$(out_forms...)]
        $(inps) = Array{Any,1}(undef, $(length(inp_inits)))
        $(outs) = Array{Any,1}(undef, $(length(out_inits)))
        $(shapes) = Shape[]
        function $(Symbol("__func_$name"))()
          nothing
        end
      end
    else
      quote
        $(inp_docs) = Any[$(inp_forms...)]
        $(out_docs) = Any[$(out_forms...)]
        $(inps) = Array{Any,1}(undef, $(length(inp_inits)))
        $(outs) = Array{Any,1}(undef, $(length(out_inits)))
        $(shapes) = Shape[]
        function $(Symbol("__func_$name"))()
          $(inp_inits...)
          __result = begin
                $(forms...)
              end
          $(out_inits...)
          nothing
        end
      end
    end
  end

define_kgh_function(name::String, body::String) =
  Base.eval(Main, create_kgh_function(name, body))
#=
fn = """
a < Number("Number", "N", "Parameter N", 0.0)
b < Any()
c < Number(5)
d < String("Foo")

f(a + b)

c = a + 1

c > Any()
d > String("Bar")
"""

create_kgh_function("foo", fn)
=#

#=
fn = """
a < Number("Number", "N", "Parameter N", 0.0)
b < Any()
c < Number(5)
d < String("Foo")
e > String("Bar")
_ > Number()

e = f(a + b)
g(b*c)
"""

create_kgh_function("foo", fn)

fn2 = """
e = f(a + b)
g(b*c)
"""

create_kgh_function("foo", fn2)
using InteractiveUtils
InteractiveUtils.@code_lowered Khepri.foo()
=#
#=
fn = """
a < Number("Number", "N", "Parameter N", 0.0)
b < Any()
c < Number(5)
d < String("Foo")

f(a + b)

c = a + 1

c > Any()
d > String("Bar")
"""

eval(create_kgh_function("foo", fn))
=#

#=
create_kgh_function("foo", raw"""
path < String()
autos < Numbers()
costs < Numbers()

open(path, "w") do f
  for (a, c) in zip(autos, costs)
    println(f, "$a $c")
  end
end
""")

dump(:("foo $bar"))

=#

#=
str1 = """# v < Type("Input", "I", "Parameter I", default)
a < Number(2)
b < Number(3)
_ > Number()

sqrt(a^2 + b^2)"""

create_kgh_function("zzz", str1)
=#
