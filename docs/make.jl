using Khepri
using Documenter

makedocs(;
    modules=[Khepri],
    authors="António Menezes Leitão <antonio.menezes.leitao@gmail.com>",
    repo="https://github.com/aptmcl/Khepri.jl/blob/{commit}{path}#L{line}",
    sitename="Khepri.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aptmcl.github.io/Khepri.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aptmcl/Khepri.jl",
    devbranch="master",
)
