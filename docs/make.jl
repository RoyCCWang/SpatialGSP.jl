using Documenter
using SpatialGSP

# # local.
makedocs(
    sitename = "SpatialGSP",
    modules = [SpatialGSP],
    #format = Documenter.HTML(),
    pages = [
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: Bernstein filtering" =>
        "generated/bernstein_filtering_lit.md",
        "Demo: Axis graph" =>
        "generated/axis_lit.md",
    ],
)

# # github.
# makedocs(
#     sitename="SpatialGSP.jl",
#     modules=[SpatialGSP],
#     format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
#     pages = [
#         "Overview" => "index.md",
#         "Public API" => "api.md",
#         "Demo: chaining & remedy" => "generated/chaining.md",
#     ],
# )
# deploydocs(
#     repo = "github.com/RWAlgorithms/SpatialGSP.jl.git",
#     target = "build",
#     branch = "gh-pages",
#     versions = ["stable" => "v^", "v#.#" ],
# )