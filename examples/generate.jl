using Literate

include("helpers/gen_utils.jl")

# `dest_dir` is where the generated files will end up in. We delete all the files in that directory first.
dest_dir = "../docs/src/generated"
reset_dir(dest_dir)

# # Bernstein filtering
# fix the URL. This is generated because we're using Documenter.jl-flavoured Markdown.
postprocfunc = xx->replace(
    xx,
    "EditURL = \"bernstein_filtering.jl\"" =>
    "EditURL = \"../../../examples/bernstein_filtering.jl\"" # as ifthe pwd() is in the `dest_dir`
)

Literate.markdown(
    "bernstein_filtering.jl";
    execute = true,
    name = "bernstein_filtering_lit", # make this different ahn "bernstein_filter.jl" so it is easier to find and delete all generated files.
    postprocess = postprocfunc,
)

move_prefix_name = "bernstein_filtering_lit"
move_md(dest_dir, move_prefix_name)


# # Axis search graph construction
postprocfunc = xx->replace(
    xx,
    "EditURL = \"axis.jl\"" =>
    "EditURL = \"../../../examples/axis.jl\"" # as ifthe pwd() is in the `dest_dir`
)

Literate.markdown(
    "axis.jl";
    execute = true,
    name = "axis_lit", # make this different ahn "bernstein_filter.jl" so it is easier to find and delete all generated files.
    postprocess = postprocfunc,
)

move_prefix_name = "axis_lit"
move_md(dest_dir, move_prefix_name)

nothing

