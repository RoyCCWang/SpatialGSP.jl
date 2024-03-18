
function visualizebernsteinbasis(
    λ_max::Real,
    L::Integer,
    fig_num::Integer;
    title_string = "Bernstein polynomial basis of degree $L",
    Nq::Integer = 1000,
    show_legend = true)

    evalbern = collect(
        tt->GSP.evalbernsteinpolynomialdirect(i, L, tt/λ_max)
        for i = 0:L
    )

    ts = LinRange(0, λ_max, Nq)
    PLT.figure(fig_num)
    fig_num += 1

    for i = 0:L
        PLT.plot(ts, evalbern[i+1].(ts), label = "$i")
    end

    PLT.title(title_string)

    if show_legend
        PLT.legend()
    end

    return fig_num
end
