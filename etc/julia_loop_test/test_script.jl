using ArgParse
using Distributed

@everywhere include("/global/home/users/drewhart/test_fn.jl")

s = ArgParseSettings()

@add_arg_table! s begin
    "-x"
	help = "this is x"
	arg_type = Int
	default = 3
end

function main()
    args = parse_args(s)
    x = args["x"]

    nw = nworkers()

    pmap(sleepy_fn, repeat([x], nw))

end

main()

