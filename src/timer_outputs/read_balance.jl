using DelimitedFiles
using Printf 

"""
    parse_balance(filename)

Reads the `balance.txt` log file and returns a tuple of `(step_number, timing, active)`.
"""
function parse_balance(filename::String)

    # read the balance file
    f = open(filename)
    lines = readlines(f)
    close(f)

    # select all lines that contain step information
    sel = findall(occursin.("Step", lines))

    # split the lines at equal signs
    split_lines = split.(lines[sel], '=')

    # read step counter 
    steps = [parse(Int64, split_lines[i][2][1:9]) for i = 1:length(split_lines)]

    # find unique steps in case run has been restarted
    unique_steps = unique(steps)

    # find positions where unique steps are in total steps
    sel = get_index_list(unique_steps, steps)

    # get time of step
    timing = [parse(Float64, split_lines[i][3][1:10]) for i in sel]

    # get number of active particles
    active = [parse(Int64, split_lines[i][4][1:11]) for i in sel]

    return steps, timing, active
end

"""
    print_performance(filename)

Basic info on performance from `balance.txt` file.
"""
function print_performance(filename::String)

    steps, timing, active = parse_balance(filename)

    println("\nStatistics: \n")
    println("Total steps:            ", steps[end])
    println("Total time:             ", sum(timing) / 3600.0, " h")
    println("Longest time for step:  ", maximum(timing) / 3600.0, " h")
    println("Mean time for step:     ", mean(timing), " s")
    println("Mean active particles:  ", mean(active))
    println("Mean time per particle: ", mean(timing ./ active), " s")
    println()

end

