using GadgetIO, Test, DelimitedFiles

@testset "Objects" begin
    @test_nowarn Header()
    @test_nowarn Info_Line()
end

@testset "Read Snapshot" begin

    snap_file = joinpath(dirname(@__FILE__), "snap_050")

    @test_nowarn read_snap(snap_file)
    @test_nowarn read_header(snap_file)
    @test_nowarn read_info(snap_file)

    d = read_snap(snap_file, "POS", 0)

    ideal_file = joinpath(dirname(@__FILE__), "pos.dat")
    d_ideal = Float32.(readdlm(ideal_file))

    @test d == d_ideal
end

@testset "Write Snapshot" begin

    # read in reference file
    ref_file = joinpath(dirname(@__FILE__), "snap_050")
    head = head_to_obj(ref_file)
    x = read_snap(ref_file, "POS", 0)

    # specify output file for testing
    output_file = joinpath(dirname(@__FILE__), "write_test.dat")
    f = open(output_file, "w")

    @test_nowarn write_header(f, head)
    @test_nowarn write_block(f, x, "POS")
    @test_nowarn write_block(f, x, "", snap_format=1)

    # @test_warn "Please specify blockname!" write_block(f, x, "")

    close(f)
end