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

    # check if read to dict works
    @test_nowarn d = read_snap(snap_file)

    @test d["PartType0"]["POS"] == d_ideal
end

@testset "Snapshot utility" begin
    
    ref_file = joinpath(dirname(@__FILE__), "snap_050")

    @test_nowarn print_blocks(ref_file)

    # Check if a block is present
    present = block_present(ref_file, "POS")
    @test present == true
    present = block_present(ref_file, "BHMA")
    @test present == false

    blocks = ["POS", "VEL", "MASS"]
    blocks_checked, no_mass_block = GadgetIO.check_blocks(ref_file, blocks)

    @test blocks_checked == ["POS", "VEL"]
    @test no_mass_block == true

    @test_nowarn GadgetIO.get_block_positions(ref_file)

    filename = GadgetIO.select_file(ref_file, 0)

    @test filename == ref_file
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
    close(f)

    pos_info = Info_Line("POS", Float32, 3, [1, 0, 0, 0, 0, 0])
    x_check = read_block_by_name(output_file, "POS", info=pos_info, parttype=0)

    # check if we read the same thing we wrote
    @test x_check == x

    f = open(output_file, "w")
    @test_nowarn write_block(f, x, "", snap_format=1)

    @test_throws ErrorException("Please specify blockname!") write_block(f, x, "")
    close(f)

    rm(output_file)
end

x = 0