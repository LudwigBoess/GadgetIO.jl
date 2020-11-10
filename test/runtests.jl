using GadgetIO, Test, DelimitedFiles, HTTP


@testset "GadgetIO" begin

function filter_dummy(filename::String)
    mtop = read_subfind(filename, "MTOP")
    return findall(mtop .> 20.0)
end

@info "downloading test data..."
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_sedov", "./snap_sedov")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/pos_sedov.dat", "./pos_sedov.dat")

download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.0", "./sub_002.0")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.1", "./sub_002.1")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.2", "./sub_002.2")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.3", "./sub_002.3")

download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0", "./snap_002.0")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1", "./snap_002.1")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2", "./snap_002.2")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3", "./snap_002.3")

download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0.key", "./snap_002.0.key")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1.key", "./snap_002.1.key")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2.key", "./snap_002.2.key")
download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3.key", "./snap_002.3.key")

download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.key.index", "./snap_002.key.index")

@info "done!"

    @testset "Objects" begin
        @test_nowarn SnapshotHeader()
        @test_nowarn InfoLine()
    end

    @testset "Read Snapshot" begin

        snap_file = joinpath(dirname(@__FILE__), "snap_sedov")

        @testset "Read blocks" begin
            @test_nowarn read_snap(snap_file)
            @test_nowarn read_header(snap_file)
            @test_nowarn read_info(snap_file)

            d = read_snap(snap_file, "POS", 0)

            ideal_file = joinpath(dirname(@__FILE__), "pos_sedov.dat")
            d_ideal = Float32.(readdlm(ideal_file))

            @test d == d_ideal

            # check if read to dict works
            @test_nowarn d = read_snap(snap_file)
            @test d["PartType0"]["POS"] == d_ideal
        end

        @testset "Read particles in box" begin
            center = [6382.1562, -1780.875, -7379.2812]
            rvir   = 194.34073
            pos = read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=false, parttype=1)

            @test pos[1,:] ≈ [6202.94, -1771.11, -7363.77]
        end

        @testset "Read particles in halo" begin
            center = [6382.1562, -1780.875, -7379.2812]
            rvir   = 194.34073
            #read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=false)
        end

    end

    @testset "Read subfind" begin

        subbase = "sub_002"
        subfile = "./sub_002.0"

        @testset "standard read" begin
            # check if standard reading works
            mtop = read_subfind(subbase * ".0", "MTOP")
            @test mtop ≈ Float32[79.19597, 50.597424, 41.497887, 30.298456, 25.498701, 24.998726, 25.69869, 16.29917, 19.599, 19.499006, 16.599154, 15.49921, 16.09918, 13.199327, 13.099333, 11.999389, 11.699404]
        end

        @testset "Filter Subfind" begin
           # check if filter works
            find_mass_gt_20(M) = ( (M > 20.0) ? true : false )
            dummy = filter_subfind("sub_002", "MTOP", find_mass_gt_20)
            @test dummy[1] == HaloID(0, 1)
            @test dummy[2] == HaloID(0, 2)

            dummy2 = filter_subfind("sub_002", filter_dummy)
            @test dummy == dummy2

            # find the most massive halo in the sample subfind output
            center, rvir, haloid = find_most_massive_halo("sub_002", 4)
            @test center ≈ [6382.1562, -1780.875, -7379.2812]
            @test rvir   ≈ 194.34073
            @test haloid == HaloID(0, 1) 
        end

        @testset "Read halo props" begin
            prop, haloid = read_halo_prop_and_id("sub_002", 22, "MTOP", 4)

            @test prop ≈ 11.299424
            @test haloid == HaloID(1, 6)
        end

        @testset "Error handling" begin
            # check error handling
            @test_throws ErrorException("Block MVIR not present!") read_subfind(subfile, "MVIR")    
        end
        
    end

    @testset "Snapshot utility" begin
        
        ref_file = joinpath(dirname(@__FILE__), "snap_sedov")

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

    @testset "Peano-Hilbert" begin
        
    end

    @testset "Write Snapshot" begin

        # read in reference file
        ref_file = joinpath(dirname(@__FILE__), "snap_sedov")
        head = head_to_obj(ref_file)
        x = read_snap(ref_file, "POS", 0)

        # specify output file for testing
        output_file = joinpath(dirname(@__FILE__), "write_test.dat")
        f = open(output_file, "w")

        @test_nowarn write_header(f, head)
        @test_nowarn write_block(f, x, "POS")
        close(f)

        pos_info = InfoLine("POS", Float32, 3, [1, 0, 0, 0, 0, 0])
        x_check = read_block(output_file, "POS", info=pos_info, parttype=0)

        # check if we read the same thing we wrote
        @test x_check == x

        f = open(output_file, "w")
        @test_nowarn write_block(f, x, "", snap_format=1)

        @test_throws ErrorException("Please specify blockname!") write_block(f, x, "")
        close(f)

        rm(output_file)
    end

@info "delete test data..."
rm("snap_sedov")
rm("pos_sedov.dat")
rm("sub_002.0")
rm("sub_002.1")
rm("sub_002.2")
rm("sub_002.3")
rm("snap_002.0")
rm("snap_002.1")
rm("snap_002.2")
rm("snap_002.3")
rm("snap_002.0.key")
rm("snap_002.1.key")
rm("snap_002.2.key")
rm("snap_002.3.key")
rm("snap_002.key.index")
@info "done!"
end

