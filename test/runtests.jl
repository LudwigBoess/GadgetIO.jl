using GadgetIO, Test, DelimitedFiles, Downloads


@testset "GadgetIO" begin

function filter_dummy(filename::String)
    mtop = read_subfind(filename, "MTOP")
    return findall(mtop .> 7.0)
end

function pass_all(snap_file)
    h = read_header(snap_file)
    return collect(1:h.npart[1])
end

@info "downloading test data..."
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_sedov", "./snap_sedov")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/pos_sedov.dat", "./pos_sedov.dat")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.0", "./sub_002.0")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.1", "./sub_002.1")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.2", "./sub_002.2")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.3", "./sub_002.3")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0", "./snap_002.0")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1", "./snap_002.1")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2", "./snap_002.2")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3", "./snap_002.3")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0.key", "./snap_002.0.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1.key", "./snap_002.1.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2.key", "./snap_002.2.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3.key", "./snap_002.3.key")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.key.index", "./snap_002.key.index")

@info "done!"

    @testset "Objects" begin
        @test_nowarn SnapshotHeader()
        @test_nowarn InfoLine()
    end

    @testset "Read Snapshot" begin

        snap_file = "snap_sedov"

        @testset "Read blocks" begin
            @test_nowarn read_snap(snap_file)
            @test_nowarn read_header(snap_file)
            @test_nowarn read_info(snap_file)

            d = read_snap(snap_file, "POS", 0)

            ideal_file = joinpath(dirname(@__FILE__), "pos_sedov.dat")
            d_ideal = copy(transpose(Float32.(readdlm(ideal_file))))

            @test d == d_ideal

            # check if read to dict works
            @test_nowarn d = read_snap(snap_file)
            @test d["PartType0"]["POS"] == d_ideal

            snap_base = "snap_002"

            blocks = ["POS", "RHO"]
            data = read_blocks_over_all_files(snap_base, blocks, filter_function=pass_all, parttype=0)

            rho = read_block(snap_base, "RHO", parttype=0)

            @test data["RHO"] == rho

            pos = read_block(snap_base, "POS", parttype=0)

            @test data["POS"] == pos
        end

        @testset "Read particles in box" begin
            center = Float32[3978.9688, -95.40625, -8845.25]
            rvir   = 118.76352
            pos = read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=false, parttype=1)

            @test pos[:,1] ≈ Float32[3882.5537, -20.574343, -8768.669]

            pos = read_particles_in_box("snap_002", "POS", center .- rvir, center .+ rvir, 
                                        use_keys=false, parttype=1, verbose=false)

            @test pos[:,1] ≈ Float32[3882.5537, -20.574343, -8768.669]

            # to do: proper tests!
            # need proper box for that!
            # x0 = [0.0, 0.0, 0.0]
            # x1 = [2_000.0, 2_000.0, 2_000.0]
            # @test_nowarn read_particles_in_box("snap_002", "POS", x0, x1, parttype=1)

        end

        @testset "Read particles in geometry" begin
            center = [3978.9688, -95.40625, -8845.25]
            rvir   = 118.76352

            @testset "Cube" begin
                
                cube = GadgetCube(center .- rvir, center .+ rvir)
                pos  = read_particles_in_geometry("snap_002", "POS", cube, use_keys=false, parttype=1)
                
                @test pos["POS"][:,1] ≈ Float32[3882.5537, -20.574343, -8768.669]
            end

            @testset "Sphere" begin
                sphere = GadgetSphere(center, rvir)

                @test_nowarn read_particles_in_geometry("snap_002", "POS", sphere, use_keys=false, parttype=1)
            end

            @testset "Cylinder" begin
                cylinder = GadgetCylinder(center .- 0.5rvir, center .+ 0.5rvir,
                                        0.5rvir)

                @test_nowarn read_particles_in_geometry("snap_002", "POS", cylinder, use_keys=false, parttype=1)
            end

            # to do: use key files!
            
        end

        @testset "Read particles in halo" begin
            pos = read_particles_in_halo("snap_002", "POS", "sub_002", HaloID(0,4), use_keys=false)

            @test pos[:,1] ≈ Float32[3909.1545, -189.9392, -8845.135]

            ids = UInt32[0x000028fc, 0x00002594, 0x00002963, 0x00002681, 0x00001af4, 0x00001ff1, 0x000022d7, 0x00002267, 0x000029c0, 0x0000277b]
            pos = read_particles_by_id("snap_002", ids, "POS")

            @test pos ≈ copy(transpose(Float32[-692.6776 -5005.1025 1474.2584; -734.53326 -4894.864 1665.7646; -756.7661 -4985.657 1942.4185; -801.0376 -4920.4683 1884.446; -907.67645 -4945.71 1895.1641; -939.883 -4893.6753 1874.1469; -932.33136 -4891.3984 1109.0826; -819.5988 -5004.6147 1254.0176; -644.03674 -4939.248 1164.3943; -667.2112 -5048.75 995.14856]))
        end

        @testset "Read positions" begin
            center = Float32[3978.9688, -95.40625, -8845.25]
            rvir   = 118.76352
            ff(filename) = filter_cube(filename, center .- rvir, center .+ rvir, parttype=1) 
            read_positions = find_read_positions("snap_002", ff)

            @test read_positions["N_part"] == 87

            @test read_positions[0]["index"][1] == 2441
            @test read_positions[0]["index"][5] == 3966

            @test read_positions[0]["n_to_read"][1] == 2
            @test read_positions[0]["n_to_read"][4] == 74

            # test IO
            save_read_positions("dummy.bin", read_positions)
            loaded_read_positions = load_read_positions("dummy.bin")
            delete!(loaded_read_positions, "N_part")

            @test read_positions == loaded_read_positions

            # delete dummy file
            rm("dummy.bin")
        end

        @testset "Error Handling" begin
            @test_throws ErrorException("Please specify particle type!") read_block("snap_002.0", "POS")  
            @test_throws ErrorException("Particle Type 5 not present in simulation!") read_block("snap_002.0", "POS", parttype=5, h=SnapshotHeader()) 
            @test_throws ErrorException("Block not present!") read_block("snap_002.0", "ABCD", parttype=0)  
            @test_throws ErrorException("Requested block ABCD not present!") GadgetIO.check_block_position("snap_002.0", "ABCD")  
            @test_throws ErrorException("Please provide either a dictionary with read positions or a filter function!") read_blocks_over_all_files("snap_002", ["POS"]) 
        end
    end

    @testset "Read subfind" begin

        subbase = "sub_002"
        subfile = "./sub_002.0"

        @testset "standard read" begin
            # check if standard reading works
            mtop = read_subfind("sub_002.0", "MTOP")
            @test mtop ≈ Float32[6.793532, 6.0309854, 7.230924, 8.415218]
        end

        @testset "Filter Subfind" begin
           # check if filter works
            find_mass_gt_7(M) = ( (M > 7.0) ? true : false )
            dummy = filter_subfind("sub_002", "MTOP", find_mass_gt_7)
            @test dummy[1] == HaloID(0, 3)
            @test dummy[2] == HaloID(0, 4)

            dummy2 = filter_subfind("sub_002", filter_dummy)
            @test dummy == dummy2

            # find the most massive halo in the sample subfind output
            center, rvir, haloid = find_most_massive_halo("sub_002", 4)
            @test center ≈ Float32[3978.9688, -95.40625, -8845.25]
            @test rvir   ≈ 118.76352
            @test haloid == HaloID(0, 4) 
        end

        @testset "Read halo props" begin
            prop, haloid = read_halo_prop_and_id("sub_002", 4, "MTOP", 4)

            @test prop ≈ 5.431016
            @test haloid == HaloID(1, 1)

            prop2 = read_halo_prop("sub_002", haloid, "MTOP")

            @test prop2 == prop
        end

        @testset "Error handling" begin
            # check error handling
            @test_throws ErrorException("Block MVIR not present!") read_subfind(subfile, "MVIR")  
            @test_throws ErrorException("Halo at index 1000 does not exist!") read_halo_prop_and_id(subfile, 1000, "MTOP", verbose=false)  
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

        info = read_info(ref_file)
        d = GadgetIO.allocate_data_dict(["POS", "ID"], 10, info, true)
        @test haskey(d, "MASS")

        @test_throws ErrorException("File dummy_snap.0 not present!") GadgetIO.select_file("dummy_snap", 0)  
        
        # shift across box border
        boxsize, boxsize_half = 10, 5
        @test GadgetIO.shift_across_box_border(1, 2, boxsize, boxsize_half) == 1
        @test GadgetIO.shift_across_box_border(8, 2, boxsize, boxsize_half) == -2
        @test GadgetIO.shift_across_box_border(2, 8, boxsize, boxsize_half) == 12
    end

    @testset "Peano-Hilbert" begin
        h_key = GadgetIO.read_keyheader("snap_002.0.key")

        @test h_key.nkeys_file[1] == Int32(597)
        @test h_key.bits == 10

        low_list, high_list, file_list = GadgetIO.read_key_index("snap_002.key.index")
        @test file_list == UInt32[0x00000000, 0x00000001, 0x00000002, 0x00000003, 0x00000000, 0x00000001, 0x00000002, 0x00000003, 0x00000000, 0x00000001, 0x00000002, 0x00000003, 0x00000000]
    
        @test GadgetIO.peano_hilbert_key(h_key.bits, 0, 0, 0) == 0
        @test GadgetIO.peano_hilbert_key(h_key.bits, 0, 0, 1) == 1
        @test GadgetIO.peano_hilbert_key(h_key.bits, 0, 1, 1) == 6

        @test GadgetIO.get_int_pos( 1000.5, h_key.domain_corners[1], h_key.domain_fac ) == 1

        @testset "Get Index Bounds" begin
            low_bounds, high_bounds = [0, 4, 7, 10], [1, 5, 8, 16]

            @test GadgetIO.get_index_bounds([1, 2, 3], low_bounds, high_bounds) == [1]
            @test GadgetIO.get_index_bounds([1, 2, 3, 10], low_bounds, high_bounds) == [1, 4]
            @test GadgetIO.get_index_bounds([1, 2, 10, 17], low_bounds, high_bounds) == [1, 4]
            @test GadgetIO.get_index_bounds([4, 5, 8, 12, 15], low_bounds, high_bounds) == [2, 3, 4]
        end

        keylist_ideal = UInt64[0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004, 0x0000000000000005, 0x0000000000000006, 0x0000000000000007]

        h_key = GadgetIO.read_keyheader("snap_002.0.key")
        x0 = [-100.0, -100.0, -100.0]
        x1 = [1_000.0, 1_000.0, 1_000.0]

        keylist = GadgetIO.get_keylist(h_key, x0, x1)

        @test keylist == keylist_ideal
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

    @testset "Find Index Locations" begin

        # create find and check arrays
        list_to_find = [2, 56, 354, 254, 653, 452, 7523, 45, 42, 742, 5423, 942, 105, 425, 815, 7821]
        list_to_check = [523, 9, 254, 653, 452, 2, 923, 815, 7821, 742, 354, 543, 942, 15, 25]

        # create sorted versions for forward search
        list_to_find_sorted = sort(list_to_find)
        list_to_check_sorted = sort(list_to_check)

        # correct indices of unsorted and sorted arrays
        indices = [3, 4, 5, 6, 8, 9, 10, 11, 13]
        indices_sorted = [1, 5, 6, 7, 10, 11, 12, 14, 15]

        # check methods
        @test issetequal(GadgetIO.get_index_list_arr(list_to_find_sorted, list_to_check_sorted), indices_sorted)
        @test issetequal(GadgetIO.get_index_list_dict(list_to_find, list_to_check), indices)
        @test issetequal(GadgetIO.get_index_list_dict(list_to_find_sorted, list_to_check_sorted), indices_sorted)
        @test issetequal(GadgetIO.get_index_list_set(list_to_find, list_to_check), indices)
        @test issetequal(GadgetIO.get_index_list_set(list_to_find_sorted, list_to_check_sorted), indices_sorted)

        # check that right methods are called
        @test get_index_list(list_to_check, list_to_find) == GadgetIO.get_index_list_set(list_to_check, list_to_find)
        @test get_index_list(list_to_check_sorted, list_to_find_sorted) == GadgetIO.get_index_list_arr(list_to_check_sorted, list_to_find_sorted)
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

