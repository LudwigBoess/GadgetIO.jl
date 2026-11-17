using GadgetIO, Test, DelimitedFiles, Downloads, Aqua, Random

function filter_dummy(filename::String)
    mtop = read_subfind(filename, "MTOP")
    return findall(mtop .> 7.0)
end

function pass_all(snap_file)
    h = read_header(snap_file)
    return collect(1:h.npart[1])
end

"""
Friends-of-friends labels from testing every particle pair.
"""
function fof_brute_force(pos, linking_length)
    N = size(pos, 2)
    parent = collect(1:N)
    root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
    for i = 1:N, j = i+1:N
        d2 = (Float64(pos[1, i]) - pos[1, j])^2 + (Float64(pos[2, i]) - pos[2, j])^2 + (Float64(pos[3, i]) - pos[3, j])^2
        if d2 <= Float64(linking_length)^2
            parent[root(j)] = root(i)
        end
    end
    return root.(1:N)
end

"""
True if the labels `a` and `b` divide the particles into the same groups.
"""
function same_partition(a, b)
    ab = Dict{eltype(a),eltype(b)}()
    ba = Dict{eltype(b),eltype(a)}()
    for i ∈ eachindex(a, b)
        get!(ab, a[i], b[i]) == b[i] || return false
        get!(ba, b[i], a[i]) == a[i] || return false
    end
    return true
end

"""
`N` positions of a Hernquist sphere with scale radius `a` around `center`, cut at `100a`.
"""
function hernquist_sphere(rng, N, a, center)
    r = Float64[]
    while length(r) < N
        u = sqrt(rand(rng))
        x = a * u / (1 - u)
        x < 100a && push!(r, x)
    end
    d = randn(rng, 3, N)
    return d ./ sqrt.(sum(abs2, d, dims=1)) .* r' .+ center
end

"""
Distance between `x` and `y` in a periodic box of side `L`.
"""
periodic_distance(x, y, L) = sqrt(sum(abs2, @. mod(x - y + L / 2, L) - L / 2))

"""
Write a single-file snapshot with the blocks POS and MASS for particle types 0-3 and POT for type 1.
"""
function write_test_snapshot(filename, h, pos, mass, pot)
    info = [InfoLine("POS", Float32, 3, [1, 1, 1, 1, 0, 0]),
            InfoLine("MASS", Float32, 1, [1, 1, 1, 1, 0, 0]),
            InfoLine("POT", Float32, 1, [0, 1, 0, 0, 0, 0])]
    f = open(filename, "w")
    write_header(f, h)
    write_info_block(f, info)
    write_block(f, pos, "POS")
    write_block(f, mass, "MASS")
    write_block(f, pot, "POT")
    close(f)
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

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.0", "./snap_144.0")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.1", "./snap_144.1")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.2", "./snap_144.2")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.3", "./snap_144.3")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.0.key", "./snap_144.0.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.1.key", "./snap_144.1.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.2.key", "./snap_144.2.key")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.3.key", "./snap_144.3.key")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.key.index", "./snap_144.key.index")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_mass_144.0", "./snap_mass_144.0")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_mass_144.1", "./snap_mass_144.1")

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/balance.txt", "./balance.txt")

@info "done!"


@testset "GadgetIO" begin

    @testset "Objects" begin
        @test_nowarn SnapshotHeader()
        @test_nowarn InfoLine("POS", Float32, 3, [1, 1, 1, 1, 1, 1])
    end

    @testset "Read Snapshot" begin

        snap_file = "snap_sedov"

        @testset "Format 2" begin

            @testset "Error handling" begin
                #@test_nowarn read_snap(snap_file)
                @test_nowarn read_header(snap_file)
                @test_nowarn read_info(snap_file)
            end

            @testset "single block" begin
                d = read_block(snap_file, "POS")

                ideal_file = joinpath(dirname(@__FILE__), "pos_sedov.dat")
                d_ideal = copy(transpose(Float32.(readdlm(ideal_file))))

                @test d == d_ideal

                snap_base = "snap_002"

                blocks = ["POS", "ID", "RHO"]
                data = read_blocks_filtered(snap_base, blocks, filter_function=pass_all, parttype=0)

                rho = read_block(snap_base, "RHO", parttype=0)

                @test length(data["RHO"]) == length(rho)

                @test data["RHO"] == rho

                id = read_block(snap_base, "ID", parttype=0)

                @test sort(data["ID"]) == sort(id)

                pos = read_block(snap_base, "POS", parttype=0)

                @test data["POS"] == pos

                @test_throws ErrorException read_block(snap_base, "RHO", parttype=1)
                @test_throws ErrorException read_block(snap_base, "RHO", parttype=-1)
                @test_throws ErrorException read_block(snap_base, "XXXX", parttype=0)
                @test_throws ErrorException read_block(snap_base, "XXXX", parttype=-1)
            end

            @testset "Read to Dict" begin
                # check if read to dict works
                ideal_file = joinpath(dirname(@__FILE__), "pos_sedov.dat")
                d_ideal = copy(transpose(Float32.(readdlm(ideal_file))))
                
                d = read_snap(snap_file)
                @test d["PartType0"]["POS"] == d_ideal
            end

            @testset "full block" begin

                snap_file = "snap_002.0"
                h = read_header(snap_file)
                n_all = sum(h.npart)

                mass = zeros(n_all)
                n_read = 0
                for parttype = 0:5
                    n_to_read = h.npart[parttype+1]
                    if !iszero(n_to_read)
                        mass[n_read+1:n_read+n_to_read] = read_block(snap_file, "MASS"; parttype)
                        n_read += n_to_read
                    end
                end

                mass_full = read_block(snap_file, "MASS", parttype=-1)

                @test mass_full == mass
                @test all(!iszero, mass)


                h = read_header(snap_file)
                n_all = sum(h.npart)


                id = zeros(UInt32, n_all)
                n_read = 0
                for parttype = 0:5
                    n_to_read = h.npart[parttype+1]
                    if !iszero(n_to_read)
                        id[n_read+1:n_read+n_to_read] = read_block(snap_file, "ID"; parttype)
                        n_read += n_to_read
                    end
                end

                id_full = read_block(snap_file, "ID", parttype=-1)

                @test id_full == id


                pos = zeros(Float32, 3, n_all)
                n_read = 0
                for parttype = 0:5
                    n_to_read = h.npart[parttype+1]
                    if !iszero(n_to_read)
                        pos[:, n_read+1:n_read+n_to_read] = read_block(snap_file, "POS"; parttype)
                        n_read += n_to_read
                    end
                end

                pos_full = read_block(snap_file, "POS", parttype=-1)

                @test pos_full == pos

                # distributed files
                # ToDo!
            end

            @testset "full block mass from header" begin

                snap_base = "snap_144"
                h = read_header(snap_base)
                n_all = sum(h.nall)

                read_mass_from_header = try
                    # check if code allows reading mass from the header
                    read_block(snap_base, "MASS"; parttype=0)
                    true
                catch
                    false
                end

                # if we allow reading the mass from the header, the read values also have to be consistent when
                # reading across multiple subfiles
                if read_mass_from_header
                    mass = zeros(n_all)
                    n_read = 0
                    for parttype = 0:5
                        n_to_read = h.nall[parttype+1]
                        if !iszero(n_to_read)
                            mass[n_read+1:n_read+n_to_read] = read_block(snap_base, "MASS"; parttype)
                            n_read += n_to_read
                        end
                    end

                    mass_full = read_block(snap_base, "MASS", parttype=-1)

                    @test mass_full == mass
                    @test all(!iszero, mass)
                end

            end

        end

        @testset "Format 1" begin
            snap_file = "snap_sedov"
            # read reference data in format 2
            pos_format2 = read_block(snap_file, "POS", parttype=0)
            # write format 1 block
            h = read_header(snap_file)
            f = open("dummy_snap_format1", "w")
            write_header(f, h, snap_format=1)
            write_block(f, pos_format2, "POS", snap_format=1)
            close(f)
            # read format 1 block
            pos_format1 = read_block("dummy_snap_format1", 1, parttype=0)
            # check if the data is the same
            @test pos_format1 == pos_format2
        end

        @testset "Read particles in box" begin

            @testset "Brute-Force" begin
                center = Float32[3978.9688, -95.40625, -8845.25]
                rvir = 118.76352
                pos = read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=false, parttype=1)

                @test pos[:, 1] ≈ Float32[3882.5537, -20.574343, -8768.669]

                pos = read_particles_in_box("snap_002", "POS", center .- rvir, center .+ rvir,
                    use_keys=false, parttype=1, verbose=false)

                @test pos[:, 1] ≈ Float32[3882.5537, -20.574343, -8768.669]
            end

            @testset "Peano-Hilbert reading" begin

                center = [32572.607, 33825.95, 26314.158]
                rvir = 1138.885

                pos = read_particles_in_volume("snap_144", "POS", center, rvir, use_keys=true, parttype=1, verbose=true)

                id = read_particles_in_volume("snap_144", "ID", center, rvir, use_keys=true, parttype=1, verbose=false)

                # check if correct number of particles were read
                @test size(pos, 2) == 11981

                # check if the positions are the same as in Klaus' IDL code
                @test pos[:, 1] ≈ Float32[31486.8, 33696.7, 25475.5]
                @test pos[:, 2] ≈ Float32[32027.1, 33818.6, 25249.8]
                @test pos[:, 3] ≈ Float32[31978.6, 33792.4, 25175.6]

                @test pos[:, 1:10] ≈ Float32[31486.781 32027.072 31978.645 32176.344 32121.607 32016.475 31969.11 31932.04 32121.26 32127.13; 33696.67 33818.57 33792.375 33933.582 34023.098 33961.69 33963.83 34092.285 34075.18 34012.707; 25475.465 25249.795 25175.584 25233.756 25304.545 25254.64 25245.88 25347.342 25327.582 25469.148]
                @test pos[:, end-10:end] ≈ Float32[32736.15 32729.447 32881.92 32840.797 31663.264 31936.314 32008.166 32007.818 32062.03 32029.402 32049.965; 34642.566 34646.56 34533.527 34535.62 34514.78 34576.348 34734.242 34731.906 34826.69 34821.555 34925.05; 25311.725 25476.305 25346.922 25362.049 25434.195 25240.049 25428.146 25396.178 25425.773 25233.172 25381.234]

                # check if IDs are correct 
                #@test id[1:10] ≈ UInt64[0x00000000000d8b92, 0x00000000000c3027, 0x00000000000df58f, 0x00000000000bddac, 0x00000000000dfe19, 0x00000000000c0864, 0x00000000000c1f2c, 0x00000000000de478, 0x00000000000bd814, 0x00000000000df867]
                #@test id[end-10:end] ≈ UInt64[0x00000000000d4cb4, 0x00000000000c3ee6, 0x00000000000bebea, 0x00000000000ce9c3, 0x00000000000e5e7d, 0x00000000000e3978, 0x00000000000c331a, 0x00000000000c10ee, 0x00000000000c13c7, 0x00000000000c24dd, 0x00000000000e364e]

                @test id[11480-10:11480] ≈ UInt64[0x00000000000db2e7, 0x00000000000d9f13, 0x00000000000db89a, 0x00000000000dc9af, 0x00000000000de350, 0x00000000000dbe65, 0x00000000000ddd9e, 0x00000000000de629, 0x00000000000dbe4b, 0x00000000000dbb72, 0x00000000000de077]
                @test id[11481:11490]    ≈ UInt64[0x00000000000d9960, 0x00000000000da206, 0x00000000000be085, 0x00000000000bc6fe, 0x00000000000d858d, 0x00000000000d9688, 0x00000000000d6ec5, 0x00000000000d57fd, 0x00000000000e1456, 0x00000000000ce048]
            end

        end

        @testset "Read particles in geometry" begin
            center = [3978.9688, -95.40625, -8845.25]
            rvir = 118.76352

            @testset "Cube" begin

                cube = GadgetCube(center .- rvir, center .+ rvir)
                pos = read_particles_in_geometry("snap_002", "POS", cube, use_keys=false, parttype=1)

                @test pos["POS"][:, 1] ≈ Float32[3882.5537, -20.574343, -8768.669]
            end

            @testset "Sphere" begin
                sphere = GadgetSphere(center, rvir)

                data = read_particles_in_geometry("snap_002", "POS", sphere, use_keys=false, parttype=1)

                @test data["POS"][:, 1:3] ≈ Float32[3904.7957 3871.5486 4038.2986; -135.69522 -79.47088 -62.441578; -8837.774 -8831.329 -8873.304]
            end

            @testset "Cylinder" begin
                cylinder = GadgetCylinder(center .- 0.5rvir, center .+ 0.5rvir,
                    0.5rvir)

                data = read_particles_in_geometry("snap_002", "POS", cylinder, use_keys=false, parttype=1)

                @test data["POS"][:, 1:3] ≈ Float32[3904.7957 4049.4988 4035.499; -135.69522 -91.40538 -105.1906; -8837.774 -8828.907 -8844.973]
            end

            # to do: use key files!

        end

        @testset "Read particles in halo" begin

            pos = read_particles_in_halo("snap_002", "POS", "sub_002", HaloID(0, 4), use_keys=false)

            #@test pos[:, 1] ≈ Float32[3909.1545, -189.9392, -8845.135]
            @test pos[:, 1] ≈ [3978.563, -96.09807, -8846.737]

            ids = UInt32[0x000028fc, 0x00002594, 0x00002963, 0x00002681, 0x00001af4, 0x00001ff1, 0x000022d7, 0x00002267, 0x000029c0, 0x0000277b]
            pos = read_particles_by_id("snap_002", ids, "POS")

            @test pos ≈ copy(transpose(Float32[-692.6776 -5005.1025 1474.2584; -734.53326 -4894.864 1665.7646; -756.7661 -4985.657 1942.4185; -801.0376 -4920.4683 1884.446; -907.67645 -4945.71 1895.1641; -939.883 -4893.6753 1874.1469; -932.33136 -4891.3984 1109.0826; -819.5988 -5004.6147 1254.0176; -644.03674 -4939.248 1164.3943; -667.2112 -5048.75 995.14856]))
        end

        @testset "Read positions" begin
            center = Float32[3978.9688, -95.40625, -8845.25]
            rvir = 118.76352

            # with filter function
            ff(filename) = filter_cube(filename, center .- rvir, center .+ rvir, parttype=1)
            read_positions = find_read_positions("snap_002", ff)

            @test read_positions["N_part"] == 87

            @test read_positions[0]["index"][1] == 2441
            @test read_positions[0]["index"][5] == 3966

            @test read_positions[0]["n_to_read"][1] == 2
            @test read_positions[0]["n_to_read"][4] == 74

            # with gadget geometry
            cube = GadgetCube(center .- rvir, center .+ rvir)
            read_positions_geo = find_read_positions("snap_002", cube, parttype=1)

            @test read_positions_geo["N_part"] == 87

            @test read_positions_geo[0]["index"][1] == 2441
            @test read_positions_geo[0]["index"][5] == 3966

            @test read_positions_geo[0]["n_to_read"][1] == 2
            @test read_positions_geo[0]["n_to_read"][4] == 74

            # test IO
            save_read_positions("dummy.bin", read_positions)
            loaded_read_positions = load_read_positions("dummy.bin")
            delete!(loaded_read_positions, "N_part")

            @test read_positions == loaded_read_positions

            # delete dummy file
            rm("dummy.bin")
        end

        @testset "Read files different infolines" begin
            subbase = "./snap_mass_144"
            subfile0 = "$subbase.0"
            subfile1 = "$subbase.1"

            mass = read_block(subbase, "MASS"; parttype=4)
            mass0 = read_block(subfile0, "MASS"; parttype=4)
            mass1 = read_block(subfile1, "MASS"; parttype=4)

            # check that reading mass did not read boundary particles by accident;
            # the reason is that the first file does not have particle type 3, but
            # the second file does
            @test mass == [mass0; mass1]
        end

        @testset "Error Handling" begin
            #@test_throws ErrorException("Please specify particle type!") read_block("snap_002.0", "POS")
            @test_throws ErrorException("Particle Type 5 not present in simulation!") read_block("snap_002.0", "POS", parttype=5, h=SnapshotHeader())
            @test_throws ErrorException("Block ABCD not present!") read_block("snap_002.0", "ABCD", parttype=0)
            @test_throws ErrorException("Requested block ABCD not present!") GadgetIO.check_block_position("snap_002.0", "ABCD")
            @test_throws ErrorException("Please provide either a dictionary with read positions or a filter function!") read_blocks_filtered("snap_002", ["POS"])
        end

    end

    @testset "Read subfind" begin

        subbase = "sub_002"
        subfile = "./sub_002.0"

        @testset "converting header" begin
            # Checking that highword conversion works for large numbers of halos (example taken from Magneticum Box2b)
            h = SnapshotHeader(Int32[22, 15453, 64805080, 22, 0, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.7986652252549261, 0.25208907108833967, Int32(1), Int32(1), UInt32[0x0168391c, 0x01a81b1a, 0x8e18f082, 0x00000cf3, 0x00000000, 0x00000000], Int32(1), Int32(1024), 640000.0, 0.272, 0.728, 0.704, Int32(1), Int32(0), UInt32[0x00000000, 0x00000000, 0x00000004, 0x00000000, 0x00000000, 0x00000000], Int32(0), Int32(0), Int32(3), 0.0f0, Int32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            h_subfind = GadgetIO.convert_header(h)
            @test h_subfind.totfof == get_total_particles(h.nall[3], h.npartTotalHighWord[3])
        end

        @testset "standard read" begin
            # check if standard reading works
            mtop = read_subfind("sub_002.0", "MTOP")
            @test mtop ≈ Float32[6.793532, 6.0309854, 7.230924, 8.415218]
        end

        @testset "Filter Subfind" begin

            @testset "subfind length" begin
                @test GadgetIO.read_subfind_length("sub_002.0", "MTOP") == 4
            end

            @testset "Filtering" begin
                # check if filter works
                dummy = filter_subfind("sub_002", filter_dummy)
                @test dummy[1] == HaloID(0, 3)
                @test dummy[2] == HaloID(0, 4)

                # find the most massive halo in the sample subfind output
                center, rvir, haloid = find_most_massive_halo("sub_002", 4)
                @test center ≈ Float32[3978.9688, -95.40625, -8845.25]
                @test rvir ≈ 118.76352
                @test haloid == HaloID(0, 4)
            end

            @testset "Read Positions" begin

                # get same filtering as before
                dummy = filter_subfind("sub_002", filter_dummy)
                read_positions = GadgetIO.halo_ids_to_read_positions(dummy)

                blocks = ["GPOS", "MTOP"]
                data_from_haloids = read_halo_prop("sub_002", blocks, dummy)

                data_from_read_positions = read_blocks_filtered("sub_002", blocks; read_positions)

                @test data_from_read_positions == data_from_haloids
            end

            @testset "IO" begin
                # check if filter works
                dummy = filter_subfind("sub_002", filter_dummy)

                filename = "halo_ids.dat"

                @test_nowarn save_halo_ids(filename, dummy)

                dummy2 = load_halo_ids(filename)
                @test dummy2 == dummy
            end

        end

        @testset "Read halo props" begin
            prop, haloid = read_halo_prop_and_id("sub_002", "MTOP", 4)

            @test prop ≈ 5.431016
            @test haloid == HaloID(1, 1)

            prop2 = read_halo_prop("sub_002", "MTOP", haloid)

            @test prop2 == prop

            mtop = read_subfind("sub_002", "MTOP")
            mtop2, haloids = read_subfind("sub_002", "MTOP", return_haloid=true)

            @test mtop == mtop2
            @test read_halo_prop("sub_002", "MTOP", haloids[end-4]) == mtop[end-4]

            ids = [0, 10, length(haloids) - 3]
            @test read_subfind("sub_002", "MTOP", ids) == mtop[ids.+1]

            mtop_from_ids, haloids_from_ids = read_subfind("sub_002", "MTOP", ids; return_haloid=true)
            @test read_subfind("sub_002", "MTOP", ids) == mtop_from_ids
            @test haloids_from_ids == haloids[ids.+1]

            # test large offset by reading last couple of halos
            ids_end = [length(haloids) - 1, length(haloids) - 3]
            @test read_subfind("sub_002", "MTOP", ids_end) == mtop[ids_end.+1]

            mtop_from_ids, haloids_from_ids = read_subfind("sub_002", "MTOP", ids_end; return_haloid=true)
            @test read_subfind("sub_002", "MTOP", ids_end) == mtop_from_ids
            @test haloids_from_ids == haloids[ids_end.+1]

            pos = read_subfind("sub_002", "SPOS")
            pos2, haloids = read_subfind("sub_002", "SPOS"; return_haloid=true)
            @test pos == pos2

            # reading halo property for HaloID 
            @test read_halo_prop("sub_002", "SPOS", haloids[end-4]) == pos[:, end-4]
            
            # reading halo property for index (0-based)
            @test read_halo_prop("sub_002", "SPOS", length(haloids) - 5) == pos[:, end-4]

            ids = [0, 10, length(haloids) - 3]
            @test read_subfind("sub_002", "SPOS", ids) == pos[:, ids.+1]

            pos_from_ids, haloids_from_ids = read_subfind("sub_002", "SPOS", ids; return_haloid=true)
            @test read_subfind("sub_002", "SPOS", ids) == pos_from_ids
            @test haloids_from_ids == haloids[ids.+1]

            # load reference data
            msub_from_ids, haloids_from_ids = read_subfind("sub_002", "MSUB", ids; return_haloid=true)

            # read multiple blocks at the same time
            blocks = ["SPOS", "MSUB"]

            # for array of indices
            data_from_ids = read_halo_prop("sub_002", blocks, ids)
            @test data_from_ids["SPOS"] == pos_from_ids
            @test data_from_ids["MSUB"] == msub_from_ids

            # for single index
            data_from_ids = read_halo_prop("sub_002", blocks, ids[1])
            @test data_from_ids["SPOS"][1] == pos_from_ids[1]
            @test data_from_ids["MSUB"][1] == msub_from_ids[1]

            # for single block and array of indices
            @test read_halo_prop("sub_002", "SPOS", ids) == pos_from_ids
            @test read_halo_prop("sub_002", "MSUB", ids) == msub_from_ids

            # for array of HaloIDs
            data_from_haloids = read_halo_prop("sub_002", blocks, haloids[ids.+1])
            @test data_from_haloids["SPOS"] == pos_from_ids
            @test data_from_haloids["MSUB"] == msub_from_ids

            # for single HaloID
            data_from_haloids = read_halo_prop("sub_002", blocks, haloids[ids.+1][1])
            @test data_from_haloids["SPOS"][1] == pos_from_ids[1]
            @test data_from_haloids["MSUB"][1] == msub_from_ids[1]

            # for single block and array of HaloIDs
            @test read_halo_prop("sub_002", "SPOS", haloids[ids.+1]) == pos_from_ids
            @test read_halo_prop("sub_002", "MSUB", haloids[ids.+1]) == msub_from_ids
        end

        @testset "Error handling" begin
            # check error handling
            @test_throws ErrorException("Block MVIR not present!") read_subfind(subfile, "MVIR")
            #@test_throws ErrorException("Halo at index 1000 does not exist!") read_halo_prop_and_id(subfile, 1000, "MTOP", verbose = false)

            filename = "pos_sedov.dat"
            @test_throws ErrorException("incorrect file format encountered when reading header of $filename") read_subfind_header(filename)

            # halo prop error handling
            @test_throws ErrorException("All requested blocks must be for the same halo type. Block MTOP is not available for halo type 1 but only for halo type 0.") read_halo_prop("sub_002", ["SPOS", "MTOP"], 1:3)

            @test_logs (:warn, "The Vector of HaloIDs is not sorted for requesting the properties from Subfind, the returned properties are returned as if they were sorted, however.") read_halo_prop("sub_002", ["GPOS", "MTOP"], [HaloID(0,3), HaloID(0,2)], verbose=false)
            @test_logs (:warn, "The Vector of i_global is not sorted for requesting the properties from Subfind, the returned properties are returned as if they were sorted, however.") read_halo_prop("sub_002", ["GPOS", "MTOP"], [3,2], verbose=false)
        end

    end

    @testset "FoF halo finder" begin

        rng = Xoshiro(42)

        @testset "FoF groups" begin

            # uniform points below, at and above the percolation threshold n l^3 ≈ 0.65
            pos = rand(rng, Float32, 3, 2000)
            for nl3 ∈ (0.3, 0.65, 1.0)
                linking_length = cbrt(nl3 / 2000)
                @test same_partition(GadgetIO.fof_labels(pos, linking_length), fof_brute_force(pos, linking_length))
            end

            # clumps of different density on a uniform background
            pos = hcat([randn(rng, 3, 300) .* 0.005k .+ rand(rng, 3) for k = 1:8]..., rand(rng, 3, 500))
            for linking_length ∈ (0.002, 0.01, 0.05)
                @test same_partition(GadgetIO.fof_labels(pos, linking_length), fof_brute_force(pos, linking_length))
            end

            # isolated pairs are linked if and only if they are no further apart than the linking length
            linking_length = 1.7
            npairs = 8000
            pos = zeros(3, 2npairs)
            for k = 1:npairs
                d = randn(rng, 3)
                pos[:, 2k-1] = 10.0 .* [k % 20, (k ÷ 20) % 20, k ÷ 400]
                pos[:, 2k] = pos[:, 2k-1] .+ d ./ sqrt(sum(abs2, d)) .* linking_length .* (0.98 + 0.04rand(rng))
            end
            label = GadgetIO.fof_labels(pos, linking_length)
            @test all((label[2k-1] == label[2k]) == (sum(abs2, pos[:, 2k] .- pos[:, 2k-1]) <= linking_length^2) for k = 1:npairs)

            # chains along the diagonal of the cells: one group below the linking length, single particles above
            for spacing ∈ (0.99, 1.01)
                pos = [1.0, 1.0, 1.0] ./ √3 .* (spacing .* (0:500))'
                @test length(unique(GadgetIO.fof_labels(pos, 1.0))) == (spacing < 1 ? 1 : 501)
            end

            @test isempty(GadgetIO.fof_labels(zeros(3, 0), 1.0))
            @test GadgetIO.fof_labels(rand(rng, 3, 1), 1.0) == [1]
            @test length(unique(GadgetIO.fof_labels(0.1 .* rand(rng, 3, 100), 1.0))) == 1

            # particles across the whole box, or along a line, and an invalid linking length
            @test_throws ErrorException GadgetIO.fof_labels([0.0 1.0e6; 0.0 1.0e6; 0.0 0.0], 1.0e-2)
            @test_throws ErrorException GadgetIO.fof_labels([0.0 0.0; 0.0 0.0; 0.0 1.0e12], 1.0)
            @test_throws ErrorException GadgetIO.fof_labels(rand(rng, 3, 10), 0.0)
        end

        @testset "Subfind FoF groups" begin

            # for linking lengths of 41.8 - 42.0 the groups are those of the on-the-fly FoF of
            # the simulation, compared by the IDs of their type 1 particles
            pos = read_block("snap_002", "POS", parttype=1)
            id = read_block("snap_002", "ID", parttype=1)
            label = GadgetIO.fof_labels(pos, 41.88)
            label_of = Dict(zip(id, label))

            glen = read_subfind("sub_002", "GLEN")
            goff = read_subfind("sub_002", "GOFF")
            pid = read_subfind("sub_002", "PID")
            for g ∈ eachindex(glen)
                ids = filter(in(keys(label_of)), pid[goff[g]+1:goff[g]+glen[g]])
                @test Set(id[label.==label_of[ids[1]]]) == Set(ids)
            end
        end

        @testset "Shrinking sphere" begin

            # Hernquist sphere with an offset clump on a uniform background
            center = [1234.5, -678.9, 42.0]
            pos = hcat(hernquist_sphere(rng, 20_000, 1.0, center),
                       hernquist_sphere(rng, 5_000, 0.5, center .+ [6.0, 0.0, 0.0]),
                       center .+ 40.0 .* (rand(rng, 3, 4_000) .- 0.5))
            @test sqrt(sum(abs2, GadgetIO.shrinking_sphere_center(pos, ones(size(pos, 2))) .- center)) < 0.05

            # coincident particles end the iteration
            @test GadgetIO.shrinking_sphere_center(repeat([1.0, 2.0, 3.0], 1, 500), ones(500)) ≈ [1.0, 2.0, 3.0]
            @test_throws ErrorException GadgetIO.shrinking_sphere_center(zeros(3, 0), Float64[])
        end

        @testset "find_main_halo" begin

            # the main halo is the most massive Subfind FoF group
            gpos = Float64.(read_subfind("sub_002", "GPOS")[:, 1])
            rtop = read_subfind("sub_002", "RTOP")[1]

            center = find_main_halo("snap_002", linking_length=41.88, center=:shrinking_sphere, verbose=false)
            @test sqrt(sum(abs2, center .- gpos)) < 0.25rtop

            center = find_main_halo("snap_002", center=:shrinking_sphere, verbose=false)
            @test sqrt(sum(abs2, center .- gpos)) < 0.25rtop

            # default linking length from the mean density in the box, ρ_crit in kpc/h, 1e10 Msun/h, km/s
            mass = read_block("snap_002", "MASS", parttype=1)
            ρ_crit = 3 * 0.1^2 / (8π * 43018.7)
            @test GadgetIO.fof_linking_length("snap_002", mass) ≈ 0.2 * cbrt(sum(Float64, mass) / length(mass) / ((0.24 - 0.04) * ρ_crit)) rtol = 1.0e-3

            @test_throws ErrorException("Block POT not present, use center=:shrinking_sphere for snapshots without potential.") find_main_halo("snap_002", verbose=false)
            @test_throws ErrorException find_main_halo("snap_002", center=:mean, verbose=false)
        end

        @testset "Periodic box and potential" begin

            # single-file copy of snap_002 with positions on a grid of 1/32, so moving them across
            # the border of the box is exact in Float32
            h = read_header("snap_002")
            h.npart .= h.nall
            h.num_files = 1
            L = h.boxsize

            pos = reduce(hcat, [read_block("snap_002", "POS"; parttype) for parttype = 0:3])
            mass = reduce(vcat, [read_block("snap_002", "MASS"; parttype) for parttype = 0:3])
            pos = round.(32pos) ./ 32

            # members of the main halo from the IDs of the most massive Subfind group
            id = read_block("snap_002", "ID", parttype=1)
            pid = read_subfind("sub_002", "PID")
            main = findall(in(Set(pid[1:read_subfind("sub_002", "GLEN")[1]])), id)
            pos1 = pos[:, get_total_particles(h, 0) .+ eachindex(id)]

            # potential of type 1 with the minimum outside of the main halo
            gpos = Float64.(read_subfind("sub_002", "GPOS")[:, 1])
            pot = Float32.(vec(sum(abs2, pos1 .- gpos, dims=1)))
            i_center = main[argmin(pot[main])]
            pot[findfirst(i -> i ∉ main, eachindex(pot))] = -1.0f10

            # move the main halo onto the border of the box
            shift = [round(32 * (L / 2 - gpos[1])) / 32, 0.0, 0.0]
            pos_shifted = Float32.(mod.(pos .+ shift .+ L / 2, L) .- L / 2)

            mktempdir() do dir
                snap = joinpath(dir, "snap")
                snap_shifted = joinpath(dir, "snap_shifted")
                write_test_snapshot(snap, h, pos, mass, pot)
                write_test_snapshot(snap_shifted, h, pos_shifted, mass, pot)

                # the centre is the member with the lowest potential, not the global minimum
                center = find_main_halo(snap, linking_length=41.88, verbose=false)
                @test center ≈ pos1[:, i_center]

                center = find_main_halo(snap_shifted, linking_length=41.88, verbose=false)
                @test periodic_distance(center, pos1[:, i_center] .+ shift, L) < 1.0e-3

                # the halo is found in one piece across the border of the box
                center_ref = find_main_halo(snap, linking_length=41.88, center=:shrinking_sphere, verbose=false)
                center = find_main_halo(snap_shifted, linking_length=41.88, center=:shrinking_sphere, verbose=false)
                @test periodic_distance(center, center_ref .+ shift, L) < 1.0e-3

                @test GadgetIO.fof_linking_length(snap_shifted, read_block(snap_shifted, "MASS", parttype=1)) ≈
                      GadgetIO.fof_linking_length(snap, read_block(snap, "MASS", parttype=1))
            end
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
        blocks_checked, no_mass_block = GadgetIO.check_blocks(ref_file, blocks, 0)

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

        h = read_header(filename)
        @test get_total_particles(h, 0) == 125000
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

        @test GadgetIO.get_int_pos(1000.5, h_key.domain_corners[1], h_key.domain_fac, h_key.bits) == 1

        @testset "Get Index Bounds" begin
            low_bounds, high_bounds = [0, 4, 7, 10], [1, 5, 8, 16]

            @test GadgetIO.get_index_bounds([1, 2, 3], low_bounds, high_bounds) == [1]
            @test GadgetIO.get_index_bounds([1, 2, 3, 10], low_bounds, high_bounds) == [1, 4]
            @test GadgetIO.get_index_bounds([1, 2, 10, 17], low_bounds, high_bounds) == [1, 4]
            @test GadgetIO.get_index_bounds([4, 5, 8, 12, 15], low_bounds, high_bounds) == [2, 3, 4]
        end

        keylist_ideal = [0, 1, 2, 3, 4, 5, 6, 7, 153391689, 153391690, 153391693, 153391694, 383479221, 383479222, 421827144, 421827145, 421827150, 421827151, 651914678, 651914679, 690262601, 920350129, 920350134, 1073741816, 1073741819, 1073741820, 1073741823]

        h_key = GadgetIO.read_keyheader("snap_002.0.key")
        x0 = [-100.0, -100.0, -100.0]
        x1 = [1_000.0, 1_000.0, 1_000.0]

        keylist = GadgetIO.get_keylist(h_key, x0, x1)

        @test keylist == keylist_ideal
    end

    @testset "Write Snapshot" begin

        # read in reference file
        ref_file = joinpath(dirname(@__FILE__), "snap_sedov")
        head = head_to_struct(ref_file)
        x = read_snap(ref_file, "POS", 0)

        # specify output file for testing
        output_file = joinpath(dirname(@__FILE__), "write_test.dat")
        f = open(output_file, "w")

        @test_logs (:info, "Writing block: HEAD") write_header(f, head)
        @test_logs (:info, "Writing block: POS ") (:info, "Writing block done.") write_block(f, x, "POS")
        close(f)

        pos_info = InfoLine("POS", Float32, 3, [1, 0, 0, 0, 0, 0])
        x_check = read_block(output_file, "POS", info=pos_info, parttype=0)

        # check if we read the same thing we wrote
        @test x_check == x

        # test info block write
        @testset "Write INFO block" begin
            ref_file = joinpath(dirname(@__FILE__), "snap_sedov")
            info_snap = read_info(ref_file)

            output_file = joinpath(dirname(@__FILE__), "write_test.dat")
            f = open(output_file, "w")
            write_info_block(f, info_snap)
            close(f)

            info_file = read_info(output_file)

            for i = 1:length(info_snap)

                @test info_file[i].block_name == info_snap[i].block_name
                @test info_file[i].data_type == info_snap[i].data_type
                @test info_file[i].n_dim == info_snap[i].n_dim
                @test info_file[i].is_present == info_snap[i].is_present

            end
        end
        # snap format 1
        f = open(output_file, "w")
        @test_logs (:info, "Writing block done.") write_block(f, x, "", snap_format=1)

        @test_throws ErrorException("Please specify blockname!") write_block(f, x, "")
        close(f)

        rm(output_file)
    end

    @testset "Find Index Locations" begin

        # create find and check arrays
        list_to_find = [1, 29, 12]
        list_to_check = [12, 1, 18, 19, 29]

        # create sorted versions for forward search
        list_to_find_sorted = sort(list_to_find)
        list_to_check_sorted = sort(list_to_check)

        # correct indices of unsorted and sorted arrays
        indices = [2, 5, 1]
        indices_sorted = [1, 2, 5]

        # check methods
        @test issetequal(GadgetIO.get_index_list_arr(list_to_find_sorted, list_to_check_sorted), indices_sorted)
        @test issetequal(GadgetIO.get_index_list_dict(list_to_find, list_to_check), indices)
        @test issetequal(GadgetIO.get_index_list_dict(list_to_find_sorted, list_to_check_sorted), indices_sorted)
        @test issetequal(GadgetIO.get_index_list_set(list_to_find, list_to_check), indices_sorted)
        @test issetequal(GadgetIO.get_index_list_set(list_to_find_sorted, list_to_check_sorted), indices_sorted)

        # check that right methods are called
        @test get_index_list(list_to_check, list_to_find) == GadgetIO.get_index_list_dict(list_to_check, list_to_find)
        @test get_index_list(list_to_check_sorted, list_to_find_sorted) == GadgetIO.get_index_list_arr(list_to_check_sorted, list_to_find_sorted)
    end

    @testset "CPU log files" begin

        @testset "balance.txt" begin

            # read the balance file 
            steps, timing, active = parse_balance("balance.txt")

            # check if all arrays contain the same number of entries 
            @test length(steps) == length(timing) == length(active)

            # check if all timesteps are read 
            @test steps[end] == 44666

            # check if correct total simulation time is read [s]
            @test sum(timing) == 195442.071

            # check if maximum number of active particles is consistent with total particles in simulation 
            @test maximum(active) == 89212289

            # check if printng works 
            @test_nowarn print_performance("balance.txt")
        end

    end

    @testset "Aqua.jl (code quality)" begin
        Aqua.test_all(GadgetIO)
    end

end

@info "delete test data..."
rm("snap_sedov")
rm("pos_sedov.dat")
rm("halo_ids.dat")
rm("sub_002.0")
rm("sub_002.1")
rm("sub_002.2")
rm("sub_002.3")
rm("dummy_snap_format1")

rm("snap_002.0")
rm("snap_002.1")
rm("snap_002.2")
rm("snap_002.3")
rm("snap_002.0.key")
rm("snap_002.1.key")
rm("snap_002.2.key")
rm("snap_002.3.key")
rm("snap_002.key.index")

rm("snap_144.0")
rm("snap_144.1")
rm("snap_144.2")
rm("snap_144.3")
rm("snap_144.0.key")
rm("snap_144.1.key")
rm("snap_144.2.key")
rm("snap_144.3.key")
rm("snap_144.key.index")

rm("snap_mass_144.0")
rm("snap_mass_144.1")

rm("balance.txt")

@info "done!"
