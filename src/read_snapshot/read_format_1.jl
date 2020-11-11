

# function snap_1_d(filename::String, data::Dict{Any,Any})

#     f = open(filename)

#     seek(f, 264)

#     N = sum(data["Header"]["npart"])
#     skipsize = read(f, Int32)
#     bit_size = skipsize/(3*N)

#     if Int(bit_size) == 4

#     elseif Int(bit_size) == 8

#     else
#         println("read error! neither 32 nor 64 bits data!")
#         return -1
#     end

#     # read positions
#     for i in 1:length(data["Header"]["PartTypes"])

#         if data["Header"]["npart"][i] != 0

#             data[data["Header"]["PartTypes"][i]] = Dict()
#             N = Int64(data["Header"]["npart"][i])
#             #dummy = read(f, Float32, (3,N))

#             if Int(bit_size) == 4
#                 data[data["Header"]["PartTypes"][i]]["POS"] = copy(transpose(read!(f, Array{Float32,2}(undef,(3,N)))))
#             else
#                 data[data["Header"]["PartTypes"][i]]["POS"] = copy(transpose(read!(f, Array{Float64,2}(undef,(3,N)))))
#             end


#         end

#     end

#     p = position(f)

#     # skip identifiers
#     seek(f, p+8)

#     # read Velocities
#     for i in 1:length(data["Header"]["PartTypes"])

#         if data["Header"]["npart"][i] != 0

#             N = Int64(data["Header"]["npart"][i])

#             if Int(bit_size) == 4
#                 data[data["Header"]["PartTypes"][i]]["VEL"] = copy(transpose(read!(f, Array{Float32,2}(undef,(3,N)))))
#             else
#                 data[data["Header"]["PartTypes"][i]]["VEL"] = copy(transpose(read!(f, Array{Float64,2}(undef,(3,N)))))
#             end

#         end

#     end

#     #dummy = 0
#     #gc()

#     p = position(f)

#     # skip identifiers
#     seek(f, p+8)

#     # Read IDs
#     for i in 1:length(data["Header"]["PartTypes"])

#         if data["Header"]["npart"][i] != 0

#             N = Int64(data["Header"]["npart"][i])
#             data[data["Header"]["PartTypes"][i]]["ID"] = copy(transpose(read!(f, Array{UInt32,2}(undef,(1,N)))))

#         end

#     end

#     p = position(f)

#     # skip identifiers
#     seek(f, p+8)

#     # Read IDs
#     for i in 1:length(data["Header"]["PartTypes"])

#         if data["Header"]["npart"][i] != 0

#             if data["Header"]["massarr"][i] == 0

#                 N = Int64(data["Header"]["npart"][i])

#                 if Int(bit_size) == 4
#                     data[data["Header"]["PartTypes"][i]]["MASS"] = copy(transpose(read!(f, Array{Float32,2}(undef,(1,N)))))
#                 else
#                     data[data["Header"]["PartTypes"][i]]["MASS"] = copy(transpose(read!(f, Array{Float64,2}(undef,(1,N)))))
#                 end

#             else

#                 N = Int64(data["Header"]["npart"][i])
#                 data[data["Header"]["PartTypes"][i]]["MASS"] = zeros(N)
#                 data[data["Header"]["PartTypes"][i]]["MASS"] .= data["Header"]["massarr"][i]

#             end

#         end

#     end

#     if data["Header"]["npart"][1] != 0

#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)

#         N = Int64(data["Header"]["npart"][1])

#         # read U
#         if Int(bit_size) == 4
#             data["PartType0"]["InternalEnergy"] = copy(transpose(read!(f, Array{Float32,2}(undef,(1,N)))))
#         else
#             data["PartType0"]["InternalEnergy"] = copy(transpose(read!(f, Array{Float64,2}(undef,(1,N)))))
#         end


#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)

#         # Read Density
#         if Int(bit_size) == 4
#             data["PartType0"]["Density"] = copy(transpose(read!(f, Array{Float32,2}(undef,(1,N)))))
#         else
#             data["PartType0"]["Density"] = copy(transpose(read!(f, Array{Float64,2}(undef,(1,N)))))
#         end

#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)

#         # Read SmoothingLength
#         if Int(bit_size) == 4
#             data["PartType0"]["SmoothingLength"] = copy(transpose(read!(f, Array{Float32,2}(undef,(1,N)))))
#         else
#             data["PartType0"]["SmoothingLength"] = copy(transpose(read!(f, Array{Float64,2}(undef,(1,N)))))
#         end

#     end

#     close(f)

#     return data


# end
