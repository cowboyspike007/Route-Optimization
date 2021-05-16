using Pkg, Random
using Dates
using PyPlot, Logging
import Base.show

struct Node
    id::Int16
    x::Float32
    y::Float32
end

function show(io::IO, n::Node)
    print(string("Node:", n.id))
end

import Base.copy
    
mutable struct Firefly
    path::Array{Node}
    cost::Float32
end

function show(io::IO, n::Firefly)
    print(string("Firefly:", n.cost, " "))
end

copy(s::Firefly) = Firefly(s.path, s.cost)

function arc_distance(f1::Firefly, f2::Firefly)::Int
    """Calculates the arc distance between two firefly paths.
    Example:
        Firefly i -> 1 2 4 6 3 5 7
        Firefly j -> 1 2 4 3 5 7 6 
        1-2, 2-4, 3-5, 5-7 arcs in firefly i are also presented
        in firefly j. But 4-6, 6-3 arcs are missing in
        firefly j. so the arc distance between firefly i
        and firefly j is 2 here.
    """
    missmatched_pair_counts = 0
    f1_pairs = zip(f1.path[1:end-1], f1.path[2:end])
    f2_pairs = zip(f2.path[1:end-1], f2.path[2:end])

    for f1_pair in f1_pairs
        if f1_pair ∉ f2_pairs missmatched_pair_counts += 1 end
    end

    return missmatched_pair_counts
end

function attractiveness(source_node::Firefly, destination_node::Firefly, λ::Float32)::Float32
    """Firefly’s attractiveness is proportional to the light
    intensity seen by adjacent fireflies, we can now define the
    attractiveness β of a firefly depending on distance between
    two fireflies.
    """
    β = source_node.cost
    r, _ = hamming_distance(source_node, destination_node)
    return β * exp(-λ ^ (r^2))
end

function create_distance_matrix(nodes::Array{Node})
    """Creates and returns the matrix of distances between given nodes.
    Uses euclidean_distance function to calculate distance.
    """
    node_num = length(nodes)
    result_matrix = zeros(node_num, node_num)

    for n1 in nodes
        for n2 in nodes
            result_matrix[n1.id,n2.id] = haversine_distance(n1, n2)
        end
    end
    return result_matrix
end

function generate_mutated_fireflies(firefly::Firefly, distance_matrix, r::Int, k::Int)::Array{Firefly}
    """Generates mutated fireflies using inversion_mutation function.
    Returns the new generated firefly array.
    """
    return [inversion_mutation(copy(firefly), distance_matrix, r) for _ in 1:k]
end


function hamming_distance(firefly_1::Firefly, firefly_2::Firefly)
    """Calculates the hamming distance between two fireflies.
    """
    dist = 0
    dist_info = []
    for (i, n) in enumerate(firefly_1.path)
        if n.id != firefly_2.path[i].id
            dist += 1
            push!(dist_info, true)
        else
            push!(dist_info, false)
        end
    end
    return dist, dist_info
end


function haversine_distance(n1::Node, n2::Node)::Float32
    #converting decimal degree to radians
    lat1=n1.y
    long1=n1.x
    lat2=n2.y
    long2=n2.x
    lat1,long1,lat2,long2=map(deg2rad,[lat1,long1,lat2,long2])
    
    # Using Haversine formula for distance calculation
    dlat=lat2-lat1
    dlong=long2-long1
    a=sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlong/2)^2
    c=2 * asin(sqrt(a))
    r=6371 #radius of earth in km
    return c*r
end

function init_firefly_paths(fireflies::Array{Firefly}, distance_matrix)
    """Initializes random paths for the given firefly array.
    """
    println("inside init_firefly_paths")
    for f in fireflies
        shuffle!(f.path)
        f.cost = path_cost(f, distance_matrix)
    end
end

function inversion_mutation(firefly::Firefly, distance_matrix, r::Int)::Firefly
    """Reverses the random part of given firefly path with
    given difference.
    """
    length_of_mutation = rand(2:r);
    max_len = length(firefly.path)-length_of_mutation
    index1 = rand(1:max_len)
    index2 = index1 + length_of_mutation
    mutated_path = reverse(firefly.path, index1, index2)
    new_f = copy(firefly)
    new_f.path = mutated_path
    new_f.cost = path_cost(new_f, distance_matrix)
    return new_f
end

function light_intensity(source_node::Firefly, destination_node::Firefly, λ::Float32)::Float32
    """Calculates the light intensity for the source firefly.
    """
    I_0 = source_node.cost
    r, _ = hamming_distance(source_node, destination_node)
    return I_0 * ℯ ^ (-λ * (r^2))
end

function move_firefly(f1::Firefly, f2::Firefly, r::Int)
    """Moving f1 to f2 less than r swap operations.
    """
    number_of_swaps = Int(ceil(rand() * (r-2) ))
    d, d_info = hamming_distance(f1, f2)

    while number_of_swaps > 0
        d, d_info = hamming_distance(f1, f2)
        random_index = rand([i for i in 1:length(d_info) if d_info[i]])
        value_to_copy = f2.path[random_index]
        index_to_move = findfirst(x -> x==value_to_copy, f1.path)

        if number_of_swaps == 1 && f1.path[index_to_move] == f2.path[random_index] && f1.path[random_index] == f2.path[index_to_move]
            break
        end

        f1.path[random_index], f1.path[index_to_move] = f1.path[index_to_move], f1.path[random_index]
        # println("swapped! ", f1.path[random_index], " -- ", f1.path[index_to_move])
        if f1.path[index_to_move] == f2.path[index_to_move] number_of_swaps -= 1 end
        number_of_swaps -= 1
    end
end

function path_cost(firefly::Firefly, distance_matrix)::Float32
    """Calculates the total path cost from given distance matrix.
    """
    path = firefly.path
    path_len = length(path) - 1
    p_cost = distance_matrix[path[end].id, path[1].id]
    for i in range(1, length=path_len)
        p_cost +=  distance_matrix[path[i].id, path[i+1].id]
    end

    return p_cost
end

function read_nodes(cluster_info)::Array{Node}
    """Reading nodes from giving file path.
    """
    nodes = Node[]
    for coord in cluster_info
            id = parse(Int16, coord[1])
            x_coord = parse(Float32, coord[2])
            y_coord = parse(Float32, coord[3])
            push!(nodes, Node(id, x_coord, y_coord))
    end
    return nodes
end

function plot_firefly(best_x, best_y, x, y, x_cost, best_cost, file_name, tsp_file_name)
    
    plot(best_x, best_y, "-")
    plot(x, y, "ro", markersize=2.0)
    xlabel("X Dimension")
    ylabel("Y Dimension")
    title("Firefly Cost: $best_cost", fontsize=10)
    suptitle("File: $file_name", fontweight="bold")
    fig_name = "graphs/best_firefly_path_$file_name-cost_$best_cost.png"
    savefig(fig_name)
    #println("Saved figure:",fig_name)
    clf()
    plot(x_cost)
    grid(true)
    ylabel("Path Cost")
    xlabel("Iteration Number")
    title("Firefly Iteration Graphic", fontsize=12)
    fig_name = "graphs/iteration_$tsp_file_name-cost_$best_cost.png"
    savefig(fig_name)
    clf()
end
    
function main(cluster_data, name)
    final_route=[]
    
    for ncluster in 1:length(cluster_data)
    
        file_name=name[ncluster]
        nodes = read_nodes(cluster_data[ncluster])
        distance_matrix = create_distance_matrix(nodes)
        # Checking if logs direcitory is available!
        # If is not create one!
        if !isdir("logs") mkdir("logs"); end
        
        tsp_file_name = file_name
        io_log = open("logs/$tsp_file_name.txt", "w+")

        # Setting Logger!
        logger = SimpleLogger(io_log)
        global_logger(logger)

        # Setting hyperparameters!
        LIGHT_ABSORPTION_COEFF = 0.2
        ATTRACTION_COEFF = 1
        ITERATION_NUMBER = 100
        POPULATION_NUMBER = 150
        NUMBER_OF_MUTATION = 10

        # Main firefly algorithm
        bests = []
        pop = [Firefly(copy(nodes), -1.0) for _ in 1:POPULATION_NUMBER]
        init_firefly_paths(pop, distance_matrix)
        println(length(pop), " Sized population created!")
        for t in 1:ITERATION_NUMBER
            new_pop = []
            for f1 in pop
                for f2 in pop
                    if f1.cost < f2.cost 
                        mutated_f1s = generate_mutated_fireflies(f2, distance_matrix, arc_distance(f1, f2), NUMBER_OF_MUTATION)
                        new_pop = [new_pop ; mutated_f1s]
                        push!(new_pop, f2)
                    end
                end
            end
            [push!(pop, n) for n in new_pop]
            sorted_pop = sort!(pop, by=p->p.cost)
            while length(pop) != POPULATION_NUMBER
                pop!(pop)
            end
            push!(bests, sorted_pop[1])
            current_best_cost = copy(bests[end].cost)
            now = string(Dates.now())
            println("Step:",t, " - Best ", bests[end])
            @info("Step:$t, Population Size:$POPULATION_NUMBER, Mutation Number:$NUMBER_OF_MUTATION, Date:$now",
                  current_best_cost)
        end

        close(io_log)

        # Checking if graphs direcitory is available for saving graphs!
        # If is not create one!
        if !isdir("graphs") mkdir("graphs"); end

        # Plotting and printing the best solutions.
        best = sort(pop, by=p->p.cost)[1]
        best_x = [n.x for n in best.path]
        best_y = [n.y for n in best.path]
        best_id= [n.id for n in best.path]
        best_cost = best.cost
        x = [n.x for n in nodes]
        y = [n.y for n in nodes]
        x_cost = [b.cost for b in bests]
        plot_firefly(best_x, best_y, x, y, x_cost, best_cost, file_name, tsp_file_name)
        z=[best_id, best_x, best_y]
        push!(final_route, z)
    end
    return final_route
end