#define OMPI_SKIP_MPICXX //I recommend not using C++ bindings
#include <mpich/mpi.h>
#include <omp.h>
#include <array>
#include <memory>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <cstring>
#include <chrono>

#define parallel

std::vector<std::size_t> g_node_counts;

std::uint32_t GetNodeId()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

std::uint32_t GetNumNodes()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

bool MasterNode()
{
    return GetNodeId() == 0u;
}

void Barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

std::uint32_t GetNodeIdFromElementId(std::size_t id)
{
    std::size_t count = 0;
    for (auto i = 0u; i < g_node_counts.size(); ++i)
    {
        count += g_node_counts[i];
        if (id < count)
        {
            return i;
        }
    }
    assert(0);
    return -1;
}

std::size_t GetLocalElementId(std::size_t global_id)
{
    std::size_t result = global_id;
    for (auto i = 1u; i <= GetNodeIdFromElementId(global_id); ++i)
    {
        result -= g_node_counts[i - 1];
    }

    return result;
}

void heapify(std::uint32_t * arr, int n)
{
#pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        // if child is bigger than parent
        if (arr[i] > arr[(i - 1) / 2])
        { 
            int j = i;

            // swap child and parent until
            // parent is smaller
            while (arr[j] > arr[(j - 1) / 2])
            { 
                std::swap(arr[j], arr[(j - 1) / 2]);
                j = (j - 1) / 2;
            }
        }
    }
}

template <typename T>
T getDigit(const T &value_r, const unsigned int &digit_r)
{
    return (value_r & (0xFF << ((sizeof(T) - digit_r - 1) * __CHAR_BIT__))) >> ( (sizeof(T) - digit_r - 1) * __CHAR_BIT__ );
}

void radixSort(std::uint32_t * sourceArray, std::vector< unsigned int> & aux,  long left_r, long right_r, const unsigned int digit_r)
{
    typedef unsigned int elem_t;

    const unsigned int num_of_digits = sizeof(elem_t);
    const unsigned int max_value = 256;

    std::vector< unsigned int > count(max_value + 1);

    if (digit_r > num_of_digits - 1 || left_r > right_r)
    {
        return;
    }
    for (unsigned int j = 0; j < max_value; j++)
    {
        count[j] = 0;
    }
    for (unsigned int i = left_r; i <= right_r; i++)
    {
        unsigned int j = getDigit(sourceArray[i], digit_r);
        count[j + 2]++;
    }
    for (unsigned int r = 0; r <= max_value-1; r++)
    {
        count[r+1] += count[r];
    }
    for (unsigned int i = left_r; i <= right_r; i++)
    {
        unsigned int j = getDigit(sourceArray[i], digit_r);
        aux[count[j+1]++] = sourceArray[i];
    }
    for (unsigned int i = left_r; i <= right_r; i++)
    {
        sourceArray[i] = aux[i - left_r];
    }

    // #pragma omp parallel shared(sourceArray, aux)
    {
    for (unsigned int r = 0; r < max_value; r++)
    {
        radixSort(sourceArray, aux, left_r + count[r], left_r + count[r+1] - 1, digit_r + 1);
    }
    }
}

void sort_1(
    std::uint32_t * data,
    int local_size,
    int global_size
)
{
    // Build heap (rearrange array)
    // heapify(data, local_size);

    std::vector< std::uint32_t > aux(2*local_size);
    std::fill(aux.begin(), aux.end(), 0);
    radixSort(data, aux, 0, local_size-1, 0);
    std::reverse(&data[0], &data[local_size]);

    Barrier();

    static const std::uint32_t sorted_marker = ~0u;

    std::size_t num_tree_elements = local_size;

    // Position where sorted part of the data starts
    for (std::size_t sorted_position = global_size - 1; sorted_position > 0; --sorted_position)
    {
        // Node id where sorted data starts
        std::uint32_t sorted_node_id = GetNodeIdFromElementId(sorted_position);

        // If number of tree elements is zero, the node is sorted
        bool is_node_sorted = (num_tree_elements == 0);

        // Gather roots
        std::vector<std::uint32_t> roots(GetNumNodes());

        // If a node is sorted, send a special marker that will be ignored
        std::uint32_t send_root = is_node_sorted ? sorted_marker : data[0];
        MPI_Allgather(&send_root, 1, MPI_UNSIGNED, roots.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);

        // Find a node with max root
        std::uint32_t max_root_node_id = std::max_element(roots.begin(), roots.end(),
            [](std::uint32_t a, std::uint32_t b)
        {
            if (a == sorted_marker) return true;
            else if (b == sorted_marker) return false;
            else return a < b;
        }) - roots.begin();

        Barrier();

        // Check if we need to swap data between nodes
        if (max_root_node_id != sorted_node_id)
        {
            std::uint32_t recv_data;

            // Here we need to swap values:
            // data[0] from node with max_root_node_id
            // And data[GetLocalElementId(sorted_position)] from node with sorted_node_id

            // 1. Send data[0] from max_root_node_id to sorted_node_id
            if (GetNodeId() == max_root_node_id)
            {
                MPI_Send(&data[0], 1, MPI_UNSIGNED, sorted_node_id, 0, MPI_COMM_WORLD);
            }
            else if (GetNodeId() == sorted_node_id)
            {
                MPI_Status status;
                MPI_Recv(&recv_data, 1, MPI_UNSIGNED, max_root_node_id, 0, MPI_COMM_WORLD, &status);
            }

            // 2. Send data[GetLocalElementId(sorted_position)] form sorted_node_id to max_root_node_id
            if (GetNodeId() == sorted_node_id)
            {
                MPI_Send(&data[GetLocalElementId(sorted_position)], 1, MPI_UNSIGNED, max_root_node_id, 0, MPI_COMM_WORLD);
            }
            else if (GetNodeId() == max_root_node_id)
            {
                MPI_Status status;
                MPI_Recv(&recv_data, 1, MPI_UNSIGNED, sorted_node_id, 0, MPI_COMM_WORLD, &status);
            }

            // 3. Swap recv_data and data[GetLocalElementId(sorted_position)] on sorted_node_id
            if (GetNodeId() == sorted_node_id)
            {
                std::swap(recv_data, data[GetLocalElementId(sorted_position)]);
                // Decrease the number of tree elements
                --num_tree_elements;
                // And heapify remaining values inside the node
                heapify(data, num_tree_elements);

                // std::fill(aux.begin(), aux.end(), 0);
                // radixSort(data, aux, 0, num_tree_elements-1, 0);
                // std::reverse(&data[0], &data[num_tree_elements]);
            }

            // 4. Swap recv_data and data[0] on max_root_node_id
            if (GetNodeId() == max_root_node_id)
            {
                std::swap(recv_data, data[0]);
                // And heapify
                heapify(data, num_tree_elements);

                // std::fill(aux.begin(), aux.end(), 0);
                // radixSort(data, aux, 0, num_tree_elements-1, 0);
                // std::reverse(&data[0], &data[num_tree_elements]);
            }
        }
        else if (GetNodeId() == max_root_node_id)
        {
            // Swap the root with the last sorted position
            std::swap(data[0], data[GetLocalElementId(sorted_position)]);
            // Decrease the number of tree elements
            --num_tree_elements;
            // And heapify remaining values inside the node
            heapify(data, num_tree_elements);

            // std::fill(aux.begin(), aux.end(), 0);
            // radixSort(data, aux, 0, num_tree_elements-1, 0);
            // std::reverse(&data[0], &data[num_tree_elements]);
        }

    }

}

void sort(
    std::uint32_t * data,
    int local_size,
    int global_size
)
{
    char * filename = "sorted.txt";

    std::ofstream file;
    // MPI_File file;

    MPI_Status status;

    if (MasterNode()) {
        file.open(filename, std::ios::out | std::ios::app);
        // int result = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        // if (result) {
        //     assert("\nUnable to open file\n");
        //     return;
        // }
    }

    // Build heap (rearrange array)
    // heapify(data, local_size);

    std::vector< std::uint32_t > aux(2*local_size);
    std::fill(aux.begin(), aux.end(), 0);
    radixSort(data, aux, 0, local_size-1, 0);

    static const std::uint32_t sorted_marker = ~0u;

    std::size_t num_tree_elements = 0;

    // Barrier();

    // // Position where sorted part of the data starts
    // for (std::size_t sorted_position = global_size; sorted_position > 0; --sorted_position)
    // {
    //     // If number of tree elements is zero, the node is sorted
    //     bool is_node_sorted = (num_tree_elements == local_size);

    //     // Gather roots
    //     std::vector<std::uint32_t> roots(GetNumNodes());

    //     // If a node is sorted, send a special marker that will be ignored
    //     std::uint32_t send_root = is_node_sorted ? sorted_marker : data[num_tree_elements];
    //     MPI_Allgather(&send_root, 1, MPI_UNSIGNED, roots.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);

    //     // Find a node with max root
    //     std::uint32_t min_root_node_id = std::min_element(roots.begin(), roots.end(),
    //         [](std::uint32_t a, std::uint32_t b)
    //     {
    //         if (a == sorted_marker) return false;
    //         else if (b == sorted_marker) return true;
    //         else return a < b;
    //     }) - roots.begin();


    //     if (GetNodeId() == min_root_node_id && num_tree_elements < local_size) {

    //         ++num_tree_elements;

    //     }

    //     if (MasterNode()) {
    //         file << roots[min_root_node_id] << std::endl;
    //         // std::string value = std::to_string(roots.at(min_root_node_id)) + '\n';
    //         // MPI_File_write(file, value.c_str(), value.length(), MPI_CHARACTER, &status);
    //     }

    // }

    if (MasterNode()) {
        // MPI_File_close(&file);
        file.close();
    }

}

auto GenerateData(std::size_t data_size)
{
    auto data = std::unique_ptr< std::uint32_t[] >(new std::uint32_t[data_size]);

    srand(GetNodeId());
    std::generate_n(data.get(), data_size, []()
    {
        return rand();
    });

    return std::move(data);
}

void GenerateDataFile(std::string filename, std::size_t data_size) {
    std::ofstream file;
    std::size_t rand_max = 4294967295;

    // srand(1234);

    file.open(filename);

    for (auto i = 0u; i < data_size; ++i) {
        file << rand() % rand_max << std::endl; // 4294967295 is max for std::uint_32t
    }

    file.close();
}

auto getDataFromFile(std::string filename) {
    std::vector< std::size_t > ranges(256);

    std::ofstream outfile;

    std::string outfilename =  "node_output_" + std::to_string(GetNodeId()) + ".txt";
    outfile.open(outfilename);

    std::vector< std::uint32_t > local_data;

    auto num_nodes = GetNumNodes();
    std::size_t delimitor = 255 / num_nodes;


    for (auto i = 0u; i < 256; ++i) {
        std::size_t value = i / delimitor;
        ranges[i] = (value > num_nodes - 1 ) ? num_nodes - 1 : value;
    }

    std::ifstream file;
    file.open(filename);

    // MPI_File file;
    // if ( !MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file) ) {
    //     std::cerr << "Can't open file from: " << GetNodeId() << std::endl;

    //     return local_data;
    // }

    auto node_id = GetNodeId();
    std::uint32_t value;

    while (true) {
        if (file.eof()) {
            break;
        }

        file >> value;

        // std::cout << "Node: " << node_id << ", MSD: " << getDigit(value, 0) << ", ranges[" << getDigit(value, 0) << "]:" << ranges[getDigit(value, 0)] << std::endl;
        if (node_id == ranges[getDigit(value, 0)] && !file.eof() ) {
            local_data.push_back(value);
        }

    }

    file.close();

    for (auto i = 0u; i < ranges.size(); ++i) {
        outfile << "Ranges[" << i << "]: " << ranges[i] << std::endl;
    }
    for (auto i = 0u; i < local_data.size(); ++i) {
        outfile << "MSD[" << i << "]: " << getDigit(local_data[i], 0) << std::endl; 
    }

    std::cout << "Node: " << GetNodeId() << " got " << local_data.size() << " values" << std::endl;

    outfile.close();

    return local_data;
}

void WriteDataToFile(std::string const& filename, std::uint32_t * data, std::size_t count)
{
    std::ofstream file;

    if (MasterNode())
    {
        file.open(filename);
    }

    Barrier();

    if (!MasterNode())
    {
        file.open(filename, std::ios::out | std::ios::app);
    }

    Barrier();

    for (auto i = 0u; i < GetNumNodes(); ++i)
    {
        if (GetNodeId() == i)
        {
            std::cout << "Node: " << GetNodeId() << " started writing values to file" << std::endl;
            for (auto j = 0u; j < count; ++j)
            {
                file << data[j] << std::endl;
            }
        }

        Barrier();
    }
}

void WriteToFileMPI(std::string const& filename, std::uint32_t * data, std::size_t count) {
    MPI_File file;

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL, &file);

    // if (MasterNode()) {
    //     if (!MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL, &file)) {
    //         assert("Unable to open file");
    //         return;
    //     }
    // } else {
    //     if (!MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_APPEND  | MPI_MODE_RDWR, MPI_INFO_NULL, &file)) {
    //         assert("Unable to open file");
    //         return;
    //     }
    // }

    Barrier();
    std::cout << "Node " << GetNodeId() << " opened file" << std::endl;
    Barrier();

    MPI_Status status;

    for (std::size_t i = 0u; i < GetNumNodes(); ++i) {
        Barrier();
        if (GetNodeId() == i) {
            for (auto j = 0; j < count; ++j) {
                std::string value = std::to_string(data[j]) + '\n';
                MPI_File_write_shared(file, value.c_str(), value.length(), MPI_CHAR, &status);
            }
        }
        Barrier();
    }

    MPI_File_close(&file);
}

std::pair< bool, bool > is_sorted(const std::uint32_t * begin, const std::uint32_t * end, const MPI_Comm comm)
{
    auto mpi_type = MPI_INT;
    // Массив на узле отсортирован?
    auto inner_result = std::is_sorted(begin, end);
    auto num_nodes = GetNumNodes();

    auto outer_result = true;
    {
        // Собираем результаты о сортировках
        auto result = static_cast<std::uint32_t>(inner_result);
        auto others_results = std::unique_ptr< std::uint32_t[] >(new std::uint32_t[num_nodes]);
        MPI_Allgather(&result, 1, mpi_type, others_results.get(), 1, mpi_type, comm);
        auto others_total = std::all_of(others_results.get(), others_results.get() + num_nodes,
            [](const auto lhs)
        {
            return lhs != 0u;
        });
        outer_result &= others_total;
    }

    {
        constexpr auto mm_size = 2;
        // Минимум и максимум
        auto mm = std::array< std::uint32_t, mm_size >{*begin, *(begin + (std::distance(begin, end) - 1))};
        // Собираем данные о максимумах и минимумов со всех узлов
        auto others_mm = std::unique_ptr< std::uint32_t[] >(new std::uint32_t[num_nodes * mm_size]);
        MPI_Allgather(mm.data(), mm_size, mpi_type, others_mm.get(), mm_size, mpi_type, comm);
        // Если пары [min max] [min max] ... упорядочены,
        // то считаем данные отсортированными
        auto others_mm_total = std::is_sorted(others_mm.get(), others_mm.get() + num_nodes * mm_size);
        outer_result &= others_mm_total;
    }
    return { inner_result, outer_result };
}

int isSorted() {
	std::ifstream myfile("sorted.txt");
	if (myfile.is_open())
	{
		unsigned long long int arrSize = 0;

		std::uint32_t prev, cur = 0;
	    myfile >> prev;
		cur = prev;
		++arrSize;

		while (true)
		{
            if (myfile.eof())
				break;

	        myfile >> cur;

            if (cur < prev) {
                std::cout << "Unsorted at " << arrSize << std::endl;
                return 0;
            }

			prev = cur;
		}
    // I should have closed the file here, but as the program was ending I was lazy	
    }
	else
	{
		std::cout << "\nUnable to open file\n";
        return 0;
	}

	std::cout << "Data is sorted" << std::endl;
	return 0;
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "Please, read README before using" << std::endl;
        return 0;
    }
    auto data_size = atoi(argv[1]);
    
    GenerateDataFile("unsorted.txt", data_size);

    auto provided = 0;
    // MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Init(&argc, &argv);

    auto num_nodes = GetNumNodes();
    auto node_id = GetNodeId();

    g_node_counts.resize(num_nodes);

    if (argc != 2)
    {
        if (node_id == 0)
        {
            std::cerr << "Please, read README before using" << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    // Get total size from cmd line
    auto global_size = atoi(argv[1]);

    std::vector< uint32_t > local_data_vector = getDataFromFile("unsorted.txt");
    const std::size_t local_size = local_data_vector.size();
    std::uint32_t * data = local_data_vector.data();

    // Calculate data size per node
    // const std::size_t local_size = (global_size / num_nodes) + ((node_id != (num_nodes - 1)) ? 0 : (global_size % num_nodes));
    Barrier();

    MPI_Allgather(&local_size, sizeof(std::size_t), MPI_CHAR, g_node_counts.data(), sizeof(std::size_t), MPI_CHAR, MPI_COMM_WORLD);

    // Generate data
    // auto data = GenerateData(local_size);

    // Issue a barrier
    Barrier();

    // WriteDataToFile("unsorted.txt", data.get(), local_size);

    std::chrono::high_resolution_clock::time_point start_time =
        std::chrono::high_resolution_clock::now();

    // Sort
    sort(data, local_size, global_size);

    // Issue a barrier
    Barrier();

    if (MasterNode())
    {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start_time).count();
        std::cout << "Execution of sort took " << elapsed << " s." << std::endl;
    }

    // WriteDataToFile("sorted.txt", data, local_size);
    WriteToFileMPI("sorted.txt", data, local_size);

    if (MasterNode()) {
        isSorted();
    }

    MPI_Finalize();

    return 0;
}
