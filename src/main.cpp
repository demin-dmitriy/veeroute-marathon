#ifdef ONLINE_JUDGE
    #pragma GCC optimize ("O3")
    #define NDEBUG
#endif

#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include <ctime>
#include <cassert>

constexpr size_t MAX_WORKERS_REQUIRED = 7;
constexpr size_t MAX_LOCATION_COUNT = 2000;
constexpr double TIME_LIMIT_SECONDS = 15;

struct Timer
{
    clock_t start_time;

    Timer()
            : start_time(clock())
    { }

    double seconds_elapsed() const
    {
        return double(clock() - start_time) / CLOCKS_PER_SEC;
    }
};

const Timer timer;

struct Vector
{
    int x;
    int y;

    int l1_norm() const
    {
        return abs(x) + abs(y);
    }
};

Vector operator-(Vector a, Vector b)
{
    return Vector
    {
        .x = a.x - b.x,
        .y = a.y - b.y
    };
}

int distance(Vector a, Vector b)
{
    return (a - b).l1_norm();
}

struct TimeWindow
{
    int from;
    int to;
};

struct Location
{
    Vector point;
    int duration;
    int workers_required;
    TimeWindow time_window;

    static Location read(std::istream& in)
    {
        Location result;
        in >> result.point.x
           >> result.point.y
           >> result.duration
           >> result.workers_required
           >> result.time_window.from
           >> result.time_window.to;
        return result;
    }
};

struct Task
{
    int n;
    std::vector<Location> locations;

    explicit Task(std::istream& in)
    {
        in >> n;
        locations.reserve(n);
        for (int i = 0; i != n; ++i)
        {
            locations.emplace_back(Location::read(in));
        }
    }
};

struct Vertex;

struct Edge
{
    Vertex* to;
    int earliest_arrive_moment;
};

struct Vertex
{
    int earliest_work_start_moment;
    const Location* const location;
    Vertex* const twin;
    std::vector<Edge> edges;

    explicit Vertex(Vertex* twin, const Location* location)
        : location(location)
        , twin(twin)
    {
        edges.reserve(location->workers_required);
    }
};

struct FullVertex
{
    Vertex forward;
    Vertex backward;

    explicit FullVertex(const Location* location)
        : forward(&backward, location)
        , backward(&forward, location)
    { }
};

struct Graph
{
    enum : size_t
    {
        BACKWARD_BASE = 0,
        FORWARD_BASE = 1
    };

    std::vector<FullVertex> vertices;

    explicit Graph(const Task& task)
    {
        vertices.reserve(task.n + 1);
        vertices.emplace_back(&task.locations[0]);

        for (const Location& location : task.locations)
        {
            vertices.emplace_back(&location);
        }
    }
};

struct Processor
{
    Task task;
    Graph graph;

    explicit Processor(std::istream& input)
        : task(input)
        , graph(task)
    { }


};

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    #ifdef READ_TASK_FROM_FILE
        assert(argc == 2);
        std::ifstream input(argv[1]);
    #else
        std::ifstream& input = std::cin;
    #endif

    Processor processor(input);


    #ifndef ONLINE_JUDGE
        std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds.\n";
    #endif
}
