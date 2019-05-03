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

using Moment = int;

constexpr size_t MAX_WORKERS_REQUIRED = 7;
constexpr size_t MAX_LOCATION_COUNT = 2000;
constexpr Moment MAX_MOMENT = 1000;
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

template <class T, size_t MAX_SIZE>
struct CheapStack
{
    size_t size = 0;
    std::array<T, MAX_SIZE> elements;

    void push(T value)
    {
        elements[size] = std::move(value);
        size += 1;
        assert(size < MAX_SIZE);
    }

    T pop()
    {
        assert(size > 0);
        size -= 1;
        return std::move(elements[size]);
    }

    void clear()
    {
        size = 0;
    }
};

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

Moment twin_moment(Moment moment)
{
    return MAX_MOMENT - moment;
}

struct TimeWindow
{
    Moment from;
    Moment to;

    TimeWindow twin() const
    {
        return TimeWindow
        {
            .from = twin_moment(from),
            .to = twin_moment(to)
        };
    }
};

struct Location
{
    Vector point;
    int duration;
    int workers_required;
    TimeWindow time_window;
    int index;

    static Location read(std::istream& in, int index)
    {
        Location result;
        in >> result.point.x
           >> result.point.y
           >> result.duration
           >> result.workers_required
           >> result.time_window.from
           >> result.time_window.to;
        result.index = index;
        return result;
    }

    Location twin() const
    {
        return Location
        {
            .point = point,
            .duration = duration,
            .workers_required = workers_required,
            .time_window = time_window.twin(),
            .index = index
        };
    }
};

struct FullLocation
{
    Location forward;
    Location backward;

    static FullLocation read(std::istream& in, int index)
    {
        FullLocation result;

        result.forward = Location::read(in, index);
        result.backward = result.forward.twin();

        return result;
    }
};

int distance(const Location& a, const Location& b)
{
    return distance(a.point, b.point);
}

struct Task
{
    int n;
    std::vector<FullLocation> locations;

    explicit Task(std::istream& in)
    {
        in >> n;
        locations.reserve(n);
        for (int i = 0; i != n; ++i)
        {
            locations.emplace_back(FullLocation::read(in, i + 1));
        }
    }
};

struct Vertex;

struct Edge
{
    Vertex* to;
    Edge* twin;
    Moment earliest_arrive_moment;
};

struct FullVertex;
int distance(const Vertex& a, const Vertex& b);

struct Vertex
{
    using Version = unsigned int; // Won't overflow in 15 seconds. Could overflow, if it were running longer.

    static inline CheapStack<Vertex*, 2 * (MAX_LOCATION_COUNT + 1)> marked; // Store both forward and backward vertices
    static inline Version current_global_version = 1;

    Moment earliest_work_start_moment; // TODO: maybe end_moment is more convenient to work with
    Version version;
    const Location* const location;
    FullVertex* const parent;
    Vertex* const twin;
    std::vector<Edge> edges;

    explicit Vertex(FullVertex* parent, Vertex* twin, const Location* location)
        : earliest_work_start_moment(-1)
        , version(0)
        , location(location)
        , parent(parent)
        , twin(twin)
    {
        edges.reserve(location->workers_required);
    }

    static void unmark_all()
    {
        current_global_version += 1;
        marked.clear();
    }

    void mark()
    {
        if (version != current_global_version)
        {
            version = current_global_version;

            for (const Edge& edge : edges)
            {
                edge.to->mark();
            }

            marked.push(this);
        }
    }

    void recalculate()
    {
        earliest_work_start_moment = location->time_window.from;

        for (const Edge& twin_edge : twin->edges)
        {
            earliest_work_start_moment = std::max(earliest_work_start_moment, twin_edge.twin->earliest_arrive_moment);
        }

        Moment earliest_work_end_moment = earliest_work_start_moment + location->duration;

        assert(earliest_work_end_moment <= location->time_window.to);

        for (Edge& edge : edges)
        {
            // TODO: this computation might be needed when new edge is added
            edge.earliest_arrive_moment = earliest_work_end_moment + distance(*this, *edge.to);
        }
    }

    static void recalculate_all_marked()
    {
        while (marked.size != 0)
        {
            // TODO: early stopping ... but how exactly?
            marked.pop()->recalculate();
        }
    }
};

int distance(const Vertex& a, const Vertex& b)
{
    return distance(*a.location, *b.location);
}

struct FullVertex
{
    Vertex forward;
    Vertex backward;

    explicit FullVertex(const FullLocation* location)
        : forward(this, &backward, &location->forward)
        , backward(this, &forward, &location->backward)
    { }

    void mark()
    {
        forward.mark();
        backward.mark();
    }

    void recalculate()
    {
        forward.recalculate();
        backward.recalculate();
    }

    void insert_between(FullVertex& a, FullVertex& b)
    {
        // TODO
    }

    void connect(FullVertex& other)
    {
        // TODO
//        forward.edges.emplace_back({.to = other.forward });
    }
};

struct Graph
{
    enum : size_t
    {
        BACKWARD_BASE = 0,
        FORWARD_BASE = 1
    };

    const Task& task;
    std::vector<FullVertex> vertices;

    explicit Graph(const Task& task)
        : task(task)
    {
        vertices.reserve(task.n + 1);
        vertices.emplace_back(&task.locations[0]);

        for (const FullLocation& location : task.locations)
        {
            vertices.emplace_back(&location);
        }
    }
};

std::vector<Moment> calculate_true_last_arrive_moments(Graph& graph)
{
    assert(0 == Vertex::marked.size);

    std::vector<Moment> last_arrive_moments(graph.task.n + 1, -1);

    const auto update_edge = [&last_arrive_moments](const Vertex& from, const Edge& edge, Moment ready_moment)
    {
        const Moment arrive_moment = ready_moment + distance(from, *edge.to);
        const size_t to_index = edge.to->location->index;
        last_arrive_moments[to_index] = std::max(last_arrive_moments[to_index], arrive_moment);
    };

    Vertex& base = graph.vertices[Graph::FORWARD_BASE].forward;
    base.mark();

    // Base is special

    Vertex::marked.pop(); // Discard base
    for (const Edge& edge : base.edges)
    {
        update_edge(base, edge, twin_moment(edge.twin->earliest_arrive_moment));
    }

    // Note: we are skipping base here, hence -1.
    while (Vertex::marked.size != 0)
    {
        const Vertex& from = *Vertex::marked.pop();
        const Moment work_start_moment = last_arrive_moments[from.location->index];
        assert(work_start_moment != -1);
        const Moment work_end_moment = work_start_moment + from.location->duration;

        for (const Edge& edge : from.edges)
        {
            update_edge(from, edge, work_end_moment);
        }
    }

    Vertex::unmark_all();

    return last_arrive_moments;
}

std::ostream& operator<<(std::ostream& output, Graph& graph)
{
    std::vector<Moment> last_arrive_moment = calculate_true_last_arrive_moments(graph);
    std::vector<int> visit_count(graph.task.n + 1, 0);

    const Vertex* const base = &graph.vertices[Graph::FORWARD_BASE].forward;
    const Vertex* const backward_base = &graph.vertices[Graph::BACKWARD_BASE].forward;

    for (const Edge& first_edge: base->edges)
    {
        Moment current_moment = twin_moment(first_edge.twin->earliest_arrive_moment);
        output << "start " << current_moment << " 1\n";

        const Vertex* current_vertex = first_edge.to;
        current_moment += distance(*base, *first_edge.to);

        while (true)
        {
            const size_t index = current_vertex->location->index;

            output << "arrive " << current_moment << " " << index << "\n";

            if (current_vertex != backward_base)
            {
                break;
            }

            const Moment start_work_moment = last_arrive_moment[index];
            const Moment end_work_moment = start_work_moment + current_vertex->location->duration;

            assert(current_moment <= start_work_moment);

            output << "work " << start_work_moment << " " << end_work_moment << " " << index << "\n";

            const Vertex* const next_vertex = current_vertex->edges[visit_count[index]].to;

            current_moment = end_work_moment + distance(*current_vertex, *next_vertex);

            visit_count[index] += 1;
            current_vertex = next_vertex;
        }

        output << "end\n";
    }

    return output;
}

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

    std::cout << processor.graph;

    #ifndef ONLINE_JUDGE
        std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds.\n";
    #endif
}
