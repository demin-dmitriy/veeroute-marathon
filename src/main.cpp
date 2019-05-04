#ifdef ONLINE_JUDGE
    #pragma GCC optimize ("O3")
    #define NDEBUG
#endif

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <cassert>
#include <ctime>
#include <cmath>

constexpr double TIME_LIMIT_SECONDS = 15;

inline namespace lib
{
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

    template <class T>
    bool contains_ptr(const std::vector<T>& xs, const T* x)
    {
        return xs.data() <= x and x < xs.data() + xs.size();
    }

    template <class Value>
    struct Table
    {
        size_t size1;
        size_t size2;
        std::vector<Value> values;

        Table(size_t size1, size_t size2)
            : size1(size1)
            , size2(size2)
            , values(size1 * size2)
        { }

        Value& at(size_t i, size_t j)
        {
            assert(i < size1);
            assert(j < size2);
            return values[i * size2 + j];
        }

        const Value& at(size_t i, size_t j) const
        {
            assert(i < size1);
            assert(j < size2);
            return values[i * size2 + j];
        }
    };

    template <class Value>
    std::vector<size_t> solve_assignment_problem(const Table<Value>& table)
    {
        // Copied with small modifications from https://e-maxx.ru/algo/assignment_hungary
        const Value INF = 1000000;

        const size_t n = table.size1;
        const size_t m = table.size2;

        std::vector<Value> u(n + 1), v(m + 1);
        std::vector<size_t> p(m + 1), way(m + 1);
        std::vector<Value> minv(m + 1, INF);
        std::vector<char> used(m + 1, false);

        for (size_t i = 1; i <= n; ++i)
        {
            p[0] = i;
            size_t j0 = 0;

            if (i > 1)
            {
                fill(begin(used), end(used), false);
                fill(begin(minv), end(minv), INF);
            }

            do
            {
                used[j0] = true;
                size_t i0 = p[j0];
                size_t j1 = 0;
                Value delta = INF;

                for (size_t j = 1; j <= m; ++j)
                {
                    if (!used[j]) {
                        const Value cur = table.at(i0 - 1, j - 1) - u[i0] - v[j];
                        if (cur < minv[j])
                        {
                            minv[j] = cur;
                            way[j] = j0;
                        }
                        if (minv[j] < delta)
                        {
                            delta = minv[j];
                            j1 = j;
                        }
                    }
                }
                for (size_t j = 0; j <= m; ++j)
                {
                    if (used[j])
                    {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    }
                    else
                    {
                        minv[j] -= delta;
                    }
                }
                j0 = j1;
            } while (p[j0] != 0);

            do
            {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }

        return p;
    }

    template <class Range>
    auto reversed(Range& range)
    {
        struct ReversedRange
        {
            Range& range;

            explicit ReversedRange(Range& range)
                : range(range)
            { }

            auto begin() { return range.rbegin(); }
            auto begin() const { return range.rbegin(); }
            auto end() { return range.rend(); }
            auto end() const { return range.rend(); }
        };

        return ReversedRange(range);
    }
}

inline namespace task
{
    using Moment = int;
    using Index = size_t;

    constexpr size_t MAX_WORKERS_REQUIRED = 7;
    constexpr size_t MAX_LOCATION_COUNT = 2000;
    constexpr Moment MAX_MOMENT = 1000;

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
                .from = twin_moment(to),
                .to = twin_moment(from)
            };
        }

        bool is_valid() const
        {
            return from <= to;
        }
    };

    struct Location
    {
        Vector point;
        int duration;
        int workers_required;
        TimeWindow time_window;
        Index index;

        static Location read(std::istream& in, Index index)
        {
            Location result;
            in >> result.point.x
               >> result.point.y
               >> result.duration
               >> result.workers_required
               >> result.time_window.from
               >> result.time_window.to;
            result.index = index;

            if (result.is_base())
            {
                result.time_window.to = MAX_MOMENT;
            }

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

        bool is_base() const
        {
            return 1 == index;
        }

        bool is_valid() const
        {
            // NOTE: incomplete implementation until more sanity checks are needed.
            return time_window.is_valid();
        }
    };

    struct FullLocation
    {
        Location forward;
        Location backward;

        static FullLocation read(std::istream& in, Index index)
        {
            FullLocation result;

            result.forward = Location::read(in, index);
            result.backward = result.forward.twin();

            assert(result.forward.is_valid());
            assert(result.backward.is_valid());

            return result;
        }
    };

    int distance(const Location& a, const Location& b)
    {
        return distance(a.point, b.point);
    }

    struct Task
    {
        size_t n;
        std::vector<FullLocation> locations;

        explicit Task(std::istream& in)
        {
            in >> n;
            assert(n > 0);
            locations.reserve(n);

            for (size_t i = 0; i != n; ++i)
            {
                locations.emplace_back(FullLocation::read(in, i + 1));
            }
        }
    };
}

inline namespace graph
{
    struct Vertex;

    struct Edge
    {
        Vertex* to;
        Edge* twin;
        Moment earliest_arrive_moment;
    };

    struct FullVertex;
    int distance(const Vertex* a, const Vertex* b);

    struct Vertex
    {
        using Version = unsigned int; // Won't overflow in 15 seconds. Could overflow, if it was running longer.

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
            if (0 == location->workers_required) // This is base
            {
                edges.reserve(MAX_WORKERS_REQUIRED * MAX_LOCATION_COUNT); // TODO: is this cache-friendly?
            }
            else
            {
                edges.reserve(location->workers_required);
            }
        }

        void fix_twin_pointers_on_realloc()
        {
            for (Edge& edge : edges)
            {
                edge.twin->twin = &edge;
            }
        }

        static void unmark_all()
        {
            current_global_version += 1;
            marked.clear();
        }

        bool is_marked() const
        {
            return version == current_global_version;
        }

        void mark_recursively()
        {
            if (not is_marked())
            {
                version = current_global_version;

                for (const Edge& edge : edges)
                {
                    edge.to->mark_recursively();
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
                edge.earliest_arrive_moment = earliest_work_end_moment + distance(this, edge.to);
            }
        }

        static void recalculate_all_marked()
        {
            while (marked.size != 0)
            {
                // TODO: early stopping ... but how exactly?
                marked.pop()->recalculate();
            }
            current_global_version += 1;
        }
    };

    Edge* link(Vertex* a, Vertex* b)
    {
        const auto add_edge = [](Vertex* v) -> Edge*
        {
            assert(v->edges.size() < v->edges.capacity());
            v->edges.emplace_back();
            return &v->edges.back();
        };

        Edge* ab = add_edge(a);
        Edge* ba = add_edge(b->twin);

        ab->to = b;
        ab->twin = ba;

        ba->to = a->twin;
        ba->twin = ab;

        return ab;
    }

    void unlink(Edge* ab)
    {
        const auto remove_edge = [](Vertex* vertex, Edge* edge)
        {
            assert(contains_ptr(vertex->edges, edge));
            const size_t i = edge - vertex->edges.data();
            vertex->edges[i] = vertex->edges.back();
            vertex->edges.pop_back();
        };

        Edge* ba = ab->twin;

        remove_edge(ba->to->twin, ab);
        remove_edge(ab->to->twin, ba);

        ab->twin = ba;
        ba->twin = ab;
    }

    void interpose(Edge* ab, Vertex* v)
    {
        Vertex* a = ab->twin->to;
        Vertex* b = ab->to;

        unlink(ab);
        link(a, v);
        link(v, b);
    }

    void cut(Vertex* v)
    {
        // TODO: solving assignment problem, and getting vector of all immediate children can be useful on its own.
        assert(v->edges.size() == v->twin->edges.size());

        const size_t m = v->edges.size();

        std::vector<Vertex*> froms;
        std::vector<Vertex*> tos;

        froms.reserve(m);
        tos.reserve(m);

        for (const Edge& edge : v->edges)
        {
            tos.emplace_back(edge.to);
        }

        for (const Edge& edge : v->twin->edges)
        {
            froms.emplace_back(edge.to->twin);
        }

        Table<int> distances(m, m);

        for (size_t i = 0; i != m; ++i)
        {
            for (size_t j = 0; j != m; ++j)
            {
                distances.at(i, j) = distance(froms[i], tos[j]);
            }
        }

        std::vector<size_t> assignment = solve_assignment_problem(distances);

        for (Vertex* u : { v, v->twin })
        {
            for (Edge& edge : u->edges)
            {
                unlink(&edge);
            }
        }

        for (size_t i = 1; i <= m; ++i)
        {
            link(froms[assignment[i] - 1], tos[i - 1]);
        }
    }

    int distance(const Vertex* a, const Vertex* b)
    {
        return distance(*a->location, *b->location);
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
            forward.mark_recursively();
            backward.mark_recursively();
        }

        void recalculate()
        {
            forward.recalculate();
            backward.recalculate();
        }

        void copy_dynamics_from(const FullVertex& other)
        {
            forward.earliest_work_start_moment = other.forward.earliest_work_start_moment;
            backward.earliest_work_start_moment = other.backward.earliest_work_start_moment;
        }
    };

    struct Graph
    {
        enum : Index
        {
            BACKWARD_BASE = 0,
            FORWARD_BASE = 1
        };

        const Task* task;
        std::vector<FullVertex> vertices;

        explicit Graph(const Task* task)
            : task(task)
        {
            vertices.reserve(task->n + 1);
            vertices.emplace_back(&task->locations[0]);

            for (const FullLocation& location : task->locations)
            {
                vertices.emplace_back(&location);
            }
        }

        Graph(Graph&&) = default;

        Graph copy() const
        {
            Graph result(this->task);

            // Ugly, but whatever.
            for (size_t i = 0; i != vertices.size(); ++i)
            {
                const FullVertex& this_full_vertex = vertices[i];
                FullVertex& that_full_vertex = result.vertices[i];

                that_full_vertex.copy_dynamics_from(this_full_vertex);

                for (size_t j = 0; j != this_full_vertex.forward.edges.size(); ++j)
                {
                    const Edge& this_edge = this_full_vertex.forward.edges[j];
                    size_t to_index = this_edge.to->parent - vertices.data();

                    Edge* that_edge = link(&that_full_vertex.forward, &result.vertices[to_index].forward);
                    that_edge->earliest_arrive_moment = this_edge.earliest_arrive_moment;
                    that_edge->twin->earliest_arrive_moment = this_edge.twin->earliest_arrive_moment;
                }
            }

            return result;
        }

        void remove_empty_paths()
        {
            const Vertex* backward_base = &vertices[BACKWARD_BASE].forward;

            for (Edge& edge : reversed(vertices[FORWARD_BASE].forward.edges))
            {
                if (edge.to == backward_base)
                {
                    unlink(&edge);
                }
            }
        }

        Graph& operator=(Graph&&) = default;
        void operator=(const Graph&) = delete;
    };

    std::vector<Moment> calculate_true_last_arrive_moments(Graph& graph)
    {
        assert(0 == Vertex::marked.size);

        std::vector<Moment> last_arrive_moments(graph.task->n + 1, -1);

        const auto update_edge = [&last_arrive_moments](const Vertex* from, const Edge& edge, Moment ready_moment)
        {
            const Moment arrive_moment = ready_moment + distance(from, edge.to);
            const Index to_index = edge.to->location->index;
            last_arrive_moments[to_index] = std::max(last_arrive_moments[to_index], arrive_moment);
        };

        Vertex* base = &graph.vertices[Graph::FORWARD_BASE].forward;
        base->mark_recursively();

        // Base is special

        Vertex::marked.pop(); // Discard base
        for (const Edge& edge : base->edges)
        {
            update_edge(base, edge, twin_moment(edge.twin->earliest_arrive_moment));
        }

        // Note: we are skipping base here, hence -1.
        while (Vertex::marked.size != 0)
        {
            const Vertex* from = Vertex::marked.pop();
            const Moment work_start_moment = last_arrive_moments[from->location->index];
            assert(work_start_moment != -1);
            const Moment work_end_moment = work_start_moment + from->location->duration;

            for (const Edge& edge : from->edges)
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
        std::vector<size_t> visit_count(graph.task->n + 1, 0);

        const Vertex* const base = &graph.vertices[Graph::FORWARD_BASE].forward;
        const Vertex* const backward_base = &graph.vertices[Graph::BACKWARD_BASE].forward;

        for (const Edge& first_edge: base->edges)
        {
            Moment current_moment = twin_moment(first_edge.twin->earliest_arrive_moment);
            output << "start " << current_moment << " 1\n";

            const Vertex* current_vertex = first_edge.to;
            current_moment += distance(base, first_edge.to);

            while (true)
            {
                const Index index = current_vertex->location->index;

                output << "arrive " << current_moment << " " << index << "\n";

                if (current_vertex == backward_base)
                {
                    break;
                }

                const Moment start_work_moment = last_arrive_moment[index];
                const Moment end_work_moment = start_work_moment + current_vertex->location->duration;

                assert(current_moment <= start_work_moment);
                assert(current_vertex->location->time_window.from <= start_work_moment);
                assert(end_work_moment <= current_vertex->location->time_window.to);

                output << "work " << start_work_moment << " " << end_work_moment << " " << index << "\n";

                assert(visit_count[index] < current_vertex->edges.size());

                const Vertex* const next_vertex = current_vertex->edges[visit_count[index]].to;

                current_moment = end_work_moment + distance(current_vertex, next_vertex);

                visit_count[index] += 1;
                current_vertex = next_vertex;
            }

            output << "end\n";
        }

        return output;
    }

}

inline namespace solvers
{
    void generate_empty_routes(Graph& graph)
    {
        const double n = graph.task->n;
        const double average_workers_required = (1. + 8.) / 2.;
        const double average_duration = (5. + 30.) / 2.;
        const double average_gap = 100. / (1. + sqrt(n));
        const double work_time_window_length = 800. - 200.;

        const double expected_workers_required =
            (n * average_workers_required * (average_duration + average_gap) - average_gap) / work_time_window_length;

        const size_t empty_route_count = static_cast<size_t>(1.5 * expected_workers_required);

        for (size_t i = 0; i != empty_route_count; ++i)
        {
            link(&graph.vertices[Graph::FORWARD_BASE].forward, &graph.vertices[Graph::BACKWARD_BASE].forward);
        }
    }

    void greedy(Graph& graph)
    {

    }
}

struct Processor
{
    Task task;
    Graph graph;

    explicit Processor(std::istream& input)
        : task(input)
        , graph(&task)
    { }

    void run_circuit()
    {
        generate_empty_routes(graph);
        greedy(graph);
        graph.remove_empty_paths();
    }
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
    processor.run_circuit();
    std::cout << processor.graph;

    #ifndef ONLINE_JUDGE
        std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds.\n";
    #endif
}
