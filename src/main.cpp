#ifdef ONLINE_JUDGE
    #pragma GCC optimize ("O3")
    #define NDEBUG
#endif

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string_view>
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

    template <class V, class Func>
    auto map(const std::vector<V>& values, Func&& func) -> std::vector<decltype(func(std::declval<V>()))>
    {
        using U = decltype(func(std::declval<V>()));
        std::vector<U> result;
        result.reserve(values.size());

        for (const V& value : values)
        {
            result.emplace_back(func(value));
        }

        return result;
    }

    template <class T>
    bool contains_ptr(const std::vector<T>& xs, const T* x)
    {
        return xs.data() <= x and x < xs.data() + xs.size();
    }

    template <class T>
    struct Pool
    {
        size_t size;
        size_t next_free;
        std::vector<T> storage;
        std::vector<T*> free;

        explicit Pool(size_t size)
            : size(size)
            , next_free(0)
        {
            storage.resize(size);
            free.resize(size);

            for (size_t i = 0; i != size; ++i)
            {
                free[i] = storage.data() + i;
            }
        }

        T* allocate() noexcept
        {
            assert(next_free < free.size());
            return free[next_free++];
        }

        void deallocate(T* ptr) noexcept
        {
            assert(contains_ptr(storage, ptr));
            assert(next_free > 0);
            free[--next_free] = ptr;
        }
    };

    template <class T>
    struct Span
    {
        T* data;
        size_t size;

        template <class Container>
        explicit Span(Container& container)
            : data(container.data())
            , size(container.size())
        { }

        T* begin() { return data; }
        T* end() { return data + size; }
        const T* begin() const { return data; }
        const T* end() const { return data + size; }
    };

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

        bool is_empty() const
        {
            return from >= to;
        }

        int length() const
        {
            return to - from;
        }

        TimeWindow intersect(TimeWindow other) const
        {
            return TimeWindow
            {
                .from = std::max(from, other.from),
                .to = std::min(to, other.to)
            };
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
            return not time_window.is_empty() or is_base();
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

    int distance(const Location* a, const Location* b)
    {
        return distance(a->point, b->point);
    }

    struct Task
    {
        size_t n;
        std::vector<FullLocation> locations;

        size_t sum_workers_required;

        explicit Task(std::istream& in)
            : sum_workers_required(0)
        {
            in >> n;
            assert(n > 0);
            locations.reserve(n);

            for (size_t i = 0; i != n; ++i)
            {
                locations.emplace_back(FullLocation::read(in, i + 1));
                sum_workers_required += locations.back().forward.workers_required;
            }
        }
    };
}

inline namespace graph
{
    struct Vertex;

    struct Edge
    {
        Vertex* to
            #ifndef NDEBUG
                = nullptr
            #endif
            ;
        size_t index
            #ifndef NDEBUG
                = 800000
            #endif
            ;
        Edge* twin
            #ifndef NDEBUG
                = nullptr
            #endif
            ;
        Moment earliest_arrive_moment
            #ifndef NDEBUG
                = 10000
            #endif
            ;
    };

    struct FullVertex;
    int distance(const Vertex* a, const Vertex* b);

    struct Vertex
    {
        using Version = unsigned int; // Won't overflow in 15 seconds. Could overflow, if it was running longer.

        static inline CheapStack<Vertex*, 2 * (MAX_LOCATION_COUNT + 1)> marked; // Stores both forward and backward vertices
        static inline Version current_global_version = 1;

        Moment earliest_work_start_moment; // TODO: maybe end_moment is more convenient to work with
        Version version;
        const Location* const location;
        Vertex* const twin;
        std::vector<Edge*> edges;

        explicit Vertex(FullVertex* parent, Vertex* twin, const Location* location)
            : earliest_work_start_moment(-1)
            , version(0)
            , location(location)
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

                for (const Edge* edge : edges)
                {
                    edge->to->mark_recursively();
                }

                marked.push(this);
            }
        }

        Moment get_earliest_work_end_moment() const
        {
            return earliest_work_start_moment + location->duration;
        }

        void recalculate()
        {
            earliest_work_start_moment = location->time_window.from;

            for (const Edge* twin_edge : twin->edges)
            {
                earliest_work_start_moment = std::max(earliest_work_start_moment, twin_edge->twin->earliest_arrive_moment);
            }

            const Moment earliest_work_end_moment = get_earliest_work_end_moment();

            assert(earliest_work_end_moment <= location->time_window.to);

            for (Edge* edge : edges)
            {
                // TODO: this computation might be needed when new edge is added
                edge->earliest_arrive_moment = earliest_work_end_moment + distance(this, edge->to);
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

    int distance(const Vertex* a, const Vertex* b)
    {
        return distance(a->location, b->location);
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
        Pool<Edge> edge_pool;

        explicit Graph(const Task* task)
            : task(task)
            , edge_pool(4 * task->sum_workers_required)
        {
            assert(task->n >= 1);

            vertices.reserve(task->n + 1);
            vertices.emplace_back(&task->locations[0]);

            for (const FullLocation& location : task->locations)
            {
                vertices.emplace_back(&location);
            }

            vertices[FORWARD_BASE].recalculate();
            vertices[BACKWARD_BASE].recalculate();
        }

        Graph(Graph&&) = default;

        Vertex* forward_start()
        {
            return &vertices[FORWARD_BASE].forward;
        }

        Vertex* forward_finish()
        {
            return &vertices[BACKWARD_BASE].forward;
        }

        Vertex* backward_start()
        {
            return &vertices[BACKWARD_BASE].backward;
        }

        Vertex* backward_finish()
        {
            return &vertices[FORWARD_BASE].backward;
        }

        Edge* link(Vertex* a, Vertex* b)
        {
            const auto add_edge = [this](Vertex* v) -> Edge*
            {
                assert(v->edges.size() < v->edges.capacity());
                Edge* edge = edge_pool.allocate();
                edge->index = v->edges.size();
                v->edges.emplace_back(edge);
                return edge;
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
            const auto remove_edge = [this](Vertex* vertex, Edge* edge)
            {
                assert(contains_ptr(edge_pool.storage, edge));
                assert(edge->index < vertex->edges.size());
                assert(vertex->edges[edge->index] == edge);
                vertex->edges[edge->index] = vertex->edges.back();
                vertex->edges[edge->index]->index = edge->index;
                vertex->edges.pop_back();
                edge_pool.deallocate(edge);
            };

            Edge* ba = ab->twin;

            remove_edge(ba->to->twin, ab);
            remove_edge(ab->to->twin, ba);
        }

        void interpose(Edge* ab, Vertex* c)
        {
            Vertex* a = ab->twin->to->twin;
            Vertex* b = ab->to;

            unlink(ab);
            Edge* ac = link(a, c);
            Edge* cb = link(c, b);

            ac->earliest_arrive_moment = a->get_earliest_work_end_moment() + distance(a, c);
            cb->twin->earliest_arrive_moment = b->twin->get_earliest_work_end_moment() + distance(b, c);

            assert(c->edges.size() == c->twin->edges.size());
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

            for (const Edge* edge : v->edges)
            {
                tos.emplace_back(edge->to);
            }

            for (const Edge* edge : v->twin->edges)
            {
                froms.emplace_back(edge->to->twin);
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
                for (Edge* edge : u->edges)
                {
                    unlink(edge);
                }
            }

            for (size_t i = 1; i <= m; ++i)
            {
                link(froms[assignment[i] - 1], tos[i - 1]);
            }
        }

        Graph copy() const;

        void dump_graphviz(std::string_view filename = "graph.dot") const
        {
            std::ofstream out(filename.data());

            const auto print = [&out](auto&& ... args) { (out << ... << args) << "\n"; };
            const int scale = 50;

            print("digraph G {");

            print("  {");
            for (const auto& full_vertex : vertices)
            {
                const Vertex& v = full_vertex.forward;
                const Location& loc = *v.location;

                if (&vertices[BACKWARD_BASE] == &full_vertex)
                {
                    continue;
                }

                print
                (
                    "    ",
                    loc.index,
                    "[",
                    "shape=point,",
                    "fontsize=7,",
                    "xlabel=\"", loc.index, ": [", loc.time_window.from, ";", loc.time_window.to, "] ",
                        v.earliest_work_start_moment, " d=", loc.duration, " p=", loc.workers_required, "\",",
                    (loc.is_base() ? "color=\"#AA1010\"," : v.edges.size() == 0 ? "color=\"#AA1010\",width=0.15," : ""),
                    "pos=\"", scale * loc.point.x, ",", scale * loc.point.y, "\"",
                    "]"
                );
            }
            print("  }");
            print("");

            for (const auto& full_vertex : vertices)
            {
                const Vertex& a = full_vertex.forward;
                for (const Edge* edge : full_vertex.forward.edges)
                {
                    const Vertex& b = *edge->to;

                    print
                    (
                        "  ",
                        a.location->index, " -> ", b.location->index, " ",
                        "[",
                        "arrowsize=0.5,",
                        "penwidth=0.5,",
                        "fontsize=10,",
                        "color=\"#9ACEEB\""
                        "]"
                    );
                }
            }

            print("}");
        }

        Graph& operator=(Graph&&) = default;
        void operator=(const Graph&) = delete;
    };

    Graph Graph::copy() const
    {
//        Graph result(this->task);
//
//        // Ugly, but whatever.
//        for (size_t i = 0; i != vertices.size(); ++i)
//        {
//            const FullVertex& this_full_vertex = vertices[i];
//            FullVertex& that_full_vertex = result.vertices[i];
//
//            that_full_vertex.copy_dynamics_from(this_full_vertex);
//
//            for (size_t j = 0; j != this_full_vertex.forward.edges.size(); ++j)
//            {
//                const Edge& this_edge = this_full_vertex.forward.edges[j];
//                size_t to_index = this_edge.to->parent - vertices.data();
//
//                Edge* that_edge = result.link(&that_full_vertex.forward, &result.vertices[to_index].forward);
//                that_edge->earliest_arrive_moment = this_edge.earliest_arrive_moment;
//                that_edge->twin->earliest_arrive_moment = this_edge.twin->earliest_arrive_moment;
//            }
//        }
//
//        return result;
        throw std::runtime_error("Not implemented yet");
    }

    bool are_all_vertices_unmarked(const Graph& graph)
    {
        for (const FullVertex & full_vertex : graph.vertices)
        {
            if (full_vertex.forward.is_marked() or full_vertex.backward.is_marked())
            {
                return false;
            }
        }
        return true;
    }

    bool are_dynamics_up_to_date(Graph& graph)
    {
        assert(0 == Vertex::marked.size);

        graph.forward_start()->mark_recursively();
        graph.backward_start()->mark_recursively();

        // TODO: use graph copy instead of manual checking
        while (Vertex::marked.size != 0)
        {
            Vertex* v = Vertex::marked.pop();

            Moment vertex_dynamic_before = v->earliest_work_start_moment;
            std::vector<int> edge_dynamic_before = map(v->edges, [](Edge* edge) { return edge->earliest_arrive_moment; });

            v->recalculate();

            Moment vertex_dynamic_after = v->earliest_work_start_moment;
            std::vector<int> edge_dynamic_after = map(v->edges, [](Edge* edge) { return edge->earliest_arrive_moment; });

            if (vertex_dynamic_before != vertex_dynamic_after or edge_dynamic_before != edge_dynamic_after)
            {
                Vertex::unmark_all();
                return false;
            }
        }

        Vertex::unmark_all();
        return true;
    }

    std::vector<Moment> calculate_true_last_arrive_moments(Graph& graph)
    {
        assert(0 == Vertex::marked.size);

        std::vector<Moment> last_arrive_moments(graph.task->n + 1, -1);

        const auto update_edge = [&last_arrive_moments](const Vertex* from, const Edge* edge, Moment ready_moment)
        {
            const Moment arrive_moment = ready_moment + distance(from, edge->to);
            const Index to_index = edge->to->location->index;
            last_arrive_moments[to_index] = std::max(last_arrive_moments[to_index], arrive_moment);
        };

        assert(are_all_vertices_unmarked(graph));
        assert(are_dynamics_up_to_date(graph));

        Vertex* base = graph.forward_start();
        base->mark_recursively();

        // Base is special

        Vertex::marked.pop(); // Discard base
        for (const Edge* edge : base->edges)
        {
            update_edge(base, edge, twin_moment(edge->twin->earliest_arrive_moment));
        }

        while (Vertex::marked.size != 0)
        {
            const Vertex* from = Vertex::marked.pop();
            assert(last_arrive_moments[from->location->index] != -1);
            const Moment work_start_moment = std::max(from->location->time_window.from, last_arrive_moments[from->location->index]);
            const Moment work_end_moment = work_start_moment + from->location->duration;

            assert(work_end_moment <= from->location->time_window.to);

            for (const Edge* edge : from->edges)
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

        const Vertex* const start = graph.forward_start();
        const Vertex* const finish = graph.forward_finish();

        for (const Edge* first_edge: start->edges)
        {
            Moment current_moment = twin_moment(first_edge->twin->earliest_arrive_moment);
            output << "start " << current_moment << " 1\n";

            const Vertex* current_vertex = first_edge->to;
            current_moment += distance(start, first_edge->to);

            while (true)
            {
                const Index index = current_vertex->location->index;

                output << "arrive " << current_moment << " " << index << "\n";

                if (current_vertex == finish)
                {
                    break;
                }

                const Moment start_work_moment = std::max(current_vertex->location->time_window.from, last_arrive_moment[index]);
                const Moment end_work_moment = start_work_moment + current_vertex->location->duration;

                assert(current_moment <= start_work_moment);
                assert(current_vertex->location->time_window.from <= start_work_moment);
                assert(end_work_moment <= current_vertex->location->time_window.to);

                output << "work " << start_work_moment << " " << end_work_moment << " " << index << "\n";

                assert(visit_count[index] < current_vertex->edges.size());

                const Vertex* const next_vertex = current_vertex->edges[visit_count[index]]->to;

                current_moment = end_work_moment + distance(current_vertex, next_vertex);

                visit_count[index] += 1;
                current_vertex = next_vertex;
            }

            output << "end\n";
        }

        return output;
    }

}

inline namespace edit
{


    TimeWindow time_window_after_interpose(const Edge* ab, const Vertex* v)
    {
        Vertex* a = ab->twin->to->twin;
        Vertex* b = ab->to;

        return v->location->time_window.intersect
        (
            TimeWindow
            {
                .from = a->get_earliest_work_end_moment() + distance(a, v),
                .to = twin_moment(b->twin->get_earliest_work_end_moment() + distance(b, v))
            }
        );
    }
}

namespace scorers
{
    int interpose_reward(const Edge* ab, const Vertex* c)
    {
        const Vertex* a = ab->twin->to->twin;
        const Vertex* b = ab->to;

        const int d = c->location->duration;
        const int p = c->location->workers_required;

        return distance(a, b) - distance(a, c) - distance(c, b) + d * (p + 4);
    }

    int radial_cost(Vector center, Vertex* x)
    {
        const int d = x->location->duration;
        const int p = x->location->workers_required;

        return 2 * distance(center, x->location->point) - d * p * (p + 4);
    }
}

namespace comparators
{
    struct VertexPtrComparator
    {
        bool operator()(Vertex* a, Vertex* b) const
        {
            return a < b;
        }
    };
}

inline namespace solvers
{
    void remove_empty_paths(Graph& graph)
    {
        const Vertex* finish = graph.forward_finish();

        for (Edge* edge : reversed(graph.forward_start()->edges))
        {
            if (edge->to == finish)
            {
                graph.unlink(edge);
            }
        }
    }

    double expected_workers_required(const Graph& graph)
    {
        const double n = graph.task->n;
        const double average_workers_required = (1. + 8.) / 2.;
        const double average_duration = (5. + 30.) / 2.;
        const double average_gap = 100. / (1. + sqrt(n));
        const double work_time_window_length = 800. - 200.;

        return (n * average_workers_required * (average_duration + average_gap) - average_gap) / work_time_window_length;
    }

    void generate_empty_routes(Graph& graph)
    {
        const size_t empty_route_count = std::min
        (
            graph.task->sum_workers_required,
            static_cast<size_t>(2 * expected_workers_required(graph))
        );

        for (size_t i = 0; i != empty_route_count; ++i)
        {
            Edge* edge = graph.link(graph.forward_start(), graph.forward_finish());
            edge->earliest_arrive_moment = edge->twin->earliest_arrive_moment = 0;
        }
    }

    size_t get_total_edge_count(const Graph& graph)
    {
        size_t count = 0;

        for (const FullVertex& full_vertex : graph.vertices)
        {
            count += full_vertex.forward.edges.size();
        }

        return count;
    }

    namespace greedy
    {
        struct Candidate
        {
            Edge* edge;
            TimeWindow time_window;
            int reward;

            bool operator<(const Candidate& other) const noexcept { return reward > other.reward; }
        };

        std::vector<Candidate> get_candidates(Graph& graph, Vertex* vertex)
        {
            std::vector<Candidate> candidates;
            candidates.reserve(get_total_edge_count(graph));

            const int duration = vertex->location->duration;

            for (FullVertex& full_vertex : graph.vertices)
            {
                for (Edge* edge : full_vertex.forward.edges)
                {
                    const TimeWindow time_window = time_window_after_interpose(edge, vertex);
                    if (time_window.length() >= duration)
                    {
                        candidates.emplace_back(Candidate
                        {
                            .edge = edge,
                            .time_window = time_window,
                            .reward = scorers::interpose_reward(edge, vertex)
                        });
                    }
                }
            }

            std::sort(begin(candidates), end(candidates));

            return candidates;
        }

        std::vector<const Candidate*> select_from_candidates(const std::vector<Candidate>& candidates, Vertex* vertex)
        {
            assert(0 == Vertex::marked.size);

            const size_t workers_required = static_cast<size_t>(vertex->location->workers_required);
            const int duration = vertex->location->duration;
            const int max_attempt_count = 8; // TODO: parameterize

            std::vector<const Candidate*> chosen;
            chosen.reserve(workers_required);

            for (int i = 0; i != max_attempt_count; ++i)
            {
                TimeWindow current_time_window { .from = 0, .to = MAX_MOMENT };

                for (const Candidate& candidate : candidates)
                {
                    const TimeWindow new_time_window = current_time_window.intersect(candidate.time_window);

                    if (new_time_window.length() < duration)
                    {
                        continue;
                    }

                    if (candidate.edge->twin->to->twin->is_marked() or candidate.edge->to->twin->is_marked())
                    {
                        continue;
                    }

                    candidate.edge->to->mark_recursively();
                    candidate.edge->twin->to->mark_recursively();

                    current_time_window = new_time_window;
                    chosen.emplace_back(&candidate);

                    if (chosen.size() == workers_required)
                    {
                        // Keep vertices marked because they're exactly those which need to be recomputed.
                        return chosen;
                    }
                }

                Vertex::unmark_all();
                chosen.clear();
            }

            return chosen;
        }

        int insert(Graph& graph, Vertex* vertex)
        {
            assert(not vertex->location->is_base());
            assert(vertex->edges.empty());
            assert(vertex->twin->edges.empty());

            const std::vector<Candidate> candidates = get_candidates(graph, vertex);
            const std::vector<const Candidate*> chosen = select_from_candidates(candidates, vertex);

            int reward = 0;

            if (not chosen.empty())
            {
                assert(chosen.size() == static_cast<size_t>(vertex->location->workers_required));

                for (const Candidate* candidate : chosen)
                {
                    reward += candidate->reward;
                    graph.interpose(candidate->edge, vertex);
                }

                assert(vertex->edges.size() == static_cast<size_t>(vertex->location->workers_required));
                assert(vertex->twin->edges.size() == static_cast<size_t>(vertex->location->workers_required));

                vertex->mark_recursively();
                vertex->twin->mark_recursively();

                Vertex::recalculate_all_marked();
            }

            assert(0 == Vertex::marked.size);

            return reward;
        }

        std::vector<Vertex*> get_orders(Graph& graph)
        {
            std::vector<Vertex*> orders;
            orders.reserve(graph.task->n - 1);

            assert(graph.vertices.size() == orders.capacity() + 2);

            for (size_t i = 2; i != graph.vertices.size(); ++i)
            {
                orders.emplace_back(&graph.vertices[i].forward);
            }

            return orders;
        }

        // TODO: параметризовать
        int optimize(Graph& graph)
        {
            std::vector<Vertex*> orders = get_orders(graph);
            // TODO: accept comparator as a parameter
            comparators::VertexPtrComparator comparator;
            std::sort(begin(orders), end(orders), comparator);

            int total_reward = 0;

            for (Vertex* order : orders)
            {
                if (order->edges.empty())
                {
                    total_reward += insert(graph, order);
                }
            }

            return total_reward;
        }
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
        greedy::optimize(graph);
        remove_empty_paths(graph);
    }
};

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    #ifndef ONLINE_JUDGE
        const auto ensure = [](bool cond, std::string_view message="")
        {
            if (not cond)
            {
                throw std::runtime_error(std::string(message));
            }
        };

        std::unique_ptr<std::ifstream> input_file = nullptr;

        for (int i = 1; i < argc; ++i)
        {
            std::string_view arg_i(argv[i]);

            if (arg_i == "--task")
            {
                ensure(argc > i + 1);
                i += 1;
                input_file = std::make_unique<std::ifstream>(argv[i]);
            }
            else
            {
                ensure(false, "unrecognized argument");
            }
        }

        std::istream& input = (input_file == nullptr) ? std::cin : *input_file;
    #else
        std::istream& input = std::cin;
    #endif

    Processor processor(input);
    processor.run_circuit();

    #ifndef ONLINE_JUDGE
        processor.graph.dump_graphviz();
    #endif

    std::cout << processor.graph;

    #ifndef ONLINE_JUDGE
        std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds." << std::endl;
    #endif
}
