#ifdef ONLINE_JUDGE
    #pragma GCC optimize ("O3")
    #ifndef NDEBUG
        #define NDEBUG
    #endif
#endif

#include <algorithm>
#include <array>
#include <deque>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <string_view>
#include <vector>

#include <cassert>
#include <ctime>
#include <cmath>

#ifdef ONLINE_JUDGE
    #define TRACE(...)
    #define TRACELN(...)
#else
    std::unique_ptr<std::ofstream> trace_file = nullptr;
    #define TRACE(...) trace(__VA_ARGS__)
    #define TRACELN(...) traceln(__VA_ARGS__)
#endif

constexpr double TIME_LIMIT_SECONDS = 15;

#ifndef ONLINE_JUDGE
    int WORKERS_COUNT = -1;
#endif

std::default_random_engine RANDOM_ENGINE(774130);

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

    enum Bool : bool { TRUE = true, FALSE = false };

    template <class T, size_t MAX_SIZE>
    struct CheapStack
    {
        size_t size = 0;
        std::array<T, MAX_SIZE> elements;

        void push(T value)
        {
            elements[size] = std::move(value);
            size += 1;
            assert(size <= MAX_SIZE);
        }

        T pop()
        {
            assert(size > 0);
            size -= 1;
            return std::move(elements[size]);
        }

        T& operator[](size_t i)
        {
            assert(i < size);
            return elements[i];
        }

        const T& operator[](size_t i) const
        {
            assert(i < size);
            return elements[i];
        }

        void swap_remove(size_t index)
        {
            using std::swap;
            assert(index < size);
            swap(elements[index], elements[size - 1]);
            size -= 1;
        }

        auto begin() { return elements.begin(); }
        auto begin() const { return elements.begin(); }
        auto end() { return elements.begin() + size; }
        auto end() const { return elements.begin() + size; }

        bool empty() const
        {
            return 0 == size;
        }

        void clear()
        {
            size = 0;
        }
    };

    template <class Container, class Func>
    auto map(Container&& values, Func&& func) -> std::vector<decltype(func(*values.begin()))>
    {
        using U = decltype(func(*values.begin()));
        std::vector<U> result;
        result.reserve(values.end() - values.begin());

        for (auto& value : values)
        {
            result.emplace_back(func(value));
        }

        return result;
    }

    template <class T>
    std::vector<T> iota(size_t n, T start = 0)
    {
        std::vector<T> result;
        result.resize(n);
        std::iota(begin(result), end(result), start);
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
        size_t next_free;
        std::vector<T> storage;
        std::vector<T*> free;

        explicit Pool(size_t size)
            : next_free(0)
        {
            storage.resize(size);
            free.resize(size);

            for (size_t i = 0; i != size; ++i)
            {
                free[i] = storage.data() + i;
            }
        }

        size_t allocated_count() const
        {
            return next_free;
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
        T* m_data;
        size_t m_size;


        template <class Container>
        Span(Container&& container)
            : m_data(container.data())
            , m_size(container.size())
        { }

        T* data() const
        {
            return m_data;
        }

        size_t size() const
        {
            return m_size;
        }

        T& operator[](size_t i) { assert(i < m_size); return m_data[i]; }
        const T& operator[](size_t i) const { assert(i < m_size); return m_data[i]; }

        T* begin() { return m_data; }
        T* end() { return m_data + m_size; }
        const T* begin() const { return m_data; }
        const T* end() const { return m_data + m_size; }
    };

    // small_vector::max_size == MAX_WOREKRS_REQUIRED
    // Try replacing this later with more efficient class. (maybe just extend CheapStack?)
    template <class T>
    using small_vector = std::vector<T>;

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

    template <class T, class Func>
    void sort_rewards_with(std::vector<T>& values, Func&& func)
    {
        const auto rewards = map(values, std::forward<Func>(func));
        std::vector<size_t> indices = iota<size_t>(values.size());
        std::sort(begin(indices), end(indices), [&rewards](size_t a, size_t b) { return rewards[a] > rewards[b]; });
        values = map(indices, [&values](size_t i) { return values[i]; });
    }

    template <class T, class Container>
    Span<T> skip(size_t k, Container&& container)
    {
        Span<T> span(container);
        if (k >= container.size())
        {
            span.m_size = 0;
        }
        else
        {
            span.m_data += k;
            span.m_size -= k;
        }
        return span;
    }

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
            auto begin() const { return ((const Range&) range).rbegin(); }
            auto end() { return range.rend(); }
            auto end() const { return ((const Range&) range).rend(); }
        };

        return ReversedRange(range);
    }

    template <class T, class Score>
    struct TopKQueue
    {
        struct Node;

        struct ScoreInfo
        {
            Score score;
            Node* node;

            bool operator<(const ScoreInfo& other) const noexcept { return score > other.score; }
        };

        struct Node
        {
            Node* prev;
            Node* next;
            typename std::multiset<ScoreInfo>::iterator iterator;
            T data;
        };

        const size_t k;
        Pool<Node> pool;
        std::multiset<ScoreInfo> score_set;
        Node begin;
        Node end;

        explicit TopKQueue(size_t k)
            : k(k)
            , pool(k + 1)
            , begin
            {
                .prev = nullptr,
                .next = &end,
                .iterator = {},
                .data = {}
            }
            , end
            {
                .prev = &begin,
                .next = nullptr,
                .iterator = {},
                .data = {}
            }
        {
            assert(k > 0);
        }

        bool empty() const
        {
            return 0 == pool.allocated_count();
        }

        T pop()
        {
            assert(pool.allocated_count() > 0);
            assert(pool.allocated_count() == score_set.size());
            assert(begin.next != &end);
            assert(begin.next != nullptr);

            Node* first = begin.next;
            T result = std::move(first->data);
            score_set.erase(first->iterator);
            remove(first);
            assert(pool.allocated_count() == score_set.size());
            return result;
        }

        void push(T value, Score score)
        {
            assert(pool.allocated_count() <= k);
            assert(pool.allocated_count() == score_set.size());

            Node* node = insert_before(std::move(value), &end);
            node->iterator = score_set.insert(ScoreInfo{.score = score, .node = node});

            if (score_set.size() == k + 1)
            {
                auto worst_iterator = score_set.begin();
                remove(worst_iterator->node);
                score_set.erase(worst_iterator);
            }

            assert(pool.allocated_count() <= k);
            assert(pool.allocated_count() == score_set.size());
        }

    private:
        Node* insert_before(T value, Node* after)
        {
            Node* node = pool.allocate();
            node->data = std::move(value);
            node->prev = after->prev;
            node->next = after;
            node->prev->next = node;
            node->next->prev = node;
            return node;
        }

        void remove(Node* node)
        {
            node->next->prev = node->prev;
            node->prev->next = node->next;
            pool.deallocate(node);
        }
    };

    #ifndef ONLINE_JUDGE
        void trace() { }

        template <class Arg, class ... Args>
        void trace(const Arg& arg, const Args& ... args)
        {
            if (trace_file)
            {
                (*trace_file) << arg << " ";
                trace(args...);
            }
        }

        template <class ... Args>
        void traceln(const Args& ... args)
        {
            if (trace_file)
            {
                trace(args...);
                (*trace_file) << "\n";
            }
        }
    #endif
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
        size_t workers_required;
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

        void update(Vertex* source);
    };

    struct FullVertex;
    int distance(const Vertex* a, const Vertex* b);

    struct Vertex
    {
        using Version = unsigned int; // Won't overflow in 15 seconds. Could overflow, if it was running longer.

        static inline CheapStack<Vertex*, 2 * (MAX_LOCATION_COUNT + 1)> marked; // Stores both forward and backward vertices
        static inline Version current_global_version = 1;

        // TODO: maybe separate class for storing dynamics
        Moment earliest_work_end_moment;
        Version version; // TODO: лучше как-то по-другому оформить. Может быть завести уникальный индекс для Vertex?
        size_t index;
        const Location* const location;
        Vertex* const twin;
        small_vector<Edge*> edges;

        explicit Vertex(size_t index, Vertex* twin, const Location* location)
            : earliest_work_end_moment(-1)
            , version(0)
            , index(index)
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

        Moment get_earliest_work_start_moment() const
        {
            return earliest_work_end_moment - location->duration;
        }

        Moment get_earliest_work_end_moment() const
        {
            return earliest_work_end_moment;
        }

        void recalculate()
        {
            Moment earliest_work_start_moment = location->time_window.from;

            for (const Edge* twin_edge : twin->edges)
            {
                earliest_work_start_moment = std::max(earliest_work_start_moment, twin_edge->twin->earliest_arrive_moment);
            }

            earliest_work_end_moment = earliest_work_start_moment + location->duration;

            assert(earliest_work_end_moment <= location->time_window.to);

            for (Edge* edge : edges)
            {
                edge->update(this);
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

    void Edge::update(Vertex* source)
    {
        earliest_arrive_moment = source->get_earliest_work_end_moment() + distance(source, to);
    }

    int distance(const Vertex* a, const Vertex* b)
    {
        return distance(a->location, b->location);
    }

    struct FullVertex
    {
        Vertex forward;
        Vertex backward;

        explicit FullVertex(size_t index, const FullLocation* location)
            : forward(2 * index, &backward, &location->forward)
            , backward(2 * index + 1, &forward, &location->backward)
        { }

        void mark_recursively()
        {
            forward.mark_recursively();
            backward.mark_recursively();
        }

        void recalculate()
        {
            forward.recalculate();
            backward.recalculate();
        }

//        void copy_dynamics_from(const FullVertex& other)
//        {
//            forward.earliest_work_start_moment = other.forward.earliest_work_start_moment;
//            backward.earliest_work_start_moment = other.backward.earliest_work_start_moment;
//        }
    };

    static_assert(sizeof(FullVertex) == 2 * sizeof(Vertex));

    std::vector<std::pair<Vertex*, Vertex*>> cut_assignment(Vertex* v)
    {
        const std::vector<Vertex*> froms = map(v->twin->edges, [](const Edge* edge) { return edge->to->twin; });
        const std::vector<Vertex*> tos = map(v->edges, [](const Edge* edge) { return edge->to; });

        assert(froms.size() == tos.size());
        const size_t m = froms.size();

        Table<int> distances(m, m);

        for (size_t i = 0; i != m; ++i)
        {
            for (size_t j = 0; j != m; ++j)
            {
                distances.at(i, j) = distance(froms[i], tos[j]);
            }
        }

        const std::vector<size_t> assignment = solve_assignment_problem(distances);

        std::vector<std::pair<Vertex*, Vertex*>> result;

        for (size_t i = 1; i <= m; ++i)
        {
            result.emplace_back(froms[assignment[i] - 1], tos[i - 1]);
        }

        return result;
    }

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
            vertices.emplace_back(vertices.size(), &task->locations[0]);

            for (const FullLocation& location : task->locations)
            {
                vertices.emplace_back(vertices.size(), &location);
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

        Vertex* get(size_t vertex_index)
        {
            assert(vertex_index < 2 * vertices.size());
            return &vertices[0].forward + vertex_index;
        }

        const Vertex* get(size_t vertex_index) const
        {
            assert(vertex_index < 2 * vertices.size());
            return &vertices[0].forward + vertex_index;
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

            ac->update(a);
            cb->twin->update(b->twin);

            assert(c->edges.size() == c->twin->edges.size());
        }

        void cut(Vertex* v)
        {
            const auto assignment = cut_assignment(v);

            for (Vertex* u : { v, v->twin })
            {
                for (Edge* edge : u->edges)
                {
                    unlink(edge);
                }
            }

            for (auto [from, to] : assignment)
            {
                Edge* edge = link(from, to);
                edge->update(from);
                edge->twin->update(to->twin);
            }
        }

        void two_opt_swap(Edge* a, Edge* b)
        {
            std::swap
            (
                a->to->twin->edges[a->twin->index],
                b->to->twin->edges[b->twin->index]
            );

            std::swap(a->twin->index, b->twin->index);
            std::swap(a->to, b->to);


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
                        v.get_earliest_work_start_moment(), " d=", loc.duration, " p=", loc.workers_required, "\",",
                    (loc.is_base() ? "color=\"#AA1010\"," : v.edges.empty() ? "color=\"#AA1010\",width=0.15," : ""),
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

            Moment vertex_dynamic_before = v->earliest_work_end_moment;
            std::vector<int> edge_dynamic_before = map(v->edges, [](Edge* edge) { return edge->earliest_arrive_moment; });

            v->recalculate();

            Moment vertex_dynamic_after = v->earliest_work_end_moment;
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

    struct ReachabilityTable
    {
        using Row = std::vector<Vertex::Version>;
        std::vector<Vertex::Version> last_version;
        std::vector<Row> table;

        explicit ReachabilityTable(Graph& graph)
            : last_version(2 * graph.vertices.size(), 0)
        {
            table.resize(2 * graph.vertices.size());

            for (Row& row : table)
            {
                // Note: row[i] stores version for FullVertex, but last_version[j] stores version for Vertex.
                row.resize(graph.vertices.size());
            }
        }

        bool is_reachable(Vertex* from, Vertex* to) const
        {
            assert(is_current(from->index));
            return table[from->index][to->index / 2] == Vertex::current_global_version;
        }

        bool mark_recursively(const Edge* edge)
        {
            return mark_recursively(edge->to) | mark_recursively(edge->twin->to);
        }

        bool mark_recursively(const Vertex* vertex)
        {
            const size_t index = vertex->index;

            assert(index < table.size());

            if (not is_current(index))
            {
                last_version[index] = Vertex::current_global_version;
                mark_recursively(table[index], vertex);
                return true;
            }

            return false;
        }

    private:
        bool is_current(size_t index) const
        {
            return last_version[index] == Vertex::current_global_version;
        }

        static void mark_recursively(Row& row, const Vertex* vertex)
        {
            const size_t index = vertex->index / 2;
            if (Vertex::current_global_version != row[index])
            {
                row[index] = Vertex::current_global_version;

                for (const Edge* edge : vertex->edges)
                {
                    mark_recursively(row, edge->to);
                }
            }
        }
    };
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

    int radial_cost(Vector center, const Vertex* x)
    {
        const int d = x->location->duration;
        const int p = x->location->workers_required;

        return 2 * distance(center, x->location->point) - d * p * (p + 4);
    }

    int removal_reward(Vertex* vertex)
    {
        // TODO: distance не играет по факту роли, если нужно и так и так ждать.
        assert(vertex->edges.size() == vertex->location->workers_required);

        const int d = vertex->location->duration;
        const int p = vertex->location->workers_required;

        int reward_now = d * p * (p + 4);

        // TODO: попробовать сделать версию этого, которая основывается на основе значениях динамики
        for (const Vertex* v : { vertex, vertex->twin })
        {
            for (const Edge* edge : v->edges)
            {
                reward_now -= distance(vertex, edge->to);
            }
        }

        int reward_after = 0;

        for (auto [from, to] : cut_assignment(vertex))
        {
            reward_after -= distance(from, to);
        }

        return reward_after - reward_now;
    }
}

inline namespace solvers
{
    std::vector<Vertex*> get_all_orders(Graph& graph)
    {
        return map(skip<FullVertex>(2, graph.vertices), [](FullVertex& full_vertex) { return &full_vertex.forward; });
    }

    std::vector<Vertex*> get_unresolved_orders(Graph& graph)
    {
        std::vector<Vertex*> orders;
        orders.reserve(graph.task->n);

        for (FullVertex& full_vertex : skip<FullVertex>(2, graph.vertices))
        {
            if (full_vertex.forward.edges.empty())
            {
                orders.emplace_back(&full_vertex.forward);
            }
        }

        return orders;
    }

    template <class Scorer>
    std::vector<Vertex*> get_sorted_unresolved_orders(Graph& graph, Scorer&& scorer)
    {
        std::vector<Vertex*> orders = get_unresolved_orders(graph);
        sort_rewards_with(orders, std::forward<Scorer>(scorer));
        return orders;
    }

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
        const double a = 0.2358594030965598;
        const double b = 54.57148230518328;
        const double sigma = 13.717673452304112;

        return a * n + b + sigma;
    }

    void generate_empty_routes(Graph& graph)
    {
        size_t empty_route_count = std::min(graph.task->sum_workers_required, static_cast<size_t>(expected_workers_required(graph)));

        #ifndef ONLINE_JUDGE
            if (WORKERS_COUNT != -1)
            {
                assert(WORKERS_COUNT > 0);
                empty_route_count = std::min(graph.task->sum_workers_required, (size_t) WORKERS_COUNT);
            }
        #endif

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
        struct EdgeProfile
        {
            Edge* edge;
            TimeWindow time_window;
            int reward;

            bool operator<(const EdgeProfile& other) const noexcept { return reward > other.reward; }
        };

        std::vector<EdgeProfile> get_edge_profiles(Graph &graph, Vertex *vertex)
        {
            std::vector<EdgeProfile> profiles;
            profiles.reserve(get_total_edge_count(graph));

            const int duration = vertex->location->duration;

            for (FullVertex& full_vertex : graph.vertices)
            {
                for (Edge* edge : full_vertex.forward.edges)
                {
                    const TimeWindow time_window = time_window_after_interpose(edge, vertex);
                    if (time_window.length() >= duration)
                    {
                        profiles.emplace_back(EdgeProfile
                        {
                            .edge = edge,
                            .time_window = time_window,
                            .reward = scorers::interpose_reward(edge, vertex)
                        });
                    }
                }
            }

            std::sort(begin(profiles), end(profiles));

            return profiles;
        }

        using EdgesSelection = CheapStack<size_t, MAX_WORKERS_REQUIRED>;

        struct EdgesSelector
        {
            static inline std::unique_ptr<ReachabilityTable> reachability_table = nullptr;

            int min_reward_threshold;
            double reward_percentile_cutoff;
            size_t width;
            int max_eval_points;
            size_t max_selections_to_consider;

            static void setup_reachability_table(Graph& graph)
            {
                if (nullptr == reachability_table)
                {
                    reachability_table = std::make_unique<ReachabilityTable>(graph);
                }

                Vertex::current_global_version += 1;
            }

            EdgesSelection operator()(Graph& graph, std::vector<EdgeProfile>& profiles, Vertex* vertex)
            {
                assert(0 == Vertex::marked.size);
                assert(max_selections_to_consider > 0);

                setup_reachability_table(graph);

                std::vector<Searcher::State> selections;
                selections.reserve(max_selections_to_consider);

                Searcher::Listener listener
                (
                    [&](const Searcher::State& state) -> Action
                    {
                        selections.push_back(state);

                        TRACE("candidate:");
                        for (size_t i : state.selection)
                        {
                            TRACE(i);
                        }
                        TRACELN("(reward", state.reward, ")");

                        if (selections.size() == max_selections_to_consider)
                        {
                            return STOP;
                        }

                        return ADD_STATE;
                    }
                );

                Searcher searcher(profiles, vertex, max_eval_points, width, listener);
                searcher.search();

                if (not selections.empty())
                {
                    int highest_reward = min_reward_threshold;

                    for (const Searcher::State& state : selections)
                    {
                        highest_reward = std::max(highest_reward, state.reward);
                    }

                    const int threshold = std::max(min_reward_threshold, int(highest_reward * reward_percentile_cutoff));

                    std::vector<const Searcher::State*> filtered_selections;
                    filtered_selections.reserve(max_selections_to_consider);

                    TRACE("filtered:");
                    for (const Searcher::State& state : selections)
                    {
                        if (state.reward >= threshold)
                        {
                            TRACE(&state - &selections[0]);
                            filtered_selections.push_back(&state);
                        }
                    }
                    TRACELN();

                    if (filtered_selections.size() > 0)
                    {
                        std::uniform_int_distribution<size_t> random(0, filtered_selections.size() - 1);
                        return filtered_selections[random(RANDOM_ENGINE)]->selection;
                    }
                }

                return {};
            }

            enum Action
            {
                ADD_STATE,
                SKIP_STATE,
                STOP
            };

            struct Searcher
            {
                struct State
                {
                    EdgesSelection selection;
                    size_t index;
                    int reward;
                };

                struct ListenerI
                {
                    virtual ~ListenerI() = default;
                    virtual Action on_found_selection(const State&) = 0;
                };

                template <class Lambda>
                struct Listener final : Searcher::ListenerI
                {
                    Lambda lambda;
                    explicit Listener(Lambda lambda) : lambda(std::move(lambda)) { }
                    Action on_found_selection(const State& state) override { return lambda(state); }
                };

                static inline const int QUEUE_POP_POINTS = 1;
                static inline const int MARK_POINTS = 40;

                const int duration;
                const size_t workers_required;
                int eval_points_left;
                const std::vector<EdgeProfile>& profiles;
                TopKQueue<State, int> queue;
                std::vector<Bool> chosen;
                ListenerI& listener;

                Searcher
                (
                    const std::vector<EdgeProfile>& profiles,
                    Vertex* vertex,
                    int eval_points_left,
                    size_t width,
                    ListenerI& listener
                )
                    : duration(vertex->location->duration)
                    , workers_required(vertex->location->workers_required)
                    , eval_points_left(eval_points_left)
                    , profiles(profiles)
                    , queue(width)
                    , chosen(profiles.size(), FALSE)
                    , listener(listener)
                { }

                int get_reward(const State& state)
                {
                    int reward = 0;

                    for (size_t i : state.selection)
                    {
                        reward += profiles[i].reward;
                    }

                    return reward;
                }

                void search()
                {
                    State start;
                    start.index = 0;

                    switch (greedy_add_to_selection(start))
                    {
                        case ADD_STATE:
                            queue.push(start, start.reward);
                            break;
                        case SKIP_STATE:
                            break;
                        case STOP:
                            return;
                    }

                    while (not queue.empty())
                    {
                        State state = queue.pop();
                        eval_points_left -= QUEUE_POP_POINTS;

                        for (size_t i = 0; i != state.selection.size; ++i)
                        {
                            size_t pos = state.selection[i];

                            if (eval_points_left <= 0)
                            {
                                return;
                            }

                            if (pos >= state.index)
                            {
                                State new_state
                                {
                                    .selection = state.selection,
                                    .index = pos + 1,
                                    .reward = {}
                                };

                                new_state.selection.swap_remove(i);

                                switch (greedy_add_to_selection(new_state))
                                {
                                    case ADD_STATE:
                                        queue.push(new_state, new_state.reward);
                                        break;
                                    case SKIP_STATE:
                                        break;
                                    case STOP:
                                        return;
                                }
                            }
                        }
                    }
                }

                bool is_marked(const EdgesSelection& selection, Edge* edge)
                {
                    Vertex* from = edge->twin->to->twin;
                    Vertex* to = edge->to->twin;

                    for (size_t i : selection)
                    {
                        const Edge* selected_edge = profiles[i].edge;

                        if (reachability_table->is_reachable(selected_edge->to, from)
                            or reachability_table->is_reachable(selected_edge->twin->to, to))
                        {
                            return true;
                        }
                    }

                    return false;
                }

                // TODO: needs some tracing to understand if this effective or not.
                Action greedy_add_to_selection(State& state)
                {
                    EdgesSelection& selection = state.selection;

                    assert(selection.size < workers_required);

                    if (state.index + workers_required - selection.size > profiles.size())
                    {
                        return SKIP_STATE;
                    }

                    TimeWindow current_time_window { .from = 0, .to = MAX_MOMENT };

                    for (size_t i : selection)
                    {
                        assert(i < profiles.size());
                        chosen[i] = TRUE;
                        current_time_window = current_time_window.intersect(profiles[i].time_window);
                    }

                    assert(current_time_window.length() >= duration);

                    for (size_t i = state.index; i < profiles.size(); ++i)
                    {
                        if (FALSE == chosen[i])
                        {
                            const EdgeProfile& profile = profiles[i];
                            TimeWindow new_time_window = current_time_window.intersect(profile.time_window);

                            if (new_time_window.length() < duration)
                            {
                                continue;
                            }

                            if (is_marked(selection, profile.edge))
                            {
                                continue;
                            }

                            if (reachability_table->mark_recursively(profile.edge))
                            {
                                eval_points_left -= MARK_POINTS;
                            }

                            current_time_window = new_time_window;
                            selection.push(i);

                            if (selection.size == workers_required)
                            {
                                break;

                            }
                        }
                    }

                    for (size_t i : selection)
                    {
                        chosen[i] = FALSE;
                    }

                    state.reward = get_reward(state);

                    if (selection.size == workers_required)
                    {
                        return listener.on_found_selection(state);
                    }

                    return ADD_STATE;
                }
            };
        };

        template <class Strategies>
        struct Inserter
        {
            size_t max_profiles;
            Strategies strategies;

            void operator()(Graph& graph, Vertex* vertex)
            {
                assert(not vertex->location->is_base());
                assert(vertex->edges.empty());
                assert(vertex->twin->edges.empty());

                // NOTE: this could be a customization point, but probably it isn't needed.
                std::vector<EdgeProfile> profiles = get_edge_profiles(graph, vertex);

                if (profiles.size() > max_profiles)
                {
                    profiles.resize(max_profiles);
                }

                const EdgesSelection selection = strategies.select_edges(graph, profiles, vertex);

                if (not selection.empty())
                {
                    assert(selection.size == vertex->location->workers_required);

                    for (size_t i : selection)
                    {
                        assert(i < profiles.size());
                        graph.interpose(profiles[i].edge, vertex);
                    }

                    assert(vertex->edges.size() == vertex->location->workers_required);
                    assert(vertex->twin->edges.size() == vertex->location->workers_required);

                    vertex->mark_recursively();
                    vertex->twin->mark_recursively();

                    Vertex::recalculate_all_marked();
                }

                assert(0 == Vertex::marked.size);
            }
        };

        template <class Strategies>
        struct Greedy
        {
            Strategies strategies;

            void operator()(Graph& graph)
            {
                for (Vertex* order : strategies.choose_orders(graph))
                {
                    assert(order->edges.empty());
                    strategies.insert(graph, order);
                }
            }
        };

        struct InserterStrategies
        {
            EdgesSelector select_edges;
        };

        struct GreedyStrategies
        {
            static inline const auto choose_orders = [](Graph& graph) -> std::vector<Vertex*>
            {
                const Vector center = graph.forward_start()->location->point;
                return get_sorted_unresolved_orders(graph, [center](const Vertex* x) { return distance(center, x->location->point); });
            };

            Inserter<InserterStrategies> insert;
        };
    }

    namespace ruin
    {
        int remove_unprofitable_orders(Graph& graph)
        {
            assert(0 == Vertex::marked.size);

            int total_reward = 0;
            bool changed = true;

            while (changed)
            {
                changed = false;

                for (FullVertex& full_vertex : skip<FullVertex>(2, graph.vertices))
                {
                    Vertex* vertex = &full_vertex.forward;

                    if (not vertex->edges.empty())
                    {
                        assert(vertex->edges.size() == vertex->location->workers_required);
                        const int reward = scorers::removal_reward(vertex); // NOTE: this reward is customization point

                        // NOTE: this condition is customization point (e.g. choose temperature).
                        if (reward >= 0)
                        {
                            full_vertex.mark_recursively();
                            total_reward += reward;
                            graph.cut(vertex); // NOTE: recalculating assignment twice here, but it doesn't really matter
                            changed = true;
                            assert(vertex->edges.empty());
                        }
                    }
                }
            }

            Vertex::recalculate_all_marked();

            return total_reward;
        }
    }

    int edge_edge_optimizer(Graph&)
    {

        return 0;
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
        greedy::GreedyStrategies strategies {
            .insert = {
                .max_profiles = 1000,
                .strategies = {
                    .select_edges = {
                        .min_reward_threshold = -10000,
                        .reward_percentile_cutoff = 1,
                        .width = 16,
                        .max_eval_points = 140000,
                        .max_selections_to_consider = 5
                    }
                }
            }
        };

        generate_empty_routes(graph);
        greedy::Greedy<greedy::GreedyStrategies>{strategies}(graph);
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
        bool report = false;
        bool graphviz = false;

        for (int i = 1; i < argc; ++i)
        {
            std::string_view arg_i(argv[i]);

            if (arg_i == "--task")
            {
                ensure(argc > i + 1);
                i += 1;
                input_file = std::make_unique<std::ifstream>(argv[i]);
            }
            else if (arg_i == "--report")
            {
                report = true;
            }
            else if (arg_i == "--graphviz")
            {
                graphviz = true;
            }
            else if (arg_i == "--workers-count")
            {
                ensure(argc > i + 1);
                i += 1;
                WORKERS_COUNT = std::stoi(argv[i]);
            }
            else if (arg_i == "--trace")
            {
                trace_file = std::make_unique<std::ofstream>("trace.log", std::ios::app);
                (*trace_file) << "\n--- --- ---\n";
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
        if (graphviz)
        {
            processor.graph.dump_graphviz();
        }
    #endif

    std::cout << processor.graph;

    #ifndef ONLINE_JUDGE
        if (report)
        {
            std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds." << std::endl;
        }
    #endif
}
