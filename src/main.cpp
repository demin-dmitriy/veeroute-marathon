#ifdef ONLINE_JUDGE
    #pragma GCC optimize ("O3")
    #define NDEBUG
#endif


#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include <ctime>
#include <cassert>


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
    int p;
    TimeWindow time_window;

    static Location read(std::istream& in)
    {
        Location result;
        in >> result.point.x
           >> result.point.y
           >> result.duration
           >> result.p
           >> result.time_window.from
           >> result.time_window.to;
        return result;
    }
};

struct Task
{
    int n;
    std::vector<Location> locations;

    void read(std::istream& in)
    {
        assert(locations.empty());

        in >> n;
        locations.reserve(n);
        for (int i = 0; i != n; ++i)
        {
            locations.emplace_back(Location::read(in));
        }
    }
};

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



Task task;
const Timer timer;
constexpr double TIME_LIMIT_SECONDS = 15;

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    #ifdef READ_TASK_FROM_FILE
        assert(argc == 2);
        {
            std::ifstream input(argv[1]);
            task.read(input);
        }
    #else
        task.read(std::cin);
    #endif

    std::cerr << "Elapsed: " << std::setprecision(2) << std::fixed << timer.seconds_elapsed() << " seconds.\n";
}
