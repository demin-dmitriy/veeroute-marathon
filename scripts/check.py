#!/usr/bin/env python3

import sys
import re

from argparse import ArgumentParser


LINE_REGEX = re.compile(r'''
    ^
    (?:
        (?P<start> start \s+ (?P<start_at_moment> \d+) \s+ (?P<start_at_location> \d+) )
        |
        (?P<arrive> arrive \s+ (?P<arrive_at_moment> \d+) \s+ (?P<arrive_at_location> \d+) )
        |
        (?P<work> work \s+ (?P<work_start_moment> \d+) \s+ (?P<work_end_moment> \d+) \s+ (?P<work_at_location> \d+) )
        |
        (?P<end> end )
    )
    \s*$
''', re.VERBOSE)


class ScheduleError(Exception):
    pass


class Location:
    def __init__(self, line):
        x, y, d, p, l, h = [ int(value) for value in line.split() ]
        self.point = (x, y)
        self.duration = d
        self.p = p
        self.time_window = (l, h)

        self.work_start_moment = None
        self.workers_count = 0

    def __str__(self):
        return f'{{point={self.point}, duration={self.duration}, p={self.p}, time_window={self.time_window}}}'

    def __repr__(self):
        return self.__str__()


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument('--task', required=True)
    return parser.parse_args(argv)


def segment(start, end):
    return range(start, end + 1)


def read_locations(f):
    lines = f.readlines()
    n = int(lines[0])
    locations = [ None ] + [ Location(line) for line in lines[1:] ]
    assert(len(locations) == n + 1)
    return locations


def dist(location1, location2):
    return abs(location1.point[0] - location2.point[0]) + abs(location1.point[1] - location2.point[1])


def ensure(condition, message=''):
    if not condition:
        raise ScheduleError(message)


def main(argv):
    args = parse_args(argv)

    with open(args.task) as f:
        locations = read_locations(f)

    current_location = None
    current_time = None
    total_reward = 0
    total_worker_count = 0
    total_distance = 0

    for line in sys.stdin:
        m = LINE_REGEX.fullmatch(line)
        assert m
        kind = m.lastgroup

        if kind == 'start':
            at_moment = int(m.group('start_at_moment'))
            at_location = int(m.group('start_at_location'))

            ensure(current_location is None)
            assert current_time is None
            ensure(at_location == 1)
            ensure(0 <= at_moment <= 1000)

            current_location = at_location
            current_time = at_moment
            total_reward -= 240
            total_reward += at_moment # At the 'end' we will subtract current_time
            total_worker_count += 1

        elif kind == 'arrive':
            at_moment = int(m.group('arrive_at_moment'))
            at_location = int(m.group('arrive_at_location'))

            distance = dist(locations[at_location], locations[current_location])

            ensure(at_moment >= current_time + distance)

            current_location = at_location
            current_time = at_moment
            total_distance += distance

        elif kind == 'work':
            start_moment = int(m.group('work_start_moment'))
            end_moment = int(m.group('work_end_moment'))
            at_location = int(m.group('work_at_location'))

            location = locations[at_location]

            ensure(current_location == at_location)
            ensure(start_moment >= current_time)
            ensure(end_moment - start_moment == location.duration)
            ensure(location.time_window[0] <= start_moment <= end_moment <= location.time_window[1])

            if location.work_start_moment is not None:
                ensure(start_moment == location.work_start_moment)
            else:
                location.work_start_moment = start_moment

            ensure(location.workers_count < location.p)

            location.workers_count += 1

            current_time = end_moment

        elif kind == 'end':
            ensure(current_location == 1)
            ensure(0 <= current_time <= 1000)

            total_reward -= current_time # See 'start'

            current_location = None
            current_time = None

        else:
            assert False, 'Invalid command in schedule'

    ensure(current_location is None)
    assert current_time is None

    location_stats = { i: 0 for i in segment(1, 7) }

    for location in locations[1:]:
        if location.work_start_moment is not None:
            ensure(location.p == location.workers_count)
            total_reward += location.duration * location.p * (location.p + 5)
            location_stats[location.p] += 1

    print({
        'reward': total_reward,
        'workers': total_worker_count,
        'locations': sum(location_stats.values()),
        'distance': total_distance,
        'locations_by_p': location_stats,
    })


if __name__ == '__main__':
    main(sys.argv[1:])
