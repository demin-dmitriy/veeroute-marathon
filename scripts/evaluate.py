#!/usr/bin/env python3

import sys
import json
import time

from math import sqrt
from pathlib import Path
from argparse import ArgumentParser
from plumbum import local
from statistics import mean, median, stdev


SCRIPT_DIR = Path(__file__).absolute().parent
EVAL_DIR = SCRIPT_DIR.parent / 'eval'

generate = local[ str(SCRIPT_DIR / 'generate.py') ]
check = local[ str(SCRIPT_DIR / 'check.py') ]


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('--dataset', required=False)
    parser.add_argument('--count', required=False, default=10, type=int)
    parser.add_argument('--no-log', required=False, action='store_true')
    return parser.parse_args(argv)


def stats(values):
    return {
        'min': min(values),
        'median': median(values),
        'max': max(values),
        'mean': round(mean(values), 2),
        'σ': round(stdev(values), 2) / sqrt(len(values)) if len(values) > 1 else 0
    }

def format_stats(d):
    return f'{d["min"]:<15.2f} {d["max"]:<15.2f} {d["mean"]:<15.2f} {d["σ"]:<15.2f}'


class Result:


    def __init__(self):
        self.elapsed = []
        self.reward = []
        self.workers = []
        self.distance = []
        self.locations = []

    def add(self, time, checker_result):
        self.elapsed.append(time)
        self.reward.append(checker_result['reward'])
        self.workers.append(checker_result['workers'])
        self.distance.append(checker_result['distance'])
        self.locations.append(checker_result['locations'])

    def print(self, file=sys.stdout):
        print('            {:<15} {:<15} {:<15} {:<15}'.format('min', 'max', 'mean', 'σ'), file=file)
        print('elapsed   |', format_stats(stats(self.elapsed)), file=file)
        print('reward    |', format_stats(stats(self.reward)), file=file)
        print('workers   |', format_stats(stats(self.workers)), file=file)
        print('distance  |', format_stats(stats(self.distance)), file=file)
        print('locations |', format_stats(stats(self.locations)), file=file)

    def clear_previous_output(self):
        clear_last_lines(6)


def clear_last_lines(count):
    assert count > 0
    print('\r', end='')
    for i in range(count):
        print('\033[1A\033[K', end='')


def main(argv):
    args = parse_args(argv)
    solver = local[args.executable]

    EVAL_DIR.mkdir(exist_ok=True)

    result = Result()

    with local.cwd(str(EVAL_DIR)):
        for i in range(args.count):
            (generate > 'test.in')()
            start = time.time()
            checker_result_str = (solver['--task', 'test.in'] | check['--task', 'test.in'])()
            end = time.time()
            result.add(end - start, json.loads(checker_result_str))

            if i > 0:
                result.clear_previous_output()
            result.print()

    if not args.no_log:
        with (EVAL_DIR / 'log.txt').open('a') as log:
            log.write(f'\n{time.strftime("%Y-%m-%d %H:%M")} count={args.count}\n')
            result.print(file=log)


if __name__ == '__main__':
    main(sys.argv[1:])
