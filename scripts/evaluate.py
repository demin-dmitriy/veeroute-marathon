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
    parser.add_argument('--message', '-m', required=False)
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


def run_solver(solver, task, result):
    start = time.time()
    checker_result_str = (solver['--task', task] | check['--task', task])()
    end = time.time()
    result.add(end - start, json.loads(checker_result_str))


def print_result_console(i, item, n, result):
    if i > 0:
        result.clear_previous_output()
        clear_last_lines(1)

    print(f'[{item} / {n}]')
    result.print()


def main(argv):
    args = parse_args(argv)
    solver = local[args.executable]

    EVAL_DIR.mkdir(exist_ok=True)

    result = Result()

    with local.cwd(str(EVAL_DIR)):
        if args.dataset is None:
            for i in range(args.count):
                (generate > 'task.in')()
                run_solver(solver, 'task.in', result)
                print_result_console(i, i, args.count, result)
        else:
            tasks = sorted(Path(args.dataset).glob('*.in'))
            for i, task in enumerate(tasks):
                run_solver(solver, task, result)
                print_result_console(i, str(task), len(tasks), result)


    if not args.no_log:
        with (EVAL_DIR / 'log.txt').open('a') as log:
            comment = f' # {args.message}' if args.message else ''
            count = f' count={args.count}' if not args.dataset else ''
            dataset = f' dataset={args.dataset}' if args.dataset else ''

            log.write(f'\n{time.strftime("%Y-%m-%d %H:%M")} {count}{dataset}{comment}\n')
            result.print(file=log)


if __name__ == '__main__':
    main(sys.argv[1:])
