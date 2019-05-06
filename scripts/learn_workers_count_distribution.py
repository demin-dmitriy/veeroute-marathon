#!/usr/bin/env python3

import sys
import json

import numpy as np

from math import sqrt
from pathlib import Path
from plumbum import local
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).absolute().parent
EVAL_DIR = SCRIPT_DIR.parent / 'eval'

generate = local[ str(SCRIPT_DIR / 'generate.py') ]
check = local[ str(SCRIPT_DIR / 'check.py') ]
solver = None


def segment(start, end):
    return range(start, end + 1)


def magic_formula(n):
    average_workers_required = (1. + 8.) / 2.
    average_duration = (5. + 30.) / 2.
    average_gap = 100. / (1. + np.sqrt(n))
    work_time_window_length = 800. - 200.
    return (n * average_workers_required * (average_duration + average_gap) - average_gap) / work_time_window_length;


def main(argv):
	assert len(argv) == 1
	global solver
	solver = local[argv[0]]

	try:
		table = [0] * 2001

		with local.cwd(str(EVAL_DIR)):
			start_search = 120

			for n in segment(500, 2000):
				(generate['--n', n] > 'task.in')()

				best_reward = -1
				best_workers_count = None
				no_improvement_streak = 0

				for workers_count in segment(start_search, 900):
					result = json.loads((solver['--task', 'task.in', '--workers-count', workers_count] | check['--task', 'task.in'])())
					reward = result['reward']

					if reward > best_reward:
						best_reward = reward
						best_workers_count = workers_count
						no_improvement_streak = 0
					else:
						no_improvement_streak += 1
						if no_improvement_streak > 15:
							break

				if best_workers_count == start_search or best_workers_count == 900:
					print('warning: maybe search segment is too narrow')
				print(f'n={n}, best={best_workers_count}, searched from {start_search}')
				start_search = best_workers_count - 20
				table[n] = best_workers_count

	finally:
		with (EVAL_DIR / 'best_workers_count.txt').open('a') as result_file:
			result_file.write(repr(table) + '\n')
		print(table)
		baseline = magic_formula(np.arange(2001))
		plt.plot(table)
		plt.plot(1.9 * baseline)
		plt.plot(baseline)
		plt.show()


if __name__ == '__main__':
	main(sys.argv[1:])
