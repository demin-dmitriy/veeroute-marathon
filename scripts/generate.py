#!/usr/bin/env python3

import sys

from argparse import ArgumentParser
from random import randint, choice

# See https://codeforces.com/contest/1160/problem/A1 on how tests are generated.

def segment(start, end):
	return range(start, end + 1)


def generate_all_possible_timewindows():
	return [
		(l, h)

		for l in segment(200, 800)
		for h in segment(l + 60, min(800, l + 300))
	]


POSSIBLE_TIMEWINDOWS = generate_all_possible_timewindows()


def rand_point(locations):
	assert len(locations) < 101 ** 2 / 2, r'There is must be more than 50% chance to find new location'

	while True:
		p = randint(0, 100), randint(0, 100)
		if p not in locations:
			return p


def parse_args(argv):
	parser = ArgumentParser()
	parser.add_argument('--n', required=False)
	return parser.parse_args(argv)


def main(argv):
	args = parse_args(argv)

	if args.n is None:
		n = randint(500, 2000)
	else:
		n = int(args.n)
		assert n > 0

	print(n)

	locations = set()

	# Generate base
	x, y = rand_point(locations)
	locations.add((x, y))
	print(x, y, 0, 0, 0, 0)

	# Generate all other points
	for i in range(n - 1):
		x, y = rand_point(locations)
		locations.add((x, y))
		d = randint(5, 30) # work duration
		p = randint(1, 7)  # worker count
		l, h = choice(POSSIBLE_TIMEWINDOWS)
		print(x, y, d, p, l, h)


if __name__ == '__main__':
	main(sys.argv[1:])
