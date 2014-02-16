from random import Random
from time import time
import inspyred

def main(prng=None, display=False):
    if prng is None:
        prng = Random()
        prng.seed(time())

    problem = inspyred.benchmarks.Binary(inspyred.benchmarks.Schwefel(2),
                                         dimension_bits=30)

    ea = inspyred.ec.GA(prng)
    ea.terminator = inspyred.ec.terminators.evaluation_termination
    final_pop = ea.evolve(generator=problem.generator,
                          evaluator=problem.evaluator,
                          pop_size=100,
                          bounder=problem.bounder,
                          maximize=problem.maximize,
                          max_evaluations=30000)

    if display:
        best = max(final_pop)
        print 'Best Solution: \n{}'.format(str(best))
    return ea

if __name__ == '__main__':
    main(display=True)