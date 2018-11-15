#!/usr/bin/env python3

# based on example code from https://pymotw.com/2/multiprocessing/basics.html
import multiprocessing
import random
import time

def worker(num):
    """A job that runs for a random amount of time between 5 and 10 seconds."""
    time.sleep(random.randrange(5,11))
    print('Worker:' + str(num) + ' finished')
    return

if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = multiprocessing.Process(target=worker, args=(i,))
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()

    print('*** All jobs finished ***')
