#!/usr/bin/env python3

# based on example code from https://pymotw.com/2/multiprocessing/basics.html
import multiprocessing
import random
import time

"""
a method able to create a "worker" which in this case is a method that gets a number and a random timer 
once the timer is done running it will print that its finished.
"""
def worker(num):
    """A job that runs for a random amount of time between 5 and 10 seconds."""
    time.sleep(random.randrange(5,11))
    print('Worker:' + str(num) + ' finished')
    return

if __name__ == '__main__':
    jobs = []
    #creates 5 workers
    for i in range(5):
        #creates a procces that runs the worker method with the argument given to args.
        p = multiprocessing.Process(target=worker, args=(i,))
        #adds the new process p to the list of jobs.
        jobs.append(p)
        #stats the process
        p.start()
    #joins the jobs to make the script wait till all jobs are done before going to the printstatement below.
    for job in jobs:
        job.join()

    print('*** All jobs finished ***')
