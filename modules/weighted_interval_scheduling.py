
## Modified code from  https://github.com/alanmc-zz/python-interval-scheduling/blob/master/python-weighted-interval-scheduling


import collections 
import bisect

"""
This algorithm returns the weight of the maximally weighted set of independent intervals
in an WeightedIntervalProblem graph.  Commonly used for scheduling optimization problems.

The algorithm works as follows:

    Imagine this scenario:  You;ve travelling to SV to raise money for you fledging startup,
    and need to schedule investor meetings during your stay to maximize the money you will potentially raise.
    Every potential investor meeting you can take has a start and end time, and you can only be at
    one meeting at a time.  In addition to start and end times, each meeting has value associated with it
    that corresponds to the expected value of money you will raise by taking the meeting.  

    Here's a trivial case (notice that the weight of the interval is independent of its length.):


      *----- 500K -----*    *--- 600K ---*       *-------- 100K --------*
    12pm              1pm  2pm         2:30pm   4pm                    8pm


    In this trival case, the maximall weighting is simply the sum of all weights, since all intervals are independent
    and can be scheduled together.

    In non-trivial cases, the intervals will overlap and tough decisions need to be made.
      
The algorithm's runtime is O(n log n) in the worst case
"""

class Interval(object):
    """docstring for Interval"""
    def __init__(self, start,end,weight,node):
        super(Interval, self).__init__()
        self.start = start
        self.end = end
        self.weight = weight
        self.node = node
    def __str__(self):
        return str(str(self.node)+', weight:'+str(self.weight)+', coord: ('+str(self.start)+','+str(self.end) +')')
        

class WeightedIntervalProblem(object):
    """docstring for WeightedIntervalProblem"""
    def __init__(self,startnode):
        super(WeightedIntervalProblem, self).__init__()
        self.startnode = startnode
        self.score = None
        self.optimal_path = []
        self.intervals = []

    def add_interval(self,interval_object):
        self.intervals.append(interval_object) 

    def weighted_interval_scheduling(self,overlap_parameter):
        '''
        Input a graph G, whose structure is a list of Interval instances.

        --- Doctest ---

        >>> G = [Interval(43,70,27,'a'),Interval(3,18,24,'b'),Interval(65,99,45,'c'),Interval(20,39,26,'d'),Interval(45,74,26,'e'),Interval(10,28,20,'f'),Interval(78,97,23,'g'),Interval(0,9,22,'h')]
        >>> i = WeightedIntervalProblem('startnode')
        >>> for x in G: i.add_interval(x)
        >>> i.weighted_interval_scheduling()
        100
        >>> G = [Interval(1,5,10,'s1'),Interval(6,10,12,'s2'),Interval(1,10,15,'s3')]
        >>> i = WeightedIntervalProblem('startnode')
        >>> for x in G: i.add_interval(x)
        >>> i.weighted_interval_scheduling()
        22
        '''
        # print([str(i) for i in self.intervals])
        G = self.intervals
        S = collections.defaultdict(int)
        # ASSUME ALREADY SORTED INPUT
        G.sort(key = lambda x: x.end)
        print([str(i) for i in G])
        start = [x.start for x in G]  
        end = [x.end for x in G]

        S[0] = 0
        S[1] = G[0].weight if G[0].start >= - overlap_parameter else 0
        for i in range(2, len(G)+1):
            S[i] = max(S[i-1], G[i-1].weight + S[ bisect.bisect_left(end , start[i-1] + overlap_parameter) ])

        self.score = S[len(G)]
        print("optimum score is:", self.score, S)
        def find_optimal_path(j):
            """ Recursive function to track down the path that gave rise to the optimal solution
            """
            if j == 0:
                return
            else:
                # print("here", j)
                if G[j-1].weight + S[ bisect.bisect_left(end, start[j-1] + overlap_parameter) ] >= S[j-1]:
                    self.optimal_path.append(self.intervals[j-1])
                    find_optimal_path(bisect.bisect_left(end, start[j-1] + overlap_parameter))
                    print("here", j, str(self.optimal_path[0]))
                else:
                    return find_optimal_path(j-1)

            return [ (i.start, i.end, i.weight, i.node) for i in self.optimal_path]

        opt = find_optimal_path(len(G))
        print(opt)
        return opt

if __name__ == "__main__":
    import doctest
    doctest.testmod()