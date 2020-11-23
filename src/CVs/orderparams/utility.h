#ifndef UTILITY_H
#define UTILITY_H

#include<vector>

// return vector containing integers ordered in range [start,end).

inline std::vector<int> range(int start, int end)
{
   std::vector<int> ret;
   ret.reserve(end - start);
   for (int i = start; i != end; ++i) {
      ret.push_back(i);
   }
   
   return ret;
}

// reindex vector in integers.  e.g. if indx1 = [0, 3, 5], it becomes
// [indx2[0], indx2[3], indx2[5].  This is used for getting the xtal
// particle indexes in absolute terms, rather than them being indexes
// into a list of xtal particles.

inline void reindex(std::vector<int>& indx1, const std::vector<int>& indx2)
{
   for (std::vector<int>::size_type i = 0; i != indx1.size(); ++i) {
      indx1[i] = indx2[indx1[i]];
   }
}

#endif
