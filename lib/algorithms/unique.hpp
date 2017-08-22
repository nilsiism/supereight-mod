#ifndef UNIQUE_HPP
#define UNIQUE_HPP

namespace algorithms {
  template <typename T>
    inline int unique(T* keys, int num_keys){
      int end = 1;
      if(num_keys < 2) return end;
      for (int i = 1; i < num_keys; ++i){
        if(keys[i] != keys[i-1]){
          keys[end] = keys[i];
          ++end;
        }
      }
      return end;
    }
}
#endif
