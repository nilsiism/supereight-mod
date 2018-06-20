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

  template <typename KeyT>
    inline int unique_multiscale(KeyT* keys, int num_keys,
      const KeyT level_mask){
      int end = 1;
      const KeyT key_mask = ~level_mask; 
      if(num_keys < 2) return end;
      KeyT prev_key = keys[0] & key_mask;
      KeyT prev_level = keys[0] & level_mask;
      for (int i = 1; i < num_keys; ++i){
        const KeyT key = keys[i] & key_mask;
        const KeyT level = keys[i] & level_mask;
        if(key != prev_key){
          keys[end] = keys[i];
          ++end;
        } else if(level > prev_level) { 
          /* end does not advance but previous entry is overwritten */
          keys[end-1] = keys[i];
        }
        prev_key = key;
        prev_level = level;
      }
      return end;
    }

  template <typename KeyT>
    inline int unique_multiscale(KeyT* keys, int num_keys,
        const KeyT level_mask, const unsigned current_level){
      int end = 1;
      const KeyT key_mask = ~level_mask; 
      if(num_keys < 2) return end;
      KeyT prev_key = keys[0] & key_mask;
      KeyT prev_level = keys[0] & level_mask;
      for (int i = 1; i < num_keys; ++i){
        const KeyT key = keys[i] & key_mask;
        const KeyT level = keys[i] & level_mask;
        // std::cout << "Key = " << key << " level = " << level << " prev_key = "
        //   << prev_key << " prev_level = " << prev_level << std::endl;
        if(level >= current_level) {
          if(key != prev_key){
             keys[end++] = keys[i];
          } else if(level > prev_level) { 
            /* end does not advance but previous entry is overwritten */
             keys[end-1] = keys[i];
          }
          prev_key = key;
          prev_level = level;
        }
      }
      return end;
    }
}
#endif
