#ifndef MINMAX_H
#define MINMAX_H

template <typename t> t min(t a, t b){
        return a < b ? a : b;
}

template <typename t> t max(t a, t b){
        return a > b ? a : b;
}

#endif