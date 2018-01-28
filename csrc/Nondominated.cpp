#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <memory>

using namespace std;

//
// Range
//
template <class V>
struct range {
    V left;
    V right;

    range(V l, V r): left(l), right(r) {}

    template <V>
    static range<V> thin(V value) {
        return range(value, value);
    }

    bool contains(V value) {
        return left <= value && value <= right;
    }

    bool operator==(range o){
        return left == o.left && right == o.right;
    }

};

//
// Utilities
//
template <class T>
T findMedian(vector<T>& list) {
    if (list.size() < 25) {
        sort(list.begin(), list.end());
        return list.get(list.size() / 2);
    }

    auto result = vector<T>();
    for (int i = 0; i + 5 < list.size(); i += 5) {
        list.sort(list.begin() + i, list.begin() + i + 5);
        result.push_back(list[i + 2]);
    }
    return findMedian(result);
}

template <class V>
vector<V> removeAdjacentDuplicates(vector<V>& list) {
    vector<V> unique = vector<V>();
    if (list.size() > 0) {
        unique.add(list.get(0));
        for (int i = 1; i < list.size(); i++) {
            if (list.get(i - 1) != list.get(i)) {
                unique.add(list.get(i));
            }
        }
    }
    return unique;
}

template <class V>
struct split_t {
    vector<V> L;
    vector<V> M;
    vector<V> R;

    split_t(vector<V>& l, vector<V>& m, vector<V>& r): L(l), M(m), R(r) {}
    split_t() {};
};

template <class V>
split_t<V> split(vector<V> list, V med, bool less(V, V)) {
    split_t<V> result = split_t<V>();
    for (V value : list) {
        if (less(value, med)) {
            result.L.push_back(value);
        } else if (less(med, value)) {
            result.R.push_back(value);
        } else {
            result.M.push_back(value);
        }
    }
    return result;
}

template <class V>
std::ostream &operator<<(std::ostream &os, split_t<V> const &split) {
    os << "L: ";
    for (auto it = split.L.begin(); it < split.L.end(); it++)
        os << *it << ", ";
    os << endl << "M: ";
    for (auto it = split.M.begin(); it < split.M.end(); it++)
        os << *it << ", ";
    os << endl << "R: ";
    for (auto it = split.R.begin(); it < split.R.end(); it++)
        os << *it << ", ";
    os << endl;

    return os;
}

template <class V>
struct cached {
    bool evaluated = false;
    V cache;
    function<V()> evaluate;

    cached(const function<V()>& eval): evaluate(eval) {};

    V operator()() {
        if (evaluated) {
            return cache;
        } else {
            auto res = evaluate();
            cache = res;
            evaluated = true;
            return res;
        }
    }
};


struct timer {
    clock_t start;

    timer(): start(clock()) {}

    double sec() {
        return (clock() - start) / (double) CLOCKS_PER_SEC;
    }

};

ostream& operator<<(ostream &os, timer t) {
    os << t.sec() << "s";
    return os;
}

int main(){
    timer timer;
    for (int i = 0; i < 1000000; i++) {
        for (int j = 0; j < 1000; j++);
    }
    cout << timer << endl;
    return 0;
}

