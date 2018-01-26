import com.sun.istack.internal.Nullable;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;

/**
 * Created by Martoon.
 */
public class Nondominated {
    public static void main(String[] args) {

    }

    static <T> T findMedian(List<T> list) {
        if (list.size() < 25) {
            list.sort(null);
            return list.get(list.size() / 2);
        }

        ArrayList<T> result = new ArrayList<>();
        for (int i = 0; i + 5 < list.size(); i += 5) {
            List<T> sub = list.subList(i, i + 5);
            sub.sort(null);
            result.add(sub.get(2));
        }
        return findMedian(result);
    }

    /**
     * Segment tree.
     */
    static class SegmentTree<K extends Comparable<K>, V, M> {
        private final Function<V, K> getKey;
        private final Function<V, M> getSummary;
        private final BinaryOperator<M> addSummaries;
        private Tree tree;

        public static <K extends Comparable<K>, V, M> SegmentTree<K, V, M>
        mkSegmentTree(Function<V, K> getKey, Function<V, M> getSummary, BinaryOperator<M> addSummaries, List<V> values) {
            List<K> keys = new ArrayList<>();
            for (V value : values) {
                keys.add(getKey.apply(value));
            }
            keys.sort(null);
            List<K> uniqKeys = Util.removeAdjacentDuplicates(keys);
            uniqKeys.add(keys.get(0));
            for (int i = 1; i < keys.size(); i++) {
                if (keys.get(i - 1) != keys.get(i)) {
                    uniqKeys.add(keys.get(i));
                }
            }
            return new SegmentTree<>(getKey, getSummary, addSummaries).buildTree(uniqKeys);
        }

        private SegmentTree<K, V, M> buildTree(List<K> sortedKeys) {
            tree = buildTreeImpl(sortedKeys);
            return this;
        }

        private Tree buildTreeImpl(List<K> keys) {
            assert keys.size() != 0;
            if (keys.size() == 1) {
                return new Leaf(Range.thin(keys.get(0)));
            } else {
                Range<K> range = new Range<>(keys.get(0), keys.get(keys.size() - 1));
                int m = keys.size() / 2;
                K med = keys.get(m);
                List<K> leftKeys = keys.subList(0, m);
                List<K> rightKeys = keys.subList(m, keys.size());
                return new Node(buildTreeImpl(leftKeys), buildTreeImpl(rightKeys), k -> k.compareTo(med) < 0, range);
            }
        }

        private SegmentTree(Function<V, K> getKey, Function<V, M> getSummary, BinaryOperator<M> addSummaries) {
            this.getKey = getKey;
            this.getSummary = getSummary;
            this.addSummaries = addSummaries;
        }

        private abstract class Tree {
            public Range<K> range;
            @Nullable
            public M summary;

            public Tree(Range<K> range) {
                this.range = range;
            }

            public abstract M request(Range<K> reqRange);

            public abstract void modify(K key, Consumer<Set<V>> modifier);
        }

        private class Leaf extends Tree {
            private final Set<V> values = new TreeSet<>();
            public Leaf(Range<K> range) {
                super(range);
            }

            @Override
            public @Nullable
            M request(Range<K> reqRange) {
                if (reqRange.contains(range.left)) {
                    return summary;
                } else {
                    return null;
                }
            }

            @Override
            public void modify(K key, Consumer<Set<V>> modifier) {
                assert this.range.contains(key);
                modifier.accept(this.values);
                updateSummary();
            }

            private void updateSummary() {
                this.summary = null;
                for (V value : values) {
                    M summary2 = getSummary.apply(value);
                    this.summary = this.summary == null
                            ? summary2
                            : addSummaries.apply(this.summary, summary2);
                }
            }
        }

        private class Node extends Tree {
            public final Tree left;
            public final Tree right;
            public final Predicate<K> isLeft;

            public Node(Tree left, Tree right, Predicate<K> isLeft, Range<K> range) {
                super(range);
                this.left = left;
                this.right = right;
                this.isLeft = isLeft;
            }

            @Override
            public M request(Range<K> reqRange) {
                if (reqRange.equals(range)) {  // totally matches
                    return summary;
                } else if (reqRange.right.compareTo(right.range.left) < 0) {  // go left
                    return left.request(reqRange);
                } else if (reqRange.left.compareTo(left.range.right) > 0) {  // go right
                    return right.request(reqRange);
                } else {
                    M leftRes = left.request(new Range<>(reqRange.left, this.range.right));
                    M rightRes = right.request(new Range<>(this.range.left, reqRange.right));
                    return addSummaries.apply(leftRes, rightRes);
                }
            }

            @Override
            public void modify(K key, Consumer<Set<V>> modifier) {
                if (left.range.contains(key)) {
                    left.modify(key, modifier);
                } else {
                    right.modify(key, modifier);
                }
                this.summary = addSummaries.apply(left.summary, right.summary);
            }
        }
    }

    static class Range<V extends Comparable<V>> {
        public final V left;
        public final V right;

        public Range(V left, V right) {
            this.left = left;
            this.right = right;
        }

        public static <V extends Comparable<V>> Range<V> thin(V value) {
            return new Range<>(value, value);
        }

        public boolean contains(V value) {
            return left.compareTo(value) < 0 && value.compareTo(right) < 0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Range<?> range = (Range<?>) o;

            return left.equals(range.left) && right.equals(range.right);
        }
    }

    /**
     * Helper in benchmarking.
     */
    static class Timer implements AutoCloseable {
        private final String name;
        private final long startTime;

        public Timer(String name) {
            this.name = name;
            this.startTime = System.currentTimeMillis();
        }

        public long getMillis() {
            return System.currentTimeMillis() - startTime;
        }

        public double getSeconds() {
            return getMillis() / 1000D;
        }

        public void printTime() {
            System.out.printf("%s passed for %s seconds%n", name, getSeconds());
        }

        @Override
        public void close() {
            printTime();
        }
    }

    /**
     * Class with utility functions.
     */
    static class Util {
        public static <V> List<V> removeAdjacentDuplicates(List<V> list) {
            ArrayList<V> unique = new ArrayList<>();
            unique.add(list.get(0));
            for (int i = 1; i < list.size(); i++) {
                if (list.get(i - 1) != list.get(i)) {
                    unique.add(list.get(i));
                }
            }
            return unique;
        }
    }
}
