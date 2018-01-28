import com.sun.istack.internal.Nullable;

import java.io.*;
import java.util.*;
import java.util.function.*;
import java.util.stream.Collectors;

/**
 * Created by Martoon.
 */
public class Nondominated {

    public static void main(String[] args) throws IOException {
        SacSolver.optimizationEdge = 50;

//        interactAndSolve(new File("sorting.in"), new File("sorting.out"));

//        TestSort.seriesTestSolver(4, 25, 100000);
//        TestSort.testSolver(19, 3, 3);

        List<Point> points = TestSort.genPoints(234, 4, 100000);
        try {
            Thread.sleep(10000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("Starting");
        new Timer("Sort").invoke(() -> new SacSolver(points).solve());

//        List<Point> points = TestSort.genPoints(234, 2, 3);
//        for (int i = 0; i < 1000; i += 10) {
//            System.out.println("For " + i + ":");
//            SacSolver.optimizationEdge = i;
//            new Timer("Sort").invoke(5, () -> new SacSolver(points).solve());
//        }

//        for (int i = 1; i <= 500; i++) {
//            List<Point> points = TestSort.genPoints(234, 4, i);
//            System.out.println(i + " points");
//            new Timer("Dumb").invoke(1000, () -> new DumbSolver(points).solve());
//            new Timer("Sac").invoke(1000, () -> new SacSolver(points).solve());
//        }

    }

    public static void interactAndSolve(File input, File output) throws IOException {
        List<Point> points = readPoints(input);

        Solver solver = new SacSolver(points);
        List<Integer> ranks = solver.solve();

        try (PrintStream out = new PrintStream(output)) {
            for (Integer rank : ranks) {
                out.print(rank + " ");
            }
        }
    }

    public static void dumpInput(File file, List<Point> points) throws FileNotFoundException {
        try (PrintStream out = new PrintStream(file)) {
            out.print(points.size());
            out.print(" ");
            out.print(points.get(0).coords.length);
            out.println();

            for (Point point : points) {
                for (Integer coord : point.coords) {
                    out.print(coord);
                    out.print(" ");
                }
                out.println();
            }
        }
    }

    public static List<Point> readPoints(File file) throws IOException {
        ArrayList<Point> points = new ArrayList<>();
        try (BufferedReader in = new BufferedReader(new FileReader(file))) {
            String[] header = in.readLine().split("\\s+");
            int n = Integer.parseInt(header[0]);
            int d = Integer.parseInt(header[1]);

            String line;
            while ((line = in.readLine()) != null) {
                ArrayList<Integer> coord = new ArrayList<>();
                for (String s : line.split("\\s+")) {
                    coord.add(Integer.parseInt(s));
                }
                points.add(new Point(coord));
            }
        }
        return points;
    }

    /**
     * Type aliases
     */
    static class Point {
        private static final int COORD_X = 0;
        private static final int COORD_Y = 1;

        private final Integer[] coords;

        public Point(List<Integer> coords) {
            this.coords = coords.toArray(new Integer[coords.size()]);
        }

        public int getDimensions() {
            return coords.length;
        }

        public boolean dominates(Point o) {
            boolean meetStrict = false;
            for (int i = 0; i < this.coords.length; i++) {
                int cmp = this.coords[i].compareTo(o.coords[i]);
                if (cmp > 0) {
                    return false;
                } else if (cmp < 0) {
                    meetStrict = true;
                }
            }
            return meetStrict;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Point point = (Point) o;

            return Arrays.equals(coords, point.coords);
        }

        @Override
        public String toString() {
            return Arrays.stream(coords)
                    .map(Object::toString)
                    .collect(Collectors.joining(", ", "(", ")"));
        }
    }

    /**
     * Solution
     */

    static abstract class Solver {
        protected final Point[] points;
        protected final Integer[] ranks;

        public Solver(List<Point> points) {
            this(points.toArray(new Point[0]), Collections.nCopies(points.size(), 0).toArray(new Integer[0]));
        }

        public Solver(Point[] points, Integer[] ranks) {
            this.points = points;
            this.ranks = ranks;

            if (points.length == 0)
                throw new IllegalArgumentException("Can't solve for no points");
        }

        public Integer getRank(int i) {
            return ranks[i];
        }

        public void updateRank(int i, int newRank) {
            ranks[i] = Math.max(ranks[i], newRank);
        }

        protected abstract void solveForIds(List<Integer> ids);

        public List<Integer> solve() {
            ArrayList<Integer> ids = new ArrayList<>();
            for (int i = 0; i < points.length; i++) {
                ids.add(i);
            }
            solveForIds(ids);

            return Arrays.asList(ranks);
        }
    }

    @SuppressWarnings("StatementWithEmptyBody")
    static class SacSolver extends Solver {
        private static int optimizationEdge = 1;

        public SacSolver(List<Point> points) {
            super(points);
        }

        @Override
        protected void solveForIds(List<Integer> ids) {
            ids.sort(lexSortComparator);
            sortSAC(points[0].getDimensions(), ids);
        }

        private void sortSAC(int dimension, List<Integer> ids) {
//            System.out.println("sortSAC " + dimension + " " + ids.size());
            assert optimizationEdge > 0;
            if (dimension <= 1) {
                throw new IllegalArgumentException("Can't solve for so few dimensions");
            } else if (ids.size() <= optimizationEdge) {
                new DumbSolver(points, ranks).sortDumb(ids);
            } else if (dimension == 2) {
                sortSweep(ids);
            } else {
                // Helpers
                Function<Integer, Integer> getPointCoord = id -> points[id].coords[dimension - 1];

                // Get median
                ArrayList<Integer> zs = new ArrayList<>();
                for (Integer id : ids) {
                    zs.add(getPointCoord.apply(id));
                }
                Integer median = Util.findMedian(zs);

                // Split
                Function<Integer, Integer> comparingToMedian = id -> getPointCoord.apply(id).compareTo(median);
                Util.SplitResult<Integer> split = Util.split(ids, comparingToMedian);

                // Call recursively
                sortSAC(dimension, split.L);
                updateSAC(dimension, split.L, split.M);
                sortSAC(dimension - 1, split.M);
                updateSAC(dimension, split.L, split.R);
                updateSAC(dimension, split.M, split.R);
                sortSAC(dimension, split.R);
            }

        }

        private void updateSAC(int dimension, List<Integer> known, List<Integer> request) {
//            System.out.println("updateSAC " + dimension + " " + known.size() + "/" + request.size());
            assert optimizationEdge >= 0;
            if (dimension <= 1) {
                throw new IllegalArgumentException("Can't solve for so few dimensions");
            } else if (request.isEmpty() || known.isEmpty()) {
            } else if (known.size() + request.size() <= optimizationEdge) {
                new DumbSolver(points, ranks).updateDumb(known, request);
            } else if (dimension == 2) {
                updateSweep(known, request);
            } else {

                // Helpers
                Function<Integer, Integer> getPointCoord = id -> points[id].coords[dimension - 1];

                // Get median
                ArrayList<Integer> zs = new ArrayList<>();
                for (Integer req : request) {
                    zs.add(getPointCoord.apply(req));
                }
                Integer median = Util.findMedian(zs);

                // Split
                Function<Integer, Integer> comparingToMedian = id -> getPointCoord.apply(id).compareTo(median);
                Util.SplitResult<Integer> knownSplit = Util.split(known, comparingToMedian);
                Util.SplitResult<Integer> requestSplit = Util.split(request, comparingToMedian);

                // Recursive calls
                updateSAC(dimension, knownSplit.L, requestSplit.L);
                updateSAC(dimension - 1, knownSplit.L, requestSplit.M);
                updateSAC(dimension - 1, knownSplit.M, requestSplit.M);
                updateSAC(dimension - 1, knownSplit.L, requestSplit.R);
                updateSAC(dimension - 1, knownSplit.M, requestSplit.R);
                updateSAC(dimension, knownSplit.R, requestSplit.R);
            }

        }

        private void sortSweep(List<Integer> ids) {
//            System.out.println("sortSweep " + ids.size());
            if (ids.isEmpty())
                return;

            // line :: (Y, Rank) -> Id
            TreeMap<Integer, Integer> line = new TreeMap<>();
            // tree :: SegmentTree<K = Y, V = Id, M = max Rank>
            SegmentTree<Integer, Integer, Integer> ranksTree = makeRanksTree(ids);

            Integer prevId = null;
            for (Integer id : ids) {
                if (prevId != null && points[prevId].equals(points[id])) {
                    updateRank(id, getRank(prevId));
                } else {
                    Integer key = getY(id);
                    Map.Entry<Integer, Integer> higherEntry = line.higherEntry(key);
                    Integer higherPoint = higherEntry == null ? null : higherEntry.getValue();

                    int newRank = getNextRankAfter(ranksTree, id);

                    if (higherPoint != null && newRank == getRank(higherPoint)) {
                        removePoint(line, ranksTree, null, higherPoint);
                    }

                    updateRank(id, newRank);
                    insertPoint(line, ranksTree, null, id);
                }
                prevId = id;
            }
        }

        private void updateSweep(List<Integer> knowns, List<Integer> requests) {
//            System.out.println("updateSweep " + knowns.size() + "/" + requests.size());
            if (knowns.isEmpty() || requests.isEmpty())
                return;

            // line :: (Y, Rank) -> Id
            TreeMap<Integer, Integer> line = new TreeMap<>();
            // tree :: SegmentTree<K = Y, V = Id, M = max Rank>
            SegmentTree<Integer, Integer, Integer> ranksTree = makeRanksTree(knowns);
            // front :: Rank -> @Nullable Id
            Integer[] front = new Integer[points.length];

            int iKnown = 0;
            for (Integer req : requests) {
                while (iKnown < knowns.size() && lexSortComparator.compare(knowns.get(iKnown), req) <= 0) {
                    Integer known = knowns.get(iKnown);
                    iKnown++;

                    // processing known point
                    int rank = getRank(known);
                    Integer oldKnown = front[rank];
                    if (oldKnown == null) {
                        insertPoint(line, ranksTree, front, known);
                    } else if (getY(oldKnown) > getY(known)) {
                        removePoint(line, ranksTree, front, oldKnown);
                        insertPoint(line, ranksTree, front, known);
                    }
                }

                // processing request point
                int newRank = getNextRankAfter(ranksTree, req);
                updateRank(req, newRank);
            }
        }

        private SegmentTree<Integer, Integer, Integer> makeRanksTree(List<Integer> ids) {
            return SegmentTree.make(this::getY, this::getRank, Math::max, ids);
        }

        private int getNextRankAfter(SegmentTree<Integer, Integer, Integer> ranksTree, Integer id) {
            Integer maxRank = ranksTree.request(new Range<>(Integer.MIN_VALUE, getY(id)));
            return maxRank == null ? 0 : maxRank + 1;
        }

        private void insertPoint(TreeMap<Integer, Integer> line, SegmentTree<Integer, Integer, Integer> ranksTree, Integer[] front, Integer id) {
            if (line != null) {
                line.put(getY(id), id);
            }
            if (ranksTree != null) {
                ranksTree.put(id);
            }
            if (front != null) {
                front[getRank(id)] = id;
            }
        }

        private void removePoint(TreeMap<Integer, Integer> line, SegmentTree<Integer, Integer, Integer> ranksTree, Integer[] front, Integer id) {
            if (line != null) {
                line.remove(getY(id));
            }
            if (ranksTree != null) {
                ranksTree.remove(id);
            }
            if (front != null) {
                front[getRank(id)] = null;
            }
        }

        private Integer getY(Integer id) {
            return points[id].coords[Point.COORD_Y];
        }

        private final Comparator<Integer> lexSortComparator = (id1, id2) -> {
            Point p1 = points[id1];
            Point p2 = points[id2];
            int x1 = p1.coords[Point.COORD_X];
            int y1 = p1.coords[Point.COORD_Y];
            int x2 = p2.coords[Point.COORD_X];
            int y2 = p2.coords[Point.COORD_Y];
            return x1 == x2 ? y1 - y2 : x1 - x2;
        };

    }

    static class DumbSolver extends Solver {
        public DumbSolver(List<Point> points) {
            super(points);
        }

        public DumbSolver(Point[] points, Integer[] ranks) {
            super(points, ranks);
        }

        @Override
        protected void solveForIds(List<Integer> ids) {
            sortDumb(ids);
        }

        public void sortDumb(List<Integer> ids) {
            if (ids.isEmpty())
                return;

            Util.Cached<Integer>[] rankEvals = new Util.Cached[points.length];
            for (Integer req : ids) {
                Util.Cached<Integer> eval = new Util.Cached<>(() -> {
                    int[] maxRank = {-1};
                    for (Integer id : ids) {
                        if (points[id].dominates(points[req])) {
                            maxRank[0] = Math.max(maxRank[0], rankEvals[id].get());
                        }
                    }
                    return maxRank[0] + 1;
                });
                rankEvals[req] = eval;
            }

            for (Integer id : ids) {
                Util.Cached<Integer> rankEval = rankEvals[id];
                updateRank(id, rankEval.get());
            }
        }

        public void updateDumb(List<Integer> known, List<Integer> request) {
            if (known.isEmpty())
                return;

            for (Integer req : request) {
                for (Integer kn : known) {
                    if (points[kn].dominates(points[req])) {
                        ranks[req] = Math.max(ranks[req], ranks[kn] + 1);
                    }
                }
            }
        }
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
        make(Function<V, K> getKey, Function<V, M> getSummary, BinaryOperator<M> addSummaries, List<V> values) {
            if (values.isEmpty())
                throw new IllegalArgumentException("Can't build empty segment tree");

            List<K> keys = new ArrayList<>();
            for (V value : values) {
                keys.add(getKey.apply(value));
            }
            keys.sort(null);
            List<K> uniqKeys = Util.removeAdjacentDuplicates(keys);
            BinaryOperator<M> addSummariesSafe = (a, b) ->
                    a == null ? b
                            : b == null ? a
                            : addSummaries.apply(a, b);
            return new SegmentTree<>(getKey, getSummary, addSummariesSafe).buildTree(uniqKeys);
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

        public M request(Range<K> range) {
            return tree.request(range);
        }

        public void put(V value) {
            tree.modify(getKey.apply(value), s -> s.add(value));
        }

        public void remove(V value) {
            tree.modify(getKey.apply(value), s -> s.remove(value));
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

            public abstract void forEach(Consumer<V> consumer);
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

            @Override
            public void forEach(Consumer<V> consumer) {
                values.forEach(consumer);
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
                } else if (reqRange.left.compareTo(range.right) > 0 || range.left.compareTo(reqRange.right) > 0) {  // totally in nowhere
                    return null;
                } else if (reqRange.right.compareTo(right.range.left) < 0) {  // go left
                    return left.request(reqRange);
                } else if (reqRange.left.compareTo(left.range.right) > 0) {  // go right
                    return right.request(reqRange);
                } else {  // go in both directions
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

            @Override
            public void forEach(Consumer<V> consumer) {
                left.forEach(consumer);
                right.forEach(consumer);
            }
        }

        @Override
        public String toString() {
            ArrayList<V> values = new ArrayList<>();
            tree.forEach(values::add);
            return values.toString();
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
            return left.compareTo(value) <= 0 && value.compareTo(right) <= 0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Range<?> range = (Range<?>) o;

            return left.equals(range.left) && right.equals(range.right);
        }

        @Override
        public String toString() {
            return String.format("(%s, %s)", left, right);
        }
    }

    /**
     * Helper in benchmarking.
     */
    @SuppressWarnings("unused")
    static class Timer {
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

        public void invoke(Runnable run) {
            invoke(1, run);
        }

        public void invoke(int times, Runnable run) {
            for (int i = 0; i < times; i++) {
                run.run();
            }
            long cur = System.currentTimeMillis();
            double duration = (double) (cur - startTime) / times;
            System.out.printf("%s passed for %f seconds%n", name, duration / 1000D);
        }
    }

    /**
     * Class with utility functions.
     */
    static class Util {
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

        public static <V> List<V> removeAdjacentDuplicates(List<V> list) {
            ArrayList<V> unique = new ArrayList<>();
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

        public static <V> List<V> concat(List<V> list1, List<V> list2) {
            ArrayList<V> result = new ArrayList<>(list1);
            result.addAll(list2);
            return result;
        }

        public static class SplitResult<V> {
            public final List<V> L = new ArrayList<>();
            public final List<V> M = new ArrayList<>();
            public final List<V> R = new ArrayList<>();

            @Override
            public String toString() {
                return String.format("\n  {%s,\n  %s,\n  %s\n  }", L, M, R);
            }
        }

        public static <V> SplitResult<V> split(List<V> list, Function<V, Integer> comparator) {
            SplitResult<V> result = new SplitResult<>();
            for (V value : list) {
                int cmp = comparator.apply(value);
                List<V> toAdd = cmp < 0 ? result.L
                        : cmp == 0 ? result.M
                        : result.R;
                toAdd.add(value);
            }
            return result;
        }

        public static class Cached<V> implements Supplier<V> {
            private V cache = null;
            private final Supplier<V> evaluate;

            public Cached(Supplier<V> evaluate) {
                this.evaluate = evaluate;
            }

            @Override
            public V get() {
                return cache != null ? cache : (cache = evaluate.get());
            }
        }

        public static <V> List<V> merge(List<V> list1, List<V> list2, Comparator<V> cmp) {
            ArrayList<V> result = new ArrayList<>(list1.size() + list2.size());
            int i = 0;
            int j = 0;
            while (true) {
                if (i >= list1.size()) {
                    while (j < list2.size())
                        result.add(list2.get(j));
                    break;
                }
                if (j >= list2.size()) {
                    while (i < list1.size())
                        result.add(list1.get(i));
                    break;
                }
                if (cmp.compare(list1.get(i), list2.get(j)) <= 0) {
                    result.add(list1.get(i));
                    i++;
                } else {
                    result.add(list2.get(j));
                    j++;
                }
            }
            return result;
        }
    }

    static class Gen {
        protected final Random random;

        public Gen() {
            this.random = new Random();
        }

        public Gen(int seed) {
            this.random = new Random(seed);
        }

        protected Supplier<Integer> integer(int bound) {
            return () -> random.nextInt(2 * bound + 1) - bound;
        }

        public static <V> Supplier<List<V>> vectorOf(int k, Supplier<V> gen) {
            return () -> {
                ArrayList<V> res = new ArrayList<>();
                for (int i = 0; i < k; i++) {
                    res.add(gen.get());
                }
                return res;
            };
        }

        protected <V> Supplier<List<V>> sublistOf(Supplier<List<V>> gen) {
            return () -> {
                List<V> list = gen.get();
                ArrayList<V> res = new ArrayList<>();
                for (V v : list) {
                    if (random.nextBoolean()) {
                        res.add(v);
                    }
                }
                return res;
            };
        }

        protected <V extends Comparable<V>> Supplier<Range<V>> rangeOf(Supplier<V> gen) {
            return () -> new Range<>(gen.get(), gen.get());
        }

        protected static <A, B> Supplier<B> map(Function<A, B> mapper, Supplier<A> sup) {
            return () -> mapper.apply(sup.get());
        }
    }

    @SuppressWarnings("unused")
    static class TestSegmentTree {
        public static void testSum(int seed, int n) {
            new Gen(seed) {{
                int inputBound = Math.max(10, n * 3);
                List<Integer> inputUnsafe = Gen.vectorOf(n, integer(inputBound)).get();
                List<Integer> input = inputUnsafe.stream().distinct().collect(Collectors.toList());
                List<Integer> chosen = sublistOf(() -> input).get();
                List<Range<Integer>> requests = vectorOf(1, rangeOf(integer(inputBound + 1))).get();

                SegmentTree<Integer, Integer, Integer> tree =
                        SegmentTree.make(Function.identity(), Function.identity(), Integer::sum, input);
                for (Integer v : chosen) {
                    tree.put(v);
                }
                ArrayList<Integer> segTreeResults = new ArrayList<>();
                for (Range<Integer> request : requests) {
                    segTreeResults.add(tree.request(request));
                }

                ArrayList<Integer> sortedChosen = new ArrayList<>(chosen);
                Collections.sort(sortedChosen);
                ArrayList<Integer> dumbResults = new ArrayList<>();
                for (Range<Integer> request : requests) {
                    Integer result = null;
                    for (Integer value : sortedChosen) {
                        if (request.contains(value)) {
                            result = result == null ? 0 : result;
                            result += value;
                        }
                    }
                    dumbResults.add(result);
                }

                assert dumbResults.size() == segTreeResults.size();
                if (dumbResults.equals(segTreeResults)) {
                    System.out.println("OK!");
                } else {
                    for (int i = 0; i < dumbResults.size(); i++) {
                        Integer expected = dumbResults.get(i);
                        Integer got = segTreeResults.get(i);
                        if (!expected.equals(got)) {
                            System.out.println("Test failure");
                            System.out.println("For all points: " + input);
                            System.out.println("For chosen: " + chosen);
                            System.out.println("For request: " + requests.get(i));
                            System.out.println("Expected " + expected + ", got " + got);
                            throw new AssertionError("Test failed");
                        }
                    }
                    throw new IllegalStateException("What failed eventually??");
                }
            }};
        }

        public static void seriesTestSum(int times, int n) {
            for (int i = 1; i <= times; i++) {
                System.out.println("Iteration #" + i);
                TestSegmentTree.testSum(i, n);
            }
        }

        public static void seriesTestSum() {
            seriesTestSum(100000, 10000);
        }

    }

    public static class TestSort {
        public static List<Point> genPoints(int seed, int d, int n) {
            return new Gen(seed) {
                List<Point> points = vectorOf(n, map(Point::new, vectorOf(d, integer(100)))).get();
            }.points;
        }

        public static void testSolver(int seed, int d, int n) {
            List<Point> points = genPoints(seed, d, n);

            Solver dumbSolver = new DumbSolver(points);
            List<Integer> dumbAns = dumbSolver.solve();

            Solver sacSolver = new SacSolver(points);
            List<Integer> sacAns = sacSolver.solve();

            if (dumbAns.equals(sacAns)) {
                System.out.println("Solved OK!");
            } else {
                for (int i = 0; i < dumbAns.size(); i++) {
                    if (!Objects.equals(dumbAns.get(i), sacAns.get(i))) {
                        System.out.println("Test failure");
                        System.out.println("For input: " + points);
                        System.out.println("All ranks expected: " + dumbAns);
                        System.out.println("All ranks got     : " + sacAns);
                        System.out.printf("For %dth value expected %s, but got %s%n", i, dumbAns.get(i), sacAns.get(i));

                        try {
                            Thread.sleep(10);
                        } catch (InterruptedException ignored) {
                        }
                        throw new AssertionError("Test failed");
                    }
                }

                throw new IllegalStateException("And where is bad answer??");
            }

        }

        public static void seriesTestSolver(int d, int n, int times) {
            for (int i = 0; i < times; i++) {
                System.out.println("Iteration #" + i);
                testSolver(i, d, n);
            }
        }
    }

}
