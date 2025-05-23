/**
 * This file is part of TreeCmp, a tool for comparing phylogenetic trees using
 * the Matching Split distance and other metrics. Copyright (C) 2011, Damian
 * Bogdanowicz
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
package treecmp.metric;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import pal.tree.Node;
import pal.tree.NodeFactory;
import pal.tree.SimpleTree;
import pal.tree.Tree;
import treecmp.common.LapSolver;
import treecmp.common.TreeCmpUtils;

/**
 * UMAST metric.
 * Implementation of Procedure 1
 * Farach, Martin and Thorup, Mikkel; Fast comparison of evolutionary trees.
 */
public class UmastMetric extends BaseMetric implements Metric {
    @Override
    public boolean isRooted() {
        return false;
    }

    @Override
    public double getDistance(Tree t1, Tree t2) {
        Tree tmpT1 = t1.getCopy();
        Tree tmpT2 = t2.getCopy();
        int result = umast(new SimpleUnrootedTreePreservingNodeNumbers(tmpT1.getRoot()), new SimpleUnrootedTreePreservingNodeNumbers(tmpT2.getRoot()));

        return Math.max(tmpT1.getExternalNodeCount(), tmpT2.getExternalNodeCount()) - result;
    }

    private static int umast(final SimpleTree tree1, final SimpleTree tree2) {
        if (tree1.getExternalNodeCount() <= 3 || tree2.getExternalNodeCount() <= 3) {
            return getLeafLabelsIntersectionSize(tree1, tree2);
        }
        final int n1 = tree1.getExternalNodeCount();
        final int n2 = tree2.getExternalNodeCount();
        final int n = tree1.getExternalNodeCount() + tree2.getExternalNodeCount();
        final int k = Math.max((int)Math.pow(2.8, Math.sqrt(Math.log(n))), 2);
        final int t1MaxPartSize = tree1.getExternalNodeCount() / k;
        final int t2MaxPartSize = tree2.getExternalNodeCount() / k;
        final CoreTree core1 = new CoreTree(tree1, t1MaxPartSize);
        final CoreTree core2 = new CoreTree(tree2, t2MaxPartSize);

        //partition side trees into balanced side forests of sizes between n/2k and n/k
        final List<List<String>> sideForests1 = getBalancedSideForests(core1, t1MaxPartSize/2, t1MaxPartSize);
        final List<List<String>> sideForests2 = getBalancedSideForests(core2, t2MaxPartSize/2, t2MaxPartSize);

        int result = 0;

        final TreeRestricter r1 = new TreeRestricter(tree1);
        final TreeRestricter r2 = new TreeRestricter(tree2);

        //for all pairs (f1, f2) of opposing side forests: b = a(f1) + a(f2), umast = max(umast, umast(t1|b, t2|b)
        for (final List<String> sideForest1 : sideForests1) {
            for (final List<String> sideForest2 : sideForests2) {
                final Set<String> union = new LinkedHashSet<String>(sideForest1);
                union.addAll(sideForest2);

                final SimpleTree t1Restricted = r1.getRestrictedTo(union);
                final SimpleTree t2Restricted = r2.getRestrictedTo(union);

                result = Math.max(result, umast(t1Restricted, t2Restricted));
            }
        }

        //for all pairs (l1,l2) of opposing core leaves: compute CRMAST(t1^l1,t2^l2)
        final CRMASTSet crmastSet = new CRMASTSet(tree1, tree2);
        final int c1ExternalNodes = core1.getExternalNodeCount();
        final int c2ExternalNodes = core2.getExternalNodeCount();
        for (int i=0; i<c1ExternalNodes; i++) {
            tree1.reroot(tree1.getInternalNode(core1.getExternalNodeIdx(i)));
            for (int j=0; j<c2ExternalNodes; j++) {
                tree2.reroot(tree2.getInternalNode(core2.getExternalNodeIdx(j)));
                crmastSet.include(RMASTMetric.crmast(tree1, tree2), tree1, tree2);
            }
        }

        //for all pairs (v1,v2) of opposing core vertices: umast = max(umast, umast(v1,v2))
        final int c1Nodes = core1.getNodeCount();
        final int c2Nodes = core2.getNodeCount();
        for (int i=0; i<c1Nodes; i++) {
            final Node v1 = tree1.getInternalNode(core1.getNode(i));
            for (int j=0; j<c2Nodes; j++) {
                final Node v2 = tree2.getInternalNode(core2.getNode(j));
                result = Math.max(result, match(v1, v2, crmastSet));
            }
        }
        return result;
    }

    private static int getLeafLabelsIntersectionSize(Tree tree1, Tree tree2) {
        final Set<String> t1Labels = getLeafLabels(tree1);
        final Set<String> t2Labels = getLeafLabels(tree2);
        t1Labels.retainAll(t2Labels);
        return t1Labels.size();
    }

    private static Set<String> getLeafLabels(Tree tree) {
        final Set<String> labels = new HashSet<String>();
        for (int i=0; i<tree.getExternalNodeCount(); i++) {
            labels.add(tree.getExternalNode(i).getIdentifier().getName());
        }
        return labels;
    }

    private static List<List<String>> getBalancedSideForests(CoreTree core, int minSize, int maxSize) {
        final List<List<String>> forests = new ArrayList<List<String>>();
        List<String> smallForest = new ArrayList<String>();
        for (int i=0; i<core.getSideTreeCount(); i++) {
            final List<String> sideTree = core.getSideTreeLeafLabels(i);
            if (sideTree.size() >= minSize) {
                forests.add(sideTree);
            } else {
                smallForest.addAll(sideTree);
                if (smallForest.size() >= minSize) {
                    forests.add(smallForest);
                    smallForest = new ArrayList<String>();
                }
            }
        }
        if (!smallForest.isEmpty()) {
            forests.add(smallForest);
        }
        return forests;
    }

    //TODO verify proper complexity of match. For two trees all matchings should take O((kn)^1.5 log n + n^2).
    private static int match(Node v1, Node v2, CRMASTSet crmastSet) {
        final Node[] v1Neighbors = TreeCmpUtils.getNeighboringNodes(v1);
        final Node[] v2Neighbors = TreeCmpUtils.getNeighboringNodes(v2);
        final int size = Math.max(v1Neighbors.length, v2Neighbors.length);

        final int w[][] = new int[size][];
        for (int i=0; i<size; i++) {
            w[i] = new int[size];
        }

        for (int i=0; i<v1Neighbors.length; i++) {
            for (int j=0; j<v2Neighbors.length; j++) {
                w[i][j] = -crmastSet.getForEdgePair(v1, v1Neighbors[i], v2, v2Neighbors[j]);
            }
        }

        final int[] rowSol = new int[size];
        final int[] colSol = new int[size];
        final int[] u = new int[size];
        final int[] v = new int[size];
        return -LapSolver.lap(size, w, rowSol, colSol, u, v);
    }

    private static final class CoreTree {
        private List<List<String>> sideTrees = new ArrayList<List<String>>();
        private List<Integer> externalNodes = new ArrayList<Integer>();
        private List<Integer> internalNodes = new ArrayList<Integer>();

        public CoreTree(Tree tree, int maxSideTreeSize) {
            final boolean[] isInternalNodeInCore = getCoreNodes(tree, maxSideTreeSize);

            for (int i=0; i<isInternalNodeInCore.length; i++) {
                if (isInternalNodeInCore[i]) {
                    final Node v = tree.getInternalNode(i);
                    final Node[] neighbors = TreeCmpUtils.getNeighboringNodes(v);

                    int numberOfCoreNeighbors = 0;
                    for (final Node neighbor : neighbors) {
                        if (!neighbor.isLeaf() && isInternalNodeInCore[neighbor.getNumber()]) {
                            numberOfCoreNeighbors += 1;
                        } else {
                            sideTrees.add(subtreeLeafs(v, neighbor, tree.getInternalNodeCount()));
                        }
                    }
                    if (numberOfCoreNeighbors > 1) {
                        internalNodes.add(i);
                    } else {
                        externalNodes.add(i);
                    }
                }
            }
        }

        public List<String> getSideTreeLeafLabels(int i) {
            return sideTrees.get(i);
        }

        public int getSideTreeCount() {
            return sideTrees.size();
        }

        public int getNodeCount() {
            return externalNodes.size() + internalNodes.size();
        }

        public int getNode(int i) {
            final int externalNodeCount = getExternalNodeCount();
            return i < externalNodeCount
                    ? externalNodes.get(i)
                    : internalNodes.get(i - externalNodeCount);
        }

        public int getExternalNodeIdx(int i) {
            return externalNodes.get(i);
        }

        public int getExternalNodeCount() {
            return externalNodes.size();
        }

        private static boolean[] getCoreNodes(Tree tree, int maxSideTreeSize) {
            final int internalNodes = tree.getInternalNodeCount();

            final int externalNodes = tree.getExternalNodeCount();
            final boolean[] isCoreNode = new boolean[internalNodes];
            final boolean[] isCoreBelow = new boolean[internalNodes];
            final int[] subTreeSizes = new int[internalNodes];
            for (final Node node : TreeCmpUtils.getNodesInPostOrder(tree)) {
                final int subTreeSize = node.isLeaf()
                        ? 1
                        : subTreeSizes[node.getNumber()];
                //determine if the node is in core tree
                if (!node.isLeaf()
                        && subTreeSize >= maxSideTreeSize) {
                    //core edge going to parent node
                    if (externalNodes - subTreeSize >= maxSideTreeSize) {
                        isCoreNode[node.getNumber()] = true;
                        if (!node.isRoot()) {
                            isCoreNode[node.getParent().getNumber()] = true;
                        }
                        //or the core tree consists of a single node
                    } else if (!isCoreBelow[node.getNumber()]) {
                        isCoreNode[node.getNumber()] = true;
                    }
                }

                //update information about the parent node
                if (!node.isRoot()) {
                    subTreeSizes[node.getParent().getNumber()] += subTreeSize;
                    if (!node.isLeaf()) {
                        isCoreBelow[node.getParent().getNumber()] |= isCoreBelow[node.getNumber()] | isCoreNode[node.getNumber()];
                    }
                }
            }

            return isCoreNode;
        }

        private List<String> subtreeLeafs(Node cameFrom, Node start, int numberOfInternalNodes) {
            if (cameFrom.isLeaf()) {
                throw new IllegalArgumentException("Tree traversal has to start from an internal node.");
            }

            if (start.isLeaf()) {
                return Arrays.asList(start.getIdentifier().getName());
            }

            final boolean[] visited = new boolean[numberOfInternalNodes];
            final List<String> leafs = new ArrayList<String>();
            final Stack<Node> stack = new Stack<Node>();
            stack.add(start);
            visited[start.getNumber()] = visited[cameFrom.getNumber()] = true;

            while (!stack.isEmpty()) {
                final Node v = stack.pop();
                for (final Node neighbor : TreeCmpUtils.getNeighboringNodes(v)) {
                    if (neighbor.isLeaf()) {
                        leafs.add(neighbor.getIdentifier().getName());
                    } else if (!visited[neighbor.getNumber()]) {
                        stack.push(neighbor);
                        visited[neighbor.getNumber()] = true;
                    }
                }
            }
            return leafs;
        }
    }

    private static final class CRMASTSet {
        private final Map<Long, Integer> values;
        private final int t1Leafs;
        private final int t2Leafs;
        private final int t1Nodes;
        private final int t2Nodes;

        private final long u1Dim;
        private final long u2Dim;
        private final long v1Dim;

        public CRMASTSet(Tree t1, Tree t2) {
            t1Leafs = t1.getExternalNodeCount();
            t2Leafs = t2.getExternalNodeCount();
            t1Nodes = t1Leafs + t1.getInternalNodeCount();
            t2Nodes = t2Leafs + t2.getInternalNodeCount();
            values = new HashMap<Long, Integer>();
            u1Dim = t1Nodes * t1Nodes * t2Nodes;
            u2Dim = t1Nodes * t2Nodes;
            v1Dim = t2Nodes;
        }

        public int getForEdgePair(Node v1, Node w1, Node v2, Node w2) {
            return values.get(u1Dim * getT1Index(v1) + u2Dim * getT1Index(w1) + v1Dim * getT2Index(v2) + getT2Index(w2));
        }

        public void include(CRMAST crmast, Tree tree1, Tree tree2) {
            for (final Node v1 : TreeCmpUtils.getAllNodes(tree1)) {
                if (!v1.isRoot()) {
                    final Node u1 = v1.getParent();
                    final int u1Idx = getT1Index(u1);
                    final int v1Idx = getT1Index(v1);
                    for (final Node v2 : TreeCmpUtils.getAllNodes(tree2)) {
                        if (!v2.isRoot()) {
                            final Node u2 = v2.getParent();
                            final int u2Idx = getT2Index(u2);
                            final int v2Idx = getT2Index(v2);
                            values.put(u1Dim * u1Idx + u2Dim * v1Idx + v1Dim * u2Idx + v2Idx, crmast.getRMAST(v1, v2));
                        }
                    }
                }
            }
        }

        private int getT1Index(Node v) {
            return v.isLeaf()
                    ? v.getNumber()
                    : t1Leafs + v.getNumber();
        }

        private int getT2Index(Node v) {
            return v.isLeaf()
                    ? v.getNumber()
                    : t2Leafs + v.getNumber();
        }
    }

    /**
     * Computes LCA of two leafs in time O(1) after O(n) preprocessing.
     * Bender, A. and Farach-Colton M.; The LCA Problem Revisited
     */
    private static final class LCACalculator {
        private final Node[] order;
        private final int[] representatives;
        private final RangeMinQuery rmqSolver;

        public LCACalculator(Tree tree) {
            order = getEulerPathVisitOrder(tree);
            int[] depths = getDepths(order, tree);
            rmqSolver = new RangeMinQuery(depths);
            representatives = getRepresentatives(order, tree.getExternalNodeCount());
        }

        public Node forLeafs(Node u, Node v) {
            final int uIdx = representatives[getLeafIndex(u)];
            final int vIdx = representatives[getLeafIndex(v)];
            final int minIdx = uIdx < vIdx
                    ? rmqSolver.getMinElementIndex(uIdx, vIdx)
                    : rmqSolver.getMinElementIndex(vIdx, uIdx);
            return order[minIdx];
        }

        private static Node[] getEulerPathVisitOrder(Tree tree) {
            final int nodes = tree.getInternalNodeCount() + tree.getExternalNodeCount();
            final int edges = nodes-1;
            final Node[] order = new Node[2*edges+1];
            int orderIdx=0;

            final Stack<Node> s = new Stack<Node>();
            final Stack<Node> pathToRoot = new Stack<Node>();

            final Node root = tree.getRoot();
            pathToRoot.push(root);
            order[orderIdx++] = root;
            pushChildren(s, root);
            while (!s.empty()) {
                final Node node = s.pop();
                final Node parent = node.getParent();
                if (pathToRoot.peek() != parent) {
                    pathToRoot.pop();
                    //return towards the root until finding parent of the node
                    Node returnToNode;
                    do {
                        returnToNode = pathToRoot.pop();
                        order[orderIdx++] = returnToNode;
                    } while (returnToNode != parent);
                    pathToRoot.push(returnToNode);
                }
                //go down to the next node
                order[orderIdx++] = node;
                pathToRoot.push(node);
                pushChildren(s, node);
            }
            pathToRoot.pop();
            while (!pathToRoot.isEmpty()) {
                order[orderIdx++] = pathToRoot.pop();
            }
            return order;
        }

        private static void pushChildren(final Stack<Node> s, final Node root) {
            for (int i=root.getChildCount()-1; i>=0; --i) {
                s.push(root.getChild(i));
            }
        }

        private static int[] getDepths(Node[] nodes, Tree tree) {
            final int[] leafDepths = new int[tree.getExternalNodeCount()];
            final int[] internalDepths = new int[tree.getInternalNodeCount()];
            final Node root = tree.getRoot();
            final Stack<Node> s = new Stack<Node>();
            if (!root.isLeaf()) {
                s.push(root);
            }
            while (!s.empty()) {
                final Node node = s.pop();
                final int childDepth = internalDepths[node.getNumber()] + 1;
                for (int i=node.getChildCount()-1; i>=0; --i) {
                    final Node child = node.getChild(i);
                    if (child.isLeaf()) {
                        leafDepths[child.getNumber()] = childDepth;
                    } else {
                        internalDepths[child.getNumber()] = childDepth;
                        s.push(child);
                    }
                }
            }

            final int[] depths = new int[nodes.length];
            for (int i=nodes.length-1; i>=0; --i) {
                final Node node = nodes[i];
                if (node.isLeaf()) {
                    depths[i] = leafDepths[node.getNumber()];
                } else {
                    depths[i] = internalDepths[node.getNumber()];
                }
            }
            return depths;
        }

        private static int[] getRepresentatives(Node[] order, int numberOfLeafs) {
            final int[] representatives = new int[numberOfLeafs];
            for (int i=order.length-1; i>=0; --i) {
                final Node node = order[i];
                if (node.isLeaf()) {
                    representatives[getLeafIndex(node)] = i;
                }
            }
            return representatives;
        }

        private static int getLeafIndex(Node u) {
            return u.getNumber();
        }
    }

    private final static class RangeMinQuery {
        private final int[] minimumInBlock;
        private final int[] minimumIndexInBlock;
        private final int[] blockHashes;
        private final Map<Integer, int[][]> blockResultCache = new HashMap<Integer, int[][]>();
        private final int[][] m;
        private final int blockSize;
        private int[] t;

        public RangeMinQuery(int[] t) {
            //logN = ceil(lg(a.length))
            this.t = t;
            final int logN = Integer.SIZE-Integer.numberOfLeadingZeros(t.length-1);
            blockSize = (1+logN)/2;
            final int numberOfBlocks = (t.length + blockSize - 1) / blockSize;

            this.blockHashes = new int[numberOfBlocks];
            this.minimumInBlock = new int[numberOfBlocks];
            this.minimumIndexInBlock = new int[numberOfBlocks];
            for (int i=0, j=0; i<t.length; i+=blockSize, j++) {
                //init block minima
                minimumInBlock[j] = t[i];
                minimumIndexInBlock[j] = i;
                for (int k=Math.min(i+blockSize, t.length)-1; k>=i; --k) {
                    if (t[k] < minimumInBlock[j]) {
                        minimumInBlock[j] = t[k];
                        minimumIndexInBlock[j] = k;
                    }
                }
                //calculate block hash
                blockHashes[j] = hash(j);
            }

            //init sparse table m
            final int logSize = Integer.SIZE - Integer.numberOfLeadingZeros(numberOfBlocks-1);
            m = new int[logSize+1][];
            m[0] = new int[minimumInBlock.length];
            for (int i=0; i<minimumInBlock.length; ++i) {
                m[0][i] = i;
            }
            for (int i=1, step=1; i<=logSize; i++, step+=step) {
                m[i] = new int[minimumInBlock.length];
                for (int j=0; j<minimumInBlock.length; ++j) {
                    int leftMinIdx = m[i-1][j];
                    int rightMinIdx = j+step < minimumInBlock.length
                            ? m[i-1][j+step]
                            : leftMinIdx;
                    m[i][j] = minimumInBlock[rightMinIdx] < minimumInBlock[leftMinIdx]
                            ? rightMinIdx
                            : leftMinIdx;
                }
            }
        }

        private int hash(int blockNum) {
            int hash = 0;
            final int blockStart = blockNum * blockSize;
            final int hashLength = Math.min(blockSize, t.length - blockStart)-1;
            for (int i=0; i<hashLength; ++i) {
                hash |= (t[blockStart+i+1] > t[blockStart+i] ? 0 : 1) << i;
            }
            return hash;
        }

        private int[] expand(int blockHash) {
            final int[] block = new int[blockSize];
            for (int i=1; i<blockSize; i++) {
                final boolean isAscending = (blockHash & (1<<(i-1))) == 0;
                block[i] = block[i-1] + (isAscending ? 1 : -1);
            }
            return block;
        }

        public int getMinElementIndex(int l, int r) {
            final int lBlock = l/blockSize;
            final int rBlock = r/blockSize;

            if (rBlock - lBlock > 0) {
                final int minElementInL = minElementIndexInBlock(lBlock, l%blockSize, blockSize-1) + lBlock * blockSize;
                final int minElementInR = minElementIndexInBlock(rBlock, 0, r%blockSize) + rBlock * blockSize;
                final int minInSideBlocks = t[minElementInR] < t[minElementInL]
                        ? minElementInR
                        : minElementInL;
                if (rBlock - lBlock > 1) {
                    final int minInInternalBlocks = minElementIndex(lBlock, rBlock);
                    return t[minInSideBlocks] < t[minInInternalBlocks]
                            ? minInSideBlocks
                            : minInInternalBlocks;
                } else {
                    return minInSideBlocks;
                }
            } else {
                return minElementIndexInBlock(lBlock, l%blockSize, r%blockSize) + lBlock * blockSize;
            }
        }

        private int minElementIndex(int lBlock, int rBlock) {
            //k = floor(lg(r-l))
            final int k = Integer.SIZE - Integer.numberOfLeadingZeros(rBlock-lBlock) - 1;
            final int twoToK = 1<<k;
            int leftMinIdx = m[k][lBlock];
            int rightMinIdx = m[k][rBlock-twoToK+1];
            return minimumInBlock[leftMinIdx] < minimumInBlock[rightMinIdx]
                    ? minimumIndexInBlock[leftMinIdx]
                    : minimumIndexInBlock[rightMinIdx];
        }

        private int minElementIndexInBlock(int blockNum, int l, int r) {
            final int blockHash = blockHashes[blockNum];
            int[][] resultsForBlock = blockResultCache.get(blockHash);
            if (resultsForBlock == null) {
                resultsForBlock = calculateBlock(blockHash);
                blockResultCache.put(blockHash, resultsForBlock);
            }
            return resultsForBlock[r][l];
        }

        private int[][] calculateBlock(int blockHash) {
            final int[][] blockResults = new int[blockSize][];
            final int[] block = expand(blockHash);
            for (int r=0; r<blockSize; ++r) {
                blockResults[r] = new int[r+1];
                blockResults[r][r] = r;
                int minValue = block[r];
                for (int l=r-1; l>=0; --l) {
                    final int currentValue = block[l];
                    if (currentValue < minValue) {
                        minValue = currentValue;
                        blockResults[r][l] = l;
                    } else {
                        blockResults[r][l] = blockResults[r][l+1];
                    }
                }
            }
            return blockResults;
        }
    }

    private static final class SimpleUnrootedTreePreservingNodeNumbers extends SimpleTree {
        private static final long serialVersionUID = -4057750013155408632L;
        private boolean initialized = false;

        public SimpleUnrootedTreePreservingNodeNumbers(Node root) {
            super(contractDegree2Root(root));
            initialized = true;
        }

        @Override
        public void createNodeList() {
            if (!initialized) {
                super.createNodeList();
            }
        }

        private static Node contractDegree2Root(Node root) {
            if (!root.isRoot()) {
                throw new IllegalArgumentException("Node has to be root.");
            }
            final int rootChildCount = root.getChildCount();
            if (rootChildCount >= 3 || rootChildCount == 0) {
                return root;
            }
            Node biggestChild = root.getChild(0);
            for (int i=1; i<rootChildCount; i++) {
                final Node child = root.getChild(i);
                if (biggestChild.getChildCount() < child.getChildCount()) {
                    biggestChild = child;
                }
            }
            for (int i=0; i<rootChildCount; i++) {
                final Node child = root.getChild(i);
                if (child != biggestChild) {
                    biggestChild.addChild(child);
                }
            }
            biggestChild.setParent(null);
            return biggestChild;
        }
    }

    private static final class TreeRestricter {
        private final LCACalculator lca;
        private final int[] internalNodeDepths;
        private final Node[] leafOrder;

        public TreeRestricter(Tree tree) {
            lca = new LCACalculator(tree);

            internalNodeDepths = new int[tree.getInternalNodeCount()];
            leafOrder = new Node[tree.getExternalNodeCount()];

            int leafs = 0;
            for (final Node node : TreeCmpUtils.getNodesInPreOrder(tree)) {
                if (node.isLeaf()) {
                    leafOrder[leafs++] = node;
                } else if (!node.isRoot()) {	//is internal
                    internalNodeDepths[node.getNumber()]
                            = internalNodeDepths[node.getParent().getNumber()] + 1;
                }
            }
        }

        public SimpleTree getRestrictedTo(Set<String> leafNames) {
            final Node[] leafs = getOrderedLeafs(leafNames);
            final int[] lcaDepths = new int[leafs.length - 1];
            for (int i=1; i<leafs.length; i++) {
                final Node u = lca.forLeafs(leafs[i-1], leafs[i]);
                lcaDepths[i-1] = internalNodeDepths[u.getNumber()];
            }

            final Node root = construct(leafs, lcaDepths, 0, leafs.length, 0);
            return new SimpleUnrootedTreePreservingNodeNumbers(root);
        }

        private Node construct(Node[] leafs, int[] lcaDepths, int l, int r, int minDepth) {
            final Node u = NodeFactory.createNode(leafs[l].getIdentifier());
            Node parent = null;
            while (l<r-1 && lcaDepths[l] >= minDepth) {
                int i=l+1;
                while (i < r-1 && lcaDepths[i] > lcaDepths[l]) {
                    i++;
                }
                if (parent == null) {
                    parent = NodeFactory.createNode();
                    parent.addChild(u);
                }
                final Node subTree = construct(leafs, lcaDepths, l+1, i+1, lcaDepths[l]);
                parent.addChild(subTree);
                l = i;
            }
            return parent != null
                    ? parent
                    : u;
        }

        private Node[] getOrderedLeafs(Set<String> leafNames) {
            final Node[] leafs = new Node[leafNames.size()];
            int i=0;
            for (final Node node : leafOrder) {
                if (leafNames.contains(node.getIdentifier().getName())) {
                    leafs[i++] = node;
                }
            }
            return leafs;
        }
    }
}