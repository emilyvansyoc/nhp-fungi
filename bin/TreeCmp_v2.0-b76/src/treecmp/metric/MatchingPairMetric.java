/** This file is part of TreeCmp, a tool for comparing phylogenetic trees
    using the Matching Split distance and other metrics.
    Copyright (C) 2014,  Damian Bogdanowicz

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */
package treecmp.metric;

import pal.misc.IdGroup;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.LapSolver;
import treecmp.common.TreeCmpUtils;

public class MatchingPairMetric extends BaseMetric implements Metric {

    protected int[] rowsol;
    protected int[] colsol;
    protected int[][] assigncost;

    public MatchingPairMetric() {
        super();
        this.rooted = true;
    }

    @Override
    public double getDistance(Tree t1, Tree t2) {
        
        if (t1.getExternalNodeCount() <= 2){
            return 0.0;
        }
        
        IdGroup id1 = TreeUtils.getLeafIdGroup(t1);
        int[][] lcaMatrix1 = TreeCmpUtils.calcLcaMatrix(t1, null);
        int[][] lcaMatrix2 = TreeCmpUtils.calcLcaMatrix(t2, id1);

        int intT1Num = t1.getInternalNodeCount();
        int intT2Num = t2.getInternalNodeCount();

        Node[] postOrderT1 = TreeCmpUtils.getNodesInPostOrder(t1);
        Node[] postOrderT2 = TreeCmpUtils.getNodesInPostOrder(t2);

        short[] cSize1 = new short[intT1Num];
        short[] cSize2 = new short[intT2Num];

        TreeCmpUtils.calcCladeSizes(t1, postOrderT1, cSize1);
        TreeCmpUtils.calcCladeSizes(t2, postOrderT2, cSize2);

        int N = t1.getExternalNodeCount();

        int size = Math.max(intT1Num, intT2Num);
        if (size <= 0) {
            return 0;
        }

        assigncost = new int[size][size];
        rowsol = new int[size];
        colsol = new int[size];
        int[] u = new int[size];
        int[] v = new int[size];
        
        //iterate by all possible pairs of leaves
        //and fill assigncont with the value of intersection size
        for (int i = 0; i < N; i++){
            for (int j = i+1; j < N; j++){
               int int1 = lcaMatrix1[i][j];
               int int2 = lcaMatrix2[i][j];
               assigncost[int1][int2]++;
            }
        }
        //count LCA pairs for t1
        int[] t1IntPairCount = new int[intT1Num];
        for (int i = 0; i < intT1Num; i++){
            //Node n = t1.getInternalNode(alias1[i]);
            Node n = t1.getInternalNode(i);
            t1IntPairCount[i] =  coutChildrenPairs(n, cSize1);
        }
        //count LCA pairs for t2
        int[] t2IntPairCount = new int[intT2Num];
        for (int i = 0; i < intT2Num; i++){
            //Node n = t2.getInternalNode(alias2[i]);
            Node n = t2.getInternalNode(i);
            t2IntPairCount[i] =  coutChildrenPairs(n, cSize2);
        }

        //calc xor valuses of pairs sets and store it in assigncost matrix
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                if (i < intT1Num && j < intT2Num){
                    assigncost[i][j] = t1IntPairCount[i]+t2IntPairCount[j] - (assigncost[i][j] << 1);
                } else if (i >= intT1Num && j < intT2Num){
                    assigncost[i][j] = t2IntPairCount[j];
                } else if (i < intT1Num && j >= intT2Num){
                    assigncost[i][j] = t1IntPairCount[i];
                }else {
                    //normally should not happen
                    assigncost[i][j] = 0;
                }    
            }
        }   
        int metric = LapSolver.lap(size, assigncost, rowsol, colsol, u, v);
        return (0.5 * (double) metric);
    }
    
    
    int coutChildrenPairs(Node n, short[] clustSizeTab) {
        int chCount = n.getChildCount();
        int[] cSize = new int[chCount];

        for (int i = 0; i < chCount; i++) {
            Node chNode = n.getChild(i);
            if (chNode.isLeaf()) {
                cSize[i] = 1;
            } else {
                cSize[i] = clustSizeTab[chNode.getNumber()];
            }
        }
        int pairCount = 0;
        for (int i = 0; i < cSize.length; i++) {
            for (int j = i + 1; j < cSize.length; j++) {
                pairCount += (cSize[i] * cSize[j]);
            }
        }
        return pairCount;
    }
}

