/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package treecmp.metric.weighted;

import pal.tree.Tree;
import pal.tree.TreeTool;
import treecmp.metric.BaseMetric;
import treecmp.metric.Metric;

/**
 *
 * @author Damian
 */
public class GeoUnrootedMetric extends BaseMetric implements Metric {

    private GeoMetricWrapper geoMetricWrapper = new GeoMetricWrapper();

    @Override
    public boolean isRooted() {
        return false;
    }

    @Override
    public double getDistance(Tree t1, Tree t2) {
        int extT1Num = t1.getExternalNodeCount();
        Tree t1u, t2u;
        if (extT1Num <= 1) {
            return 0.0;
        }
        else if (extT1Num > 2) {
            t1u = TreeTool.getUnrooted(t1);
            t2u = TreeTool.getUnrooted(t2);
        }
        else {
            t1u = t1;
            t2u = t2;
        }
        double dist = geoMetricWrapper.getDistance(t1u, t2u, false, null);
        return dist;
    }
}
