/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.spr;

import treecmp.metric.MatchingClusterMetricO3;
import treecmp.metric.Metric;
import treecmp.metric.RFClusterMetric;

/**
 *
 * @author Damian
 */
public class SprHeuristicMCMetric extends SprHeuristicBaseMetric{

 @Override
protected Metric getMetric(){
    return new MatchingClusterMetricO3();
 }
}
