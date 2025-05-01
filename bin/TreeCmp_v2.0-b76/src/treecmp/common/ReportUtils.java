/** This file is part of TreeCmp, a tool for comparing phylogenetic trees
    using the Matching Split distance and other metrics.
    Copyright (C) 2011,  Damian Bogdanowicz

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

package treecmp.common;

import java.util.Locale;
import treecmp.config.IOSettings;


public class ReportUtils {

    private final static int ROW_PRECISION = 4;
    public final static String ROW_DATA_FORMAT="%1$."+ROW_PRECISION+"f";
    private static IOSettings ioSet = IOSettings.getIOSettings();
    private static boolean pruneTrees = ioSet.isPruneTrees();
    private static boolean randomComparison = ioSet.isRandomComparison();

    public final static String NUM_COLUMN = "No";
    public final static String T1_COLUMN = "Tree1";
    public final static String T2_COLUMN = "Tree2";
    public final static String T_COLUMN = "Tree";
    public final static String RT_COLUMN = "RefTree";
    public final static String T1_TAXA= "Tree1_taxa";
    public final static String T2_TAXA= "Tree2_taxa";
    public final static String T_TAXA= "Tree_taxa";
    public final static String RT_TAXA= "RefTree_taxa";
    public final static String COMMON_TAXA= "Common_taxa";
    public final static String YULE_FRAC= "_toYuleAvg";
    public final static String UNIF_FRAC= "_toUnifAvg";
    public final static String NA_FRAC= "N/A";
    //it t2 ==-1 do not print t2

    private static int rowCount;

    public static void setRowCount(int metricsCount) {
        int rc = metricsCount;
        if (randomComparison) rc *= 3;
        if (pruneTrees) rc += 3;
        rc += 3;
        rowCount = rc;
    }

     public static Object[] getResultRow(int rowNum, int t1, int t2, StatCalculator[] stats){

         Object[] row = new Object[rowCount];
         double dist, distToYuleAvg, distToUnifAvg;
         String distStr;
         int cellNum = 0;
         row[cellNum++] = rowNum;
         row[cellNum++] = t1;
         row[cellNum++] = t2;

        //to do if prune enabled
        if (pruneTrees && stats.length > 0){
            row[cellNum++] = stats[0].getT1TaxaNum();
            row[cellNum++] = stats[0].getT2TaxaNum();
            row[cellNum++] = stats[0].getCommonTaxaNum();
        }

        for (int i=0; i< stats.length; i++){

             dist = stats[i].getLastDist();
             distStr = String.format(Locale.US,ROW_DATA_FORMAT, dist);
             row[cellNum++] = distStr;
             if (randomComparison){
                 distToYuleAvg = stats[i].getLastDistToYuleAvg();
                 if (distToYuleAvg != Double.NEGATIVE_INFINITY)
                    distStr = String.format(Locale.US,ROW_DATA_FORMAT, distToYuleAvg);
                 else
                     distStr = NA_FRAC;
                 row[cellNum++] = distStr;

                 distToUnifAvg = stats[i].getLastDistToUnifAvg();
                 if (distToUnifAvg != Double.NEGATIVE_INFINITY)
                    distStr = String.format(Locale.US,ROW_DATA_FORMAT, distToUnifAvg);
                 else
                     distStr = NA_FRAC;
                 row[cellNum++] = distStr;
             }
        }

        return row;
    }

    public static Object[] getHeaderRow(StatCalculator[] stats){

        return getHeaderRow(stats, false);
    }

/*     public static String getHeaderRow(StatCalculator[] stats){
         
         return getHeaderRow(stats, false);
     }*/

     /*public static String getHeaderRow(StatCalculator[] stats, boolean ifReefTreeMode){
        StringBuilder sb = new StringBuilder();

        String metricName;

         if (!ifReefTreeMode){
            sb.append(NUM_COLUMN);
            sb.append(sep);
            sb.append(T1_COLUMN);
        }else{
            sb.append(NUM_COLUMN);
            sb.append(sep);
            sb.append(RT_COLUMN);
        }
       
        //used in reference tree mode
        if (!ifReefTreeMode){
            sb.append(sep);
            sb.append(T2_COLUMN);
        } else {
            sb.append(sep);
            sb.append(T_COLUMN);
        }
        //to do if prune enabled
        if (pruneTrees && stats.length > 0) {
            sb.append(sep);

            if (!ifReefTreeMode) {
                sb.append(T1_TAXA);
            } else {
                sb.append(RT_TAXA);
            }

            sb.append(sep);

            if (!ifReefTreeMode) {
                sb.append(T2_TAXA);
            } else {
                sb.append(T_TAXA);
            }

            sb.append(sep);
            sb.append(COMMON_TAXA);
        }

        for (int i = 0; i < stats.length; i++) {
            sb.append(sep);
            metricName = stats[i].getName();
            sb.append(metricName);
            if (randomComparison){
                sb.append(sep);             
                sb.append(metricName);
                sb.append(YULE_FRAC);
                sb.append(sep);
                sb.append(metricName);
                sb.append(UNIF_FRAC);
            }
        }
        //sb.append("\n");
        return sb.toString();
    }*/

    public static Object[] getHeaderRow(StatCalculator[] stats, boolean ifReefTreeMode){

        Object[] row = new Object[rowCount];
        String metricName;
        int cellNum = 0;
        if (!ifReefTreeMode){
            row[cellNum++] = NUM_COLUMN;
            row[cellNum++] = T1_COLUMN;
        }else{
            row[cellNum++] = NUM_COLUMN;
            //row[cellNum++] = sep;
            row[cellNum++] = RT_COLUMN;
        }

        //used in reference tree mode
        if (!ifReefTreeMode){
            //row[cellNum++] = sep;
            row[cellNum++] = T2_COLUMN;
        } else {
            //row[cellNum++] = sep;
            row[cellNum++] = T_COLUMN;
        }
        //to do if prune enabled
        if (pruneTrees && stats.length > 0) {
            //row[cellNum++] = sep;

            if (!ifReefTreeMode) {
                row[cellNum++] = T1_TAXA;
            } else {
                row[cellNum++] = RT_TAXA;
            }

            if (!ifReefTreeMode) {
                row[cellNum++] = T2_TAXA;
            } else {
                row[cellNum++] = T_TAXA;
            }

            row[cellNum++] = COMMON_TAXA;
        }

        for (int i = 0; i < stats.length; i++) {
            metricName = stats[i].getName();
            row[cellNum++] = metricName;
            if (randomComparison){
                row[cellNum++] = metricName + YULE_FRAC;
                row[cellNum++] = metricName + UNIF_FRAC;
            }
        }
        //row[cellNum++] = "\n";
        return row;
    }

    public static void update() {
        IOSettings ioSet = IOSettings.getIOSettings();
        String sep = ioSet.getSSep();
        pruneTrees = ioSet.isPruneTrees();
        randomComparison = ioSet.isRandomComparison();
    }
}
