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

import treecmp.io.ResultWriter;
import treecmp.config.IOSettings;
import treecmp.metric.Metric;


public class SummaryStatCalculator {


    private Metric met;
    private int count;
    private double sum;
    private double sq_sum;
    private double min;
    private double max;

    public SummaryStatCalculator() {

        this.count=0;
        this.sum=0.0;
        this.sq_sum=0.0;

        this.min=Double.MAX_VALUE;
        this.max=-Double.MAX_VALUE;
    }


    /**
     *
     * @param _met
     */
    public SummaryStatCalculator(Metric _met) {

        //call non-parameter constructor
        this();
        this.met=_met;
    }

    public void addMetric(Metric _met)
    {
        this.met=_met;
        this.clear();
    }


   public void clear()
  {
        this.count=0;
        this.sum=0.0;
        this.sq_sum=0.0;

        this.min=Double.MAX_VALUE;
        this.max=-Double.MAX_VALUE;

    }


  public double getMax()
  {
   return this.max;
  }

  public double getMin()
  {
    return this.min;
  }

  public double getAvg()
  {
      double avg=Double.POSITIVE_INFINITY;
      if (count>0)
          avg=this.sum/(double)count;

      return avg;
  }

  public double getVariance()
  {
      double var=Double.POSITIVE_INFINITY;
      double avg;
      if (count>0)
      {
          avg=this.getAvg();
          var=this.sq_sum/(double)count-avg*avg;
      }
      return var;
  }

  public double getStd()
  {
        double std=Double.POSITIVE_INFINITY;
        double var;

        if(count>0)
        {
            var=this.getVariance();
            std=Math.sqrt(var);
        }
        return std;
  }


  public int getCount()
  {
      return this.count;
  }

  public void insertValue(double dist)
  {
        sum+=dist;
        count++;
        sq_sum+=dist*dist;

        if(dist<min) min=dist;
        if(dist>max) max=dist;

  }

    public String getName() {
        return this.met.getName();
    }

     public String getCommandLineName() {
        return this.met.getCommandLineName();
    }


    public static void printSummary(ResultWriter out, SummaryStatCalculator[] sStatCalc) {

        if (IOSettings.getIOSettings().isGenSummary()) {
            int size = sStatCalc.length;
            String separator = IOSettings.getIOSettings().getSSep();
            Object[] line = new Object[6];

            line[0] = "---------";
            out.setRow(line);
            out.write();
            line[0] = "Summary:";

            out.setRow(line);
            out.write();

            //name-avg-std-min-max-count
            line[0] = "Name";
            line[1] = "Avg";
            line[2] = "Std";
            line[3] = "Min";
            line[4] = "Max";
            line[5] = "Count";

            out.setRow(line);
            out.write();

            for (int i = 0; i < size; i++) {
                line[0] = sStatCalc[i].getName();
                line[1] = sStatCalc[i].getAvg();
                line[2] = sStatCalc[i].getStd();
                line[3] = sStatCalc[i].getMin();
                line[4] = sStatCalc[i].getMax();
                line[5] = sStatCalc[i].getCount();
                out.setRow(line);
                out.write();
            }
        }
    }
}
