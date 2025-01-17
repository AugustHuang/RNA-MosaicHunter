/*
 * The MIT License
 *
 * Copyright (c) 2016 Center for Bioinformatics, Peking University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package cn.edu.pku.cbi.mosaichunter.math;

import java.util.Arrays;

import org.apache.commons.math3.distribution.NormalDistribution;

public class WilcoxonRankSumTest {
    
    private static double EPS = 1e-9;
    
    private WilcoxonRankSumTest() {
        // private constructor 
    }
    
    public static double twoSided(double[] x, double[] y) {
        int nx = x.length;
        int ny = y.length;
        double[] v = new double[nx + ny];
        Arrays.sort(x);
        System.arraycopy(x, 0, v, 0, nx);
        System.arraycopy(y, 0, v, nx, ny);
        Arrays.sort(v);
        
        double[] u = new double[v.length];
        double[] rank = new double[v.length];
        int[] cnt = new int[v.length];
        int n = 0;
        for (int i = 0; i < v.length; ++i) {
            if (i == 0 || v[i] - v[i - 1] > EPS) {
                u[n] = v[i];
                n++;
            }
            cnt[n - 1]++;
            rank[n - 1] += i + 1;
        }
        for (int i = 0; i < n; ++i) {
            rank[i] /= cnt[i];
        }
        
        double stats = 0;
        int j = 0;
        for (int i = 0; i < nx; ++i) {
            while(x[i] > u[j] + EPS) {
                j++;
            }
            stats += rank[j];
        }
        stats -= nx * (nx + 1) / 2;
        double z = stats - x.length * y.length / 2;

        double tmp = 0;
        for (int i = 0; i < n; ++i) {
            tmp += (double) cnt[i] * cnt[i] * cnt[i] - cnt[i];
        }
        double sigma = Math.sqrt(
                ((double) nx * ny / 12) * (nx + ny + 1 - tmp / (nx + ny) / (nx + ny - 1)));
        if (z > EPS) {
            z -= 0.5;
        } else if (z < -EPS) {
            z += 0.5;
        }
        z /= sigma;
        
        double p = new NormalDistribution().cumulativeProbability(z);
        double pValue = 2 * Math.min(p, 1 - p); 
        if (Double.isNaN(pValue)) {
            return 1;
        }
        return pValue;
    }
}
