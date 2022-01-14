/*
 * Copyright (c) 2015 Google, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.+ *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package com.google.earthengine.examples.experimental;

import com.google.common.annotations.UsedReflectively;
import com.google.common.annotations.VisibleForTesting;
import com.google.earthengine.api.array.EEArray;
import com.google.earthengine.api.base.AlgorithmBase;
import com.google.earthengine.api.base.BB;
import com.google.earthengine.api.base.UserInput;
import com.google.earthengine.api.image.Band;
import com.google.earthengine.api.image.BandFunctor;
import com.google.earthengine.api.image.Image;
import com.google.earthengine.api.image.ImageConstructor;
import com.google.earthengine.api.image.OutputBand;
import com.google.earthengine.api.image.PixelFunctionCollection;
import com.google.earthengine.api.image.TypedImageCollection;
import com.google.earthengine.api.pixel.PixelFunction;
import com.google.earthengine.api.pixel.PixelIterator;
import com.google.earthengine.api.pixel.PixelMask;
import com.google.earthengine.api.pixel.PixelStore;
import com.google.earthengine.api.pixel.PixelType;
import com.google.earthengine.api.proj.Proj;

/**
 * Wrapper for the Verdet Solver.
 */
public class VerdetWrapper extends ImageConstructor<VerdetWrapper.Args> {
  public static final VerdetWrapper INSTANCE = new VerdetWrapper();

  @VisibleForTesting
  static final String BAND_NAME = "score";

  @VisibleForTesting
  static final PixelType BAND_TYPE = PixelType.DOUBLE.withDimensions(1);

  /**
   *
   */
  @AlgorithmBase.FirstArgs
  public static class Args extends Verdet.Args {
    @Doc(help = "Collection from which to extract VeRDET scores. This collection is "
        + "expected to contain 1 image for each year, sorted temporally.")
    @Required
    TypedImageCollection timeSeries;
  }

  @Override
  public String doc() {
    return "Vegetation Regeneration and Disturbance Estimates through Time, forest change"
        + " detection algorithm. This algorithm generates a yearly clear-sky composite from"
        + " satellite imagery, calculates a spectral vegetation index for each pixel in"
        + " that composite, spatially segments the vegetation index image into patches,"
        + " temporally divides the time series into differently sloped segments, and then"
        + " labels those segments as disturbed, stable, or regenerating. Segmentation at both"
        + " the spatial and temporal steps are performed using total variation regularization.\n"
        + "The output consists of a 1D array per pixel containing the slope of fitted trend lines."
        + " Negative values indicate disturbance and positive values regeneration.\n"
        + "See: Hughes, M.J., Kaylor, S.D. and Hayes, D.J., 2017. "
        + " Patch-based forest change detection from Landsat time series. Forests, 8(5), p.166.";
  }

  @Override
  protected Image apply(Image output, Args args) {
    new VerdetFunctor(output.newBand(BAND_NAME), args);
    return output;
  }

  /**
   * Band functor wrapping the Verdet solver.
   */
  public static final class VerdetFunctor extends BandFunctor {
    private final Args args;

    public VerdetFunctor(Band output, Args args) {
      this.args = args;
      if (args.timeSeries.bandNames().size() > 1) {
        throw UserInput.error("VeRDET only operates on collections of images "
            + "of one band.");
      }
      addInputCollection(args.timeSeries, args.timeSeries.bandNames().get(0));
      addOutputBand(output, BAND_TYPE, Proj.WGS84, BB.UNBOUNDED);
    }

    @UsedReflectively
    private void computeOutput(
        OutputBand output,
        PixelIterator iter,
        final PixelFunctionCollection[] inputs) {
      // Create a solver
      final Verdet solver = new Verdet(args);

      int size = inputs[0].allImages().length;
      double[] a = new double[size];
      int[] lengths = { size };
      PixelStore ps = output.getPixelStore();

      for (; iter.isValid(); iter.next()) {
        int count = 0;
        // Collect the inputs.
        for (int i = 0; i < size; i++) {
          PixelFunction pf = inputs[0].getPixelFunction(i);
          if (pf.mask().getMask(iter.x(), iter.y()) == PixelMask.INVALID) {
            break;
          }
          a[i] = inputs[0].getPixelFunction(i).getFloatOrMissingDataValue(iter.x(), iter.y());
          count++;
        }
        // In case there were any masked pixels, abort.
        if (count == size) {
          EEArray result = EEArray.wrap(lengths, solver.getResult(a));
          iter.setArray(ps, result);
        }
      }
    }
  }
}