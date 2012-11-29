/*
Copyright (C) 2012 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ImageGraphCut_H
#define ImageGraphCut_H

// Custom
#include "PixelDifference.h"

// Submodules
#include "Mask/ForegroundBackgroundSegmentMask.h"

// ITK
#include "itkImage.h"
#include "itkSampleToHistogramFilter.h"
#include "itkHistogram.h"
#include "itkListSample.h"

// STL
#include <vector>

// Kolmogorov's code
#include "Kolmogorov/graph.h"
typedef Graph GraphType;

template <typename TImage, typename TPixelDifferenceFunctor = RGBPixelDifference<typename TImage::PixelType> >
class ImageGraphCut
{
public:

  ImageGraphCut(){}

  ImageGraphCut(TPixelDifferenceFunctor pixelDifferenceFunctor) :
    PixelDifferenceFunctor(pixelDifferenceFunctor){}

  TPixelDifferenceFunctor PixelDifferenceFunctor;

  /** This is a special type to keep track of the graph node labels. */
  typedef itk::Image<void*, 2> NodeImageType;

  /** The type of the histograms. */
  typedef itk::Statistics::Histogram< float,
          itk::Statistics::DenseFrequencyContainer2 > HistogramType;

  /** The type of a list of pixels/indexes. */
  typedef std::vector<itk::Index<2> > IndexContainer;

  /** Several initializations are done here. */
  void SetImage(TImage* const image);

  /** Get the image that we are segmenting. */
  TImage* GetImage();

  /** Create and cut the graph (The main driver function). */
  void PerformSegmentation();

  /** Return a list of the selected (via scribbling) pixels. */
  IndexContainer GetSources();
  IndexContainer GetSinks();

  /** Set the selected pixels. */
  void SetSources(const IndexContainer& sources);
  void SetSinks(const IndexContainer& sinks);

  /** Get the output of the segmentation. */
  ForegroundBackgroundSegmentMask* GetSegmentMask();

  /** Set the weight between the regional and boundary terms. */
  void SetLambda(const float);

  /** Set the number of bins per dimension of the foreground and background histograms. */
  void SetNumberOfHistogramBins(const int);

protected:

  /** A Kolmogorov graph object */
  GraphType* Graph;

  /** The output segmentation */
  ForegroundBackgroundSegmentMask::Pointer ResultingSegments;

  /** User specified foreground points */
  IndexContainer Sources;

  /** User specified background points */
  IndexContainer Sinks;

  /** The weighting between unary and binary terms */
  float Lambda = 0.01f;

  /** The number of bins per dimension of the foreground and background histograms */
  int NumberOfHistogramBins = 10;

  /** An image which keeps tracks of the mapping between pixel index and graph node id */
  NodeImageType::Pointer NodeImage;

  // Typedefs
  typedef typename TImage::PixelType PixelType;
  typedef itk::Statistics::ListSample<PixelType> SampleType;
  typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;

  /** Create the histograms from the users selections */
  void CreateSamples();

  /** Estimate the "camera noise" */
  double ComputeNoise();

  /** Create a Kolmogorov graph structure from the image and selections */
  void CreateGraph();

  /** Perform the s-t min cut */
  void CutGraph();

  /** The ITK data structure for storing the values that we will compute the histogram of. */
  typename SampleType::Pointer ForegroundSample;
  typename SampleType::Pointer BackgroundSample;

  /** The histograms. */
  const HistogramType* ForegroundHistogram = nullptr;
  const HistogramType* BackgroundHistogram = nullptr;

  /** ITK filters to create histograms. */
  typename SampleToHistogramFilterType::Pointer ForegroundHistogramFilter;
  typename SampleToHistogramFilterType::Pointer BackgroundHistogramFilter;

  /** The image to be segmented */
  typename TImage::Pointer Image;

};

#include "ImageGraphCut.hpp"

#endif
