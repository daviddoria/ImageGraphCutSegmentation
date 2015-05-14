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

// Boost
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

/** Perform graph cut based segmentation on an image. Image pixels can contain any
  * number of components (i.e. grayscale, RGB, RGBA, RGBD, etc.).
  * This is an implementation of the technique described here:
  * http://www.eecs.berkeley.edu/~efros/courses/AP06/Papers/boykov-iccv-01.pdf
  */
template <typename TImage, typename TPixelDifferenceFunctor = RGBPixelDifference<typename TImage::PixelType> >
class ImageGraphCut
{
public:
  // Typedefs

  /** This is a special type to keep track of the graph node labels. */
  typedef itk::Image<unsigned int, 2> NodeImageType;

  /** The type of the histograms. */
  typedef itk::Statistics::Histogram< float,
          itk::Statistics::DenseFrequencyContainer2 > HistogramType;

  /** The type of a list of pixels/indexes. */
  typedef std::vector<itk::Index<2> > IndexContainer;

  typedef typename TImage::PixelType PixelType;
  typedef itk::Statistics::ListSample<PixelType> SampleType;
  typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;

  /** If nothing else is provided, this is the default background likelihood function. */
  float InternalForegroundLikelihood(const PixelType& pixel);

  /** If nothing else is provided, this is the default background likelihood function. */
  float InternalBackgroundLikelihood(const PixelType& pixel);

  ImageGraphCut()
  {
      this->ForegroundLikelihood =
              boost::bind(
                  &ImageGraphCut::
                  InternalForegroundLikelihood, this, _1);

      this->BackgroundLikelihood =
              boost::bind(
                  &ImageGraphCut::
                  InternalBackgroundLikelihood, this, _1);
  }

  ImageGraphCut(TPixelDifferenceFunctor pixelDifferenceFunctor) :
    PixelDifferenceFunctor(pixelDifferenceFunctor){}

  TPixelDifferenceFunctor PixelDifferenceFunctor;

  /** Provide the image to segment. */
  void SetImage(TImage* const image);

  /** Several initializations are done here. */
  void Initialize();

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

  void SetForegroundLikelihoodFunction(boost::function<float (const PixelType& pixel)> f)
  {
    this->ForegroundLikelihood = f;
    this->CustomLikelihood = true;
  }

  void SetBackgroundLikelihoodFunction(boost::function<float (const PixelType& pixel)> f)
  {
    this->BackgroundLikelihood = f;
    this->CustomLikelihood = true;
  }

protected:

  /** A graph object for Kolmogorov*/
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
      boost::no_property,
      boost::property<boost::edge_index_t, std::size_t> > GraphType;

  typedef boost::graph_traits<GraphType>::vertex_descriptor VertexDescriptor;
  typedef boost::graph_traits<GraphType>::edge_descriptor EdgeDescriptor;
  typedef boost::graph_traits<GraphType>::vertices_size_type VertexIndex;
  typedef boost::graph_traits<GraphType>::edges_size_type EdgeIndex;

  /** Store the list of edges and their corresponding reverse edges. */
  std::vector<EdgeDescriptor> ReverseEdges;

  /** Create an edge on the graph. */
  unsigned int AddBidirectionalEdge(unsigned int numberOfEdges, const unsigned int source,
                            const unsigned int target,
                            const float weight);

  /** The main graph object. */
  GraphType Graph;

  /** Maintain a list of all of the edge weights. */
  std::vector<float> EdgeWeights;

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

  /** Create the histograms from the users selections */
  void CreateSamples();

  /** Estimate the "camera noise" */
  double ComputeNoise();

  /** Create a Kolmogorov graph structure from the image and selections */
  void CreateGraph();

  /** Create the edges between pixels and neighboring pixels (the grid). */
  void CreateNEdges();

  /** Create the edges between pixels and the terminals (source and sink). */
  void CreateTEdges();

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

  /** The node id of the foreground terminal. */
  unsigned int SourceNodeId;

  /** The node id of the background terminal. */
  unsigned int SinkNodeId;

  /** The function pointer that gets called to determine the likelihood that the pixel belongs to the foreground. */
  boost::function<float (const PixelType& pixel)> ForegroundLikelihood;

  /** The function pointer that gets called to determine the likelihood that the pixel belongs to the background. */
  boost::function<float (const PixelType& pixel)> BackgroundLikelihood;

  /** If a custom likelihood function is set, we don't need to compute histograms internally. */
  bool CustomLikelihood = false;
};

#include "ImageGraphCut.hpp"

#endif
