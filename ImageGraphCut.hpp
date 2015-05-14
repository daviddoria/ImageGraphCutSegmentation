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

#ifndef ImageGraphCut_HPP
#define ImageGraphCut_HPP

#include "ImageGraphCut.h"

// Submodules
#include "Mask/ITKHelpers/Helpers/Helpers.h"
#include "Mask/ITKHelpers/ITKHelpers.h"

// ITK
#include "itkImageRegionIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkMaskImageFilter.h"

// STL
#include <cmath>
#include <algorithm>

// Boost
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetImage(TImage* const image)
{
  this->Image = TImage::New();
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::Initialize()
{
    // Setup the output (mask) image
    this->ResultingSegments = ForegroundBackgroundSegmentMask::New();
    this->ResultingSegments->SetRegions(this->Image->GetLargestPossibleRegion());
    this->ResultingSegments->Allocate();

    // Setup the image to store the node ids
    this->NodeImage = NodeImageType::New();
    this->NodeImage->SetRegions(this->Image->GetLargestPossibleRegion());
    this->NodeImage->Allocate();

    // Initializations
    this->ForegroundSample = SampleType::New();
    this->BackgroundSample = SampleType::New();

    this->ForegroundHistogramFilter = SampleToHistogramFilterType::New();
    this->BackgroundHistogramFilter = SampleToHistogramFilterType::New();

    // Blank the NodeImage
    ITKHelpers::SetImageToConstant(this->NodeImage.GetPointer(), 0);

    // Blank the output image
    ITKHelpers::SetImageToConstant(this->ResultingSegments.GetPointer(),
                                   ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND);

}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CutGraph()
{
  // Compute mininum cut
  std::cout << "CutGraph()..." << std::endl;
  boost::graph_traits<GraphType>::vertex_descriptor s = vertex(this->SourceNodeId, this->Graph);
  boost::graph_traits<GraphType>::vertex_descriptor t = vertex(this->SinkNodeId, this->Graph);

  std::vector<int> groups(num_vertices(this->Graph));
  std::vector<float> residual_capacity(num_edges(this->Graph)); //this needs to be initialized to 0

  boykov_kolmogorov_max_flow(this->Graph,
          boost::make_iterator_property_map(&EdgeWeights[0], get(boost::edge_index, this->Graph)),
          boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, this->Graph)),
          boost::make_iterator_property_map(&ReverseEdges[0], get(boost::edge_index, this->Graph)),
          boost::make_iterator_property_map(&groups[0], get(boost::vertex_index, this->Graph)),
          get(boost::vertex_index, this->Graph),
          s,
          t);

  std::cout << "Finished max_flow()." << std::endl;

  // Iterate over the node image, querying the graph object for the association of each pixel and storing them as the output mask
  itk::ImageRegionConstIterator<NodeImageType>
      nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  // From the documentation:
  // http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/boykov_kolmogorov_max_flow.html
  // If the color of a vertex after running the algorithm is black the vertex belongs to
  // the source tree else it belongs to the sink-tree (used for minimum cuts).
  while(!nodeImageIterator.IsAtEnd())
  {
    if(groups[nodeImageIterator.Get()] == groups[this->SourceNodeId])
    {
      this->ResultingSegments->SetPixel(nodeImageIterator.GetIndex(),
                                        ForegroundBackgroundSegmentMaskPixelTypeEnum::FOREGROUND);
    }
    else
    {
      this->ResultingSegments->SetPixel(nodeImageIterator.GetIndex(),
                                        ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND);
    }
    ++nodeImageIterator;
  }

  std::cout << "Finished CutGraph()." << std::endl;
}

// This function assumes that the ReverseEdges and EdgeWeights members are already large enough to accept the
// new data (otherwise we would have to push_back/resize millions of times).
template <typename TImage, typename TPixelDifferenceFunctor>
unsigned int ImageGraphCut<TImage, TPixelDifferenceFunctor>::
AddBidirectionalEdge(unsigned int numberOfEdges, const unsigned int source, const unsigned int target, const float weight)
{
    // skip "edge already exists check" if the node ids are too high (because it segfaults)
    if(source < num_vertices(this->Graph)+1 && target < num_vertices(this->Graph)+1)
    {
        // If the edge already exists, do nothing
        //auto edgeCheck = boost::edge(source, target, this->Graph);
//        std::cout << "Checking edge..." << std::endl;
        std::pair<EdgeDescriptor, bool> edgeCheck = boost::edge(source, target, this->Graph);
//        std::cout << "Edge check complete." << std::endl;
        if(edgeCheck.second == true)
        {
            return numberOfEdges;
        }
    }

    // Add edges between grid vertices. We have to create the edge and the reverse edge,
    // then add the reverseEdge as the corresponding reverse edge to 'edge', and then add 'edge'
    // as the corresponding reverse edge to 'reverseEdge'
//    int nextEdgeId = num_edges(this->Graph); // Calling this every time is VERY slow
    int nextEdgeId = numberOfEdges;

    EdgeDescriptor edge;
    bool inserted;

    boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, this->Graph);
    if(!inserted)
    {
        std::cerr << "Not inserted!" << std::endl;
        return numberOfEdges;
    }

    EdgeDescriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, this->Graph).first;

    this->ReverseEdges[nextEdgeId] = reverseEdge;
    this->ReverseEdges[nextEdgeId + 1] = edge;
    this->EdgeWeights[nextEdgeId] = weight;
    this->EdgeWeights[nextEdgeId + 1] = weight;

    return numberOfEdges + 2;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::PerformSegmentation()
{
  // This function performs some initializations and then creates and cuts the graph

  this->Initialize();

  this->CreateGraph();
  this->CutGraph();
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateSamples()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms
  std::cout << "CreateSamples()" << std::endl;

  // Ensure at least one pixel has been specified for both the foreground and background
  std::cout << "Currently there are " << this->Sources.size() << " sources and "
            << this->Sinks.size() << " sinks." << std::endl;
  if((this->Sources.size() <= 0) || (this->Sinks.size() <= 0))
  {
    std::cerr << "At least one source (foreground) pixel and one sink (background) "
                 "pixel must be specified!" << std::endl;
    return;
  }

  unsigned int numberOfComponentsPerPixel =
      this->Image->GetNumberOfComponentsPerPixel();

  // We want the histogram bins to take values from 0 to 255 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(numberOfComponentsPerPixel);
  HistogramType::MeasurementVectorType binMaximum(numberOfComponentsPerPixel);
  for(unsigned int i = 0; i < numberOfComponentsPerPixel; i++)
  {
    // If one channel (often, the alpha channel) is all 255, for example, then however
    // the ITK histogram constructs the bins with range (0,255) does not include the 255 pixels,
    // and therefore none of the pixels are included in the histogram at all!
    // Using slightly less than 0 (-.5) and slightly more than 255 (255.5) as the bin limits
    // fixes this problem
//    binMinimum[i] = 0;
//    binMaximum[i] = 255;
    binMinimum[i] = -0.5f;
    binMaximum[i] = 255.5f;
  }

  // Setup the histogram size
//  std::cout << "Image components per pixel: "
//            << numberOfComponentsPerPixel << std::endl;
  typename SampleToHistogramFilterType::HistogramSizeType
      histogramSize(numberOfComponentsPerPixel);
  histogramSize.Fill(this->NumberOfHistogramBins);

  // Create foreground samples and histogram
  this->ForegroundSample->Clear();
  this->ForegroundSample->SetMeasurementVectorSize(numberOfComponentsPerPixel);
  //std::cout << "Measurement vector size: " << this->ForegroundSample->GetMeasurementVectorSize() << std::endl;
  //std::cout << "Pixel size: " << this->Image->GetPixel(this->Sources[0]).GetNumberOfElements() << std::endl;
  
  for(unsigned int i = 0; i < this->Sources.size(); i++)
  {
    this->ForegroundSample->PushBack(this->Image->GetPixel(this->Sources[i]));
  }

  this->ForegroundHistogramFilter->SetHistogramSize(histogramSize);
  this->ForegroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  this->ForegroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  this->ForegroundHistogramFilter->SetAutoMinimumMaximum(false);
  this->ForegroundHistogramFilter->SetInput(this->ForegroundSample);
  this->ForegroundHistogramFilter->Modified();
  this->ForegroundHistogramFilter->Update();

  this->ForegroundHistogram = this->ForegroundHistogramFilter->GetOutput();

  // Create background samples and histogram
  this->BackgroundSample->Clear();
  this->BackgroundSample->SetMeasurementVectorSize(numberOfComponentsPerPixel);
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
  {
    this->BackgroundSample->PushBack(this->Image->GetPixel(this->Sinks[i]));
  }

  this->BackgroundHistogramFilter->SetHistogramSize(histogramSize);
  this->BackgroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  this->BackgroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  this->BackgroundHistogramFilter->SetAutoMinimumMaximum(false);
  this->BackgroundHistogramFilter->SetInput(this->BackgroundSample);
  this->BackgroundHistogramFilter->Modified();
  this->BackgroundHistogramFilter->Update();

  this->BackgroundHistogram = BackgroundHistogramFilter->GetOutput();

  std::cout << "Finished CreateSamples()" << std::endl;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateNEdges()
{
  std::cout << "CreateNEdges()" << std::endl;
  // Create n-edges and set n-edge weights (links between image nodes)

  itk::Size<2> imageSize = this->Image->GetLargestPossibleRegion().GetSize();

  unsigned int expectedNumberOfNEdges = 2*(
                                          imageSize[0] * (imageSize[1] - 1) + // vertical edges
                                          (imageSize[0]-1) * imageSize[1]  // horizontal edges
                                          );
  std::cout << "Resizing for " << expectedNumberOfNEdges << " N-edges." << std::endl;
  this->EdgeWeights.resize(num_edges(this->Graph) + expectedNumberOfNEdges);
  this->ReverseEdges.resize(num_edges(this->Graph) + expectedNumberOfNEdges);

  // We are only using a 4-connected structure,
  // so the kernel (iteration neighborhood) must only be
  // 3x3 (specified by a radius of 1)
  itk::Size<2> radius;
  radius.Fill(1);

  typedef itk::ShapedNeighborhoodIterator<TImage> IteratorType;

  // Traverse the image adding an edge between the current pixel
  // and the pixel below it and the current pixel and the pixel to the right of it.
  // This prevents duplicate edges (i.e. we cannot add an edge to
  // all 4-connected neighbors of every pixel or almost every edge would be duplicated.
  std::vector<typename IteratorType::OffsetType> neighbors;
  typename IteratorType::OffsetType bottom = {{0,1}};
  neighbors.push_back(bottom);
  typename IteratorType::OffsetType right = {{1,0}};
  neighbors.push_back(right);

  typename IteratorType::OffsetType center = {{0,0}};

  IteratorType iterator(radius, this->Image, this->Image->GetLargestPossibleRegion());
  iterator.ClearActiveList();
  iterator.ActivateOffset(bottom);
  iterator.ActivateOffset(right);
  iterator.ActivateOffset(center);

  // Estimate the "camera noise"
  double sigma = this->ComputeNoise();

  unsigned int currentNumberOfEdges = num_edges(this->Graph);
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
  {
    PixelType centerPixel = iterator.GetPixel(center);

    for(unsigned int i = 0; i < neighbors.size(); i++)
    {
      bool valid;
      iterator.GetPixel(neighbors[i], valid);

      // If the current neighbor is outside the image, skip it
      if(!valid)
      {
        continue;
      }

      PixelType neighborPixel = iterator.GetPixel(neighbors[i]);

      // Compute the Euclidean distance between the pixel intensities
      float pixelDifference = PixelDifferenceFunctor.Difference(centerPixel, neighborPixel);

      // Compute the edge weight
      float weight = exp(-pow(pixelDifference,2)/(2.0*sigma*sigma));
      assert(weight >= 0);

      // Add the edge to the graph
      unsigned int node1 = this->NodeImage->GetPixel(iterator.GetIndex(center));
      unsigned int node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));

      currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, node1, node2, weight);
    }
  }

  std::cout << "Finished CreateNEdges()" << std::endl;
}

template <typename TImage, typename TPixelDifferenceFunctor>
float ImageGraphCut<TImage, TPixelDifferenceFunctor>::InternalForegroundLikelihood(const PixelType& pixel)
{
    HistogramType::MeasurementVectorType measurementVector(pixel.Size());
    for(unsigned int i = 0; i < pixel.Size(); i++)
    {
      measurementVector[i] = pixel[i];
    }

    HistogramType::IndexType foregroundIndex;
    this->ForegroundHistogram->GetIndex(measurementVector, foregroundIndex);
    float sourceHistogramValue =
        this->ForegroundHistogram->GetFrequency(foregroundIndex);

    sourceHistogramValue /= this->ForegroundHistogram->GetTotalFrequency();

    return sourceHistogramValue;
}

template <typename TImage, typename TPixelDifferenceFunctor>
float ImageGraphCut<TImage, TPixelDifferenceFunctor>::InternalBackgroundLikelihood(const PixelType& pixel)
{
    HistogramType::MeasurementVectorType measurementVector(pixel.Size());
    for(unsigned int i = 0; i < pixel.Size(); i++)
    {
      measurementVector[i] = pixel[i];
    }

    HistogramType::IndexType backgroundIndex;
    this->BackgroundHistogram->GetIndex(measurementVector, backgroundIndex);
    float sinkHistogramValue =
        this->BackgroundHistogram->GetFrequency(backgroundIndex);

    sinkHistogramValue /= this->BackgroundHistogram->GetTotalFrequency();

    return sinkHistogramValue;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateTEdges()
{
  // Setup links from pixel nodes to terminal nodes
  std::cout << "CreateTEdges()" << std::endl;

  itk::Size<2> imageSize = this->Image->GetLargestPossibleRegion().GetSize();

  unsigned int expectedNumberOfTEdges = 2*2*(imageSize[0] * imageSize[1]);

  this->EdgeWeights.resize(num_edges(this->Graph) + expectedNumberOfTEdges);
  this->ReverseEdges.resize(num_edges(this->Graph) + expectedNumberOfTEdges);

  // Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node)

  // Compute the histograms of the selected foreground and background pixels
  if(!this->CustomLikelihood)
  {
    CreateSamples();
  }

  itk::ImageRegionIteratorWithIndex<TImage>
      imageIterator(this->Image,
                    this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NodeImageType>
      nodeIterator(this->NodeImage,
                   this->NodeImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  nodeIterator.GoToBegin();

  // Since the t-weight function takes the log of the histogram value,
  // we must handle bins with frequency = 0 specially (because log(0) = -inf)
  // For empty histogram bins we use tinyValue instead of 0.
  float tinyValue = 1e-10;

  unsigned int currentNumberOfEdges = num_edges(this->Graph);

  while(!imageIterator.IsAtEnd())
  {
    itk::Index<2> currentIndex = imageIterator.GetIndex();

    // Skip the computation and edge creation if the current pixel already has a fixed assignment
    if(Helpers::Contains(this->Sinks, currentIndex) || Helpers::Contains(this->Sources, currentIndex))
    {
        ++imageIterator;
        ++nodeIterator;
        continue;
    }
    PixelType pixel = imageIterator.Get();
    //std::cout << "Pixels have size: " << pixel.Size() << std::endl;

    HistogramType::MeasurementVectorType measurementVector(pixel.Size());
    for(unsigned int i = 0; i < pixel.Size(); i++)
    {
      measurementVector[i] = pixel[i];
    }

    float sourceLikelihood = ForegroundLikelihood(pixel);
    float sinkLikelihood = BackgroundLikelihood(pixel);

    if(sourceLikelihood <= 0)
    {
      sourceLikelihood = tinyValue;
    }

    if(sinkLikelihood <= 0)
    {
      sinkLikelihood = tinyValue;
    }

    // Add the edge to the graph and set its weight
    // log() is the natural log
    currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, nodeIterator.Get(),
                                                this->SinkNodeId, -this->Lambda*log(sourceLikelihood));
    currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, nodeIterator.Get(),
                                                this->SourceNodeId, -this->Lambda*log(sinkLikelihood));

    ++imageIterator;
    ++nodeIterator;
  }

  // Set very high source weights for the pixels that were
  // selected as foreground by the user
  for(unsigned int i = 0; i < this->Sources.size(); i++)
  {
    currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sources[i]),
                                                this->SourceNodeId,  std::numeric_limits<float>::max());

    currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sources[i]),
                                                this->SinkNodeId, 0);
  }

  // Set very high sink weights for the pixels that
  // were selected as background by the user
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
  {
      currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sinks[i]),
                                                  this->SourceNodeId, 0);

      currentNumberOfEdges = AddBidirectionalEdge(currentNumberOfEdges, this->NodeImage->GetPixel(this->Sinks[i]),
                                                  this->SinkNodeId, std::numeric_limits<float>::max());
  }

  std::cout << "Finished CreateTEdges()" << std::endl;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateGraph()
{
  std::cout << "CreateGraph()" << std::endl;

  // Add all of the nodes to the graph and store their IDs in a "node image"
  itk::ImageRegionIterator<NodeImageType> nodeImageIterator(this->NodeImage,
                                                            this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  unsigned int nodeId = 0;
  while(!nodeImageIterator.IsAtEnd())
  {
    nodeImageIterator.Set(nodeId);
    nodeId++;
    ++nodeImageIterator;
  }

  // Set the sink and source ids to be the two numbers immediately following the number of vertices in the grid
  this->SinkNodeId = nodeId;
  nodeId++;
  this->SourceNodeId = nodeId;

  CreateNEdges();
  CreateTEdges();

  {
  itk::Size<2> imageSize = this->NodeImage->GetLargestPossibleRegion().GetSize();

  std::cout << "Number of edges " << num_edges(this->Graph) << std::endl;
  int expectedEdges = imageSize[0]*imageSize[1] * 2 * 2 + // one '2' is because there is an edge to both the source and sink from each pixel, and the other '2' is because they are double edges (bidirectional)
                      2*(imageSize[0]-1)*imageSize[1] + 2*imageSize[0]*(imageSize[1]-1); // both '2's are for the double edges (this is the number of horizontal edges + the number of vertical edges)
  std::cout << "(Should be " << expectedEdges << " edges.)" << std::endl;
  }
}

template <typename TImage, typename TPixelDifferenceFunctor>
double ImageGraphCut<TImage, TPixelDifferenceFunctor>::ComputeNoise()
{
  // Compute an estimate of the "camera noise". This is used in the N-weight function.

  // Since we use a 4-connected neighborhood, the kernel must be 3x3 (a rectangular radius of 1 creates a kernel side length of 3)
  itk::Size<2> radius;
  radius[0] = 1;
  radius[1] = 1;

  typedef itk::ShapedNeighborhoodIterator<TImage> IteratorType;

  std::vector<typename IteratorType::OffsetType> neighbors;
  typename IteratorType::OffsetType bottom = {{0,1}};
  neighbors.push_back(bottom);
  typename IteratorType::OffsetType right = {{1,0}};
  neighbors.push_back(right);

  typename IteratorType::OffsetType center = {{0,0}};

  IteratorType iterator(radius, this->Image, this->Image->GetLargestPossibleRegion());
  iterator.ClearActiveList();
  iterator.ActivateOffset(bottom);
  iterator.ActivateOffset(right);
  iterator.ActivateOffset(center);

  double sigma = 0.0;
  int numberOfEdges = 0;

  // Traverse the image collecting the differences between neighboring pixel intensities
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
  {
    PixelType centerPixel = iterator.GetPixel(center);

    for(unsigned int i = 0; i < neighbors.size(); i++)
    {
      bool valid;
      iterator.GetPixel(neighbors[i], valid);
      if(!valid)
      {
        continue;
      }

      PixelType neighborPixel = iterator.GetPixel(neighbors[i]);

      float colorDifference = PixelDifferenceFunctor.Difference(centerPixel, neighborPixel);
      sigma += colorDifference;
      numberOfEdges++;
    }
  }

  // Normalize
  sigma /= static_cast<double>(numberOfEdges);

  return sigma;
}

template <typename TImage, typename TPixelDifferenceFunctor>
typename ImageGraphCut<TImage, TPixelDifferenceFunctor>::IndexContainer ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSources()
{
  return this->Sources;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetLambda(const float lambda)
{
  this->Lambda = lambda;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetNumberOfHistogramBins(int bins)
{
  this->NumberOfHistogramBins = bins;
}

template <typename TImage, typename TPixelDifferenceFunctor>
ForegroundBackgroundSegmentMask* ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSegmentMask()
{
  return this->ResultingSegments;
}

template <typename TImage, typename TPixelDifferenceFunctor>
typename ImageGraphCut<TImage, TPixelDifferenceFunctor>::IndexContainer ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetSinks()
{
  return this->Sinks;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetSources(const IndexContainer& sources)
{
  this->Sources = sources;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetSinks(const IndexContainer& sinks)
{
  this->Sinks = sinks;
}

template <typename TImage, typename TPixelDifferenceFunctor>
TImage* ImageGraphCut<TImage, TPixelDifferenceFunctor>::GetImage()
{
  return this->Image;
}

#endif
