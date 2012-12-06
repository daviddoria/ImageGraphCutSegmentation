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

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::SetImage(TImage* const image)
{
  this->Image = TImage::New();
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());

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
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CutGraph()
{
  // Compute max-flow
  this->Graph->maxflow();

  // Iterate over the node image, querying the Kolmorogov graph object for the association of each pixel and storing them as the output mask
  itk::ImageRegionConstIterator<NodeImageType>
      nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
  {
    if(this->Graph->what_segment(nodeImageIterator.Get()) == GraphType::SOURCE)
    {
      this->ResultingSegments->SetPixel(nodeImageIterator.GetIndex(), ForegroundBackgroundSegmentMaskPixelTypeEnum::FOREGROUND);
    }
    else if(this->Graph->what_segment(nodeImageIterator.Get()) == GraphType::SINK)
    {
      this->ResultingSegments->SetPixel(nodeImageIterator.GetIndex(), ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND);
    }
    ++nodeImageIterator;
  }

  delete this->Graph;
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::PerformSegmentation()
{
  // This function performs some initializations and then creates and cuts the graph

  // Ensure at least one pixel has been specified for both the foreground and background
  if((this->Sources.size() <= 0) || (this->Sinks.size() <= 0))
  {
    std::cerr << "At least one source (foreground) pixel and one sink (background) "
                 "pixel must be specified!" << std::endl;
    std::cerr << "Currently there are " << this->Sources.size()
              << " and " << this->Sinks.size() << " sinks." << std::endl;
    return;
  }

  // Blank the NodeImage
  itk::ImageRegionIterator<NodeImageType>
      nodeImageIterator(this->NodeImage,
                        this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
  {
    nodeImageIterator.Set(nullptr);
    ++nodeImageIterator;
  }

  // Blank the output image
  ITKHelpers::SetImageToConstant(this->ResultingSegments.GetPointer(),
                                 ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND);

  this->CreateGraph();
  this->CutGraph();
}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateSamples()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms

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
  std::cout << "Image components per pixel: "
            << numberOfComponentsPerPixel << std::endl;
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

}

template <typename TImage, typename TPixelDifferenceFunctor>
void ImageGraphCut<TImage, TPixelDifferenceFunctor>::CreateGraph()
{
  // Form the graph
  this->Graph = new GraphType;

  // Add all of the nodes to the graph and store their IDs in a "node image"
  itk::ImageRegionIterator<NodeImageType> nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
    {
    nodeImageIterator.Set(this->Graph->add_node());
    ++nodeImageIterator;
    }

  // Estimate the "camera noise"
  double sigma = this->ComputeNoise();

  ////////// Create n-edges and set n-edge weights
  ////////// (links between image nodes)

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
      void* node1 = this->NodeImage->GetPixel(iterator.GetIndex(center));
      void* node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
      this->Graph->add_edge(node1, node2, weight, weight);
      }
    }

  ////////// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node) //////////

  // Compute the histograms of the selected foreground and background pixels
  CreateSamples();

  itk::ImageRegionIterator<TImage>
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

  while(!imageIterator.IsAtEnd())
  {
    PixelType pixel = imageIterator.Get();
    //std::cout << "Pixels have size: " << pixel.Size() << std::endl;

    HistogramType::MeasurementVectorType measurementVector(pixel.Size());
    for(unsigned int i = 0; i < pixel.Size(); i++)
    {
      measurementVector[i] = pixel[i];
    }

    HistogramType::IndexType backgroundIndex;
    this->BackgroundHistogram->GetIndex(measurementVector, backgroundIndex);
    float sinkHistogramValue =
        this->BackgroundHistogram->GetFrequency(backgroundIndex);

    HistogramType::IndexType foregroundIndex;
    this->ForegroundHistogram->GetIndex(measurementVector, foregroundIndex);
    float sourceHistogramValue =
        this->ForegroundHistogram->GetFrequency(foregroundIndex);

    // Conver the histogram value/frequency to make it as if it came from a normalized histogram
    if(this->BackgroundHistogram->GetTotalFrequency() == 0 ||
       this->ForegroundHistogram->GetTotalFrequency() == 0)
    {
      throw std::runtime_error("The foreground or background histogram TotalFrequency is 0!");
    }

    sinkHistogramValue /= this->BackgroundHistogram->GetTotalFrequency();
    sourceHistogramValue /= this->ForegroundHistogram->GetTotalFrequency();

    if(sinkHistogramValue <= 0)
    {
      sinkHistogramValue = tinyValue;
    }
    if(sourceHistogramValue <= 0)
    {
      sourceHistogramValue = tinyValue;
    }

    // Add the edge to the graph and set its weight
    // log() is the natural log
    this->Graph->add_tweights(nodeIterator.Get(),
                              -this->Lambda*log(sinkHistogramValue),
                              -this->Lambda*log(sourceHistogramValue));
    ++imageIterator;
    ++nodeIterator;
  }

  // Set very high source weights for the pixels that were
  // selected as foreground by the user
  for(unsigned int i = 0; i < this->Sources.size(); i++)
  {
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sources[i]),
                              this->Lambda * std::numeric_limits<float>::max(),0);
  }

  // Set very high sink weights for the pixels that
  // were selected as background by the user
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
  {
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sinks[i]),
                              0,this->Lambda * std::numeric_limits<float>::max());
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
